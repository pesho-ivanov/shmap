#pragma once

#include <set>
#include <map>

#include "index.h"
#include "io.h"
#include "sketch.h"
#include "utils.h"
#include "handler.h"
#include "mapper.h"
#include "cmath"

#include "edlib.h"
#include <ankerl/unordered_dense.h>
#include "doctest.h"

namespace sweepmap {

enum class Metric {
	JACCARD,
	CONTAINMENT_INDEX
};

class SHMapper : public Mapper {
	const SketchIndex &tidx;
	Handler *H;

public:
	typedef vector<Seed> Seeds;
	using Matches = vector<Match>;
	using Bucket = pair<bucket_t, qpos_t>;
	//using Buckets = std::unordered_map<bucket_t, qpos_t>;
	//using Buckets = gtl::flat_hash_map<bucket_t, qpos_t>;
	using Buckets = ankerl::unordered_dense::map<bucket_t, qpos_t>;
	//using Hist = unordered_map<hash_t, qpos_t>;
	//using Hist = gtl::flat_hash_map<hash_t, qpos_t>;
	using Hist = ankerl::unordered_dense::map<hash_t, qpos_t>;

	Seeds select_kmers(sketch_t& p, int &nonzero) {
		H->T.start("seeding");
		H->T.start("collect_kmer_info");
			Seeds kmers;
			sort(p.begin(), p.end(), [](const Kmer &a, const Kmer &b) { return a.h < b.h; });
			qpos_t strike = 0;
			for (size_t ppos = 0; ppos < p.size(); ++ppos) {
				++strike;
				if (ppos == p.size()-1 || p[ppos].h != p[ppos+1].h) {
					rpos_t hits_in_t = tidx.count(p[ppos].h);
					if (hits_in_t > 0) {
						nonzero += strike;
					}
						if (H->params.max_matches == -1 || hits_in_t <= H->params.max_matches) {
							Seed el(p[ppos], hits_in_t, kmers.size());
							//strike = 1; // comment out for Weighted metric
							el.occs_in_p = strike;
							kmers.push_back(el);
							if (H->params.max_seeds != -1 && (rpos_t)kmers.size() >= H->params.max_seeds)  // TODO maybe account for occs_in_p
								break;
						}
					//}
					strike = 0;
				}
			}
		H->T.stop("collect_kmer_info");

		H->T.start("sort_kmers");
			//pdqsort_branchless(kmers.begin(), kmers.end(), [](const Seed &a, const Seed &b) {
			sort(kmers.begin(), kmers.end(), [](const Seed &a, const Seed &b) {
				return a.hits_in_T < b.hits_in_T;
			});
		H->T.stop("sort_kmers");
		H->T.stop("seeding");

		return kmers;
	}

private:

	bool do_overlap(const Mapping &a, const Mapping &b) {
		if (a.segm_id != b.segm_id)
			return false;
		return a.T_r >= b.T_l && a.T_l <= b.T_r;
	}

	double overlap(const Mapping &a, const Mapping &b) {
		if (a.segm_id != b.segm_id)
			return 0.0;
		int cap = std::max(0, std::min(a.T_r, b.T_r) - std::max(a.T_l, b.T_l));
		int cup = std::max(a.T_r, b.T_r) - std::min(a.T_l, b.T_l);
		assert(cup >= 0 && cap >=0 && cup >= cap);
		return 1.0 * cap / cup;
	}

	double ContainmentIndex(qpos_t intersection, qpos_t m) {
		assert(intersection <= s_sz);
		assert(intersection <= m);
		assert(1.0*intersection / m <= 1.0);
		return 1.0*intersection / m;
	}

	double Jaccard(qpos_t intersection, qpos_t m, qpos_t s_sz) {
		assert(intersection <= s_sz);
		assert(intersection <= m);
		assert(1.0*intersection / (m + s_sz - intersection) <= 1.0);
		return 1.0*intersection / (m + s_sz - intersection);
	}

	Mapping bestIncludedJaccard(const vector<Match> &M, const qpos_t P_sz, qpos_t lmin, qpos_t lmax, const qpos_t m, const bucket_t &bucket, Hist diff_hist) {
		qpos_t intersection = 0;
		qpos_t same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise
		Mapping best;
		
		for (int i=1; i<(int)M.size(); i++)
			assert(M[i-1].hit.tpos < M[i].hit.tpos);
		assert(M.size() == 0 || (int)M.size() <= M.back().hit.tpos - M.front().hit.tpos + 1);

		auto l = M.begin(), r = M.begin();
		for(; l != M.end(); ++l) {
			// move r to the left until l+lmin<r
			for(; l + 1 < r && l->hit.tpos + lmin < r->hit.tpos; ) {
				--r;
				if (++diff_hist[r->seed.kmer.h] >= 1) {
					--intersection;
					same_strand_seeds -= r->is_same_strand() ? +1 : -1;
				}
			}

			// move r to the right until l+lmax<=r
			for(;  r != M.end() && r->hit.tpos <= l->hit.tpos + lmax; ++r) {
				assert(diff_hist.contains(r->seed.kmer.h));
				if (--diff_hist[r->seed.kmer.h] >= 0) {
					++intersection;
					same_strand_seeds += r->is_same_strand() ? +1 : -1;
				}
				
				assert (l->hit.r <= r->hit.r);
				if (l < r) {
					int s_sz = r->hit.tpos - l->hit.tpos + 1;
					double J = Jaccard(intersection, m, s_sz);
					if (J > best.J)
						best = Mapping(H->params.k, P_sz, m, l->hit.r, r->hit.r, l->hit.segm_id, intersection, J, same_strand_seeds, l, r, bucket);
				}
			}

			assert(diff_hist.contains(l->seed.kmer.h));
			if (++diff_hist[l->seed.kmer.h] >= 1) {
				--intersection;
				same_strand_seeds -= l->is_same_strand() ? +1 : -1;
			}

			assert(intersection >= 0);
		}
		assert(l == M.end());
		assert(r == M.end());
		assert(intersection == 0);
		
		return best;
	}

	Mapping bestFixedLength(const vector<Match> &M, const qpos_t P_sz, qpos_t lmax, const qpos_t m, const bucket_t &bucket, Hist &diff_hist, Metric metric) {
		qpos_t intersection = 0;
		qpos_t same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise
		Mapping best;
		
		for (int i=1; i<(int)M.size(); i++)
			assert(M[i-1].hit.tpos < M[i].hit.tpos);
		assert(M.size() == 0 || (int)M.size() <= M.back().hit.tpos - M.front().hit.tpos + 1);

		auto l = M.begin(), r = M.begin();
		for(; l != M.end(); ++l) {
			for(;  r != M.end()  && r->hit.tpos <= l->hit.tpos + lmax; ++r) {
				//&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				//&& r->hit.tpos + H->params.k <= l->hit.tpos + P_sz
				//&& r->hit.r + H->params.k <= l->hit.r + P_sz
				assert(diff_hist.contains(r->seed.kmer.h));
				if (--diff_hist[r->seed.kmer.h] >= 0) {
					++intersection;
					same_strand_seeds += r->is_same_strand() ? +1 : -1;
				}
				
				assert (l->hit.r <= r->hit.r);
			}

			double J;
			if (metric == Metric::JACCARD) {
				int s_sz = prev(r)->hit.tpos - l->hit.tpos + 1;
				J = Jaccard(intersection, m, s_sz);
			} else if (metric == Metric::CONTAINMENT_INDEX) {
				J = ContainmentIndex(intersection, m);
			} else {
				J = -1.0;
				assert(false);
			}

			assert(J >= -0.0);
			assert(J <= 1.0);
			if (J > best.J)
				best = Mapping(H->params.k, P_sz, m, l->hit.r, prev(r)->hit.r, l->hit.segm_id, intersection, J, same_strand_seeds, l, prev(r), bucket);

			assert(diff_hist.contains(l->seed.kmer.h));
			if (++diff_hist[l->seed.kmer.h] >= 1) {
				--intersection;
				same_strand_seeds -= l->is_same_strand() ? +1 : -1;
			}

			assert(intersection >= 0);
		}
		assert(l == M.end());
		assert(r == M.end());
		assert(intersection == 0);
		
		return best;
	}

  public:
	SHMapper(const SketchIndex &tidx, Handler *H) : tidx(tidx), H(H) {
        H->C.inc("seeds_limit_reached", 0);
        H->C.inc("mapped_reads", 0);
        if (H->params.theta < 0.0 || H->params.theta > 1.0) {
            cerr << "tThres = " << H->params.theta << " outside of [0,1]." << endl;
            exit(1);
        }
    }
	virtual ~SHMapper() = default;

	double hseed(qpos_t p, qpos_t seeds, qpos_t matches) {
		return 1.0 - double(seeds - matches) / p;
	}

	bool seed_heuristic_pass(const vector<Mapping> &maps, const Seeds &kmers, qpos_t bucket_l, qpos_t m, bucket_t b, rpos_t max_matches, qpos_t i, qpos_t seeds,
			int best_idx, double *lowest_sh, double thr_init) {
		if (H->params.no_bucket_pruning)
			return true;

		*lowest_sh = 1.0; // should only get lower

		double thr2 = best_idx == -1 ? thr_init : maps[best_idx].J-0.015;

		if (hseed(m, seeds, max_matches) < thr2)
			return false;

		H->T.start("seed_heuristic");
		bool ret = true;
		for (; i < (qpos_t)kmers.size(); i++) {
			seeds += kmers[i].occs_in_p;
			if (tidx.matches_in_bucket(kmers[i], b, bucket_l))
				max_matches += kmers[i].occs_in_p;
			double sh = hseed(m, seeds, max_matches);
			assert(sh <= *lowest_sh);
			*lowest_sh = min(*lowest_sh, sh);	
			if (sh < thr2) {
				ret = false;
				break;
			}
		}
		H->T.stop("seed_heuristic");
		return ret;
	}
	
	double sigmas_diff(qpos_t X, qpos_t Y) {
		if (X==-1) X=0;
		if (Y==-1) Y=0;
		return std::abs(X - Y) / std::sqrt(X + Y);
	}

	int mapq(Mapping m, double new_theta) {
		// minimap2: mapQ = 40 (1-f2/f1) min(1, m/10) log f1, where m is #anchors on primary chain
		assert(m.J >= 0.0);
		if (m.J2 < new_theta)
			m.J2 = new_theta;
		if (m.J - m.J2 > 0.015) {
			if (abs(m.same_strand_seeds) < m.intersection/2)
				return 5;
			return 60;
		} else {
			return 0;
		}
	}

	Matches collect_matches(const bucket_t &b, qpos_t bucket_l, const SketchIndex &tidx, const unordered_map<hash_t, Seed> &p_ht) {
		Matches M;
		for (rpos_t i = b.b*bucket_l; i < std::min((b.b+2)*bucket_l, (rpos_t)tidx.T[b.segm_id].kmers.size()); i++) {
			assert(b.segm_id >= 0 && b.segm_id < (rpos_t)tidx.T.size());
			const auto &kmer = tidx.T[b.segm_id].kmers[i];
			const auto seed_it = p_ht.find(kmer.h);
			if (seed_it != p_ht.end()) {
				M.push_back(Match(seed_it->second, Hit(kmer, i, b.segm_id)));
				assert(M[ M.size()-1 ].hit.tpos == M[ M.size()-1 ].hit.tpos);
			}
		}
		return M;
	}

	void match_rest(qpos_t seed_matches, qpos_t P_sz, qpos_t lmax, qpos_t bucket_l, qpos_t m, const Seeds &kmers, const vector<Bucket> &B_vec, Hist &diff_hist, int seeds, int first_kmer_after_seeds, const unordered_map<hash_t, Seed> &p_ht,
			vector<Mapping> &maps, int &total_matches, int &best_idx, int &final_buckets, double thr_init, int forbidden_idx, const string &query_id, int *lost_on_pruning) {
		double best_J = -1.0;
		total_matches = seed_matches;

		int lost_on_seeding = (0);
		H->C.inc("lost_on_seeding", lost_on_seeding);

		final_buckets = 0;
		H->T.start("seed_heuristic"); H->T.stop("seed_heuristic");  // init
		H->T.start("match_collect"); H->T.stop("match_collect");
		H->T.start("sweep"); H->T.stop("sweep");
		for (auto &[b, seed_matches]: B_vec) {
			double lowest_sh;
			if (seed_heuristic_pass(maps, kmers, bucket_l, m, b, seed_matches, first_kmer_after_seeds, seeds, best_idx, &lowest_sh, thr_init)) {
				if (do_overlap(query_id, b, bucket_l, tidx))
					*lost_on_pruning = 0;
				H->T.start("match_collect");
					Matches M = collect_matches(b, bucket_l, tidx, p_ht);
				H->T.stop("match_collect");
				total_matches += M.size();
				++final_buckets;

				H->T.start("sweep");
					auto bucket_best = bestFixedLength(M, P_sz, lmax, m, b, diff_hist, Metric::CONTAINMENT_INDEX);
					assert(bucket_best.J <= lowest_sh + 1e-7);
					if (bucket_best.J >= H->params.theta) {
						maps.push_back(bucket_best);
						if (bucket_best.J > best_J) {
							if (forbidden_idx == -1 || overlap(maps.back(), maps[forbidden_idx]) < 0.5) {
								best_J = max(best_J, bucket_best.J);
								best_idx = maps.size()-1;
							}
						}
					}
				H->T.stop("sweep");
			}
		}
		//assert(best_idx != -1 || maps.empty());

		H->C.inc("final_buckets", final_buckets);
		H->C.inc("final_mappings", maps.size());
		H->C.inc("total_matches", total_matches);
	}

    bool do_overlap(const string& query_id, const bucket_t& b, int bucket_l, const SketchIndex &tidx) {
        auto parsed = ParsedQueryId::parse(query_id);
        if (!parsed.valid) return false;
		//cerr << "name: " << tidx.T[b.first.segm_id].name<< ", parsed.segm_id: " << parsed.segm_id << endl;
		assert(b.segm_id >= 0 && b.segm_id < (rpos_t)tidx.T.size());
		const auto &segm = tidx.T[b.segm_id];
		if (segm.name == parsed.segm_id) {
			rpos_t bucket_start_t = b.b * bucket_l;
			rpos_t bucket_end_t = bucket_start_t + 2 * bucket_l;
			assert(bucket_start_t >= 0 && bucket_start_t < (rpos_t)segm.kmers.size());
			rpos_t bucket_start_T = segm.kmers[bucket_start_t].r;
			//assert(bucket_end_t >= 0 && bucket_end_t < (rpos_t)segm.kmers.size());
			rpos_t bucket_end_T = bucket_end_t < (rpos_t)segm.kmers.size() ? segm.kmers[bucket_end_t].r : segm.sz;

//			cerr << "bucket_start_t: " << bucket_start_t << ", bucket_end_t: " << bucket_end_t << ", bucket_start_T: " << bucket_start_T << ", bucket_end_T: " << bucket_end_T << endl;
			if (bucket_start_T < parsed.end_pos && bucket_end_T >= parsed.start_pos) {  // if overlap
				return true;
			}
		}
		return false;
	}

    bool lost_correct_mapping(const string& query_id, const Buckets& B, int bucket_l, const SketchIndex &tidx) {
		bool res = true;
        for (const auto& b : B)
			if (do_overlap(query_id, b.first, bucket_l, tidx)) {
				res = false;
				break;
			}
        return res;
    }

 	void gt_Jaccard(const string& query_id, int P_sz, Hist diff_hist, int m, unordered_map<hash_t, Seed> p_ht, const SketchIndex &tidx, int lmax, int bucket_l, Mapping *best_mapping) {
		auto parsed = ParsedQueryId::parse(query_id);
		assert(parsed.valid);
		
		// Find index i where tidx.T[i].seq matches parsed.segm_id
		int segm_id = -1;
		for (int i = 0; i < (int)tidx.T.size(); i++) {
			if (tidx.T[i].name == parsed.segm_id) {
				segm_id = i;
				break;
			}
		}
		assert(segm_id >= 0 && segm_id < (rpos_t)tidx.T.size());
		
		const auto &segm = tidx.T[segm_id];

		qpos_t start = lower_bound(segm.kmers.begin(), segm.kmers.end(), parsed.start_pos, [](const auto &kmer, const auto &pos) { return kmer.r < pos; }) - segm.kmers.begin();
		qpos_t end  = lower_bound(segm.kmers.begin(), segm.kmers.end(), parsed.end_pos, [](const auto &kmer, const auto &pos) { return kmer.r < pos; }) - segm.kmers.begin();
		int gt_lmin = end - start - 1;
		int gt_lmax = end - start + 1;
		//qpos_t lmax = segm.kmers.size();
		//qpos_t lmin = lmax;

		//int bucket_l = lmax;
		assert(bucket_l > gt_lmax);
		bucket_t b(segm_id, start/bucket_l);
		auto M = collect_matches(b, bucket_l, tidx, p_ht);
		//cerr << "lmin: " << lmin << ", lmax: " << lmax << ", start: " << start << ", end: " << end << ", segm.kmers[start].r: " << segm.kmers[start].r << ", segm.kmers[end].r: " << segm.kmers[end].r << ", parsed.start_pos: " << parsed.start_pos << ", parsed.end_pos: " << parsed.end_pos << ", M.front(): " << M.front().hit.tpos << ", M.back(): " << M.back().hit.tpos << endl;	
		
		best_mapping->gt_J 			= bestIncludedJaccard(M, P_sz, gt_lmin, gt_lmax, m, b, diff_hist).J;
		best_mapping->gt_C 			= bestFixedLength(M, P_sz, lmax, m, b, diff_hist, Metric::CONTAINMENT_INDEX).J;
		best_mapping->gt_C_bucket 	= bestFixedLength(M, P_sz, bucket_l, m, b, diff_hist, Metric::CONTAINMENT_INDEX).J;

		//return mapping_J;
	}

	std::tuple<int, int, double, double> calc_FDR(vector<Mapping> &maps, double theta, qpos_t lmax, qpos_t bucket_l, const SketchIndex &tidx, const unordered_map<hash_t, Seed> &p_ht, qpos_t P_sz, qpos_t lmin, qpos_t m, const Hist &diff_hist) {
		int FP = 0;
		for (auto &mapping: maps) {
			auto M = collect_matches(mapping.bucket, bucket_l, tidx, p_ht);
			auto mapping_best_J = bestIncludedJaccard(M, P_sz, lmin, lmax, m, mapping.bucket, diff_hist);
			if (mapping_best_J.J < H->params.theta)
				++FP;
		}

		int PP = maps.size();
		double FDR = 1.0 * FP / PP;
		double FPTP = (PP - FP > 0) ?  1.0 * FP / (PP - FP) : 0.0;
		return {PP, FP, FDR, FPTP};
	}

	void map(const string &pFile) {
		cerr << "Mapping reads using SHmap..." << endl;

		H->C.inc("kmers", 0);
		H->C.inc("kmers_notmatched", 0);
		H->C.inc("seeds", 0);
		H->C.inc("matches", 0);
		H->C.inc("seed_matches", 0);
		H->C.inc("matches_freq", 0);
		H->C.inc("spurious_matches", 0);
		H->C.inc("mappings", 0);
		H->C.inc("J_best", 0);
		H->C.inc("sketched_kmers", 0);
		H->C.inc("total_edit_distance", 0);
		H->C.inc("intersection_diff", 0);
		H->C.inc("mapq60", 0);
		H->C.inc("matches_in_reported_mappings", 0);
		H->C.inc("lost_on_seeding", 0);
		H->C.inc("lost_on_pruning", 0);

		H->T.start("mapping");
		H->T.start("query_reading");

		read_fasta_klib(pFile, [this](const string &query_id, const string &P) {
			H->C.inc("reads");
			H->T.stop("query_reading");
			H->T.start("query_mapping");
				H->T.start("sketching");
					//const char *P = seq->seq.s;
					sketch_t p = H->sketcher.sketch(P);
					H->T.stop("sketching");

				H->T.start("prepare");
					//char *query_id = seq->name.s;
					//qpos_t P_sz = (qpos_t)seq->seq.l;
					qpos_t P_sz = P.size();

					H->C.inc("read_len", P_sz);

					int nonzero = 0;
					Seeds kmers = select_kmers(p, nonzero);
					qpos_t m = 0;
					for (const auto kmer: kmers)
						m += kmer.occs_in_p;
					assert(m <= (int)p.size());
					H->C.inc("kmers_notmatched", m - nonzero);
					//cerr << "notmatched: " << m - nonzero << endl;

					unordered_map<hash_t, Seed> p_ht;
					Hist diff_hist;
					rpos_t possible_matches(0);
					for (const auto kmer: kmers) {
						p_ht.insert(make_pair(kmer.kmer.h, kmer));
						//diff_hist[kmer.kmer.h] = 1;
						diff_hist[kmer.kmer.h] = kmer.occs_in_p;
						possible_matches += kmer.hits_in_T;
					}
					H->C.inc("possible_matches", possible_matches);

					//qpos_t lmax = qpos_t(m / H->params.theta);					// maximum length of a similar mapping
					//qpos_t lmin = qpos_t(ceil(p.size() * H->params.theta));					// maximum length of a similar mapping
					qpos_t lmax = qpos_t(p.size() / H->params.theta);					// maximum length of a similar mapping
					qpos_t bucket_l = lmax;

					double coef = 1.0;// * nonzero / p.size();
					//cerr << "coef: " << coef << endl;
					double new_theta = coef * H->params.theta;
					//cerr << "theta: " << H->params.theta << " -> " << new_theta << endl;

					qpos_t S = qpos_t((1.0 - new_theta) * m) + 1;			// any similar mapping includes at least 1 seed match
					Buckets B;  			// B[segment][b] -- #matched kmers[0...i] in [bl, (b+2)l)
					qpos_t seeds = 0;

					H->C.inc("kmers_sketched", p.size());
					H->C.inc("kmers", m);
					H->C.inc("kmers_unique", kmers.size());
					H->C.inc("kmers_seeds", S);
				H->T.stop("prepare");

				H->T.start("match_seeds");
					int matches_in_B = 0;
					rpos_t seed_matches(0), max_seed_matches(0), seeded_buckets(0);  // stats
					qpos_t i = 0;
					for (; i < (qpos_t)kmers.size() && seeds < S; i++) {
						Seed seed = kmers[i];
						if (seed.hits_in_T > 0) {
							seed_matches += seed.hits_in_T;
							max_seed_matches = max(max_seed_matches, seed.hits_in_T);
							Matches seed_matches;
							tidx.get_matches(&seed_matches, seed);

							//auto B2 = B;	
							//////////////// correct B
							Buckets b2m;
							for (auto &m: seed_matches) {
								rpos_t b(m.hit.tpos / bucket_l);
								//if (b > 0) ++b2m[ bucket_t(m.hit.segm_id, b-1) ];
								//++b2m[ bucket_t(m.hit.segm_id, b) ];
								b2m[ bucket_t(m.hit.segm_id, b) ] = seed.occs_in_p;
								if (b > 0) b2m[ bucket_t(m.hit.segm_id, b-1) ] = seed.occs_in_p;
							}
							for (const auto &[b, matches]: b2m) {
								//B[b] += min(seed.occs_in_p, matches);
								B[b] += matches;
								matches_in_B += matches;
							}

							//bucket_t prev_b(-1, -1);
							//int matches_in_prev_bucket = 0;
							//for (int i=0; i < (int)seed_matches.size();) {
							//	bucket_t curr_b(seed_matches[i].hit.segm_id, seed_matches[i].hit.tpos / bucket_l);
							//	
							//	if (prev_b.segm_id != -1 && (curr_b.segm_id != prev_b.segm_id || curr_b.b != prev_b.b + 1)) {
							//		B2[prev_b] += min(matches_in_prev_bucket, seed.occs_in_p);
							//		matches_in_prev_bucket = 0;
							//	}

							//	int matches_in_curr_bucket = 1;
							//	while (++i < (int)seed_matches.size()
							//			&& seed_matches[i].hit.tpos / lmax == curr_b.b
							//			&& seed_matches[i].hit.segm_id == curr_b.segm_id) {
							//		++matches_in_curr_bucket;
							//	}
							//	prev_b = curr_b;
							//	if (--curr_b.b >= 0) {
							//		B2[curr_b] += min(matches_in_prev_bucket + matches_in_curr_bucket, seed.occs_in_p);
							//	}

							//	matches_in_prev_bucket = matches_in_curr_bucket;
							//}
							//if (prev_b.segm_id != -1)
							//	B2[prev_b] += min(matches_in_prev_bucket, seed.occs_in_p);
						}
						seeds += seed.occs_in_p;
					}
					seeded_buckets = B.size();
					H->C.inc("seeded_buckets", seeded_buckets);
					H->C.inc("seed_matches", seed_matches);
				H->T.stop("match_seeds");

				int lost_on_seeding = lost_correct_mapping(query_id, B, bucket_l, tidx);
				H->C.inc("lost_on_seeding", lost_on_seeding);

				int lost_on_pruning = 1;
				H->T.start("match_rest");
					vector<Mapping> maps;
					vector<Bucket> B_vec(B.begin(), B.end());
					sort(B_vec.begin(), B_vec.end(), [](const Bucket &a, const Bucket &b) { return a.second > b.second; });  // TODO: sort intervals by decreasing number of matches
					int total_matches, best_idx=-1, best2_idx=-1, final_buckets;
					match_rest(seed_matches, P_sz, lmax, bucket_l, m, kmers, B_vec, diff_hist, seeds, i, p_ht, maps, total_matches, best_idx, final_buckets, H->params.theta, -1, query_id, &lost_on_pruning);
					auto [FP, PP, FDR, FPTP] = tuple(-1, -1, -1, -1);
					//auto [FP, PP, FDR, FPTP] = calc_FDR(maps, H->params.theta, lmax, bucket_l, tidx, p_ht, P_sz, lmin, m, diff_hist);
					H->C.inc("FP", FP);
					H->C.inc("PP", PP);
					H->C.inc("FDR", FDR);
					H->C.inc("FPTP", FPTP);

					if (best_idx != -1) {
						// find second best mapping for mapq computation
						match_rest(seed_matches, P_sz, lmax, bucket_l, m, kmers, B_vec, diff_hist, seeds, i, p_ht, maps, total_matches, best2_idx, final_buckets, maps[best_idx].J, best_idx, query_id, &lost_on_pruning);
					}
				H->T.stop("match_rest");
				H->C.inc("lost_on_pruning", lost_on_pruning);
			H->T.stop("query_mapping");

			{
				//if (best_idx == -1) {
				//	maps.push_back(Mapping());
				//	best_idx = maps.size() - 1;
				//}
				// TODO: remove this
				//if (best_idx != -1)
				//	gt_Jaccard(query_id, P_sz, diff_hist, m, p_ht, tidx, lmax, bucket_l, &maps[best_idx]);
			}

			H->T.start("output");
				if (H->params.onlybest && maps.size() >= 1) {
					H->C.inc("mapped_reads");
					if (best_idx != -1) 
					//if (maps[best_idx].J >= H->params.theta)
					{
						assert(0 <= best_idx && best_idx < (int)maps.size());
						Mapping &m = maps[best_idx];

						const auto &segm = tidx.T[m.segm_id];
						m.seeds = S;
						m.total_matches = total_matches;
						m.match_inefficiency = 1.0 * m.total_matches / m.intersection;
						m.max_seed_matches = max_seed_matches;
						m.seed_matches = seed_matches;
						m.seeded_buckets = seeded_buckets;
						m.final_buckets = final_buckets;
						m.FPTP = FPTP;
						m.query_id = query_id.c_str();
						m.segm_name = segm.name.c_str();
						m.segm_sz = segm.sz;
						m.map_time = H->T.secs("query_mapping");

						if (best2_idx != -1) {
							//cerr << "best_idx: " << best_idx << ", best2_idx: " << best2_idx << endl;
							m.J2 = maps[best2_idx].J;
							m.bucket2 = maps[best2_idx].bucket;
							m.intersection2 = maps[best2_idx].intersection;
							m.sigmas_diff = sigmas_diff(m.intersection2, m.intersection);
							//assert(m.bucket.segm_id != m.bucket2.segm_id || abs(m.bucket.b - m.bucket2.b) > 1);
							//m.ed2 = edit_distance(m.bucket, P, P_sz, m, kmers);
						}

						m.mapq = mapq(m, new_theta);
						m.print_paf();

						if (m.mapq == 60) H->C.inc("mapq60");
						H->C.inc("matches_in_reported_mappings", m.intersection);
						H->C.inc("J_best", rpos_t(10000.0*m.J));
						H->C.inc("mappings");
					}
				}
			H->T.stop("output");
			H->T.start("query_reading");
		});
		H->T.stop("query_reading");
		H->T.stop("mapping");

		print_stats();
		print_time_stats();
	}
	
	void print_stats() {
		cerr << std::fixed << std::setprecision(1);
		//cerr << "Mapping:" << endl;
		cerr << " | Total reads:           " << H->C.count("reads") << " (~" << 1.0*H->C.count("read_len") / H->C.count("reads") << " nb p/ read)" << endl;
		cerr << " |  | lost on seeding:      " << H->C.count("lost_on_seeding") << " (" << H->C.perc("lost_on_seeding", "reads") << "%)" << endl;
		cerr << " |  | lost on pruning:      " << H->C.count("lost_on_pruning") << " (" << H->C.perc("lost_on_pruning", "reads") << "%)" << endl;
		cerr << " |  | mapped:               " << H->C.count("mapped_reads") << " (" << H->C.perc("mapped_reads", "reads") << "%)" << endl;
//		cerr << " |  |  | intersect. diff:     " << H->C.frac("intersection_diff", "mapped_reads") << " p/ mapped read" << endl;
		cerr << " | Kmers:                 " << H->C.frac("kmers", "reads") << " p/ read" << endl;
		cerr << " |  | sketched:               " << H->C.frac("kmers_sketched", "reads") << " (" << H->C.perc("kmers_sketched", "kmers") << "%)" << endl;
		cerr << " |  | not matched:            " << H->C.frac("kmers_notmatched", "reads") << " (" << H->C.perc("kmers_notmatched", "kmers") << "%)" << endl;
		cerr << " |  | unique:                 " << H->C.frac("kmers_unique", "reads") << " (" << H->C.perc("kmers_unique", "kmers") << "%)" << endl;
		cerr << " |  | seeds:                  " << H->C.frac("kmers_seeds", "reads") << " (" << H->C.perc("kmers_seeds", "kmers") << "%)" << endl;
		cerr << " | Matches:               " << H->C.frac("total_matches", "reads") << " p/ read" << endl;
		cerr << " |  | seed matches:           " << H->C.frac("seed_matches", "reads") << " (" << H->C.perc("seed_matches", "total_matches") << "%)" << endl;
		cerr << " |  | in reported mappings:   " << H->C.frac("matches_in_reported_mappings", "reads") << " (match inefficiency: " << H->C.frac("total_matches", "matches_in_reported_mappings") << "x)" << endl;
		cerr << " |  | possible matches:       " << H->C.frac("possible_matches", "reads") << " (" <<H->C.frac("possible_matches", "total_matches") << "x)" << endl;
//		cerr << " |  | frequent:               " << H->C.count("matches_freq") << " (" << H->C.perc("matches_freq", "total_matches") << "%)" << endl;
//		cerr << " |  | Seed h. reduction:      " << H->C.frac("possible_matches", "seed_matches") << "x" << endl;
		//cerr << " | Seed limit reached:    " << H->C.count("seeds_limit_reached") << " (" << H->C.perc("seeds_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Matches limit reached: " << H->C.count("matches_limit_reached") << " (" << H->C.perc("matches_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Spurious matches:      " << H->C.count("spurious_matches") << " (" << H->C.perc("spurious_matches", "matches") << "%)" << endl;
		//cerr << " | Discarded seeds:       " << H->C.count("discarded_seeds") << " (" << H->C.perc("discarded_seeds", "collected_seeds") << "%)" << endl;
		cerr << " | Mappings:              " << H->C.count("mappings") << " (" << H->C.perc("mappings", "reads") << "\% of reads)" << endl;
		cerr << " | | Seeded buckets:          " << H->C.frac("seeded_buckets", "mappings") << " /mapping" << endl;
		cerr << " | | Final buckes:            " << H->C.frac("final_buckets", "mappings") << " /mapping" << endl;
		cerr << " | | Final mappings:          " << H->C.frac("final_mappings", "mappings") << " /mapping" << endl;
		cerr << " | | Average best sim.:       " << std::fixed << std::setprecision(3) << H->C.frac("J_best", "mappings") / 10000.0 << endl;
		cerr << " | | mapq=60:                 " << H->C.count("mapq60") << " (" << H->C.perc("mapq60", "mappings") << "\% of mappings)" << endl;
		cerr << " | | FDR:                     " << H->C.perc("FP", "PP") << "%" << " = " << H->C.frac("FP", "reads") << " / " << H->C.frac("PP", "reads") << " per reads" << endl;

		//cerr << " | Average edit dist:     " << H->C.frac("total_edit_distance", "mappings") << endl;
		//print_time_stats();
        //printMemoryUsage();
	}

    void print_time_stats() {
        cerr << std::fixed << std::setprecision(1);
        cerr << " | Runtime:                "    << setw(5) << right << H->T.secs("mapping")       << " sec, " << 1.0 * H->C.count("reads") / H->T.secs("mapping")  << " reads/sec (" << setw(5) << right << H->T.range_ratio("query_mapping") << "x)" << endl; //setw(4) << right << H->C.count("reads") / H->T.secs("total") << " reads p/ sec)" << endl;
        cerr << " |  | load reads:             " << setw(5) << right << H->T.secs("query_reading")     << " (" << setw(4) << right << H->T.perc("query_reading", "mapping")       << "\%, " << setw(6) << right << H->T.range_ratio("query_reading") << "x)" << endl;
//        cerr << " |  | prepare:                "     << setw(5) << right << H->T.secs("prepare")           << " (" << setw(4) << right << H->T.perc("prepare", "mapping")             << "\%, " << setw(6) << right << H->T.range_ratio("prepare") << "x)" << endl;
        cerr << " |  | sketch reads:           " << setw(5) << right << H->T.secs("sketching")         << " (" << setw(4) << right << H->T.perc("sketching", "mapping")           << "\%, " << setw(6) << right << H->T.range_ratio("sketching") << "x)" << endl;
        cerr << " |  | seed:                   " << setw(5) << right << H->T.secs("seeding")           << " (" << setw(4) << right << H->T.perc("seeding", "mapping")             << "\%, " << setw(6) << right << H->T.range_ratio("seeding") << "x)" << endl;
        //cerr << " |  |  | collect seed info:       " << setw(5) << right << H->T.secs("collect_kmer_info") << " (" << setw(4) << right << H->T.perc("collect_kmer_info", "seeding")   << "\%, " << setw(6) << right << H->T.range_ratio("collect_kmer_info") << "x)" << endl;
        //cerr << " |  |  | sort seeds:              " << setw(5) << right << H->T.secs("sort_kmers")        << " (" << setw(4) << right << H->T.perc("sort_kmers", "seeding")          << "\%, " << setw(6) << right << H->T.range_ratio("sort_kmers") << "x)" << endl;
//        cerr << " |  |  | thin sketch:             " << setw(5) << right << H->T.secs("thin_sketch")       << " (" << setw(4) << right << H->T.perc("thin_sketch", "seeding")         << "\%, " << setw(6) << right << H->T.range_ratio("thin_sketch") << "x)" << endl;
//        cerr << " |  |  | unique seeds:            " << setw(5) << right << H->T.secs("unique_seeds")      << " (" << setw(4) << right << H->T.perc("unique_seeds", "seeding")        << "\%, " << setw(6) << right << H->T.range_ratio("unique_seeds") << "x)" << endl;
        cerr << " |  | match seeds:            " << setw(5) << right << H->T.secs("match_seeds")  << " (" << setw(4) << right << H->T.perc("match_seeds", "mapping")   << "\%, " << setw(6) << right << H->T.range_ratio("match_seeds") << "x)" << endl;
        cerr << " |  | match rest:             " << setw(5) << right << H->T.secs("match_rest")   << " (" << setw(4) << right << H->T.perc("match_rest", "mapping")     << "\%, " << setw(6) << right << H->T.range_ratio("match_rest") << "x)" << endl;
        cerr << " |  |  | seed heuristic:         " << setw(5) << right << H->T.secs("seed_heuristic")    << " (" << setw(4) << right << H->T.perc("seed_heuristic", "match_rest")     << "\%, " << setw(6) << right << H->T.range_ratio("seed_heuristic") << "x)" << endl;
        cerr << " |  |  | matches collect:        " << setw(5) << right << H->T.secs("match_collect")     << " (" << setw(4) << right << H->T.perc("match_collect", "match_rest")      << "\%, " << setw(6) << right << H->T.range_ratio("match_collect") << "x)" << endl;
        cerr << " |  |  | sweep:                  " << setw(5) << right << H->T.secs("sweep")             << " (" << setw(4) << right << H->T.perc("sweep", "match_rest")              << "\%, " << setw(6) << right << H->T.range_ratio("sweep") << "x)" << endl;

//        cerr << " |  |  | get intervals:           " << setw(5) << right << H->T.secs("get_intervals")     << " (" << setw(4) << right << H->T.perc("get_intervals", "mapping")      << "\%, " << setw(5) << right << H->T.range_ratio("get_intervals") << "x)" << endl;
//        cerr << " |  |  | filter promising buckets:" << setw(5) << right << H->T.secs("filter_promising_buckets")    << " (" << setw(4) << right << H->T.perc("filter_promising_buckets", "mapping")     << "\%, " << setw(5) << right << H->T.range_ratio("filter_promising_buckets") << "x)" << endl;
//        cerr << " |  | edit distance:          "     << setw(5) << right << H->T.secs("edit_distance")     << " (" << setw(4) << right << H->T.perc("edit_distance", "mapping")       << "\%, " << setw(6) << right << H->T.range_ratio("edit_distance") << "x)" << endl;
        cerr << " |  | output:                 " << setw(5) << right << H->T.secs("output")            << " (" << setw(4) << right << H->T.perc("output", "mapping")              << "\%, " << setw(6) << right << H->T.range_ratio("output") << "x)" << endl;
//        cerr << " |  | post proc:              "     << setw(5) << right << H->T.secs("postproc")          << " (" << setw(4) << right << H->T.perc("postproc", "mapping")            << "\%, " << setw(6) << right << H->T.range_ratio("postproc") << "x)" << endl;
    }
};

TEST_CASE("testing the seeding") {
	Handler* H;
	SketchIndex* tidx;
	SHMapper* mapper;

	params_t params;
	params.k = 25;
	params.hFrac = 0.05;
	params.theta = 0.7;

	H = new Handler(params);
	tidx = new SketchIndex(H);
	mapper = new SHMapper(*tidx, H);

	sketch_t p;
	p.emplace_back(10, 111111, false);
	p.emplace_back(20, 222222, false);
	p.emplace_back(30, 111111, false);
	p.emplace_back(40, 444444, false);
	p.emplace_back(50, 555555, false);

	int nonzero = 0;
	SHMapper::Seeds kmers = mapper->select_kmers(p, nonzero);

	CHECK(kmers.size() == 4);
	return;
	
	CHECK(kmers[0].kmer.h == 1);
	CHECK(kmers[0].occs_in_p == 2);
	
	CHECK(kmers[1].kmer.h == 2);
	CHECK(kmers[1].occs_in_p == 1);
	
	CHECK(kmers[2].kmer.h == 3);
	CHECK(kmers[2].occs_in_p == 2);

	for (size_t i = 1; i < kmers.size(); i++)
		CHECK(kmers[i-1].hits_in_T == kmers[i].hits_in_T);
}

}  // namespace sweepmap
