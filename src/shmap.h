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
//#include "phmap.hpp"
#include <ankerl/unordered_dense.h>

namespace sweepmap {
	
class SHMapper : public Mapper {
	const SketchIndex &tidx;
	Handler *H;

	using Seeds = vector<Seed>;
	using Matches = vector<Match>;
	using Bucket = pair<bucket_t, qpos_t>;
	//using Buckets = std::unordered_map<bucket_t, qpos_t>;
	//using Buckets = gtl::flat_hash_map<bucket_t, qpos_t>;
	using Buckets = ankerl::unordered_dense::map<bucket_t, qpos_t>;
	//using Hist = unordered_map<hash_t, qpos_t>;
	//using Hist = gtl::flat_hash_map<hash_t, qpos_t>;
	using Hist = ankerl::unordered_dense::map<hash_t, qpos_t>;

	Seeds select_kmers(sketch_t& p) {
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
						if (H->params.max_matches == -1 || hits_in_t <= H->params.max_matches) {
							Seed el(p[ppos], hits_in_t, kmers.size());
							//strike = 1; // comment out for Weighted metric
							el.occs_in_p = strike;
							kmers.push_back(el);
							if (H->params.max_seeds != -1 && (rpos_t)kmers.size() >= H->params.max_seeds)  // TODO maybe account for occs_in_p
								break;
						}
					}
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

	Mapping sweep(const vector<Match> &M, const qpos_t P_sz, qpos_t lmax, const qpos_t m, const Seeds &kmers, const bucket_t &bucket, Hist &diff_hist) {
		qpos_t intersection = 0;
		qpos_t same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise
		Mapping best;
		
		for (int i=1; i<(int)M.size(); i++)
			assert(M[i-1].hit.tpos < M[i].hit.tpos);
		assert(M.size() == 0 || (int)M.size() <= M.back().hit.tpos - M.front().hit.tpos + 1);
		
		// Increase the left point end of the window [l,r) one by one. O(matches)
		for(auto l = M.begin(), r = M.begin(); l != M.end(); ++l) {
			// Increase the right end of the window [l,r) until it gets out.
			for(;  r != M.end()
				//&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				//&& r->hit.tpos + H->params.k <= l->hit.tpos + P_sz
				//&& r->hit.r + H->params.k <= l->hit.r + P_sz
				&& r->hit.tpos <= l->hit.tpos + lmax 
				; ++r) {
				assert(diff_hist.contains(r->seed.kmer.h));
				if (--diff_hist[r->seed.kmer.h] >= 0) {
					++intersection;
					same_strand_seeds += r->is_same_strand() ? +1 : -1;
				}
				
				assert (l->hit.r <= r->hit.r);
			}

			//double J = 1.0*intersection / kmers.size();
			double J = 1.0*intersection / m;
			//int s_sz = prev(r)->hit.tpos - l->hit.tpos + 1;
			//assert(s_sz >= 1);
			//assert(intersection <= s_sz);
			//assert(intersection <= m);
			//double J = 1.0*intersection / (m + s_sz - intersection);

			//if (fabs(J-J_) > 0.1)
			//	cerr << J << " " << J_ << endl;
			//cerr << std::fixed << std::setprecision(3) << J << " " << 1.0*intersection / m << endl;
			//double J = 1.0*intersection / m; //kmers.size();
			//cerr << "l=" << l->hit.r << ", r=" << prev(r)->hit.r << ", intersection=" << intersection << ", P_sz=" << P_sz << ", J= " << J << endl;
			//assert(m > 0);
			//double J = 1.0*intersection / m;
			//J = 1.0*intersection / p_sz;
			//J = 1.0*intersection / std::min(p_sz, s_sz);

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
		assert(intersection == 0);
		assert(same_strand_seeds == 0);
		
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

	double hseed(qpos_t p, qpos_t seeds, qpos_t matches) {
		return 1.0 - double(seeds - matches) / p;
	}

	bool seed_heuristic_pass(const vector<Mapping> &maps, const Seeds &kmers, qpos_t lmax, qpos_t m, bucket_t b, rpos_t max_matches, qpos_t i, qpos_t seeds,
			int best_idx, int best2_idx, int bests_idx[4], double *lowest_sh, double thr_init) {
		if (H->params.no_bucket_pruning)
			return true;

		*lowest_sh = 1.0; // should only get lower

		double thr2 = thr_init; //H->params.theta;  // safe threshold for the second best
		for (int i=3; i>=1; i--)
			if (bests_idx[i] != -1) {
				thr2 = maps[bests_idx[i]].J;
				break;
			}
//		double thr2 = best_idx == -1 ? H->params.theta : maps[best_idx].J;

		if (hseed(m, seeds, max_matches) < thr2)
			return false;

		H->T.start("seed_heuristic");
		bool ret = true;
		for (; i < (qpos_t)kmers.size(); i++) {
			seeds += kmers[i].occs_in_p;
			if (tidx.matches_in_bucket(kmers[i], b, lmax))
				max_matches += kmers[i].occs_in_p;
			//double thr2 = best2_idx == -1 ? H->params.theta : maps[best2_idx].J; 
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
		return std::abs(X - Y) / std::sqrt(X + Y);
	}

	int mapq(Mapping m) {
		// minimap2: mapQ = 40 (1-f2/f1) min(1, m/10) log f1, where m is #anchors on primary chain
		assert(m.J >= 0.0);
		if (m.J2 < H->params.theta)
			m.J2 = H->params.theta;
		if (m.J - m.J2 > 0.015) {
			if (abs(m.same_strand_seeds) < 200)
				return 5;
			return 60;
		} else {
			return 0;
		}
	}

	void match_rest(qpos_t seed_matches, qpos_t P_sz, qpos_t lmax, qpos_t m, const Seeds &kmers, const vector<Bucket> &B_vec, Hist &diff_hist, int seeds, int first_kmer_after_seeds, const unordered_map<hash_t, Seed> &p_ht,
			vector<Mapping> &maps, int &total_matches, int &best_idx, int &best2_idx, int &final_buckets, double thr_init, int forbidden_idx) {
		double best_J = -1.0, best_J2 = -1.0;
		total_matches = seed_matches;

		int lost_on_seeding = (0);
		H->C.inc("lost_on_seeding", lost_on_seeding);

		int bests_idx[4] = {-1, -1, -1, -1};  // best_idx[0] -- best bucket, best_idx[1] -- best bucket size, the rest are the 3 next best (at least one of them will not be adjacent to best_idx[0])
		int maps_idx = 0;
		final_buckets = 0;
		H->T.start("seed_heuristic"); H->T.stop("seed_heuristic");  // init
		H->T.start("match_collect"); H->T.stop("match_collect");
		H->T.start("sweep"); H->T.stop("sweep");
		for (auto &[b, seed_matches]: B_vec) {
			if (forbidden_idx != -1 && b.segm_id == maps[forbidden_idx].bucket.segm_id && abs(b.b - maps[forbidden_idx].bucket.b) <= 1)
				continue;
			double lowest_sh;
			if (seed_heuristic_pass(maps, kmers, lmax, m, b, seed_matches, first_kmer_after_seeds, seeds, best_idx, best2_idx, bests_idx, &lowest_sh, thr_init)) {
				H->T.start("match_collect");
					Matches M;
					for (rpos_t i = b.b*lmax; i < std::min((b.b+2)*lmax, (rpos_t)tidx.T[b.segm_id].kmers.size()); i++) {
						assert(b.segm_id >= 0 && b.segm_id < (rpos_t)tidx.T.size());
						const auto &kmer = tidx.T[b.segm_id].kmers[i];
						const auto seed_it = p_ht.find(kmer.h);
						if (seed_it != p_ht.end()) {
							M.push_back(Match(seed_it->second, Hit(kmer, i, b.segm_id)));
							assert(M[ M.size()-1 ].hit.tpos == M[ M.size()-1 ].hit.tpos);
						}
					}
				H->T.stop("match_collect");
				total_matches += M.size();
				++final_buckets;

				H->T.start("sweep");
					auto bucket_best = sweep(M, P_sz, lmax, m, kmers, b, diff_hist);
					assert(bucket_best.J <= lowest_sh + 1e-7);
					if (bucket_best.J >= H->params.theta)
						maps.push_back(bucket_best);
				H->T.stop("sweep");
				
				for (; maps_idx < static_cast<int>(maps.size()); ++maps_idx) {
					for (int i = 0; i < 4; ++i) {
						if (bests_idx[i] == -1 || maps[bests_idx[i]].J < maps[maps_idx].J) {
							for (int j=3; j>=i+1; --j)
								bests_idx[j] = bests_idx[j-1];
							bests_idx[i] = maps_idx;

							best_idx = bests_idx[0];
							int j;
							for (j=1; j<4; ++j) {
								assert(bests_idx[j] == -1 || maps[ bests_idx[j-1] ].J >= maps[ bests_idx[j] ].J);
								if (bests_idx[j] == -1 || overlap(maps[bests_idx[j]], maps[best_idx]) < 0.5) {
									best2_idx = bests_idx[j];
									break;
								}
							}
							assert(j < 4);
							break;
						}
					}
				}
				for (int i=0; i<4; ++i) {
					assert(bests_idx[i] == -1 || maps[bests_idx[i]].J >= 0);
					assert(i==0 || bests_idx[i] == -1 || maps[bests_idx[i-1]].J >= maps[bests_idx[i]].J);
				}
				assert(best_idx == -1 || maps[best_idx].J >= bucket_best.J);

				best_J = max(best_J, bucket_best.J);
				if (best2_idx != -1) best_J2 = max(best_J2, maps[best2_idx].J);
				assert(best_idx == -1 || fabs(maps[best_idx].J - best_J) < 1e-7);
			}
		}
		assert(best_idx != -1 || maps.empty());
		
		int lost_on_pruning = (best_idx == -1);
		H->C.inc("final_buckets", final_buckets);
		H->C.inc("final_mappings", maps.size());
		H->C.inc("lost_on_pruning", lost_on_pruning);
		H->C.inc("total_matches", total_matches);
	}

	void map(const string &pFile) {
		cerr << "Mapping reads using SHmap..." << endl;

		H->C.inc("kmers", 0);
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

					Seeds kmers = select_kmers(p);
					qpos_t m = 0;
					for (const auto kmer: kmers)
						m += kmer.occs_in_p;

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

					qpos_t lmax = m/0.8;
					//qpos_t lmax = qpos_t(m / H->params.theta);					// maximum length of a similar mapping
					qpos_t S = qpos_t((1.0 - H->params.theta) * m) + 1;			// any similar mapping includes at least 1 seed match
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
								rpos_t b(m.hit.tpos / lmax);
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
							//	bucket_t curr_b(seed_matches[i].hit.segm_id, seed_matches[i].hit.tpos / lmax);
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
					//seed_matches *= 2;  // since we match each seed to 2 buckets
					seeded_buckets = B.size();
					H->C.inc("seeded_buckets", seeded_buckets);
					H->C.inc("seed_matches", seed_matches);
				H->T.stop("match_seeds");

				//cerr << "seed matches=" << seed_matches << ", matches_in_B=" << matches_in_B << ", max_seed_matches=" << max_seed_matches << ", seeds=" << seeds << ", S=" << S << ", i=" << i << endl;

				H->T.start("match_rest");
					vector<Mapping> maps;
					vector<Bucket> B_vec(B.begin(), B.end());
					sort(B_vec.begin(), B_vec.end(), [](const Bucket &a, const Bucket &b) { return a.second > b.second; });  // TODO: sort intervals by decreasing number of matches
					int total_matches, best_idx=-1, best2_idx=-1, final_buckets;
					match_rest(seed_matches, P_sz, lmax, m, kmers, B_vec, diff_hist, seeds, i, p_ht, maps, total_matches, best_idx, best2_idx, final_buckets, H->params.theta, -1);
				H->T.stop("match_rest");
			H->T.stop("query_mapping");

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

						m.mapq = mapq(m);
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

}  // namespace sweepmap
