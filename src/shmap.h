#pragma once

#include <csignal>
#include <optional>

#include "analyse_simulated.h"
#include "buckets.h"
#include "handler.h"
#include "index.h"
#include "io.h"
#include "mapper.h"
//#include "refine.h"
#include "sketch.h"
#include "types.h"
#include "utils.h"

namespace sweepmap {
	
template <bool no_bucket_pruning, bool one_sweep, bool only_sh, bool abs_pos>
class SHMapper : public Mapper {
	using BucketsType = Buckets<abs_pos>;

	const SketchIndex &tidx;
	Handler *H;
	//Matcher matcher;

	Counters C;   // local counters: cleared for each read
	Timers T;     // local timers: cleared for each read

private:
	static SHMapper* current_instance;

	static void handle_signal(int signal) {
		if (current_instance) {
			std::cout << "Signal received: " << signal << "\n";
			current_instance->print_stats();
			current_instance->print_time_stats();
		}
		std::exit(0);
	}

public:
	//using Buckets = std::unordered_map<bucket_t, qpos_t>;
	//using Buckets = gtl::flat_hash_map<bucket_t, qpos_t>;
	//using Buckets = ankerl::unordered_dense::map<bucket_t, qpos_t>;
	//using h2cnt = unordered_map<hash_t, qpos_t>;
	//using h2cnt = gtl::flat_hash_map<hash_t, qpos_t>;

	// Returns unique kmers with at least one match in T
	Seeds select_kmers(sketch_t& p) {
		T.start("group_kmers");
			sort(p.begin(), p.end(), [](const Kmer &a, const Kmer &b) { return a.h < b.h; });
		T.stop("group_kmers");

		T.start("collect_kmer_info");
			Seeds p_unique;
			qpos_t strike = 0;
			rpos_t nonzero = 0;
			for (size_t ppos = 0; ppos < p.size(); ++ppos) {
				++strike;
				if (ppos == p.size()-1 || p[ppos].h != p[ppos+1].h) {
					rpos_t hits_in_t = tidx.count(p[ppos].h);  // TODO: optimixze to 1 query; take the average strand among all hits
					if (hits_in_t > 0)
						nonzero += strike;
					p_unique.emplace_back(p[ppos], hits_in_t, strike, p_unique.size());
					strike = 0;
				}
			}
		T.stop("collect_kmer_info");

		T.start("sort_kmers");
			//pdqsort_branchless(p_unique.begin(), p_unique.end(), [](const Seed &a, const Seed &b) {
			sort(p_unique.begin(), p_unique.end(), [](const Seed &a, const Seed &b) {
				return a.hits_in_T < b.hits_in_T;
			});
		T.stop("sort_kmers");

		C.inc("kmers_notmatched", p.size() - nonzero);

		return p_unique;
	}

	void match_seeds(const Seeds &p_unique, BucketsType &B, qpos_t S) {
		rpos_t seed_matches(0), seeded_buckets(0);  // stats
		for (; B.i < (qpos_t)p_unique.size() && B.seeds < S; B.i++) {
			Seed seed = p_unique[B.i];
			if (seed.hits_in_T > 0) {
				seed_matches += seed.hits_in_T;

				C["max_seed_matches"] = max(C["max_seed_matches"], Counter(seed.hits_in_T));
				//tidx.get_matches(&seed_matches, seed);
				
				if (seed.hits_in_T == 1) {
					//T.start("match_seeds_single");
						const auto &hit = tidx.h2single.at(seed.kmer.h);
						BucketContent content(1, 0, hit.strand == seed.kmer.strand ? 1 : -1, hit.r, hit.r);
						B.add_to_pos(hit, content);
					//T.stop("match_seeds_single");
				//}
				//else if (seed.occs_in_p == 1) {
				//	T.start("match_seeds_single");
				//		for (const auto &hit: tidx.h2multi.at(seed.kmer.h)) {
				//			BucketContent content(1, 0, hit.strand == seed.kmer.strand ? 1 : -1, hit.r, hit.r);
				//			B.add_to_bucket(BucketLoc(hit.segm_id, hit.r/B.get_bucket_len(), &B), content);
				//		}
				//	T.stop("match_seeds_single");
				} else {
					BucketsType b2m(B.get_bucket_len());
					for (const auto &hit: tidx.h2multi.at(seed.kmer.h)) {
						BucketContent content(1, 0, hit.strand == seed.kmer.strand ? 1 : -1, hit.r, hit.r);
						b2m.add_to_pos(hit, content);
					}
					for (auto it = b2m.buckets.begin(); it != b2m.buckets.end(); ++it) {
						BucketContent content(min(it->second.matches, seed.occs_in_p), 0, it->second.codirection, it->second.r_min, it->second.r_max);
						B.add_to_bucket(it->first, content);
						//matches_in_B += it->second.matches;
					}
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
			B.seeds += seed.occs_in_p;
		}
		seeded_buckets = B.size();
		C.inc("seeded_buckets", seeded_buckets);
		C.inc("seed_matches", seed_matches);
		C.inc("total_matches", seed_matches);
	}

  public:
	SHMapper(const SketchIndex &tidx, Handler *H)
			: tidx(tidx), H(H) { //, matcher(tidx) {
        if (H->params.theta < 0.0 || H->params.theta > 1.0) {
            cerr << "tThres = " << H->params.theta << " outside of [0,1]." << endl;
            exit(1);
        }
        current_instance = this;  // Set current instance in constructor
    }
	virtual ~SHMapper() {
		if (current_instance == this) {
			current_instance = nullptr;
		}
	}

	double hseed(qpos_t p, qpos_t seeds, qpos_t matches) {
		assert(seeds >= matches);
		return 1.0 - double(seeds - matches) / p;
	}

	void matches_in_bucket(const BucketsType &B, const BucketLoc &b, BucketContent *bucket, const Seed &s) const {
		bucket->seeds += s.occs_in_p;
		if (s.hits_in_T == 0) {
			return;
		} else if (s.hits_in_T == 1) {
			const auto &hit = tidx.h2single.at(s.kmer.h);
			if (abs_pos ? B.begin(b) <= hit.r    && hit.r    < B.end(b)
			 		    : B.begin(b) <= hit.tpos && hit.tpos < B.end(b)) {
				bucket->matches += 1;
				bucket->codirection += hit.strand == s.kmer.strand ? 1 : -1;
				bucket->r_min = min(bucket->r_min, hit.r);
				bucket->r_max = max(bucket->r_max, hit.r);
			}
		} else if (s.hits_in_T > 1) {
			const auto &hits = tidx.h2multi.at(s.kmer.h);
			auto it = lower_bound(hits.begin(), hits.end(), B.begin(b), [&b](const Hit &hit, rpos_t pos) {
				if (hit.segm_id != b.segm_id)
					return hit.segm_id < b.segm_id;
				if constexpr (abs_pos) {
					return hit.r < pos;
				} else {
					return hit.tpos < pos;
				}
			});

			rpos_t matches = 0;
			for (; it != hits.end(); ++it) {
				if constexpr (abs_pos) {
					if (!(it->segm_id == b.segm_id && it->r < B.end(b))) break;					
				} else {
					if (!(it->segm_id == b.segm_id && it->tpos < B.end(b))) break;
				}
				matches += 1;
				bucket->codirection += it->strand == s.kmer.strand ? 1 : -1;
				bucket->r_min = min(bucket->r_min, it->r);
				bucket->r_max = max(bucket->r_max, it->r);
			}
			bucket->matches += min(matches, s.occs_in_p);
		}
	}

	bool seed_heuristic_pass(const BucketsType &B, const Seeds &p_unique, qpos_t m, const BucketLoc &b, BucketContent *bucket, double *sh, double thr) {
		//T.start("seed_heuristic");
		if constexpr (!no_bucket_pruning) {
			while (1) {
				*sh = hseed(m, bucket->seeds, bucket->matches);
				if (*sh < thr)
					return false;
				if (bucket->i >= (qpos_t)p_unique.size())
					break;
				matches_in_bucket(B, b, bucket, p_unique[bucket->i]);
				bucket->i++;
			}
		}
		//T.stop("seed_heuristic");
		return true;
	}
	
	Mapping bestFixedLength(const RefSegment &segm, const BucketsType &B, const BucketLoc &b, const h2seed_t &p_ht, h2cnt &diff_hist, const qpos_t P_sz, const qpos_t m) {
		qpos_t intersection = 0;
		qpos_t same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise
		auto t = segm.kmers;

		assert(B.begin(b) < t.size());
		auto l = B.begin(b);
		auto r = l;
		auto end = min(rpos_t(t.size()), B.end(b));
		assert(l < t.size());

		Mapping best;
		for(; l < end; ++l) {
			//for(;  r != t.begin() + end  && l-> hit.tpos <= l->hit.tpos + m; ++r) {
				//&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				//&& r->hit.tpos + H->params.k <= l->hit.tpos + P_sz
				//&& r->hit.r + H->params.k <= l->hit.r + P_sz
			for(;  r < end; ++r) {
				if constexpr (abs_pos) {
					if (!(t[r].r < t[l].r + m)) break;
				} else {
					if (!(r < l + m)) break;
				}

				if (auto it = p_ht.find(t[r].h); it != p_ht.end()) {
					same_strand_seeds += codirection(it->second.kmer, t[r]);
					assert(diff_hist.contains(t[r].h));
					if (--diff_hist[t[r].h] >= 0)
						++intersection;
				}
				assert (l <= r);
			}

			double C = 1.0 * intersection / m;
			assert(C >= -0.0 && C <= 1.0);
			if (l < r && C > best.score())
				best.update(0, P_sz-1, t[l].r, t[r-1].r, segm, intersection, C, same_strand_seeds, t[r-1].r-t[l].r+1);

			if (auto it = p_ht.find(t[l].h); it != p_ht.end()) {
				same_strand_seeds -= codirection(it->second.kmer, t[l]);
				assert(diff_hist.contains(t[l].h));
				if (++diff_hist[t[l].h] >= 1)
					--intersection;
			}

			assert(intersection >= 0);
		}
		assert(intersection == 0);
		
		return best;
	}
	
	std::optional<Mapping> match_rest(qpos_t P_sz, qpos_t m, const Seeds &p_unique, const BucketsType &B, vector<typename BucketsType::bucket_map_t::value_type> &sorted_buckets, h2cnt &diff_hist, const h2seed_t &p_ht,
			double thr, std::optional<Mapping> forbidden, const string &query_id, int *lost_on_pruning, double max_overlap) {
		int lost_on_seeding = (0);
		std::optional<Mapping> best;
		C.inc("lost_on_seeding", lost_on_seeding);

		for (auto &[b, content]: sorted_buckets) {
			double sh = 1.0;
			if (seed_heuristic_pass(B, p_unique, m, b, &content, &sh, thr)) {
				T.start("sweep");
				//cerr << forbidden_idx << ", passed SH:" << b_it->first << " " << b_it->second << endl;
				//if (matcher.do_overlap(query_id, b_it->first))
				//	*lost_on_pruning = 0;
				//T.start("match_collect");
				//	Matches M = matcher.collect_matches(b_it->first, p_ht);
				//  C.inc("total_matches", M.size());
				//T.stop("match_collect");
					//auto best_in_bucket = matcher.bestFixedLength(M, P_sz-H->params.k, lmax, m, Metric::CONTAINMENT_INDEX);
				
				//cerr << b_it->first << ", " << b_it->second << "| sh: " << sh << " best_score: " << best_score << endl;

				Mapping best_in_bucket;
				if constexpr (only_sh)
					best_in_bucket.update(0, P_sz-1, content.r_min, content.r_max, tidx.T[b.segm_id], content.matches, sh, content.codirection, content.r_max - content.r_min); // TODO: should it be prev(r) instead?
				else
					best_in_bucket = bestFixedLength(tidx.T[b.segm_id], B, b, p_ht, diff_hist, P_sz-H->params.k, m);
				best_in_bucket.set_bucket(b);
				best_in_bucket.set_sh(sh);

				if (best_in_bucket.score() > thr) {
					C.inc("final_buckets");
					if (!forbidden || Mapping::overlap(best_in_bucket, forbidden.value()) < max_overlap) {
						best = best_in_bucket;
						thr = best_in_bucket.score();
						// TODO (optimize): if (forbidden_idx != -1) break;
					}
				}
				T.stop("sweep");
			}
		}

		return best;
	}


	inline void map_read(const params_t &params, const FracMinHash &sketcher, ofstream &paulout, ofstream &unmapped_out, const string &query_id, const string &P) {
		C.clear();
		T.clear();

		C.init("seeds_limit_reached", "mapped_reads", "kmers", "kmers_notmatched", "seeds", "matches",
			   "seed_matches", "max_seed_matches", "matches_freq", "spurious_matches", "mappings", "J_best",
			   "sketched_kmers", "total_edit_distance", "intersection_diff", "mapq60", "matches_in_reported_mappings",
			   "lost_on_seeding", "lost_on_pruning", "final_buckets");

		T.init("seed_heuristic", "match_collect", "sweep");

		T.start("query_mapping");
			T.start("sketching");
				//const char *P = seq->seq.s;
				sketch_t p = sketcher.sketch(P);
				qpos_t m = p.size();
				T.stop("sketching");

			T.start("prepare");
				//char *query_id = seq->name.s;
				//qpos_t P_sz = (qpos_t)seq->seq.l;
				C.inc("read_len", P.size());

				T.start("seeding");
					Seeds p_unique = select_kmers(p);
				T.stop("seeding");

				h2seed_t p_ht;
				h2cnt diff_hist;
				rpos_t possible_matches(0);
				for (const auto kmer: p_unique) {
					p_ht.insert(make_pair(kmer.kmer.h, kmer));
					diff_hist[kmer.kmer.h] = kmer.occs_in_p;
					//if (params.verbose >= 2) {
					//	matcher.update(diff_hist);
					//}
					possible_matches += kmer.hits_in_T;
				}
				C.inc("possible_matches", possible_matches);

				//qpos_t lmax = qpos_t(m / params.theta);					// maximum length of a similar mapping
				//qpos_t lmin = qpos_t(ceil(p.size() * params.theta));					// maximum length of a similar mapping
				//qpos_t lmax = qpos_t(p.size() / params.theta);					// maximum length of a similar mapping
				//qpos_t lmax = p.size();					// maximum length of a similar mapping
				//qpos_t bucket_l = lmax;

				//double coef = 1.0;// * nonzero / p.size();
				//cerr << "coef: " << coef << endl;
				//double new_theta = params.theta;
				//cerr << "theta: " << params.theta << " -> " << new_theta << endl;

				//double theta2 = params.theta + params.min_diff;
				double theta2 = params.theta - params.min_diff;
				qpos_t S = qpos_t((1.0 - theta2) * m) + 1;			// any similar mapping includes at least 1 seed match

				BucketsType B;
				if constexpr (abs_pos)
					B.set_bucket_halflen(P.size());
				else
					B.set_bucket_halflen(m/params.theta);

				C.inc("kmers_sketched", m);
				C.inc("kmers", m);
				C.inc("kmers_unique", p_unique.size());
				C.inc("kmers_seeds", S);
			T.stop("prepare");

			T.start("match_seeds");
				match_seeds(p_unique, B, S);
			T.stop("match_seeds");

			B.propagate_seeds_to_buckets();

			//int lost_on_seeding = matcher.lost_correct_mapping(query_id, B);
			//C.inc("lost_on_seeding", lost_on_seeding);
			
			auto sorted_buckets = B.get_sorted_buckets();

			int lost_on_pruning = 1;
			std::optional<Mapping> best, best2;
			T.start("match_rest");
				T.start("match_rest_for_best");
					best = match_rest(P.size(), m, p_unique, B, sorted_buckets, diff_hist, p_ht, params.theta, std::nullopt, query_id, &lost_on_pruning, params.max_overlap);
				T.stop("match_rest_for_best");
				T.start("match_rest_for_best2");
					if (best) {
						double second_best_thr = best->score() * (1.0 - params.min_diff);
						best2 = match_rest(P.size(), m, p_unique, B, sorted_buckets, diff_hist, p_ht, second_best_thr, best, query_id, &lost_on_pruning, params.max_overlap);
					}
				T.stop("match_rest_for_best2");
			T.stop("match_rest");
			
			auto [FP, PP, FDR, FPTP] = tuple(-1, -1, -1, -1);
			//auto [FP, PP, FDR, FPTP] = calc_FDR(maps, params.theta, lmax, bucket_l, tidx, p_ht, P_sz, lmin, m, diff_hist);
			//C.inc("FP", FP);
			//C.inc("PP", PP);
			//C.inc("FDR", FDR);
			//C.inc("FPTP", FPTP);

			C.inc("lost_on_pruning", lost_on_pruning);
		T.stop("query_mapping");

		T.start("output");
			std::ostream* os;
			if (best) {
				if (best2)
					best->set_second_best(best2.value());

				if (best->mapq() == 60) C.inc("mapq60");
				C.inc("matches_in_reported_mappings", best->intersection());
				C.inc("J_best", rpos_t(10000.0*best->score()));
				C.inc("mappings");
				C.inc("mapped_reads");
				os = &cout;
			} else {
				os = &unmapped_out;
			}
			best->set_global_stats(params.theta, params.min_diff, m, query_id.c_str(), P.size(), params.k, S, C.count("total_matches"), C.count("max_seed_matches"), C.count("seed_matches"), C.count("seeded_buckets"), C.count("final_buckets"), FPTP, T.secs("query_mapping"));  // TODO: disable by flag
			best->print_paf(*os);

			if (params.verbose >= 2) {
				AnalyseSimulatedReads<abs_pos> gt(query_id, P, P.size(), diff_hist, m, p_ht, tidx, B, params.theta);
				if (best) {
					gt.print_paf(cout);
					auto gt_overlap = Mapping::overlap(gt.gt_mapping, best.value());
					*os << "\tgt_mapping_len:i:" << (gt.gt_mapping.paf.T_r - gt.gt_mapping.paf.T_l)
					   << "\treported_mapping_len:i:" << (best->paf.T_r - best->paf.T_l)
					   << "\tgt_overlap:f:" << std::fixed << std::setprecision(3) << gt_overlap;
				}
				if (!params.paramsFile.empty())
					gt.print_tsv(paulout);
			}

			*os << endl;
		T.stop("output");
	}

	void map_reads(const string &pFile) {
		//std::signal(SIGINT, handle_signal);

		cerr << "Mapping reads using SHmap..." << endl;
		ProgressBar progress_bar("Mapping");

		string pauls_fn = H->params.paramsFile + ".paul.tsv";
		string unmapped_fn = H->params.paramsFile + ".unmapped.paf";
		ofstream paulout;
		ofstream unmapped_out;
		if (!H->params.paramsFile.empty()) {
			cerr << "Unmapped reads to " << unmapped_fn << endl;
			unmapped_out = ofstream(unmapped_fn);
			if (H->params.verbose >= 2) {
				cerr << "Paul's experiment to " << pauls_fn << endl;
				paulout = ofstream(pauls_fn);
			}
		}

		H->T.start("mapping");
		H->T.start("query_reading");
		read_fasta_klib(pFile, [this, &paulout, &unmapped_out, &progress_bar](const string &query_id, const string &P, float mapping_progress) {
			H->C.inc("reads");
			H->T.stop("query_reading");
				map_read(H->params, H->sketcher, paulout, unmapped_out, query_id, P);
			H->T.start("query_reading");
			H->C += C;
			H->T += T;
			if (H->C.count("mapped_reads") % 100 == 0)
				progress_bar.update(mapping_progress);
		});
		H->T.stop("query_reading");
		H->T.stop("mapping");
		std::cerr << std::endl;  // for the progress bar

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
		cerr << " | Buckets:              " << endl;
		cerr << " | | Seeded buckets:          " << H->C.frac("seeded_buckets", "reads") << " /read" << endl;
		cerr << " | Mappings:              " << H->C.count("mappings") << " (" << H->C.perc("mappings", "reads") << "\% of reads)" << endl;
		cerr << " | | Final buckes:            " << H->C.frac("final_buckets", "mappings") << " /mapping" << endl;
//		cerr << " | | Final mappings:          " << H->C.frac("final_mappings", "mappings") << " /mapping" << endl;
		cerr << " | | Average best sim.:       " << std::fixed << std::setprecision(3) << H->C.frac("J_best", "mappings") / 10000.0 << endl;
		cerr << " | | mapq=60:                 " << H->C.count("mapq60") << " (" << H->C.perc("mapq60", "mappings") << "\% of mappings)" << endl;
		//cerr << " | | FDR:                     " << H->C.perc("FP", "PP") << "%" << " = " << H->C.frac("FP", "reads") << " / " << H->C.frac("PP", "reads") << " per reads" << endl;

		//cerr << " | Average edit dist:     " << H->C.frac("total_edit_distance", "mappings") << endl;
		//print_time_stats();
        //printMemoryUsage();
	}

    void print_time_stats() {
        cerr << std::fixed << std::setprecision(1);
        cerr << " | Runtime:                "    << setw(5) << right << H->T.secs("mapping")       << " sec, " << 1.0 * H->C.count("reads") / H->T.secs("mapping")  << " reads/sec (" << setw(5) << right << H->T.range_ratio("query_mapping") << "x)" << endl; //setw(4) << right << H->C.count("reads") / H->T.secs("total") << " reads p/ sec)" << endl;
        cerr << " |  | load reads:             " << setw(5) << right << H->T.secs("query_reading")     << " (" << setw(4) << right << H->T.perc("query_reading", "mapping")       << "\%, " << setw(6) << right << H->T.range_ratio("query_reading") << "x)" << endl;
        cerr << " |  | sketch reads:           " << setw(5) << right << H->T.secs("sketching")         << " (" << setw(4) << right << H->T.perc("sketching", "mapping")           << "\%, " << setw(6) << right << H->T.range_ratio("sketching") << "x)" << endl;
//        cerr << " |  | prepare:                " << setw(5) << right << H->T.secs("prepare")           << " (" << setw(4) << right << H->T.perc("prepare", "mapping")             << "\%, " << setw(6) << right << H->T.range_ratio("prepare") << "x)" << endl;
        cerr << " |  | seed:                   " << setw(5) << right << H->T.secs("seeding")           << " (" << setw(4) << right << H->T.perc("seeding", "mapping")             << "\%, " << setw(6) << right << H->T.range_ratio("seeding") << "x)" << endl;
		cerr << " |  |  | group_kmers:            " << setw(5) << right << H->T.secs("group_kmers")       << " (" << setw(4) << right << H->T.perc("group_kmers", "seeding")      << "\%, " << setw(6) << right << H->T.range_ratio("group_kmers") << "x)" << endl;
		cerr << " |  |  | collect_kmer_info:      " << setw(5) << right << H->T.secs("collect_kmer_info") << " (" << setw(4) << right << H->T.perc("collect_kmer_info", "seeding")   << "\%, " << setw(6) << right << H->T.range_ratio("collect_kmer_info") << "x)" << endl;
		cerr << " |  |  | sort_kmers:             " << setw(5) << right << H->T.secs("sort_kmers")        << " (" << setw(4) << right << H->T.perc("sort_kmers", "seeding")          << "\%, " << setw(6) << right << H->T.range_ratio("sort_kmers") << "x)" << endl;
        //cerr << " |  |  | collect seed info:       " << setw(5) << right << H->T.secs("collect_kmer_info") << " (" << setw(4) << right << H->T.perc("collect_kmer_info", "seeding")   << "\%, " << setw(6) << right << H->T.range_ratio("collect_kmer_info") << "x)" << endl;
        //cerr << " |  |  | sort seeds:              " << setw(5) << right << H->T.secs("sort_kmers")        << " (" << setw(4) << right << H->T.perc("sort_kmers", "seeding")          << "\%, " << setw(6) << right << H->T.range_ratio("sort_kmers") << "x)" << endl;
//        cerr << " |  |  | thin sketch:             " << setw(5) << right << H->T.secs("thin_sketch")       << " (" << setw(4) << right << H->T.perc("thin_sketch", "seeding")         << "\%, " << setw(6) << right << H->T.range_ratio("thin_sketch") << "x)" << endl;
//        cerr << " |  |  | unique seeds:            " << setw(5) << right << H->T.secs("unique_seeds")      << " (" << setw(4) << right << H->T.perc("unique_seeds", "seeding")        << "\%, " << setw(6) << right << H->T.range_ratio("unique_seeds") << "x)" << endl;
        cerr << " |  | match seeds:            " << setw(5) << right << H->T.secs("match_seeds")  << " (" << setw(4) << right << H->T.perc("match_seeds", "mapping")   << "\%, " << setw(6) << right << H->T.range_ratio("match_seeds") << "x)" << endl;
//		cerr << " |  |  | match seeds single:     " << setw(5) << right << H->T.secs("match_seeds_single")  << " (" << setw(4) << right << H->T.perc("match_seeds_single", "match_seeds")   << "\%, " << setw(6) << right << H->T.range_ratio("match_seeds_single") << "x)" << endl;
        cerr << " |  | match rest:             " << setw(5) << right << H->T.secs("match_rest")   << " (" << setw(4) << right << H->T.perc("match_rest", "mapping")     << "\%, " << setw(6) << right << H->T.range_ratio("match_rest") << "x): " << H->T.perc("match_rest_for_best2", "match_rest") << "% for second best" << endl;
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

// Initialize static member
template <bool no_bucket_pruning, bool one_sweep, bool only_sh, bool abs_pos>
SHMapper<no_bucket_pruning, one_sweep, only_sh, abs_pos>* 
	SHMapper<no_bucket_pruning, one_sweep, only_sh, abs_pos>::current_instance = nullptr;

}  // namespace sweepmap

//std::tuple<int, int, double, double> calc_FDR(vector<Mapping> &maps, double theta, qpos_t lmax, const SketchIndex &tidx, const h2seed_t &p_ht, qpos_t P_sz, qpos_t lmin, qpos_t m, const h2cnt &diff_hist) {
//    int FP = 0;
//    for (auto &mapping: maps) {
//        auto M = matcher.collect_matches(mapping.bucket(), p_ht);
//        auto mapping_best_J = matcher.bestIncludedJaccard(M, P_sz, lmin, lmax, m);
//        mapping_best_J.set_bucket(mapping.bucket());
//        if (mapping_best_J.score() < H->params.theta)
//            ++FP;
//    }
//
//    int PP = maps.size();
//    double FDR = 1.0 * FP / PP;
//    double FPTP = (PP - FP > 0) ?  1.0 * FP / (PP - FP) : 0.0;
//    return {PP, FP, FDR, FPTP};
//}
//