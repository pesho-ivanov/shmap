#pragma once

#include <set>
#include <map>

#include "analyse_simulated.h"
#include "index.h"
#include "io.h"
#include "sketch.h"
#include "utils.h"
#include "handler.h"
#include "mapper.h"
#include "buckets.h"
#include "refine.h"

//#include "indicators.hpp"
//using namespace indicators;

namespace sweepmap {

class SHMapper : public Mapper {
	const SketchIndex &tidx;
	Handler *H;
	Matcher matcher;

public:
	//typedef pair<bucket_t, qpos_t> Bucket;
	//using Buckets = std::unordered_map<bucket_t, qpos_t>;
	//using Buckets = gtl::flat_hash_map<bucket_t, qpos_t>;
	//using Buckets = ankerl::unordered_dense::map<bucket_t, qpos_t>;
	//using Hist = unordered_map<hash_t, qpos_t>;
	//using Hist = gtl::flat_hash_map<hash_t, qpos_t>;

	// Returns unique kmers with at least one match in T
	Seeds select_kmers(sketch_t& p_all, int &nonzero) {
		H->T.start("seeding");
		H->T.start("collect_kmer_info");
			Seeds kmers;
			sort(p_all.begin(), p_all.end(), [](const Kmer &a, const Kmer &b) { return a.h < b.h; });
			qpos_t strike = 0;
			for (size_t ppos = 0; ppos < p_all.size(); ++ppos) {
				++strike;
				if (ppos == p_all.size()-1 || p_all[ppos].h != p_all[ppos+1].h) {
					rpos_t hits_in_t = tidx.count(p_all[ppos].h);
					if (hits_in_t > 0) {
						nonzero += strike;
					}
						if (H->params.max_matches == -1 || hits_in_t <= H->params.max_matches) {
							Seed el(p_all[ppos], hits_in_t, kmers.size());
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

  public:
	SHMapper(const SketchIndex &tidx, Handler *H)
			: tidx(tidx), H(H), matcher(tidx) {
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

	bool seed_heuristic_pass(const Seeds &kmers, qpos_t m, Buckets::Bucket b, rpos_t seed_matches, qpos_t i, qpos_t seeds, int best_idx, double *lowest_sh, double thr) {
		if (H->params.no_bucket_pruning)
			return true;

		*lowest_sh = 1.0; // should only get lower

		if (hseed(m, seeds, seed_matches) < thr)
			return false;

		H->T.start("seed_heuristic");
		bool ret = true;
		for (; i < (qpos_t)kmers.size(); i++) {
			seeds += kmers[i].occs_in_p;
			if (tidx.matches_in_bucket(kmers[i], b))
				seed_matches += kmers[i].occs_in_p;
			double sh = hseed(m, seeds, seed_matches);
			assert(sh <= *lowest_sh);
			*lowest_sh = min(*lowest_sh, sh);	
			if (sh < thr) {
				ret = false;
				break;
			}
		}
		H->T.stop("seed_heuristic");
		return ret;
	}
	
	void match_rest(qpos_t P_sz, qpos_t lmax, qpos_t m, const Seeds &kmers, Buckets &B, Hist &diff_hist, int seeds, int first_kmer_after_seeds, const unordered_map<hash_t, Seed> &p_ht,
			vector<Mapping> &maps, int &total_matches, int &best_idx, int &final_buckets, double thr, double min_diff, int forbidden_idx, const string &query_id, int *lost_on_pruning, double max_overlap) {
		double best_J = -1.0;

		int lost_on_seeding = (0);
		H->C.inc("lost_on_seeding", lost_on_seeding);

		final_buckets = 0;
		H->T.start("seed_heuristic"); H->T.stop("seed_heuristic");  // init
		H->T.start("match_collect"); H->T.stop("match_collect");
		H->T.start("sweep"); H->T.stop("sweep");
		for (auto b_it = B.ordered_begin(); b_it != B.ordered_end(); ++b_it) {
			double lowest_sh;
			if (seed_heuristic_pass(kmers, m, b_it->first, b_it->second, first_kmer_after_seeds, seeds, best_idx, &lowest_sh, thr)) {
				//cerr << forbidden_idx << ", passed SH:" << b_it->first << " " << b_it->second << endl;
				//if (matcher.do_overlap(query_id, b_it->first))
				//	*lost_on_pruning = 0;
				H->T.start("match_collect");
					Matches M = matcher.collect_matches(b_it->first, p_ht);
				H->T.stop("match_collect");
				total_matches += M.size();
				++final_buckets;

				H->T.start("sweep");
					auto best_in_bucket = matcher.bestFixedLength(M, P_sz, lmax, m, Metric::CONTAINMENT_INDEX);
					best_in_bucket.set_bucket(b_it->first);
					assert(best_in_bucket.score() <= lowest_sh + 1e-7);
					if (best_in_bucket.score() < thr - min_diff)
						B.delete_bucket(b_it->first);
					if (best_in_bucket.score() >= thr) {
						maps.emplace_back(best_in_bucket);
						if (best_in_bucket.score() > best_J) {
							if (forbidden_idx == -1 || Mapping::overlap(maps.back(), maps[forbidden_idx]) < max_overlap) {
								best_J = max(best_J, best_in_bucket.score());
								best_idx = maps.size()-1;
								thr = best_in_bucket.score();
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

	std::tuple<int, int, double, double> calc_FDR(vector<Mapping> &maps, double theta, qpos_t lmax, const SketchIndex &tidx, const unordered_map<hash_t, Seed> &p_ht, qpos_t P_sz, qpos_t lmin, qpos_t m, const Hist &diff_hist) {
		int FP = 0;
		for (auto &mapping: maps) {
			auto M = matcher.collect_matches(mapping.bucket(), p_ht);
			auto mapping_best_J = matcher.bestIncludedJaccard(M, P_sz, lmin, lmax, m);
			mapping_best_J.set_bucket(mapping.bucket());
			if (mapping_best_J.score() < H->params.theta)
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

		read_fasta_klib(pFile, [this, &paulout, &unmapped_out](const string &query_id, const string &P, float progress) {
			H->C.inc("reads");
			H->T.stop("query_reading");
			H->T.start("query_mapping");
				H->T.start("sketching");
					//const char *P = seq->seq.s;
					sketch_t p_all = H->sketcher.sketch(P);
					H->T.stop("sketching");

				H->T.start("prepare");
					//char *query_id = seq->name.s;
					//qpos_t P_sz = (qpos_t)seq->seq.l;
					qpos_t P_sz = P.size();

					H->C.inc("read_len", P_sz);

					int nonzero = 0;
					Seeds kmers = select_kmers(p_all, nonzero);
					qpos_t m = 0;
					for (const auto kmer: kmers)
						m += kmer.occs_in_p;
					assert(m <= (int)p_all.size());
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
					matcher.update(diff_hist);

					//qpos_t lmax = qpos_t(m / H->params.theta);					// maximum length of a similar mapping
					//qpos_t lmin = qpos_t(ceil(p.size() * H->params.theta));					// maximum length of a similar mapping
					qpos_t lmax = qpos_t(p_all.size() / H->params.theta);					// maximum length of a similar mapping
					//qpos_t lmax = p_all.size();					// maximum length of a similar mapping
					//qpos_t bucket_l = lmax;

					//double coef = 1.0;// * nonzero / p.size();
					//cerr << "coef: " << coef << endl;
					double new_theta = H->params.theta;
					//cerr << "theta: " << H->params.theta << " -> " << new_theta << endl;

					qpos_t S = qpos_t((1.0 - new_theta) * m) + 1;			// any similar mapping includes at least 1 seed match
					Buckets B(lmax);  			// B[segment][b] -- #matched kmers[0...i] in [bl, (b+2)l)
					qpos_t seeds = 0;

					H->C.inc("kmers_sketched", p_all.size());
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
							Buckets b2m(B.get_bucket_len());
							for (auto &m: seed_matches) {
								//rpos_t b(m.hit.tpos / B.get_bucket_len());
								b2m.assign_to_pos(Buckets::Pos(m.hit.segm_id, m.hit.tpos), seed.occs_in_p);
								//if (b > 0) b2m.add_to_bucket(bucket_t(m.hit.segm_id, b-1), seed.occs_in_p);
							}
							for (auto it = b2m.unordered_begin(); it != b2m.unordered_end(); ++it) {
								//B[b] += min(seed.occs_in_p, matches);
								B.add_to_bucket(it->first, it->second);
								matches_in_B += it->second;
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

				//int lost_on_seeding = matcher.lost_correct_mapping(query_id, B);
				//H->C.inc("lost_on_seeding", lost_on_seeding);

				//	cerr << "seeded buckets:" << b_it->first << " " << b_it->second << endl;
				//}

				int lost_on_pruning = 1;
				double min_diff = H->params.min_diff;
				H->T.start("match_rest");
				H->T.start("match_rest_for_best");
					vector<Mapping> maps;
					//vector<Bucket> B_vec(B.begin(), B.end());
					//sort(B_vec.begin(), B_vec.end(), [](const Bucket &a, const Bucket &b) { return a.second > b.second; });  // TODO: sort intervals by decreasing number of matches
					int total_matches=seed_matches, best_idx=-1, best2_idx=-1, final_buckets;
					double thr_best = H->params.theta;
					match_rest(P_sz, p_all.size(), m, kmers, B, diff_hist, seeds, i, p_ht, maps, total_matches, best_idx, final_buckets, thr_best, min_diff, -1, query_id, &lost_on_pruning, H->params.max_overlap);
				H->T.stop("match_rest_for_best");

					auto [FP, PP, FDR, FPTP] = tuple(-1, -1, -1, -1);
					//auto [FP, PP, FDR, FPTP] = calc_FDR(maps, H->params.theta, lmax, bucket_l, tidx, p_ht, P_sz, lmin, m, diff_hist);
					H->C.inc("FP", FP);
					H->C.inc("PP", PP);
					H->C.inc("FDR", FDR);
					H->C.inc("FPTP", FPTP);

					//cerr << "best_idx: " << best_idx << " " << maps[best_idx] << endl;

				H->T.start("match_rest_for_best2");
					if (best_idx != -1) {
						// find second best mapping for mapq computation
						double second_best_thr = maps[best_idx].score() - H->params.min_diff;
						match_rest(P_sz, p_all.size(), m, kmers, B, diff_hist, seeds, i, p_ht, maps, total_matches, best2_idx, final_buckets, second_best_thr, min_diff, best_idx, query_id, &lost_on_pruning, H->params.max_overlap);
					}
				H->T.stop("match_rest_for_best2");
				H->T.stop("match_rest");
				H->C.inc("lost_on_pruning", lost_on_pruning);
			H->T.stop("query_mapping");

			H->T.start("output");
				std::ostream* os;
				Mapping best;
				if (best_idx != -1) {
					best = maps[best_idx];
					//const auto &segm = tidx.T[m.segm_id()];
					//m.set_local_stats(segm);
					if (best2_idx != -1)
						best.set_second_best(maps[best2_idx]);

					if (best.mapq() == 60) H->C.inc("mapq60");
					H->C.inc("matches_in_reported_mappings", best.intersection());
					H->C.inc("J_best", rpos_t(10000.0*best.score()));
					H->C.inc("mappings");
					H->C.inc("mapped_reads");
					os = &cout;
				} else {
					os = &unmapped_out;
				}
				best.set_global_stats(H->params.theta, H->params.min_diff, p_all.size(), query_id.c_str(), P_sz, H->params.k, S, total_matches, max_seed_matches, seed_matches, seeded_buckets, final_buckets, FPTP, H->T.secs("query_mapping"));  // TODO: disable by flag
				best.print_paf(*os);

				if (H->params.verbose >= 2) {
					AnalyseSimulatedReads gt(query_id, P, P_sz, diff_hist, m, p_ht, tidx, B, H->params.theta);
					if (best_idx != -1) {
						gt.print_paf(cout);
						auto gt_overlap = Mapping::overlap(gt.gt_mapping, maps[best_idx]);
						*os << "\tgt_mapping_len:" << (gt.gt_mapping.paf.T_r-gt.gt_mapping.paf.T_l)
						   << "\treported_mapping_len:" << (maps[best_idx].paf.T_r-maps[best_idx].paf.T_l)
						   << "\tgt_overlap:" << std::fixed << std::setprecision(3) << gt_overlap;
					}
					if (!H->params.paramsFile.empty())
						gt.print_tsv(paulout);
				}

				*os << endl;
			H->T.stop("output");

			if (H->C.count("mapped_reads") % 100 == 0) {
//				std::ostringstream msg_stream;
//				msg_stream
//					<< std::fixed << std::setprecision(1)
//					<< H->C.count("mapped_reads") << " reads: "
//					<< H->C.perc("mapped_reads", "reads") << "% mapped, ";
//					//<< H->C.perc("mapq60", "mapped_reads") << "% with Q60";
//				string msg = msg_stream.str();
				printProgress(std::cerr, progress, "Mapping");

//				cerr << endl;
//				print_stats();
//				print_time_stats();
//				cerr << "\033[32A\033[G";
//				cerr.flush();
			}

			H->T.start("query_reading");
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
		cerr << " | | Seeded buckets:          " << H->C.frac("seeded_buckets", "reads") << " /reads" << endl;
		cerr << " | Mappings:              " << H->C.count("mappings") << " (" << H->C.perc("mappings", "reads") << "\% of reads)" << endl;
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

}  // namespace sweepmap
