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

	bool overlap(const Mapping &a, const Mapping &b) {
		if (a.segm_id != b.segm_id)
			return false;
		return a.T_r >= b.T_l && a.T_l <= b.T_r;
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
				same_strand_seeds += r->is_same_strand() ? +1 : -1;
				assert(diff_hist.contains(r->seed.kmer.h));
				if (--diff_hist[r->seed.kmer.h] >= 0)
					++intersection;
				
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

			same_strand_seeds -= l->is_same_strand() ? +1 : -1;
			assert(diff_hist.contains(l->seed.kmer.h));
			if (++diff_hist[l->seed.kmer.h] >= 1)
				--intersection;

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

	bool seed_heuristic_pass(const vector<Mapping> &maps, const Seeds &kmers, qpos_t lmax, qpos_t m, bucket_t b, rpos_t max_matches, qpos_t i, qpos_t seeds, int best_idx, int best2_idx, int bests_idx[4], double *lowest_sh) {
		if (H->params.no_bucket_pruning)
			return true;

		*lowest_sh = 1.0; // should only get lower

		double thr2 = H->params.theta;  // safe threshold for the second best
		for (int i=3; i>=1; i--)
			if (bests_idx[i] != -1) {
				thr2 = maps[bests_idx[i]].J;
				break;
			}

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
		//if (abs(m.intersection - m.intersection2) == 0)
		//	return 0;
		return m.J - m.J2 >= 0.015 ? 60 : 0;
		//return sigmas_diff(m.intersection, m.intersection2) >= 0.2 ? 60 : 0;

		//if (m.J < H->params.theta)
		//	return 5;
		//double bound = m.J * 0.9;
		//double r = max(m.J2 - bound, 0.0) / (m.J - bound);  // low is good
		//double J_fl = 60.0 * (1.0 - 1.0 * r);  // high is good
		//return int(J_fl/10.0) * 10;
		//return 60.0 * (1.0 - 1.0 * pow(J_second / J_best, 0.5) );
		//return  60.0 * (1.0 - 1.0 * pow(r, 2.0));
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

		H->T.start("mapping");
		H->T.start("query_reading");

		read_fasta_klib(pFile, [this](kseq_t *seq) {
			H->C.inc("reads");
			H->T.stop("query_reading");
			H->T.start("query_mapping");
				H->T.start("sketching");
					const char *P = seq->seq.s;
					sketch_t p = H->sketcher.sketch(P);
					H->T.stop("sketching");

				H->T.start("prepare");
					char *query_id = seq->name.s;
					qpos_t P_sz = (qpos_t)seq->seq.l;

					H->C.inc("read_len", P_sz);
					Timer read_mapping_time;  // TODO: change to H->T.start
					read_mapping_time.start();

					Seeds kmers = select_kmers(p);
					qpos_t m = 0;
					for (const auto kmer: kmers)
						m += kmer.occs_in_p;

					unordered_map<hash_t, Seed> p_ht;
					Hist diff_hist;
					rpos_t potential_matches(0);
					for (const auto kmer: kmers) {
						p_ht.insert(make_pair(kmer.kmer.h, kmer));
						//diff_hist[kmer.kmer.h] = 1;
						diff_hist[kmer.kmer.h] = kmer.occs_in_p;
						potential_matches += kmer.hits_in_T;
					}
					H->C.inc("potential_matches", potential_matches);

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
					seeded_buckets = B.size();
					H->C.inc("seeded_buckets", seeded_buckets);
					H->C.inc("seed_matches", seed_matches);
				H->T.stop("match_seeds");

				//cerr << "seed matches=" << seed_matches << ", matches_in_B=" << matches_in_B << ", max_seed_matches=" << max_seed_matches << ", seeds=" << seeds << ", S=" << S << ", i=" << i << endl;

				H->T.start("match_rest");
					double best_J = -1.0, best_J2 = -1.0;
					vector<Mapping> maps;
					qpos_t total_matches = seed_matches;
					vector<Bucket> B_vec(B.begin(), B.end());
					sort(B_vec.begin(), B_vec.end(), [](const Bucket &a, const Bucket &b) { return a.second > b.second; });  // TODO: sort intervals by decreasing number of matches

					int lost_on_seeding = (0);
					H->C.inc("lost_on_seeding", lost_on_seeding);

					int best_idx(-1), best2_idx(-1);
					int bests_idx[4] = {-1, -1, -1, -1};  // best_idx[0] -- best bucket, best_idx[1] -- best bucket size, the rest are the 3 next best (at least one of them will not be adjacent to best_idx[0])
					int maps_idx = 0;
					vector<Bucket> final_buckets;
					H->T.start("seed_heuristic"); H->T.stop("seed_heuristic");  // init
					H->T.start("match_collect"); H->T.stop("match_collect");
					H->T.start("sweep"); H->T.stop("sweep");
					for (auto &[b, seed_matches]: B_vec) {
						double lowest_sh;
						if (seed_heuristic_pass(maps, kmers, lmax, m, b, seed_matches, i, seeds, best_idx, best2_idx, bests_idx, &lowest_sh)) {
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
							final_buckets.push_back( Bucket(b, seed_matches) );

							//cerr << "seed matches=" << seed_matches << ", max_matches_in_bucket=" << max_matches_in_bucket << ", M.size()=" << M.size() << endl;
							//assert((rpos_t)M.size() >= seed_matches);

							H->T.start("sweep");
								auto bucket_best = sweep(M, P_sz, lmax, m, kmers, b, diff_hist);
								//if (best.J >= H->params.theta)
								//cerr << "B=" << B.size() << ", bucket= " << bucket << ", best=" << best;
								//if (best.J > lowest_sh + 1e-7) {
								//	cerr << std::fixed << std::setprecision(5);
								//	cerr << "lowest_sh=" << lowest_sh << ", best=" << best << ", b=" << b << ", seed_matches=" << seed_matches << ", seeds=" << seeds << ", i=" << i << ", best_idx=" << best_idx << ", best2_idx=" << best2_idx << endl;
								//}
								assert(bucket_best.J <= lowest_sh + 1e-7);
								if (bucket_best.J >= H->params.theta)
									maps.push_back(bucket_best);
							H->T.stop("sweep");
							
							//cerr << "Bucket " << b << " (" << seed_matches << " matches) has " << maps.size() << " good mappings" << endl;
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
											if (bests_idx[j] == -1
												|| !overlap(maps[bests_idx[j]], maps[best_idx])) {
//												|| (maps[bests_idx[j]].bucket.segm_id != maps[best_idx].bucket.segm_id)
//												|| abs(maps[bests_idx[j]].bucket.b - maps[best_idx].bucket.b) > 1) {
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

							//if (best2_idx != -1 && fabs(maps[best2_idx].J - best_J2) >= 1e-7) {
							//	cerr << std::fixed << std::setprecision(4);
							//	cerr << "best_idx: " << best_idx << ", best2_idx: " << best2_idx << ", best_J: " << best_J << ", best_J2: " << best_J2 << endl;
							//	for (int i=0; i<4; ++i) {
							//		cerr << "bests_idx[" << i << "] = " << bests_idx[i];
							//		if (bests_idx[i] != -1) cerr << ", J=" << maps[bests_idx[i]].J << ", mapping=" << maps[bests_idx[i]];
							//		cerr << endl;
							//	}
							//	cerr << "bucket_best: " << bucket_best << endl;
							//	cerr << "maps[best2_idx]: " << maps[best2_idx] << endl;
							//}
							//assert(best2_idx == -1 || fabs(maps[best2_idx].J - best_J2) < 1e-7);

							//assert(best2_idx == -1 || abs(maps[best_idx].bucket - maps[best2_idx].bucket) > 1);
						}
					}
					assert(best_idx != -1 || maps.empty());
					//cerr << "read " << query_id << ", seeded buckets: " << B.size() << ", final buckets: " << final_buckets.size() << ", maps.size(): " << maps.size() << ", best_idx: " << best_idx << ", best2_idx: " << best2_idx << endl;
					
					int lost_on_pruning = (best_idx == -1);
					H->C.inc("final_buckets", final_buckets.size());
					H->C.inc("final_mappings", maps.size());
					H->C.inc("lost_on_pruning", lost_on_pruning);
					H->C.inc("total_matches", total_matches);
				H->T.stop("match_rest");

				H->T.start("edit_distance");
				H->T.stop("edit_distance");
				
				H->T.start("output");
					vector<Mapping> final_mappings;

					if (H->params.onlybest && maps.size() >= 1) {
						H->C.inc("mapped_reads");
						if (best_idx != -1) 
						//if (maps[best_idx].J >= H->params.theta)
						{
							assert(0 <= best_idx && best_idx < (int)maps.size());
							Mapping &m = maps[best_idx];

							const auto &segm = tidx.T[m.segm_id];
							m.total_matches = total_matches;
							m.max_seed_matches = max_seed_matches;
							m.seed_matches = seed_matches;
							m.seeded_buckets = seeded_buckets;
							m.final_buckets = final_buckets.size();
							m.query_id = query_id;
							m.segm_name = segm.name.c_str();
							m.segm_sz = segm.sz;

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
							if (m.mapq == 60)
								H->C.inc("mapq60");

							//if (m.mapq > 0)
								final_mappings.push_back(m);
						}
					}

					read_mapping_time.stop();

					for (auto &m: final_mappings) {
						m.map_time = read_mapping_time.secs() / (double)final_mappings.size();
		//				if (H->params.sam) {
		//					auto ed = m.print_sam(query_id, segm, (int)matches.size(), seq->seq.s, seq->seq.l);
		//					H->C.inc("total_edit_distance", ed);
		//				}
		//				else
							m.print_paf();
						//  H->C.inc("spurious_matches", spurious_matches(m, matches));
						assert(m.J >= 0);
						H->C.inc("J_best", rpos_t(10000.0*m.J));
						H->C.inc("mappings");
					}
		//			H->C.inc("matches", matches.size());
					//H->T.stop("postproc");

				H->T.stop("output");
				H->T.stop("query_mapping");
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
		cerr << " |  | potential_matches:      " << H->C.frac("potential_matches", "reads") << " (" <<H->C.frac("potential_matches", "total_matches") << "x)" << endl;
		cerr << " |  | seed matches:           " << H->C.frac("seed_matches", "reads") << " (" << H->C.perc("seed_matches", "total_matches") << "%)" << endl;
//		cerr << " |  | frequent:               " << H->C.count("matches_freq") << " (" << H->C.perc("matches_freq", "total_matches") << "%)" << endl;
//		cerr << " |  | Seed h. reduction:      " << H->C.frac("potential_matches", "seed_matches") << "x" << endl;
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

//	double covered_frac(rpos_t bucket_a, rpos_t bucket_b, rpos_t gt_a, rpos_t gt_b) {
//		return 1.0*max(0, min(bucket_b, gt_b) - max(bucket_a, gt_a)) / (gt_b - gt_a);
//	}
	
	// TODO: works only in original coordinates
//	bool is_safe(const string &query_id, const vector<pair<int,int>> &potential_buckets, int lmax, int *gt_a, int *gt_b) {
//		std::vector<std::string> tokens;
//		std::stringstream ss(query_id);
//		std::string token;
//
//		while (std::getline(ss, token, '!'))
//			tokens.push_back(token);
//
//		if (tokens.size() < 4)
//			return true;
//
//		*gt_a = std::stoi(tokens[2]);
//		*gt_b = std::stoi(tokens[3]) + 1;
//
//		//cerr << "safety check: " << query_id << " " << *gt_a << " " << *gt_b << " among " << potential_buckets.size() << " buckets." << endl;
//		for (const auto &[b, matches]: potential_buckets) {
//			int bucket_a = b*lmax;
//			int bucket_b = (b+2)*lmax;
//			//if (bucket_a <= *gt_a && *gt_b <= bucket_b)
//			//cerr << "safety check: " << query_id << " " << bucket_a << " " << bucket_b << " " << *gt_a << " " << *gt_b << endl;
//			if (covered_frac(bucket_a, bucket_b, *gt_a, *gt_b) >= 0.9) {
//				//cerr << "SAFE: " << query_id << " " << bucket_a << " " << bucket_b << " " << *gt_a << " " << *gt_b << endl;
//				return true;
//			}
//		}
//		return false;
//	}

//	int mapq_ed(int ed_best, int ed_second) {
//		if (ed_best == -1)
//			return 0;
//		if (ed_second == -1)
//			return 60;
//		assert(ed_best <= ed_second);
//		double bound = ed_best * 1.25;
//		double r = max((bound - ed_second) / (bound - ed_best), 0.0);  // small r is good
//		assert(r <= 1.0);
//		return int(60.0 * (1.0 - 1.0 * r));  // big score is good
//	}
	
//	int edit_distance(int b, const char *P, const pos_t P_sz, int lmax) {
//		int max_ed = -1;
//		char *s = tidx.T[0].seq.c_str() + (b*lmax);
//		int S_sz = 2*lmax;
//		auto config = edlibNewAlignConfig(max_ed, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0);
//		EdlibAlignResult result = edlibAlign(P, P_sz, s, S_sz, config);
//		assert(result.status == EDLIB_STATUS_OK);
//		int ed = result.editDistance;
//		
//		//string cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
//		//cerr << "S_sz=" << S_sz << ", P_sz=" << P_sz << ", edit distance: " << result.editDistance << ", mapping: " << *mapping << endl;
//		//cerr << cigar << endl;
//		edlibFreeAlignResult(result);
//		return ed;
//	}

	// returns cigar
//	string add_edit_distance(Mapping *mapping, const char *P, const pos_t P_sz) {
//		int max_ed = 1000; // -1 for none // TODO: make it a function of theta
//		int delta = 10000;  // TODO: remove delta
//
//		auto segm_id = mapping->segm_id;
//		int S_a = mapping->T_l, S_b = mapping->T_r;
//		char strand = mapping->strand;
//		auto &ref = tidx.T[segm_id].seq;
//		
//		// TODO: remove
//		S_a = max(0, S_a - delta);
//		S_b = min((int)ref.size()-1, S_b + delta);
//		
//		S_a -= H->params.k + 1;
//		assert(S_a >= 0 && S_a < (int)ref.size());
//		assert(S_b >= 0 && S_b < (int)ref.size());
//		
//		int S_sz = S_b - S_a;
//		const char *s;
//		string s_rev;
//		if (strand == '+')
//			s = ref.c_str() + S_a;
//		else {
//			assert(strand == '-');
//			s_rev = Mapping::reverseComplement(ref.substr(S_a, S_sz));
//			s = s_rev.c_str();
//		}
//
//		// edlib query: s
//		// edlib target: p
//		//auto config = edlibNewAlignConfig(max_ed, EDLIB_MODE_HW, EDLIB_TASK_DISTANCE, NULL, 0);
//		auto config = edlibNewAlignConfig(max_ed, EDLIB_MODE_HW, EDLIB_TASK_PATH, NULL, 0);
//		//EdlibAlignResult result = edlibAlign(s, S_sz, P, P_sz, config);
//		EdlibAlignResult result = edlibAlign(P, P_sz, s, S_sz, config);
//		assert(result.status == EDLIB_STATUS_OK);
//		mapping->ed = result.editDistance;
//		
//		string cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
//		//cerr << "S_sz=" << S_sz << ", P_sz=" << P_sz << ", edit distance: " << result.editDistance << ", mapping: " << *mapping << endl;
//		//cerr << cigar << endl;
//		edlibFreeAlignResult(result);
//		return cigar;
//	}
				//int best_ed_idx = -1, best2_ed_idx = -1; 
				//if (use_ed) {
				//	for (int i=0; i < (int)maps.size(); ++i) {
				//		auto &mapping = maps[i];
				//		add_edit_distance(&mapping, P, P_sz, m, kmers);
				//		if (mapping.ed != -1 && (best_ed_idx == -1 || mapping.ed < maps[best_ed_idx].ed)) {
				//			best2_ed_idx = best_ed_idx;
				//			best_ed_idx = i;
				//		} else if (mapping.ed != -1 && (best2_ed_idx == -1 || mapping.ed > maps[best2_ed_idx].ed)) {
				//			best2_ed_idx = i;
				//		}
				//	}
				//}

	//						if (best2_idx != -1 && final_mapping.mapq == 0) {
	//							string cigar1, cigar2;
	//							cigar1 = add_edit_distance(&maps[best_idx], P, P_sz, m, kmers);
	//							cigar2 = add_edit_distance(&maps[best2_idx], P, P_sz, m, kmers);
	//							if (maps[best_idx].ed != maps[best2_idx].ed) {
	//								cerr << endl;
	//								cerr << "mapq = 0 for read " << query_id << endl;
	//								cerr << maps[best_idx] << " " << cigar1 << endl;
	//								cerr << maps[best2_idx] << " " << cigar2 << endl;
	//							}
	//						}
	//				int gt_a, gt_b;
	//				lost_on_seeding = !is_safe(query_id, B_vec, lmax, &gt_a, &gt_b);
	//				if (!is_safe(query_id, B_vec, lmax, &gt_a, &gt_b))
	//					cerr << "Before bucket pruning, ground-truth mapping is lost: query_id=" << query_id << endl; 

//							if (use_ed) {
	//							//cerr << "best_ed_idx=" << best_ed_idx << " best2_ed_idx=" << best2_ed_idx << " final_mapping.ed=" << final_mapping.ed << endl;
	//							assert(best_ed_idx != -1);
	//							assert(final_mapping.ed != -1);
	//							final_mapping.ed2 = best2_ed_idx == -1 ? -1 : maps[best2_ed_idx].ed;
	//							final_mapping.mapq = mapq_ed(final_mapping.ed, final_mapping.ed2);
	//							//if (final_mapping.mapq > 0)
	//								final_mappings.push_back(final_mapping);
//							} else {
	//				if (maps.size() == 1) { // && maps.front().mapq > 0) {
	//					if (!is_safe(query_id, final_buckets, lmax, &gt_a, &gt_b)) {
	//						cerr << "After edit distance, ground-truth mapping is lost: query_id=" << query_id << endl;
	//						//cerr << "         mapq: " << maps.front().mapq << endl;
	//						if (best_idx != -1) {
	//							int bb = maps[best_idx].bucket;
	//							cerr << "         best: " << best_idx <<  ", bucket=" << bb << "[" << bb*lmax << ", " << (bb+2)*lmax << ")"; if (best_idx != -1) cerr << ", " << std::setprecision(5) << maps[best_idx] << endl; else cerr << endl;
	//						}
	//						if (best2_idx != -1) {
	//							int bb2 = maps[best2_idx].bucket;
	//							cerr << "  best2_idx: " << best2_idx << ", bucket=" << bb2 << "[" << bb2*lmax << ", " << (bb2+2)*lmax << ")"; if (best2_idx != -1) cerr << ", " << std::setprecision(5) << maps[best2_idx] << endl; else cerr << endl;
	//						}
	//					} else
	//						lost_on_pruning = 0;
	//				}

			
//	void no_sweep(const vector<Match> &B, const qpos_t P_sz, qpos_t lmax, const qpos_t m, const Seeds &kmers, vector<Mapping> *maps, rpos_t bucket, Hist diff_hist) {
//		int intersection = 0;
//		int same_strand_seeds = 0;
//		//Mapping best;
//		auto l = B.begin(), r = B.begin();
//		
//		for(auto it = B.begin(); it != B.end(); ++it) {
//			if (it->hit.tpos < l->hit.tpos) l = it;
//			if (it->hit.tpos > r->hit.tpos) r = it;
//
//			same_strand_seeds += it->is_same_strand() ? +1 : -1;
//			if (--diff_hist[it->seed.kmer.h] >= 0)
//				++intersection;
//		}
//
//		double J = 1.0*intersection / m;
//		assert(J <= 1.0);
//		Mapping mapping(H->params.k, P_sz, m, l->hit.r, r->hit.r, l->hit.segm_id, intersection, J, same_strand_seeds, l, r, bucket);
//		
//		//if (best.J >= H->params.theta)
//		if (mapping.J >= 0)
//			maps->push_back(mapping);
//	}


    //using Kmer = Seed;
//	using hist_t = vector<int>;
//	using Intervals = vector<pair<int, int>>;