#pragma once
#include <span>
#include <set>
#include <map>

#include "index.h"
#include "io.h"
#include "sketch.h"
#include "utils.h"
#include "handler.h"
#include "mapper.h"
#include "mapping.h"
#include "buckets.h"
#include "refine.h"
#include "types.h"

#include "analyse_simulated.h"

enum class MapResult {
    UNIQUE,
    AMBIG,
    NONE,
    _LAST // Sentinel value for getting enum size
};

inline string str(MapResult r) {
    switch(r) {
        case MapResult::UNIQUE: return "UNIQUE";
        case MapResult::AMBIG: return "AMBIG"; 
        case MapResult::NONE: return "NONE";
        default: throw std::invalid_argument("Invalid MapResult");
    }
}

namespace sweepmap {

class SHSingleReadMapper {
	const SketchIndex &tidx;
	Handler *H;
	Matcher &matcher;
	const string &query_id;
	const string &P;
	ofstream &paulout;

	Counters C;

	using SeedSpan = span<const Seed>;

public: // for testing
	// Returns (matched_kmers, m), where:
	// `matched_kmers` is the list of kmers in sketch `p` (without repetitions, sorted by increasing number of hits in `tidx`)
	// `m` is the total number of kmers in sketch `p` (including multiples)
	pair<Seeds, qpos_t> select_elements(sketch_t& p) {
		H->T.start("collect_element_info");
			int nonzero = 0;
			Seeds matched_elements;
			// After sorting, repetitive elements are next to each other.
			sort(p.begin(), p.end(), [](const El &a, const El &b) { return a.h < b.h; }); 
			qpos_t strike = 0;
			for (size_t ppos = 0; ppos < p.size(); ++ppos) {
				++strike;
				if (ppos == p.size()-1 || p[ppos].h != p[ppos+1].h) {
					rpos_t hits_in_t = tidx.count(p[ppos].h);
					if (hits_in_t > 0) {
						nonzero += strike;
					}
						if (H->params.max_matches == -1 || hits_in_t <= H->params.max_matches) {
							Seed el(p[ppos], hits_in_t, strike, matched_elements.size());
							//strike = 1; // comment out for Weighted metric
							matched_elements.push_back(el);
							if (H->params.max_seeds != -1 && (rpos_t)matched_elements.size() >= H->params.max_seeds)  // TODO maybe account for occs_in_p
								break;
						}
					//}
					strike = 0;
				}
			}

			qpos_t m = 0;
			for (const auto matched_el: matched_elements)
				m += matched_el.occs_in_p;
			assert(m <= (int)p.size());
			H->C.inc("elements_notmatched", m - nonzero);
		H->T.stop("collect_element_info");

		H->T.start("sort_elements");
			//pdqsort_branchless(kmers.begin(), kmers.end(), [](const Seed &a, const Seed &b) {
			sort(matched_elements.begin(), matched_elements.end(), [](const Seed &a, const Seed &b) {
				return a.hits_in_T < b.hits_in_T;
			});
		H->T.stop("sort_elements");

		return std::make_pair(matched_elements, m);
	}

	// Returns two spans:
	// - the first span contains the seeds such that at least 1 matches in any mapping with score >= theta
	// - the second span contains the remaining kmers
	// - the third integer is the number of seeds, including repeats
	static tuple<SeedSpan, SeedSpan, int> choose_seeds(const Seeds &all_elements, qpos_t m, double theta) {
		qpos_t S = qpos_t((1.0 - theta) * m) + 1;

		int i, seeds_with_repeats = 0;
		for (i = 0; i < (int)all_elements.size() && seeds_with_repeats < S; ++i)
			seeds_with_repeats += all_elements[i].occs_in_p;

		return {SeedSpan(all_elements.data(), i), 
		        SeedSpan(all_elements.data() + i, all_elements.size() - i),
				seeds_with_repeats};
	}

	// Takes at least `S` seeds from `kmers`. Matches each seed to the buckets `B`.  Returns the number of seeds.
	Buckets match_seeds(SeedSpan seeds, qpos_t lmax) {
		Buckets B(lmax);
		C["seed_matches"] = 0;
		C["max_seed_matches"] = 0;
		qpos_t seeds_with_repeats = 0;
		for (const auto &seed : seeds) {
				C["seed_matches"] += (int64_t)seed.hits_in_T;
				C["max_seed_matches"] = max(C["max_seed_matches"].count(), (int64_t)seed.hits_in_T);
				Matches seed_matches;
				tidx.get_matches(&seed_matches, seed);

				//auto B2 = B;	
				//////////////// correct B
				Buckets b2m(B.get_bucket_len());
				for (auto &m: seed_matches) {
					//rpos_t b(m.hit.tpos / B.get_bucket_len());
					b2m.add_to_pos(Buckets::Pos(m.hit.segm_id, m.hit.tpos));
					//if (b > 0) b2m.add_to_bucket(bucket_t(m.hit.segm_id, b-1), seed.occs_in_p);
				}
				for (auto it = b2m.unordered_begin(); it != b2m.unordered_end(); ++it) {
					//B[b] += min(seed.occs_in_p, matches);
					B.add_to_bucket(it->first, min(it->second, seed.occs_in_p));
					//matches_in_B += it->second;
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
			seeds_with_repeats += seed.occs_in_p;
		}
		C["seeded_buckets"] = B.size();
		H->C["seeded_buckets"] += C["seeded_buckets"];
		H->C["seed_matches"] += C["seed_matches"];

		return B;
	}

	double hseed(qpos_t p, qpos_t seeds, qpos_t matches) {
		return 1.0 - double(seeds - matches) / p;
	}

	bool seed_heuristic_pass(const vector<Mapping> &maps, SeedSpan remainingElemenets, qpos_t m, Buckets::Bucket b, rpos_t seed_matches,
							 qpos_t seeds, int best_idx, double *lowest_sh, double thr_init) {
		if (H->params.no_bucket_pruning)
			return true;

		H->T.start("seed_heuristic");
		*lowest_sh = 1.0; // should only get lower
		double thr2 = best_idx == -1 ? thr_init : maps[best_idx].score()-H->params.best_score_delta;  // TODO: add a parameter for 0.015

		if (hseed(m, seeds, seed_matches) < thr2) {
			H->T.stop("seed_heuristic");
			return false;
		}

		bool ret = true;
		for (const auto& el : remainingElemenets) {
			seeds += el.occs_in_p;
			auto hit_cnt = tidx.matches_in_bucket(el, b);
			if (hit_cnt > 0)
				seed_matches += min(hit_cnt, el.occs_in_p);
			double sh = hseed(m, seeds, seed_matches);
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
	
	void match_rest(qpos_t lmax, qpos_t m, SeedSpan remainingElemenets, const Buckets &B, Hist &diff_hist, int seeds_with_repeats,
					const unordered_map<hash_t, Seed> &p_ht, vector<Mapping> &maps, int &best_idx, double thr_init, int forbidden_idx) {	
		double best_J = -1.0;

		int lost_on_seeding = (0);
		H->C.inc("lost_on_seeding", lost_on_seeding);

		C["final_buckets"] = 0;
		H->T.start("seed_heuristic"); H->T.stop("seed_heuristic");  // init
		H->T.start("match_collect"); H->T.stop("match_collect");
		H->T.start("sweep"); H->T.stop("sweep");
		for (auto b_it = B.ordered_begin(); b_it != B.ordered_end(); ++b_it) {
			double lowest_sh;
			if (seed_heuristic_pass(maps, remainingElemenets, m, b_it->first, b_it->second, seeds_with_repeats, best_idx, &lowest_sh, thr_init)) {
				if (matcher.do_overlap(query_id, b_it->first))
					C["lost_on_pruning"] = 0;
				H->T.start("match_collect");
					Matches M = matcher.collect_matches(b_it->first, p_ht);
				H->T.stop("match_collect");
				C["total_matches"] += M.size();
				C.inc("final_buckets");

				H->T.start("sweep");
					auto best_in_bucket = matcher.bestFixedLength(M, P.size(), lmax, m, Metric::CONTAINMENT_INDEX);
					best_in_bucket.set_bucket(b_it->first);
					assert(best_in_bucket.score() <= lowest_sh + 1e-7);
					if (best_in_bucket.score() >= H->params.theta) {
						maps.emplace_back(best_in_bucket);
						if (best_in_bucket.score() > best_J) {
							if (forbidden_idx == -1 || Mapping::overlap(maps.back(), maps[forbidden_idx]) < 0.5) {
								best_J = max(best_J, best_in_bucket.score());
								best_idx = maps.size()-1;
							}
						}
					}
				H->T.stop("sweep");
			}
		}
		//assert(best_idx != -1 || maps.empty());

        H->C.inc("final_mappings", maps.size());
        H->C["total_matches"] += C["total_matches"];
        H->C["final_buckets"] += C["final_buckets"];
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

    tuple<vector<Mapping>, int, int, double> find_best_and_best2(qpos_t lmax, qpos_t m, SeedSpan remainingElemenets,
			const Buckets &B, Hist &diff_hist, int seeds_with_repeats, const unordered_map<hash_t, Seed> &p_ht) {
        vector<Mapping> maps;
        //vector<Bucket> B_vec(B.begin(), B.end());
        //sort(B_vec.begin(), B_vec.end(), [](const Bucket &a, const Bucket &b) { return a.second > b.second; });  // TODO: sort intervals by decreasing number of matches
        C["total_matcher"] = C["seed_matches"];
        int best_idx=-1, best2_idx=-1;
		match_rest(lmax, m, remainingElemenets, B, diff_hist, seeds_with_repeats, p_ht, maps, best_idx, H->params.theta, -1);
		auto [FP, PP, FDR, FPTP] = tuple(-1, -1, -1, -1);
		//auto [FP, PP, FDR, FPTP] = calc_FDR(maps, H->params.theta, lmax, bucket_l, tidx, p_ht, P_sz, lmin, m, diff_hist);
		H->C.inc("FP", FP);
		H->C.inc("PP", PP);
		H->C.inc("FDR", FDR);
		H->C.inc("FPTP", FPTP);

		if (best_idx != -1) {
			// find second best mapping for mapq computation
			match_rest(lmax, m, remainingElemenets, B, diff_hist, seeds_with_repeats, p_ht, maps, best2_idx, maps[best_idx].score(), best_idx);
		}
		return {maps, best_idx, best2_idx, FPTP};
	}

	MapResult output(const qpos_t seeds_with_repeats, vector<Mapping> &maps, int best_idx, int best2_idx, double FPTP) {
		if (maps.size() >= 1) {
			H->C.inc("mapped_reads");
			if (best_idx != -1) 
			//if (maps[best_idx].J >= H->params.theta)
			{
				assert(0 <= best_idx && best_idx < (int)maps.size());
				Mapping &m = maps[best_idx];

				const auto &segm = tidx.T[m.segm_id()];
				m.set_global_stats(C, query_id.c_str(), P.size(), seeds_with_repeats, FPTP, segm.name, segm.sz, H->T.secs("query_mapping"));  // TODO: disable by flag

				if (best2_idx != -1)
					m.set_second_best(maps[best2_idx]);

				m.set_mapq(H->params.theta, H->params.best_score_delta);
				m.print_paf();

				if (m.mapq() == 60) H->C.inc("mapq60");
				H->C.inc("matches_in_reported_mappings", m.intersection());
				H->C.inc("J_best", rpos_t(10000.0*m.score()));
				H->C.inc("mappings");

				return m.mapq() > 0 ? MapResult::UNIQUE : MapResult::AMBIG;
			}
		}
		return MapResult::NONE;
	}

public:
	SHSingleReadMapper(const SketchIndex &tidx, Handler *H, Matcher &matcher, const string &query_id, const string &P, ofstream &paulout)
		: tidx(tidx), H(H), matcher(matcher), query_id(query_id), P(P), paulout(paulout) {}

	// returns true if the read was mapped
	MapResult map_read(double theta) {
		assert (theta >= H->params.theta - 1e-7);
		assert (theta <= 1.0 + 1e-7);

		H->T.start("query_mapping");

		H->T.start("sketching");
			sketch_t p = H->sketcher.sketch(P);
		H->T.stop("sketching");

		H->T.start("select_elements");
			auto [elements, m] = select_elements(p);
		H->T.stop("select_elements");

		H->T.start("prepare");
			H->C.inc("read_len", P.size());
			//cerr << "notmatched: " << m - nonzero << endl;

			unordered_map<hash_t, Seed> p_ht;
			Hist diff_hist;
			rpos_t possible_matches(0);
			for (const auto el: elements) {
				p_ht.insert(make_pair(el.el.h, el));
				//diff_hist[el.el.h] = 1;
				diff_hist[el.el.h] = el.occs_in_p;
				possible_matches += el.hits_in_T;
			}
			H->C.inc("possible_matches", possible_matches);
			matcher.update(diff_hist);

			//qpos_t lmax = qpos_t(m / H->params.theta);					// maximum length of a similar mapping
			//qpos_t lmin = qpos_t(ceil(p.size() * H->params.theta));					// maximum length of a similar mapping
			qpos_t lmax = qpos_t(p.size() / theta);					// maximum length of a similar mapping
			//qpos_t bucket_l = lmax;

			//double coef = 1.0;// * nonzero / p.size();
			//cerr << "coef: " << coef << endl;
			//double new_theta = theta;
			//cerr << "theta: " << H->params.theta << " -> " << new_theta << endl;

			H->C.inc("elements_sketched", p.size());
			H->C.inc("elements", m);
			H->C.inc("elements_unique", elements.size());
		H->T.stop("prepare");

		H->T.start("choose_seeds");
			auto [seeds, remainingElemenets, seeds_with_repeats] = choose_seeds(elements, m, theta);
			H->C.inc("elements_seeds", seeds.size());
		H->T.stop("choose_seeds");

		H->T.start("match_seeds");
			Buckets B = match_seeds(seeds, lmax);
		H->T.stop("match_seeds");

		int lost_on_seeding = matcher.lost_correct_mapping(query_id, B);
		H->C.inc("lost_on_seeding", lost_on_seeding);

		C["lost_on_pruning"] = 1;
		H->T.start("match_rest");
			auto [maps, best_idx, best2_idx, FPTP] = find_best_and_best2(lmax, m, remainingElemenets, B, diff_hist, seeds_with_repeats, p_ht);
		H->T.stop("match_rest");
		H->C["lost_on_pruning"] += C["lost_on_pruning"];

		H->T.start("extra");
		//if (!H->params.paramsFile.empty()) {
		//	// todo: comment out
		//	AnalyseSimulatedReads sim(query_id, P, P_sz, diff_hist, m, p_ht, tidx, B, H->params.theta);
		//	sim.PaulsExperiment(paulout);
		//}
		H->T.stop("extra");

		H->T.stop("query_mapping");

		H->T.start("output");
			// total_matches, max_seed_matches, seed_matches, seeded_buckets
			MapResult map_result = output(seeds_with_repeats, maps, best_idx, best2_idx, FPTP);
		H->T.stop("output");

		H->C[str(map_result)] += 1;
		return map_result;
	}
};

class SHMapper : public Mapper {
	const SketchIndex &tidx;
	Handler *H;
	Matcher matcher;
	Counters C;  // per mapping counters

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

	void map_all_reads(const string &pFile) {
		cerr << "Mapping reads using SHmap..." << endl;

		H->C.inc("elements", 0);
		H->C.inc("elements_notmatched", 0);
		H->C.inc("seeds", 0);
		H->C.inc("matches", 0);
		H->C.inc("seed_matches", 0);
		H->C.inc("matches_freq", 0);
		H->C.inc("spurious_matches", 0);
		H->C.inc("mappings", 0);
		H->C.inc("J_best", 0);
		//H->C.inc("sketched_kmers", 0);
		H->C.inc("total_edit_distance", 0);
		H->C.inc("intersection_diff", 0);
		H->C.inc("mapq60", 0);
		H->C.inc("matches_in_reported_mappings", 0);
		H->C.inc("lost_on_seeding", 0);
		H->C.inc("lost_on_pruning", 0);

		H->C["repeated_mappings"] = 0;
		for (int i = 0; i < static_cast<int>(MapResult::_LAST); ++i)
			H->C.inc(str(static_cast<MapResult>(i)), 0);

		H->T.start("mapping");
		H->T.start("query_reading");

		string pauls_fn = H->params.paramsFile + ".paul.tsv";
		ofstream paulout;
		if (!H->params.paramsFile.empty()) {
			cerr << "Paul's experiment to " << pauls_fn << endl;
			paulout = ofstream(pauls_fn);
		}

		read_fasta_klib(pFile, [this, &paulout](const string &query_id, const string &P) {
			H->C.inc("reads");
			H->T.stop("query_reading");

			SHSingleReadMapper mapper(tidx, H, matcher, query_id, P, paulout);
			//mapper.map_read(H->params.theta);

			//for (double theta = 0.85; theta >= H->params.theta; theta -= 1.0)
			//	if (mapper.map_read(theta) != MapResult::NONE)
			//		break;

			if (mapper.map_read(0.8) == MapResult::NONE) {
				mapper.map_read(H->params.theta);
				H->C["repeated_mappings"] += 1;
			}

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
		cerr << " | Elements:              " << H->C.frac("elements", "reads") << " p/ read" << endl;
		cerr << " |  | sketched:               " << H->C.frac("elements_sketched", "reads") << " (" << H->C.perc("elements_sketched", "elements") << "%)" << endl;
		cerr << " |  | not matched:            " << H->C.frac("elements_notmatched", "reads") << " (" << H->C.perc("elements_notmatched", "elements") << "%)" << endl;
		cerr << " |  | unique:                 " << H->C.frac("elements_unique", "reads") << " (" << H->C.perc("elements_unique", "elements") << "%)" << endl;
		cerr << " |  | seeds:                  " << H->C.frac("elements_seeds", "reads") << " (" << H->C.perc("elements_seeds", "elements") << "%)" << endl;
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
		cerr << " | | Mapped results:          ";
		for (int i = 0; i < static_cast<int>(MapResult::_LAST); ++i)
			cerr << str(static_cast<MapResult>(i)) << ": " << H->C.count(str(static_cast<MapResult>(i))) << ", ";
		cerr << endl;
		cerr << " | | Repeated mappings:       " << H->C.count("repeated_mappings") << endl;

		//cerr << " | Average edit dist:     " << H->C.frac("total_edit_distance", "mappings") << endl;
		//print_time_stats();
        //printMemoryUsage();
	}

    void print_time_stats() {
        cerr << std::fixed << std::setprecision(1);
        cerr << " | Runtime:                "    << setw(5) << right << H->T.secs("mapping")       << " sec, " << 1.0 * H->C.count("reads") / H->T.secs("mapping")  << " reads/sec (" << setw(5) << right << H->T.range_ratio("query_mapping") << "x)" << endl; //setw(4) << right << H->C.count("reads") / H->T.secs("total") << " reads p/ sec)" << endl;
        cerr << " |  | load reads:             " << setw(5) << right << H->T.secs("query_reading")     << " (" << setw(4) << right << H->T.perc("query_reading", "mapping")       << "\%, " << setw(6) << right << H->T.range_ratio("query_reading") << "x)" << endl;
        cerr << " |  | sketch reads:           " << setw(5) << right << H->T.secs("sketching")         << " (" << setw(4) << right << H->T.perc("sketching", "mapping")           << "\%, " << setw(6) << right << H->T.range_ratio("sketching") << "x)" << endl;
		cerr << " |  | select elements:        " << setw(5) << right << H->T.secs("select_elements")   << " (" << setw(4) << right << H->T.perc("select_elements", "mapping") << "\%, " << setw(6) << right << H->T.range_ratio("select_elements") << "x)" << endl;
        //cerr << " |  |  | collect element info:    " << setw(5) << right << H->T.secs("collect_element_info") << " (" << setw(4) << right << H->T.perc("collect_element_info", "seeding")   << "\%, " << setw(6) << right << H->T.range_ratio("collect_element_info") << "x)" << endl;
        //cerr << " |  |  | sort elements:          " << setw(5) << right << H->T.secs("sort_elements")        << " (" << setw(4) << right << H->T.perc("sort_elements", "seeding")          << "\%, " << setw(6) << right << H->T.range_ratio("sort_elements") << "x)" << endl;
//        cerr << " |  |  | thin sketch:             " << setw(5) << right << H->T.secs("thin_sketch")       << " (" << setw(4) << right << H->T.perc("thin_sketch", "seeding")         << "\%, " << setw(6) << right << H->T.range_ratio("thin_sketch") << "x)" << endl;
//        cerr << " |  |  | unique seeds:            " << setw(5) << right << H->T.secs("unique_seeds")      << " (" << setw(4) << right << H->T.perc("unique_seeds", "seeding")        << "\%, " << setw(6) << right << H->T.range_ratio("unique_seeds") << "x)" << endl;
        cerr << " |  | prepare:                "     << setw(5) << right << H->T.secs("prepare")           << " (" << setw(4) << right << H->T.perc("prepare", "mapping")             << "\%, " << setw(6) << right << H->T.range_ratio("prepare") << "x)" << endl;
        cerr << " |  | match seeds:            " << setw(5) << right << H->T.secs("match_seeds")  << " (" << setw(4) << right << H->T.perc("match_seeds", "mapping")   << "\%, " << setw(6) << right << H->T.range_ratio("match_seeds") << "x)" << endl;
        cerr << " |  | match rest:             " << setw(5) << right << H->T.secs("match_rest")   << " (" << setw(4) << right << H->T.perc("match_rest", "mapping")     << "\%, " << setw(6) << right << H->T.range_ratio("match_rest") << "x)" << endl;
        cerr << " |  |  | seed heuristic:         " << setw(5) << right << H->T.secs("seed_heuristic")    << " (" << setw(4) << right << H->T.perc("seed_heuristic", "match_rest")     << "\%, " << setw(6) << right << H->T.range_ratio("seed_heuristic") << "x)" << endl;
        cerr << " |  |  | matches collect:        " << setw(5) << right << H->T.secs("match_collect")     << " (" << setw(4) << right << H->T.perc("match_collect", "match_rest")      << "\%, " << setw(6) << right << H->T.range_ratio("match_collect") << "x)" << endl;
        cerr << " |  |  | sweep:                  " << setw(5) << right << H->T.secs("sweep")             << " (" << setw(4) << right << H->T.perc("sweep", "match_rest")              << "\%, " << setw(6) << right << H->T.range_ratio("sweep") << "x)" << endl;

//        cerr << " |  |  | get intervals:           " << setw(5) << right << H->T.secs("get_intervals")     << " (" << setw(4) << right << H->T.perc("get_intervals", "mapping")      << "\%, " << setw(5) << right << H->T.range_ratio("get_intervals") << "x)" << endl;
//        cerr << " |  |  | filter promising buckets:" << setw(5) << right << H->T.secs("filter_promising_buckets")    << " (" << setw(4) << right << H->T.perc("filter_promising_buckets", "mapping")     << "\%, " << setw(5) << right << H->T.range_ratio("filter_promising_buckets") << "x)" << endl;
//        cerr << " |  | edit distance:          "     << setw(5) << right << H->T.secs("edit_distance")     << " (" << setw(4) << right << H->T.perc("edit_distance", "mapping")       << "\%, " << setw(6) << right << H->T.range_ratio("edit_distance") << "x)" << endl;
        cerr << " |  | output:                 " << setw(5) << right << H->T.secs("output")            << " (" << setw(4) << right << H->T.perc("output", "mapping")              << "\%, " << setw(6) << right << H->T.range_ratio("output") << "x)" << endl;
        cerr << " |  | extra:                  " << setw(5) << right << H->T.secs("extra")             << " (" << setw(4) << right << H->T.perc("extra", "mapping")               << "\%, " << setw(6) << right << H->T.range_ratio("extra") << "x)" << endl;
//        cerr << " |  | post proc:              "     << setw(5) << right << H->T.secs("postproc")          << " (" << setw(4) << right << H->T.perc("postproc", "mapping")            << "\%, " << setw(6) << right << H->T.range_ratio("postproc") << "x)" << endl;
    }
};

}  // namespace sweepmap

					//const char *P = seq->seq.s;
					//char *query_id = seq->name.s;
					//qpos_t P_sz = (qpos_t)seq->seq.l;