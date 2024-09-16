#pragma once

#include <algorithm>
#include <cassert>
#include <deque>
#include <iomanip>
#include <set>
#include <unordered_set>

//#include "../ext/edlib.h"
#include "../ext/pdqsort.h"

#include "index.h"
#include "io.h"
#include "sketch.h"
#include "utils.h"
#include "handler.h"
#include "mapper.h"
#include "rmq.h"
//#include "fenwick.h"

namespace sweepmap {

using std::setw;
using std::right;
using std::vector;

class RMQMapper : public Mapper {
	const SketchIndex &tidx;
	Handler *H;
	SegmentTree hist;

	using hist_t = vector<int>;
	using Seeds = vector<Seed>;
	using Matches = vector<Match>;
	using Intervals = vector<pair<int, int>>;

	Seeds select_seeds(sketch_t& p) {
		H->T.start("unique_kmers");
		pdqsort_branchless(p.begin(), p.end(), [](const Kmer &a, const Kmer &b) {
            return a.h < b.h;
		});
		p.erase(std::unique(p.begin(), p.end()), p.end());
		H->T.stop("unique_kmers");

		H->T.start("collect_seed_info");
		Seeds seeds;
		seeds.reserve(p.size());

		// TODO: limit The number of kmers in the pattern p
		for (int ppos = 0; ppos < (int)p.size(); ++ppos) {
			const auto &kmer = p[ppos];
			const auto count = tidx.count(kmer.h);
			if (count > 0)
				seeds.push_back(Seed(kmer, p[ppos].r, p[ppos].r, count, seeds.size()));
		}
		H->T.stop("collect_seed_info");
        H->C.inc("collected_seeds", seeds.size());

		H->T.start("sort_seeds");
		pdqsort_branchless(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b) {
            return a.hits_in_T < b.hits_in_T;
		});
		H->T.stop("sort_seeds");

		return seeds;
	}
	
	// Initializes the histogram of the pattern and the list of matches
	Matches match_infreq_seeds(const Seeds &seeds_infreq) {
		Matches matches_infreq;
		for (const auto seed: seeds_infreq)
			if (seed.hits_in_T > 0)
				tidx.get_matches(&matches_infreq, seed);
        sort(matches_infreq.begin(), matches_infreq.end(), [](const Match &a, const Match &b) {
			return a.hit.r < b.hit.r;
		});
		return matches_infreq;
	}

	Intervals get_intervals_of_matches(const Matches &matches_infreq) {
        Intervals intervals_infreq;
		if (matches_infreq.empty())
			return intervals_infreq;

		int l = hist.from(matches_infreq[0].hit);
		int r = hist.to(matches_infreq[0].hit);
		for (int i=1; i<(int)matches_infreq.size(); i++) {
			int curr_l = hist.from(matches_infreq[i].hit);
			int curr_r = hist.to(matches_infreq[i].hit);
			assert(curr_l <= curr_r);
			if (curr_l <= r) {
				r = curr_r;
			} else {
				intervals_infreq.push_back({l, r});
				l = curr_l;
				r = curr_r;
			}
		}
		intervals_infreq.push_back({l, r});

		return intervals_infreq;
	}

	tuple<Matches, int> match_frequent_seeds(const Seeds &seeds_freq, int t_abs, const Intervals &intervals_infreq) {
		Matches matches_freq;

        vector<Match> freq_matches;
		// match frequent seeds to all buckets with an infrequent seed
        for (auto [l, r]: intervals_infreq) {
            int max_matches = hist.query(l, r);
			int rem_seeds = (int)seeds_freq.size();
            for (const auto &seed: seeds_freq) {
				if (max_matches + rem_seeds < t_abs)
					break;
                max_matches = tidx.match_seed_around_hit(&hist, seed, l, r, &matches_freq);
				--rem_seeds;
			}
            if (max_matches > t_abs) {
                //cerr << "Better max_matches: " << t_abs << " -> " << max_matches << endl;
                t_abs = max_matches;
            }
        }

		return {matches_freq, t_abs};
	}

	const vector<Mapping> sweep(vector<Match> &M, const pos_t P_sz, const int p_sz) {
//		const int MAX_BL = 100;
		vector<Mapping> mappings;	// List of tripples <i, j, score> of matches

        unordered_multiset<int> diff_hist;
		int intersection = 0;
		Mapping best(H->params.k, P_sz, 0, -1, -1, -1, -1, -1, 0, M.end(), M.end());
		Mapping second = best;
		int same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise

		sort(M.begin(), M.end(), [](const Match &a, const Match &b) { return a.hit.r < b.hit.r; }); 

		cerr << "sweep" << endl;
		// print all matches
		//for (int m=0; m<(int)M.size(); m++) {
		//	if (m>0 && M[m].hit.r - M[m-1].hit.r>P_sz) {
		//		cerr << "-----------------" << endl;
		//	}
		//	cerr << M[m] << endl;
		//}

		// Increase the left point end of the window [l,r) one by one. O(matches)
		for(auto l = M.begin(), r = M.begin(); l != M.end(); ++l) {
			// Increase the right end of the window [l,r) until it gets out.
			for(;  r != M.end()
				&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				&& r->hit.r + H->params.k <= l->hit.r + P_sz
				; ++r) {
				same_strand_seeds += r->is_same_strand() ? +1 : -1;  // change to r inside the loop
				// If taking this kmer from T increases the intersection with P.
				// TODO: iterate following seeds
				if (!diff_hist.contains(r->seed.seed_num))
					++intersection;
				diff_hist.insert(r->seed.seed_num);
				assert (l->hit.r <= r->hit.r);
			}

			auto m = Mapping(H->params.k, P_sz, p_sz, l->hit.r, prev(r)->hit.r, l->hit.segm_id, pos_t(r-l), intersection, same_strand_seeds, l, r);

			// second best without guarantees
			// Wrong invariant:
			// best[l,r) -- a mapping best.l<=l with maximal J
			// second_best[l,r) -- a mapping second_best.l \notin [l-90%|P|; l+90%|P|] with maximal J
			//if (m.J >= H->params.tThres) {
				if (H->params.onlybest) {
					if (m.intersection > best.intersection) {  // if (intersection > best.intersection)
						if (best.T_l < m.T_l - 0.9*P_sz)
							second = best;
						best = m;
					} else if (m.intersection > second.intersection && m.T_l > best.T_l + 0.9*P_sz) {
						second = m;
					}
				} else {
					mappings.push_back(m);
				}
			//}

			// Prepare for the next step by moving `l` to the right.
			auto it = diff_hist.find(l->seed.seed_num);
			assert(it != diff_hist.end());
			diff_hist.erase(it); 
			if (!diff_hist.contains(l->seed.seed_num))
				--intersection;
			same_strand_seeds -= l->is_same_strand() ? +1 : -1;

			assert(intersection >= 0);
		}
		assert(intersection == 0);
		assert(same_strand_seeds == 0);

		if (H->params.onlybest && best.intersection != -1) { // && best.J > H->params.tThres)
			best.mapq = (best.intersection > 5 && best.J > 0.1 && best.J > 1.2*second.J) ? 60 : 0;
			best.J2 = second.J;
			mappings.push_back(best);
		}

		return mappings;
	}


  public:
	RMQMapper(const SketchIndex &tidx, Handler *H) : tidx(tidx), H(H), hist(tidx.T[0].sz) {
        H->C.inc("seeds_limit_reached", 0);
        H->C.inc("mapped_reads", 0);
        if (H->params.tThres < 0.0 || H->params.tThres > 1.0) {
            cerr << "tThres = " << H->params.tThres << " outside of [0,1]." << endl;
            exit(1);
        }
    }

	void map(const string &pFile) {
		cerr << "Mapping reads using BucketMap " << pFile << "..." << endl;

		H->C.inc("kmers", 0);
		H->C.inc("seeds", 0);
		H->C.inc("matches", 0);
		H->C.inc("matches_infreq", 0);
		H->C.inc("matches_freq", 0);
		H->C.inc("spurious_matches", 0);
		H->C.inc("J", 0);
		H->C.inc("mappings", 0);
		H->C.inc("sketched_kmers", 0);
		H->C.inc("total_edit_distance", 0);
		H->C.inc("intersection_diff", 0);

		H->T.start("mapping");
		H->T.start("query_reading");
        //fenwick_tree<int> hist(tidx.T[0].sz);

		read_fasta_klib(pFile, [this](kseq_t *seq) {
			H->C.inc("reads");
			H->T.stop("query_reading");
			H->T.start("query_mapping");
				H->T.start("sketching");
					sketch_t p = H->sketcher.sketch(seq->seq.s);
					H->T.stop("sketching");

				string query_id = seq->name.s;
				pos_t P_sz = (pos_t)seq->seq.l;
				hist.clear();
				hist.P_sz = P_sz;  // TODO: remove this hack

				H->C.inc("kmers", p.size());
				H->C.inc("read_len", P_sz);

				Timer read_mapping_time;  // TODO: change to H->T.start
				read_mapping_time.start();
				H->T.start("seeding");
					Seeds seeds = select_seeds(p);
					int t_abs = H->params.tThres * seeds.size();  // flooring is safe
					int number_of_infreq_seeds = (1.0 - H->params.tThres) * seeds.size();  // flooring is safe
					Seeds seeds_infreq(seeds.begin(), seeds.begin() + number_of_infreq_seeds);
					Seeds seeds_freq(seeds.begin() + number_of_infreq_seeds, seeds.end());
					H->T.stop("seeding");
					H->C.inc("seeds", seeds.size());

				H->T.start("matching");
					H->T.start("match_infrequent");
						Matches matches_infreq = match_infreq_seeds(seeds_infreq);
						H->T.stop("match_infrequent");

					H->T.start("get_intervals");
						Intervals intervals_infreq = get_intervals_of_matches(matches_infreq);
						//for (int i=0; i<(int)intervals_infreq.size(); i++) {
						//	if (i>0) assert(intervals_infreq[i-1].second <= intervals_infreq[i].first);
						//	cerr << "Interval: " << intervals_infreq[i].first << ", " << intervals_infreq[i].second << endl;
						//}
						H->T.stop("get_intervals");

					H->T.start("match_frequent");
						auto [matches_freq, max_t_abs] = match_frequent_seeds(seeds_freq, t_abs, intervals_infreq);
						H->T.stop("match_frequent");

					vector<Match> matches;
					matches.insert(matches.end(), matches_infreq.begin(), matches_infreq.end());
					matches.insert(matches.end(), matches_freq.begin(), matches_freq.end());

					H->C.inc("matches", matches.size());
					H->C.inc("matches_infreq", matches_infreq.size());
					H->C.inc("matches_freq", matches_freq.size());
					H->T.stop("matching");

				H->T.start("sweep");
				vector<Mapping> mappings = sweep(matches, P_sz, seeds.size());
				if (mappings.size() > 0) {
					H->C.inc("mapped_reads");
					H->C.inc("intersection_diff", t_abs - mappings[0].intersection);
				}
				H->T.stop("sweep");

				cerr << query_id << ": seeds: " << seeds.size() << ", I: " << int(H->params.tThres * seeds.size()) << " -> " << t_abs << ", matches: " << matches.size() << ", matches_infreq: " << matches_infreq.size() << ", matches_freq: " << matches_freq.size() << ", mappings: " << mappings.size();
				cerr << mappings[0] << endl;

				read_mapping_time.stop();

				for (auto &m: mappings) {
					m.map_time = read_mapping_time.secs() / (double)mappings.size();
					const auto &segm = tidx.T[m.segm_id];
	//				if (H->params.sam) {
	//					auto ed = m.print_sam(query_id, segm, (int)matches.size(), seq->seq.s, seq->seq.l);
	//					H->C.inc("total_edit_distance", ed);
	//				}
	//				else
						m.print_paf(query_id, segm, matches.size());
					//  H->C.inc("spurious_matches", spurious_matches(m, matches));
					H->C.inc("J", int(10000.0*m.J));
					H->C.inc("mappings");
				}
	//			H->C.inc("matches", matches.size());
				//H->T.stop("postproc");

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
		cerr << "Mapping:" << endl;
		cerr << " | Total reads:           " << H->C.count("reads") << " (~" << 1.0*H->C.count("read_len") / H->C.count("reads") << " nb per read)" << endl;
		cerr << " |  | mapped:               " << H->C.count("mapped_reads") << " (" << H->C.perc("mapped_reads", "reads") << "%)" << endl;
		cerr << " |  |  | intersect. diff:     " << H->C.frac("intersection_diff", "mapped_reads") << " per mapped read" << endl;
		cerr << " | Read kmers (total):    " << H->C.count("kmers") << " (" << H->C.frac("kmers", "reads") << " per read)" << endl;
		cerr << " |  | unique:                 " << H->C.count("seeds") << " (" << H->C.frac("seeds", "kmers") << ")" << endl;
		cerr << " | Matches:               " << H->C.count("matches") << " (" << H->C.frac("matches", "reads") << " per read)" << endl;
		cerr << " |  | infrequent:             " << H->C.count("matches_infreq") << " (" << H->C.perc("matches_infreq", "matches") << "%)" << endl;
		cerr << " |  | frequent:               " << H->C.count("matches_freq") << " (" << H->C.perc("matches_freq", "matches") << "%)" << endl;
		//cerr << " | Seed limit reached:    " << H->C.count("seeds_limit_reached") << " (" << H->C.perc("seeds_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Matches limit reached: " << H->C.count("matches_limit_reached") << " (" << H->C.perc("matches_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Spurious matches:      " << H->C.count("spurious_matches") << " (" << H->C.perc("spurious_matches", "matches") << "%)" << endl;
		//cerr << " | Discarded seeds:       " << H->C.count("discarded_seeds") << " (" << H->C.perc("discarded_seeds", "collected_seeds") << "%)" << endl;
		//cerr << " | Average Jaccard:       " << H->C.frac("J", "mappings") / 10000.0 << endl;
		//cerr << " | Average edit dist:     " << H->C.frac("total_edit_distance", "mappings") << endl;
		//print_time_stats();
        //printMemoryUsage();
	}

    void print_time_stats() {
        cerr << std::fixed << std::setprecision(1);
        cerr << " | Runtime:                   "     << setw(5) << right << H->T.secs("mapping")           << " (" << setw(5) << right << H->T.range_ratio("query_mapping") << "x)" << endl; //setw(4) << right << H->C.count("reads") / H->T.secs("total") << " reads per sec)" << endl;
        cerr << " |  | load queries:           "     << setw(5) << right << H->T.secs("query_reading")     << " (" << setw(4) << right << H->T.perc("query_reading", "mapping")       << "\%, " << setw(5) << right << H->T.range_ratio("query_reading") << "x)" << endl;
        cerr << " |  | sketch reads:           "     << setw(5) << right << H->T.secs("sketching")         << " (" << setw(4) << right << H->T.perc("sketching", "mapping")           << "\%, " << setw(5) << right << H->T.range_ratio("sketching") << "x)" << endl;
        cerr << " |  | seeding:                "     << setw(5) << right << H->T.secs("seeding")           << " (" << setw(4) << right << H->T.perc("seeding", "mapping")             << "\%, " << setw(5) << right << H->T.range_ratio("seeding") << "x)" << endl;
        cerr << " |  |  | collect seed info:       " << setw(5) << right << H->T.secs("collect_seed_info") << " (" << setw(4) << right << H->T.perc("collect_seed_info", "seeding")   << "\%, " << setw(5) << right << H->T.range_ratio("collect_seed_info") << "x)" << endl;
//        cerr << " |  |  | thin sketch:             " << setw(5) << right << H->T.secs("thin_sketch")       << " (" << setw(4) << right << H->T.perc("thin_sketch", "seeding")         << "\%, " << setw(5) << right << H->T.range_ratio("thin_sketch") << "x)" << endl;
//        cerr << " |  |  | unique seeds:            " << setw(5) << right << H->T.secs("unique_seeds")      << " (" << setw(4) << right << H->T.perc("unique_seeds", "seeding")        << "\%, " << setw(5) << right << H->T.range_ratio("unique_seeds") << "x)" << endl;
        cerr << " |  |  | sort seeds:              " << setw(5) << right << H->T.secs("sort_seeds")        << " (" << setw(4) << right << H->T.perc("sort_seeds", "seeding")          << "\%, " << setw(5) << right << H->T.range_ratio("sort_seeds") << "x)" << endl;
        cerr << " |  | matching seeds:         "     << setw(5) << right << H->T.secs("matching")          << " (" << setw(4) << right << H->T.perc("matching", "mapping")            << "\%, " << setw(5) << right << H->T.range_ratio("matching") << "x)" << endl;
        cerr << " |  |  | match infrequent:        " << setw(5) << right << H->T.secs("match_infrequent")  << " (" << setw(4) << right << H->T.perc("match_infrequent", "matching")   << "\%, " << setw(5) << right << H->T.range_ratio("match_infrequent") << "x)" << endl;
        cerr << " |  |  | match frequent:          " << setw(5) << right << H->T.secs("match_frequent")    << " (" << setw(4) << right << H->T.perc("match_frequent", "matching")     << "\%, " << setw(5) << right << H->T.range_ratio("match_frequent") << "x)" << endl;
//        cerr << " |  |  | filter promising buckets:" << setw(5) << right << H->T.secs("filter_promising_buckets")    << " (" << setw(4) << right << H->T.perc("filter_promising_buckets", "matching")     << "\%, " << setw(5) << right << H->T.range_ratio("filter_promising_buckets") << "x)" << endl;
        cerr << " |  | sweep:                  "     << setw(5) << right << H->T.secs("sweep")             << " (" << setw(4) << right << H->T.perc("sweep", "mapping")               << "\%, " << setw(5) << right << H->T.range_ratio("sweep") << "x)" << endl;
//        cerr << " |  | post proc:              "     << setw(5) << right << H->T.secs("postproc")          << " (" << setw(4) << right << H->T.perc("postproc", "mapping")            << "\%, " << setw(5) << right << H->T.range_ratio("postproc") << "x)" << endl;
    }
};

}  // namespace sweepmap