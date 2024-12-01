#pragma once

#include <algorithm>
#include <cassert>
#include <deque>
#include <iomanip>
#include <set>

//#include "../ext/edlib.h"
#include "../ext/pdqsort.h"

#include "index.h"
#include "io.h"
#include "sketch.h"
#include "utils.h"
#include "handler.h"
#include "mapper.h"

namespace sweepmap {

using std::setw;
using std::right;
using std::vector;

class SweepMapper : public Mapper {
	const SketchIndex &tidx;
	Handler *H;

	using hist_t = vector<int>;

	vector<Seed> select_seeds(const sketch_t& p, hist_t *hist) {
		H->T.start("collect_seed_info");
		vector<Seed> seeds;
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

		H->T.start("thin_sketch");
		// TODO: add all seeds to hist
		int total_seeds = (int)seeds.size();
		if ((int)H->params.max_seeds != -1 && (int)H->params.max_seeds < total_seeds) {
			total_seeds = (int)H->params.max_seeds;
			H->C.inc("seeds_limit_reached");
		}
        std::nth_element(seeds.begin(), seeds.begin() + total_seeds, seeds.end(), [](const Seed &a, const Seed &b) {
            return a.hits_in_T < b.hits_in_T;
        });
		H->T.stop("thin_sketch");

		H->T.start("sort_seeds");
		pdqsort_branchless(seeds.begin(), seeds.begin() + total_seeds, [](const Seed &a, const Seed &b) {
			return a.kmer.h < b.kmer.h;
		});
		H->T.stop("sort_seeds");

		H->T.start("unique_seeds");
		vector<Seed> thin_seeds;
		thin_seeds.reserve(total_seeds);
		hist->reserve(total_seeds+1);
		hist->push_back(0);
		int min_r = p.size(), max_r = -1;
		for (int i=0; i<total_seeds; i++) {
			min_r = std::min(min_r, seeds[i].r_first);
			max_r = std::max(max_r, seeds[i].r_last);
			//++(hist->back());
			hist->back() = 1;
			if (i==total_seeds-1 || seeds[i].kmer.h != seeds[i+1].kmer.h) {
				assert(min_r <= max_r);
				seeds[i].r_first = min_r;
				seeds[i].r_last = max_r;
				min_r = p.size(), max_r = -1;
				hist->push_back(0);
				thin_seeds.push_back(seeds[i]);
			}
		}

		H->T.stop("unique_seeds");
        H->C.inc("discarded_seeds", seeds.size() - thin_seeds.size());
		return thin_seeds;
	}

	// Initializes the histogram of the pattern and the list of matches
	vector<Match> match_seeds(const vector<Seed> &seeds) {
		H->T.start("collect_matches");
		vector<Match> matches;
		matches.reserve(2*(int)seeds.size());
		for (int seed_num=0; seed_num<(int)seeds.size(); seed_num++)
			if (seeds[seed_num].hits_in_T > 0)
				tidx.get_matches(&matches, seeds[seed_num]);
		H->T.stop("collect_matches");

		H->T.start("sort_matches");
		//sort
		pdqsort_branchless(matches.begin(), matches.end(), [](const Match &a, const Match &b) {
			// Preparation for sweeping: sort M by ascending positions within one reference segment.
			// NO SPEEDUP: return (int64_t(a.hit.segm_id) << 32 | a.hit.r) < (int64_t(b.hit.segm_id) << 32 | b.hit.r);
			if (a.hit.segm_id != b.hit.segm_id)
				return a.hit.segm_id < b.hit.segm_id;
			return a.hit.r < b.hit.r;
		});
		H->T.stop("sort_matches");

		return matches;
	}

	// TODO: limit diff_hist to 0s and 1s to support Jaccard (no repetitions)
	// vector<hash_t> diff_hist;  // diff_hist[kmer_hash] = #occurences in `p` - #occurences in `s`
	// vector<Match> M;   	   // for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
	const vector<Mapping> sweep(hist_t &diff_hist, const vector<Match> &M, const pos_t P_sz, const int p_sz) {
//		const int MAX_BL = 100;
		vector<Mapping> mappings;	// List of tripples <i, j, score> of matches

		int intersection = 0;
		Mapping best(H->params.k, P_sz, 0, -1, -1, -1, -1, -1, 0, M.end(), M.end());
		Mapping second = best;
		int same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise

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
				if (--diff_hist[r->seed.seed_num] >= 0)
					++intersection;
				assert (l->hit.r <= r->hit.r);
			}

			double J = 1.0*intersection / p_sz;
			auto m = Mapping(H->params.k, P_sz, p_sz, l->hit.r, prev(r)->hit.r, l->hit.segm_id, intersection, J, same_strand_seeds, l, r);

			// second best without guarantees
			// Wrong invariant:
			// best[l,r) -- a mapping best.l<=l with maximal J
			// second_best[l,r) -- a mapping second_best.l \notin [l-90%|P|; l+90%|P|] with maximal J
			if (m.J >= H->params.theta) {
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
			}

			// Prepare for the next step by moving `l` to the right.
			if (++diff_hist[l->seed.seed_num] >= 1)
				--intersection;
			same_strand_seeds -= l->is_same_strand() ? +1 : -1;

			assert(intersection >= 0);
		}
		assert(intersection == 0);
		assert(same_strand_seeds == 0);

		if (H->params.onlybest && best.intersection != -1 && best.J > H->params.theta) {
			best.mapq = (best.intersection > 5 && best.J > 0.1 && best.J > 1.2*second.J) ? 60 : 0;
			best.J2 = second.J;
			mappings.push_back(best);
		}

		return mappings;
	}

	// Return only reasonable matches (i.e. those that are not J-dominated by
	// another overlapping match). Runs in O(|all|).
	vector<Mapping> filter_reasonable(const vector<Mapping> &all, const pos_t P_len) {
		vector<Mapping> reasonable;
		std::deque<Mapping> recent;

		// Minimal separation between mappings to be considered reasonable
		pos_t sep = pos_t((1.0 - H->params.theta) * double(P_len));

		// The deque `recent' is sorted decreasingly by J
		//					  _________`recent'_________
		//                   /                          \
		// ---------------- | High J ... Mid J ... Low J | current J
		// already removed    deque.back ... deque.front    to add next
		for (const auto &next: all) {
			// 1. Prepare for adding `curr' by removing from the deque back all
			//    mappings that are too far to the left. This keeps the deque
			//    within |P| from back to front. A mapping can become reasonable
			//    only after getting removed.
			while(!recent.empty() && next.T_l - recent.back().T_l > sep) {
				// If the mapping is not marked as unreasonable (coverted by a preivous better mapping)
				if (!recent.back().unreasonable) {
					// Take the leftmost mapping.
					reasonable.push_back(recent.back());
					// Mark the next closeby mappings as not reasonable
					for (auto it=recent.rbegin(); it!=recent.rend() && it->T_l - recent.back().T_l < sep; ++it)
						it->unreasonable = true;
				}
				// Remove the mapping that is already too much behind.
				recent.pop_back();
			}
			assert(recent.empty() || (recent.back().T_r <= recent.front().T_r && recent.front().T_r <= next.T_r));

			// Now all the mappings in `recent' are close to `curr' 
			// 2. Remove from the deque front all mappings that are strictly
			//    less similar than the current J. This keeps the deque sorted
			//    descending in J from left to right
			while(!recent.empty() && recent.front().intersection < next.intersection)
				recent.pop_front();
			assert(recent.empty() || (recent.back().intersection >= recent.front().intersection && recent.front().intersection >= next.intersection));

			// 3. Add the next mapping to the front
			recent.push_front(next);	 

			// 4. If there is another mapping in the deque, it is near and better.
			if (recent.size() > 1)
				recent.front().unreasonable = true;
		}

		// 5. Add the last mapping if it is reasonable
		if (!recent.empty() && !recent.back().unreasonable)
			reasonable.push_back(recent.back());

		return reasonable;
	}

    // TODO: disable in release
    int spurious_matches(const Mapping &m, const vector<Match> &matches) {
        int included = 0;
        for (auto &match: matches)
            if (match.hit.segm_id == m.segm_id && match.hit.r >= m.T_l && match.hit.r <= m.T_r)
                included++;
        return matches.size() - included;
    }

  public:
	SweepMapper(const SketchIndex &tidx, Handler *H) : tidx(tidx), H(H) {
		H->C.inc("seeds_limit_reached", 0);
		H->C.inc("unmapped_reads", 0);
		if (H->params.theta < 0.0 || H->params.theta > 1.0) {
			cerr << "tThres = " << H->params.theta << " outside of [0,1]." << endl;
			exit(1);
		}
	}

	void map_all_reads(const string &pFile) {
		cerr << "Mapping reads using SweepMap " << pFile << "..." << endl;

		H->C.inc("kmers", 0);
		H->C.inc("seeds", 0);
		H->C.inc("matches", 0);
		H->C.inc("spurious_matches", 0);
		H->C.inc("J", 0);
		H->C.inc("mappings", 0);
		H->C.inc("total_edit_distance", 0);

		H->T.start("mapping");
		H->T.start("query_reading");
		read_fasta_klib(pFile, [this](kseq_t *seq) {
			H->T.stop("query_reading");
			H->T.start("query_mapping");
			H->T.start("sketching");
			sketch_t p = H->sketcher.sketch(seq->seq.s);
			H->T.stop("sketching");

			char *query_id = seq->name.s;
			pos_t P_sz = (pos_t)seq->seq.l;
			H->C.inc("kmers", p.size());
			H->C.inc("read_len", P_sz);

			hist_t p_hist;
			Timer read_mapping_time;
			read_mapping_time.start();
			H->T.start("seeding");
			vector<Seed> thin_seeds = select_seeds(p, &p_hist);
			H->T.stop("seeding");
			H->C.inc("seeds", thin_seeds.size());

			H->T.start("matching");
			vector<Match> matches = match_seeds(thin_seeds);
			H->T.stop("matching");

			H->T.start("sweep");
			vector<Mapping> mappings = sweep(p_hist, matches, P_sz, thin_seeds.size());
			H->T.stop("sweep");

			H->T.start("postproc");
			if (!H->params.overlaps)
				mappings = filter_reasonable(mappings, P_sz);
			read_mapping_time.stop();

			for (auto &m: mappings) {
				const auto &segm = tidx.T[m.segm_id];
				m.map_time = read_mapping_time.secs() / (double)mappings.size();
				m.query_id = query_id;
				m.segm_sz = segm.sz;
				m.segm_name = segm.name.c_str();
				m.total_matches = matches.size();
				if (H->params.sam) {
					auto ed = m.print_sam(query_id, segm, (int)matches.size(), seq->seq.s, seq->seq.l);
					H->C.inc("total_edit_distance", ed);
				}
				else m.print_paf();
                H->C.inc("spurious_matches", spurious_matches(m, matches));
				H->C.inc("J", int(10000.0*m.J));
				H->C.inc("mappings");
			}
			H->C.inc("matches", matches.size());
			H->C.inc("reads");
			if (mappings.empty())
				H->C.inc("unmapped_reads");
			H->T.stop("postproc");

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
		cerr << "Mapping stats:" << endl;
		cerr << " | Total reads:           " << H->C.count("reads") << " (~" << 1.0*H->C.count("read_len") / H->C.count("reads") << " nb per read)" << endl;
		cerr << " |  | Unmapped reads:        " << H->C.count("unmapped_reads") << " (" << H->C.perc("unmapped_reads", "reads") << "%)" << endl;
		cerr << " | Read kmers (total):    " << H->C.count("kmers") << " (" << H->C.frac("kmers", "reads") << " per read)" << endl;
		cerr << " |  | unique:                 " << H->C.count("seeds") << " (" << H->C.frac("seeds", "kmers") << ")" << endl;
		cerr << " | Matches:               " << H->C.count("matches") << " (" << H->C.frac("matches", "reads") << " per read)" << endl;
		cerr << " | Seed limit reached:    " << H->C.count("seeds_limit_reached") << " (" << H->C.perc("seeds_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Matches limit reached: " << H->C.count("matches_limit_reached") << " (" << H->C.perc("matches_limit_reached", "reads") << "%)" << endl;
		cerr << " | Spurious matches:      " << H->C.count("spurious_matches") << " (" << H->C.perc("spurious_matches", "matches") << "%)" << endl;
		cerr << " | Discarded seeds:       " << H->C.count("discarded_seeds") << " (" << H->C.perc("discarded_seeds", "collected_seeds") << "%)" << endl;
		cerr << " | Average Jaccard:       " << H->C.frac("J", "mappings") / 10000.0 << endl;
		cerr << " | Average edit dist:     " << H->C.frac("total_edit_distance", "mappings") << endl;
        //printMemoryUsage();
	}

    void print_time_stats() {
        cerr << std::fixed << std::setprecision(1);
        cerr << "Time [sec]:           "             << setw(5) << right << H->T.secs("total")             << endl;
        cerr << " | Index:                 "         << setw(5) << right << H->T.secs("indexing")          << " (" << setw(4) << right << H->T.perc("indexing", "total")              << "\%)" << endl;
        cerr << " |  | loading:                "     << setw(5) << right << H->T.secs("index_reading")     << " (" << setw(4) << right << H->T.perc("index_reading", "indexing")      << "\%)" << endl;
        cerr << " |  | sketch:                 "     << setw(5) << right << H->T.secs("index_sketching")   << " (" << setw(4) << right << H->T.perc("index_sketching", "indexing")    << "\%)" << endl;
        cerr << " |  | initialize:             "     << setw(5) << right << H->T.secs("index_initializing")<< " (" << setw(4) << right << H->T.perc("index_initializing", "indexing") << "\%)" << endl;
        cerr << " | Map:                   "         << setw(5) << right << H->T.secs("mapping")           << " (" << setw(4) << right << H->T.perc("mapping", "total")               << "\%, " << setw(5) << right << H->T.range_ratio("query_mapping") << "x, " << setw(4) << right << H->C.count("reads") / H->T.secs("total") << " reads per sec)" << endl;
        cerr << " |  | load queries:           "     << setw(5) << right << H->T.secs("query_reading")     << " (" << setw(4) << right << H->T.perc("query_reading", "mapping")       << "\%, " << setw(5) << right << H->T.range_ratio("query_reading") << "x)" << endl;
        cerr << " |  | sketch reads:           "     << setw(5) << right << H->T.secs("sketching")         << " (" << setw(4) << right << H->T.perc("sketching", "mapping")           << "\%, " << setw(5) << right << H->T.range_ratio("sketching") << "x)" << endl;
        cerr << " |  | seeding:                "     << setw(5) << right << H->T.secs("seeding")           << " (" << setw(4) << right << H->T.perc("seeding", "mapping")             << "\%, " << setw(5) << right << H->T.range_ratio("seeding") << "x)" << endl;
        cerr << " |  |  | collect seed info:       " << setw(5) << right << H->T.secs("collect_seed_info") << " (" << setw(4) << right << H->T.perc("collect_seed_info", "seeding")   << "\%, " << setw(5) << right << H->T.range_ratio("collect_seed_info") << "x)" << endl;
        cerr << " |  |  | thin sketch:             " << setw(5) << right << H->T.secs("thin_sketch")       << " (" << setw(4) << right << H->T.perc("thin_sketch", "seeding")         << "\%, " << setw(5) << right << H->T.range_ratio("thin_sketch") << "x)" << endl;
        cerr << " |  |  | sort seeds:              " << setw(5) << right << H->T.secs("sort_seeds")        << " (" << setw(4) << right << H->T.perc("sort_seeds", "seeding")          << "\%, " << setw(5) << right << H->T.range_ratio("sort_seeds") << "x)" << endl;
        cerr << " |  |  | unique seeds:            " << setw(5) << right << H->T.secs("unique_seeds")      << " (" << setw(4) << right << H->T.perc("unique_seeds", "seeding")        << "\%, " << setw(5) << right << H->T.range_ratio("unique_seeds") << "x)" << endl;
        cerr << " |  | matching seeds:         "     << setw(5) << right << H->T.secs("matching")          << " (" << setw(4) << right << H->T.perc("matching", "mapping")            << "\%, " << setw(5) << right << H->T.range_ratio("matching") << "x)" << endl;
        cerr << " |  |  | collect matches:         " << setw(5) << right << H->T.secs("collect_matches")   << " (" << setw(4) << right << H->T.perc("collect_matches", "matching")    << "\%, " << setw(5) << right << H->T.range_ratio("collect_matches") << "x)" << endl;
        cerr << " |  |  | sort matches:            " << setw(5) << right << H->T.secs("sort_matches")      << " (" << setw(4) << right << H->T.perc("sort_matches", "matching")       << "\%, " << setw(5) << right << H->T.range_ratio("sort_matches") << "x)" << endl;
        cerr << " |  | sweep:                  "     << setw(5) << right << H->T.secs("sweep")             << " (" << setw(4) << right << H->T.perc("sweep", "mapping")               << "\%, " << setw(5) << right << H->T.range_ratio("sweep") << "x)" << endl;
        cerr << " |  | post proc:              "     << setw(5) << right << H->T.secs("postproc")          << " (" << setw(4) << right << H->T.perc("postproc", "mapping")            << "\%, " << setw(5) << right << H->T.range_ratio("postproc") << "x)" << endl;
    //		cerr << "Virtual memory [MB]:  "             << setw(5) << right << C.count("total_memory_MB")  << endl;
    //		cerr << " | Index:                 "         << setw(5) << right << C.count("index_memory_MB") << " (" << setw(4) << right << C.perc("index_memory_MB", "total_memory_MB") << "\%)" << endl;
    }

};

}  // namespace sweepmap