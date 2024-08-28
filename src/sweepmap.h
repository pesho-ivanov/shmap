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

namespace sweepmap {

using std::setw;
using std::right;
using std::vector;

class SweepMap {
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
				seeds.push_back(Seed(kmer, p[ppos].r, p[ppos].r, count));
		}
		H->T.stop("collect_seed_info");
        H->C.inc("collected_seeds", seeds.size());

		H->T.start("thin_sketch");
		// TODO: add all seeds to hist
		int total_seeds = (int)seeds.size();
		if ((int)H->params.max_seeds < total_seeds) {
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
		for (int i=0; i<total_seeds-1; i++) {
			min_r = std::min(min_r, seeds[i].r_first);
			max_r = std::max(max_r, seeds[i].r_last);
			++(hist->back());
			if (seeds[i].kmer.h != seeds[i+1].kmer.h) {
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
	vector<Match> match_seeds(pos_t p_sz, const vector<Seed> &seeds) {
		H->T.start("collect_matches");
		vector<Match> matches;
		matches.reserve(2*(int)seeds.size());
		for (int seed_num=0; seed_num<(int)seeds.size(); seed_num++)
			tidx.add_matches(&matches, seeds[seed_num], seed_num);
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

	// vector<hash_t> diff_hist;  // diff_hist[kmer_hash] = #occurences in `p` - #occurences in `s`
	// vector<Match> M;   	   // for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
	const vector<Mapping> sweep(hist_t &diff_hist, const sketch_t &p, const vector<Match> &M, const pos_t P_len, const int thin_seeds_cnt) {
//		const int MAX_BL = 100;
		vector<Mapping> mappings;	// List of tripples <i, j, score> of matches

		int xmin = 0;
		Mapping best(H->params.k, P_len, 0, -1, -1, -1, -1, -1, 0, M.end(), M.end());
		Mapping second = best;
		int same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise

		// Increase the left point end of the window [l,r) one by one. O(matches)
		for(auto l = M.begin(), r = M.begin(); l != M.end(); ++l) {
			// Increase the right end of the window [l,r) until it gets out.
			for(;  r != M.end()
				&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				&& r->hit.r + H->params.k <= l->hit.r + P_len
				; ++r) {
				same_strand_seeds += r->is_same_strand() ? +1 : -1;  // change to r inside the loop
				// If taking this kmer from T increases the intersection with P.
				// TODO: iterate following seeds
				if (--diff_hist[r->seed_num] >= 0)
					++xmin;
				assert (l->hit.r <= r->hit.r);
			}

			auto m = Mapping(H->params.k, P_len, thin_seeds_cnt, l->hit.r, prev(r)->hit.r, l->hit.segm_id, pos_t(r-l), xmin, same_strand_seeds, l, r);

			// second best without guarantees
			// Wrong invariant:
			// best[l,r) -- a mapping best.l<=l with maximal J
			// second_best[l,r) -- a mapping second_best.l \notin [l-90%|P|; l+90%|P|] with maximal J
			if (H->params.onlybest) {
				if (m.xmin > best.xmin) {  // if (xmin > best.xmin)
					if (best.T_l < m.T_l - 0.9*P_len)
						second = best;
					best = m;
				} else if (m.xmin > second.xmin && m.T_l > best.T_l + 0.9*P_len) {
					second = m;
				}
			} else {
				if (m.J > H->params.tThres) {
					mappings.push_back(m);
				}
			}

			// Prepare for the next step by moving `l` to the right.
			if (++diff_hist[l->seed_num] > 0)
				--xmin;
			same_strand_seeds -= l->is_same_strand() ? +1 : -1;

			assert(xmin >= 0);
		}
		assert(xmin == 0);
		assert(same_strand_seeds == 0);

		if (H->params.onlybest && best.xmin != -1) { // && best.J > H->params.tThres)
			best.mapq = (best.xmin > 5 && best.J > 0.1 && best.J > 1.2*second.J) ? 60 : 0;
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
		pos_t sep = pos_t((1.0 - H->params.tThres) * double(P_len));

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
			while(!recent.empty() && recent.front().xmin < next.xmin)
				recent.pop_front();
			assert(recent.empty() || (recent.back().xmin >= recent.front().xmin && recent.front().xmin >= next.xmin));

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
	SweepMap(const SketchIndex &tidx, Handler *H) : tidx(tidx), H(H) {
			H->C.inc("seeds_limit_reached", 0);
			H->C.inc("unmapped_reads", 0);
			if (H->params.tThres < 0.0 || H->params.tThres > 1.0) {
				cerr << "tThres = " << H->params.tThres << " outside of [0,1]." << endl;
				exit(1);
			}
		}

	void map(const string &pFile) {
		cerr << "Mapping reads " << pFile << "..." << endl;

		H->C.inc("spurious_matches", 0);
		H->C.inc("J", 0);
		H->C.inc("mappings", 0);
		H->C.inc("sketched_kmers", 0);
		H->C.inc("total_edit_distance", 0);

		H->T.start("mapping");
		H->T.start("query_reading");
		read_fasta_klib(pFile, [this](kseq_t *seq) {
			H->T.stop("query_reading");
			H->T.start("query_mapping");
			H->T.start("sketching");
			sketch_t p = H->sketcher.sketch(seq->seq.s);
			H->T.stop("sketching");

			string query_id = seq->name.s;
			pos_t P_sz = (pos_t)seq->seq.l;

			H->C.inc("read_len", P_sz);
			hist_t p_hist;

			Timer read_mapping_time;
			read_mapping_time.start();
			H->T.start("seeding");
			vector<Seed> thin_seeds = select_seeds(p, &p_hist);
			H->T.stop("seeding");

			H->T.start("matching");
			vector<Match> matches = match_seeds(p.size(), thin_seeds);
			H->T.stop("matching");

			H->T.start("sweep");
			vector<Mapping> mappings = sweep(p_hist, p, matches, P_sz, thin_seeds.size());
			H->T.stop("sweep");

			H->T.start("postproc");
			if (!H->params.overlaps)
				mappings = filter_reasonable(mappings, P_sz);
			read_mapping_time.stop();

			for (auto &m: mappings) {
				const auto &segm = tidx.T[m.segm_id];
				m.map_time = read_mapping_time.secs() / (double)mappings.size();
				if (H->params.sam) {
					auto ed = m.print_sam(query_id, segm, (int)matches.size(), seq->seq.s, seq->seq.l);
					H->C.inc("total_edit_distance", ed);
				}
				else m.print_paf(query_id, segm, matches);
                H->C.inc("spurious_matches", spurious_matches(m, matches));
				H->C.inc("J", int(10000.0*m.J));
				H->C.inc("mappings");
				H->C.inc("sketched_kmers", m.seeds);
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
	}

	void print_stats() {
		cerr << std::fixed << std::setprecision(1);
		cerr << "Mapping stats:" << endl;
		cerr << " | Total reads:           " << H->C.count("reads") << " (~" << 1.0*H->C.count("read_len") / H->C.count("reads") << " nb per read)" << endl;
		cerr << " | Sketched read kmers:   " << H->C.count("sketched_kmers") << " (" << H->C.frac("sketched_kmers", "reads") << " per read)" << endl;
		cerr << " | Kmer matches:          " << H->C.count("matches") << " (" << H->C.frac("matches", "reads") << " per read)" << endl;
		cerr << " | Seed limit reached:    " << H->C.count("seeds_limit_reached") << " (" << H->C.perc("seeds_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Matches limit reached: " << H->C.count("matches_limit_reached") << " (" << H->C.perc("matches_limit_reached", "reads") << "%)" << endl;
		cerr << " | Spurious matches:      " << H->C.count("spurious_matches") << " (" << H->C.perc("spurious_matches", "matches") << "%)" << endl;
		cerr << " | Discarded seeds:       " << H->C.count("discarded_seeds") << " (" << H->C.perc("discarded_seeds", "collected_seeds") << "%)" << endl;
		cerr << " | Unmapped reads:        " << H->C.count("unmapped_reads") << " (" << H->C.perc("unmapped_reads", "reads") << "%)" << endl;
		cerr << " | Average Jaccard:       " << H->C.frac("J", "mappings") / 10000.0 << endl;
		cerr << " | Average edit dist:     " << H->C.frac("total_edit_distance", "mappings") << endl;
	}
};

}  // namespace sweepmap