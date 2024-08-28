#pragma once

#include <algorithm>
#include <cassert>
#include <deque>
#include <iomanip>
#include <set>

#include "../ext/edlib.h"
#include "../ext/pdqsort.h"

#include "index.h"
#include "io.h"
#include "sketch.h"

namespace sweepmap {

using std::setw;
using std::right;
using std::vector;

struct Mapping {
	int k; 	   // kmer size
	pos_t P_sz;     // pattern size |P| bp 
	pos_t seeds;     // number of seeds (subset of the sketch kmers)
	pos_t T_l;      // the position of the leftmost nucleotide of the mapping
	pos_t T_r;      // the position of the rightmost nucleotide of the mapping
	segm_t segm_id;
	pos_t s_sz;      // the position of the rightmost nucleotide of the mapping
	int xmin;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	double J, J2;     // Jaccard similarity in [0;1] for the best and for the second best mapping
	double map_time;
	int mapq;
	char strand;    // '+' or '-'
	bool unreasonable;  // reserved for filtering matches
	vector<Match>::const_iterator l, r;

    Mapping() {}
	Mapping(int k, pos_t P_sz, int seeds, pos_t T_l, pos_t T_r, segm_t segm_id, pos_t s_sz, int xmin, int same_strand_seeds, vector<Match>::const_iterator l, vector<Match>::const_iterator r)
		: k(k), P_sz(P_sz), seeds(seeds), T_l(T_l), T_r(T_r), segm_id(segm_id), s_sz(s_sz), xmin(xmin), J(double(xmin) / std::max(seeds, s_sz)), mapq(255), strand(same_strand_seeds > 0 ? '+' : '-'), unreasonable(false), l(l), r(r) {}

	// --- https://github.com/lh3/miniasm/blob/master/PAF.md ---
    void print_paf(const string &query_id, const RefSegment &segm, vector<Match> matches) const {
		int P_start = P_sz, P_end = -1;
		for (auto m = l; m != r; ++m) {
			P_start = std::min(P_start, m->seed.r_first);
			P_end = std::max(P_end, m->seed.r_last);
		}
		if (!(0 <= P_start && P_start <= P_end && P_end <= P_sz))
			std::cerr << "P_start=" << P_start << " P_end=" << P_end << " P_sz=" << P_sz << std::endl;
		assert(0 <= P_start && P_start <= P_end && P_end <= P_sz);

		auto T_l_predicted = std::max(T_l-P_start, 0);  // -P_start -- P_start too big
		auto T_r_predicted = std::min(T_r+(P_sz-P_end), segm.sz);  // +(P_sz-P_end) too big

		std::cout << query_id  			// Query sequence name
			<< "\t" << P_sz     // query sequence length
			<< "\t" << P_start   // query start (0-based; closed)
			<< "\t" << P_end  // query end (0-based; open)
			<< "\t" << strand   // '+' or '-'
			<< "\t" << segm.name // reference segment name
			<< "\t" << segm.sz // T_sz -- target sequence length
			<< "\t" << T_l_predicted  // target start on original strand (0-based)
			<< "\t" << T_r_predicted  // target start on original strand (0-based)
			<< "\t" << P_sz  // TODO: fix; Number of residue matches (number of nucleotide matches)
			<< "\t" << P_sz  // TODO: fix; Alignment block length: total number of sequence matches, mismatches, insertions and deletions in the alignment
			<< "\t" << mapq  // Mapping quality (0-255; 255 for missing)
// ----- end of required PAF fields -----
			<< "\t" << "k:i:" << k
			<< "\t" << "p:i:" << seeds  // sketches
			<< "\t" << "M:i:" << matches.size() // kmer matches in T
			<< "\t" << "s:i:" << s_sz
			<< "\t" << "I:i:" << xmin  // intersection of `p` and `s` [kmers]
			<< "\t" << "J:f:" << J   // Jaccard similarity [0; 1]
			<< "\t" << "J2:f:" << J2   // second best mapping Jaccard similarity [0; 1]
			<< "\t" << "t:f:" << map_time
			<< endl;
	}

    int print_sam(const string &query_id, const RefSegment &segm, const int matches, const char *query, const size_t query_size) const {
		int T_start = std::max(T_l-k, 0);
		int T_end = std::max(T_r, T_l-k+P_sz);
		int T_d = T_end - T_start;
		string s = segm.seq.substr(T_start, T_d);
		if (strand == '-')
		   s = reverseComplement(s);
		assert(T_d >= 0);
		auto max_edit_dist = -1; //10000;
		auto cfg = edlibNewAlignConfig(max_edit_dist, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
		EdlibAlignResult result = edlibAlign(query, query_size, s.c_str(), T_d, cfg);
		assert(result.status == EDLIB_STATUS_OK);
		char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
		//printf("query=%s, s=%s, ", query, s.c_str());
		//printf("ed=%d, cigar=%s\n", result.editDistance, cigar);

		int ed = result.editDistance;
		int flag = 0;
		if (strand == '-') flag |= 0x10;
		std::cout << query_id 				// 1 QNAME String [!-?A-~]{1,254} Query template NAME
			<< "\t" << flag 				// 2 FLAG Int [0, 216 − 1] bitwise FLAG
			<< "\t" << segm.name  			// 3 RNAME String \*|[:rname:∧*=][:rname:]* Reference sequence NAME11
			<< "\t" << T_start+1  			// 4 POS Int [0, 231 − 1] 1-based leftmost mapping POSition
			<< "\t" << mapq  				// 5 MAPQ Int [0, 28 − 1] MAPping Quality
			<< "\t" << cigar  				// 6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
			<< "\t" << "="  				// 7 RNEXT String \*|=|[:rname:∧*=][:rname:]* Reference name of the mate/next read
			<< "\t" << 0  					// 8 PNEXT Int [0, 231 − 1] Position of the mate/next read
			<< "\t" << 0  					// 9 TLEN Int [−231 + 1, 231 − 1] observed Template LENgth
			<< "\t" << query  				// 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
			<< "\t" << "!"  				// 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
			<< "\t" << "NM:i:" << ed  		// edit distance
//			<< "\t" << "q:s:" << query		// query string
//			<< "\t" << "s:s:" << s			// window string in T
			<< endl;

		free(cigar);
		edlibFreeAlignResult(result);
		return ed;
	}
};

class SweepMap {
	const SketchIndex &tidx;
	const params_t &params;
	Counters *C;
	Timers *T;
	const FracMinHash &sketcher;

	using hist_t = vector<int>;

	vector<Seed> select_seeds(const sketch_t& p, hist_t *hist) {
		T->start("collect_seed_info");
		vector<Seed> seeds;
		seeds.reserve(p.size());

		// TODO: limit The number of kmers in the pattern p
		for (int ppos = 0; ppos < (int)p.size(); ++ppos) {
			const auto &kmer = p[ppos];
			const auto count = tidx.count(kmer.h);
			if (count > 0)
				seeds.push_back(Seed(kmer, p[ppos].r, p[ppos].r, count));
		}
		T->stop("collect_seed_info");
        C->inc("collected_seeds", seeds.size());

		T->start("thin_sketch");
		// TODO: add all seeds to hist
		int total_seeds = (int)seeds.size();
		if ((int)params.max_seeds < total_seeds) {
			total_seeds = (int)params.max_seeds;
			C->inc("seeds_limit_reached");
		}
        std::nth_element(seeds.begin(), seeds.begin() + total_seeds, seeds.end(), [](const Seed &a, const Seed &b) {
            return a.hits_in_T < b.hits_in_T;
        });
		T->stop("thin_sketch");

		T->start("sort_seeds");
		pdqsort_branchless(seeds.begin(), seeds.begin() + total_seeds, [](const Seed &a, const Seed &b) {
			return a.kmer.h < b.kmer.h;
		});
		T->stop("sort_seeds");

		T->start("unique_seeds");
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

		T->stop("unique_seeds");
        C->inc("discarded_seeds", seeds.size() - thin_seeds.size());
		return thin_seeds;
	}

	// Initializes the histogram of the pattern and the list of matches
	vector<Match> match_seeds(pos_t p_sz, const vector<Seed> &seeds) {
		T->start("collect_matches");
		vector<Match> matches;
		matches.reserve(2*(int)seeds.size());
		for (int seed_num=0; seed_num<(int)seeds.size(); seed_num++)
			tidx.add_matches(&matches, seeds[seed_num], seed_num);
		T->stop("collect_matches");

		T->start("sort_matches");
		//sort
		pdqsort_branchless(matches.begin(), matches.end(), [](const Match &a, const Match &b) {
			// Preparation for sweeping: sort M by ascending positions within one reference segment.
			// NO SPEEDUP: return (int64_t(a.hit.segm_id) << 32 | a.hit.r) < (int64_t(b.hit.segm_id) << 32 | b.hit.r);
			if (a.hit.segm_id != b.hit.segm_id)
				return a.hit.segm_id < b.hit.segm_id;
			return a.hit.r < b.hit.r;
		});
		T->stop("sort_matches");

		return matches;
	}

	// vector<hash_t> diff_hist;  // diff_hist[kmer_hash] = #occurences in `p` - #occurences in `s`
	// vector<Match> M;   	   // for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
	const vector<Mapping> sweep(hist_t &diff_hist, const sketch_t &p, const vector<Match> &M, const pos_t P_len, const int thin_seeds_cnt) {
//		const int MAX_BL = 100;
		vector<Mapping> mappings;	// List of tripples <i, j, score> of matches

		int xmin = 0;
		Mapping best(params.k, P_len, 0, -1, -1, -1, -1, -1, 0, M.end(), M.end());
		Mapping second = best;
		int same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise

		// Increase the left point end of the window [l,r) one by one. O(matches)
		for(auto l = M.begin(), r = M.begin(); l != M.end(); ++l) {
			// Increase the right end of the window [l,r) until it gets out.
			for(;  r != M.end()
				&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				&& r->hit.r + params.k <= l->hit.r + P_len
				; ++r) {
				same_strand_seeds += r->is_same_strand() ? +1 : -1;  // change to r inside the loop
				// If taking this kmer from T increases the intersection with P.
				// TODO: iterate following seeds
				if (--diff_hist[r->seed_num] >= 0)
					++xmin;
				assert (l->hit.r <= r->hit.r);
			}

			auto m = Mapping(params.k, P_len, thin_seeds_cnt, l->hit.r, prev(r)->hit.r, l->hit.segm_id, pos_t(r-l), xmin, same_strand_seeds, l, r);

			// second best without guarantees
			// Wrong invariant:
			// best[l,r) -- a mapping best.l<=l with maximal J
			// second_best[l,r) -- a mapping second_best.l \notin [l-90%|P|; l+90%|P|] with maximal J
			if (params.onlybest) {
				if (m.xmin > best.xmin) {  // if (xmin > best.xmin)
					if (best.T_l < m.T_l - 0.9*P_len)
						second = best;
					best = m;
				} else if (m.xmin > second.xmin && m.T_l > best.T_l + 0.9*P_len) {
					second = m;
				}
			} else {
				if (m.J > params.tThres) {
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

		if (params.onlybest && best.xmin != -1) { // && best.J > params.tThres)
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
		pos_t sep = pos_t((1.0 - params.tThres) * double(P_len));

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
	SweepMap(const SketchIndex &tidx, const params_t &params, Counters *C, Timers *T, const FracMinHash &sketcher)
		: tidx(tidx), params(params), C(C), T(T), sketcher(sketcher) {
			C->inc("seeds_limit_reached", 0);
			C->inc("unmapped_reads", 0);
			if (params.tThres < 0.0 || params.tThres > 1.0) {
				cerr << "tThres = " << params.tThres << " outside of [0,1]." << endl;
				exit(1);
			}
		}

	void map(const string &pFile) {
		C->inc("spurious_matches", 0);
		C->inc("J", 0);
		C->inc("mappings", 0);
		C->inc("sketched_kmers", 0);
		C->inc("total_edit_distance", 0);

		T->start("mapping");
		T->start("query_reading");
		read_fasta_klib(pFile, [this](kseq_t *seq) {
			T->stop("query_reading");
			T->start("query_mapping");
			T->start("sketching");
			sketch_t p = sketcher.sketch(seq->seq.s);
			T->stop("sketching");

			string query_id = seq->name.s;
			pos_t P_sz = (pos_t)seq->seq.l;

			C->inc("read_len", P_sz);
			hist_t p_hist;

			Timer read_mapping_time;
			read_mapping_time.start();
			T->start("seeding");
			vector<Seed> thin_seeds = select_seeds(p, &p_hist);
			T->stop("seeding");

			T->start("matching");
			vector<Match> matches = match_seeds(p.size(), thin_seeds);
			T->stop("matching");

			T->start("sweep");
			vector<Mapping> mappings = sweep(p_hist, p, matches, P_sz, thin_seeds.size());
			T->stop("sweep");

			T->start("postproc");
			if (!params.overlaps)
				mappings = filter_reasonable(mappings, P_sz);
			read_mapping_time.stop();

			for (auto &m: mappings) {
				const auto &segm = tidx.T[m.segm_id];
				m.map_time = read_mapping_time.secs() / (double)mappings.size();
				if (params.sam) {
					auto ed = m.print_sam(query_id, segm, (int)matches.size(), seq->seq.s, seq->seq.l);
					C->inc("total_edit_distance", ed);
				}
				else m.print_paf(query_id, segm, matches);
                C->inc("spurious_matches", spurious_matches(m, matches));
				C->inc("J", int(10000.0*m.J));
				C->inc("mappings");
				C->inc("sketched_kmers", m.seeds);
			}
			C->inc("matches", matches.size());
			C->inc("reads");
			if (mappings.empty())
				C->inc("unmapped_reads");
			T->stop("postproc");

			T->stop("query_mapping");
			T->start("query_reading");
		});
		T->stop("query_reading");
		T->stop("mapping");

		print_stats();
	}

	void print_stats() {
		cerr << std::fixed << std::setprecision(1);
		cerr << "Mapping stats:" << endl;
		cerr << " | Total reads:           " << C->count("reads") << " (~" << 1.0*C->count("read_len") / C->count("reads") << " nb per read)" << endl;
		cerr << " | Sketched read kmers:   " << C->count("sketched_kmers") << " (" << C->frac("sketched_kmers", "reads") << " per read)" << endl;
		cerr << " | Kmer matches:          " << C->count("matches") << " (" << C->frac("matches", "reads") << " per read)" << endl;
		cerr << " | Seed limit reached:    " << C->count("seeds_limit_reached") << " (" << C->perc("seeds_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Matches limit reached: " << C->count("matches_limit_reached") << " (" << C->perc("matches_limit_reached", "reads") << "%)" << endl;
		cerr << " | Spurious matches:      " << C->count("spurious_matches") << " (" << C->perc("spurious_matches", "matches") << "%)" << endl;
		cerr << " | Discarded seeds:       " << C->count("discarded_seeds") << " (" << C->perc("discarded_seeds", "collected_seeds") << "%)" << endl;
		cerr << " | Unmapped reads:        " << C->count("unmapped_reads") << " (" << C->perc("unmapped_reads", "reads") << "%)" << endl;
		cerr << " | Average Jaccard:       " << C->frac("J", "mappings") / 10000.0 << endl;
		cerr << " | Average edit dist:     " << C->frac("total_edit_distance", "mappings") << endl;
	}
};

}  // namespace sweepmap