#ifndef SWEEP_HPP
#define SWEEP_HPP

#include <algorithm>
#include <cassert>
#include <deque>
#include <iomanip>
#include <set>

#include "io.h"
#include "sketch.h"

struct Mapping {
	int k; 	   // kmer size
	pos_t P_sz;     // pattern size |P| bp 
	pos_t p_sz;     // pattern sketch size |p| kmers
	int matches;  // L.size() -- total number of matches in `t' 
	pos_t T_l;      // the position of the leftmost nucleotide of the mapping
	pos_t T_r;      // the position of the rightmost nucleotide of the mapping
	pos_t s_sz;      // the position of the rightmost nucleotide of the mapping
	int xmin;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	double sim;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	double map_time;

    Mapping() {}
	Mapping(int k, pos_t P_sz, pos_t p_sz, int matches, pos_t T_l, pos_t T_r, pos_t s_sz, int xmin)
		: k(k), P_sz(P_sz), p_sz(p_sz), matches(matches), T_l(T_l), T_r(T_r), s_sz(s_sz), xmin(xmin), sim(double(xmin) / std::max(p_sz, s_sz)) {}

    void print_paf(const params_t &params, const string &seqID, const pos_t T_sz, const string &T_name) const {
		// --- https://github.com/lh3/miniasm/blob/master/PAF.md ---
		cout << seqID  			// Query sequence name
			<< "\t" << P_sz     // query sequence length
			<< "\t" << 0   // query start (0-based; closed)
			<< "\t" << P_sz  // query end (0-based; open)
			<< "\t" << "+"   // m.get_strand() //strand; TODO
			<< "\t" << T_name    // reference name
			<< "\t" << T_sz  // target sequence length
			<< "\t" << T_l  // target start on original strand (0-based)
			<< "\t" << T_r  // target start on original strand (0-based)
			<< "\t" << P_sz  // TODO: fix; Number of residue matches (number of nucleotide matches)
			<< "\t" << P_sz  // TODO: fix; Alignment block length: total number of sequence matches, mismatches, insertions and deletions in the alignment
			<< "\t" << 60  // Mapping quality (0-255; 255 for missing)
		// ----- end of required PAF fields -----
			<< "\t" << "k:i:" << k
			<< "\t" << "M:i:" << matches // matches of `p` in `s` [kmers]
			<< "\t" << "p:i:" << p_sz 
			<< "\t" << "s:i:" << s_sz
			<< "\t" << "I:i:" << xmin  // intersection of `p` and `s` [kmers]
			<< "\t" << "S:f:" << sim   // similarity [0; 1]
			<< "\t" << "t:f:" << map_time
			<< endl;
	}
};

class SweepMap {
	const SketchIndex &tidx;
	const params_t &params;
	Timers *T;
	Counters *C;

	using seeds_t = vector<kmer_hits_t>;

	// return if the kmer has not been seen before and construct hash2num[kmer_hash] = kmer_num, which is not more than |p| 
	inline bool add2hist(const hash_t kmer_hash, vector<int> *hist, ankerl::unordered_dense::map<hash_t, kmer_num_t> *hash2num, kmer_num_t *kmer_num) {
		auto num_it = hash2num->find(kmer_hash);
		if (num_it != hash2num->end()) {
			*kmer_num = num_it->second;
			++(*hist)[*kmer_num];
			return false;
		} else {
			*kmer_num = (kmer_num_t)hist->size();
			hash2num->insert({kmer_hash, *kmer_num});
			hist->push_back(1);
			return true;
		}
	}

	// Returns true if the window [l,r) should be extended to the right
	inline bool should_extend_right(
			const vector<int> &hist,
			const vector<Match> &L,
			const vector<Match>::const_iterator l, const vector<Match>::const_iterator r,
			const int xmin, const int Plen_nucl, const int k) {
		if (r == L.end())
			return false;
		bool extension_stays_within_window = r->T_r + k <= l->T_r + Plen_nucl;
		return extension_stays_within_window;
	}

	//void print_matches(const vector<Match> &L) const {
	//	cout << "Matches:" << endl;
	//	for (auto &m: L) {
	//		cout << "kmer_ord=" << m.kmer_ord << ", P_l=" << m.P_l << ", T_r=" << m.T_r << ", t_pos=" << m.t_pos << endl;
	//	}
	//}

	seeds_t choose_seeds(const Sketch& p, vector<int> *p_hist) {
		ankerl::unordered_dense::map<hash_t, kmer_num_t> hash2num;
		seeds_t seeds;
		seeds.reserve(p.size());

		T->start("collect_seed_info");
		for (auto &[P_l, kmer_hash] : p) {
			// TODO: limit the number of kmers in the pattern p
			kmer_num_t kmer_num;

			// Add pattern kmers from the pattern to p_hist 
			if (add2hist(kmer_hash, p_hist, &hash2num, &kmer_num)) {
				auto it = tidx.h2pos.find(kmer_hash);
				if (it != tidx.h2pos.end())
					seeds.emplace_back(P_l, kmer_num, &it->second);
			}
		}
		T->stop("collect_seed_info");

		// TODO: instead of sorting, place the seeds in a heap which retains only top-K
		// TODO: instead of sorting, find top-K value, and filter seeds below of it
		// TODO: make p_hist smaller since we don't use all seeds

		T->start("sort_seeds");
		// Sort the matching kmers by number of hits in the reference
		sort(seeds.begin(), seeds.end(), [](const kmer_hits_t &a, const kmer_hits_t &b) {
			return a.kmers_in_T->size() < b.kmers_in_T->size();
		});
		T->stop("sort_seeds");

		return seeds;
	}

	// Initializes the histogram of the pattern and the list of matches
	vector<Match> match_seeds(const Sketch& p, const seeds_t &seeds, vector<int> *p_hist) {
		T->start("collect_matches");
		// Get MAX_SEEDS of kmers with the lowest number of hits.
		vector<Match> L;
		L.reserve(P_MULTIPLICITY * p.size());

		int total_hits = 0;
		for (int seed=0; seed<(int)seeds.size(); seed++) {
			if (seed > (int)params.max_seeds) {
				C->inc("seeds_limit_reached");
				break;
			}
			const kmer_hits_t &res = seeds[seed];
			for (const auto &hit: *(res.kmers_in_T))
				L.emplace_back(res.kmer_num, res.P_l, hit.first, hit.second);
			if ((total_hits += (int)res.kmers_in_T->size()) > (int)params.max_matches) {
				C->inc("matches_limit_reached");
				break;
			}
		}
		T->stop("collect_matches");

		T->start("sort_matches");
		// Sort L by ascending positions in reference so that we can next sweep through it.
		sort(L.begin(), L.end(), [](const Match &a, const Match &b) {
			return a.T_r < b.T_r;
		});
		T->stop("sort_matches");

		return L;
	}

	// vector<int> diff_hist;  // diff_hist[kmer_hash] = #occurences in `p` - #occurences in `s`
	// vector<Match> L;   	   // for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
	const vector<Mapping> sweep(vector<int> &diff_hist, const vector<Match> &L, const pos_t p_len, const pos_t P_len) {
		vector<Mapping> mappings;	// List of tripples <i, j, score> of matches

		int xmin = 0;
		Mapping best;
		best.k = params.k;
		best.P_sz = P_len;
		best.p_sz = p_len;
		best.xmin = 0;
		best.sim = 0.0;

		// Increase the left point end of the window [l,r) one by one. O(matches)
		int i = 0, j = 0;
		for(auto l = L.begin(), r = L.begin(); l != L.end(); ++l, ++i) {
			// Increase the right end of the window [l,r) until it gets out.
			for(; should_extend_right(diff_hist, L, l, r, xmin, P_len, params.k); ++r, ++j) {
				// If taking this kmer from T increases the intersection with P.
				if (--diff_hist[r->kmer_ord] >= 0)
					++xmin;
				assert (l->T_r <= r->T_r);
			}

			auto m = Mapping(params.k, P_len, p_len, (int)L.size(), l->T_r, prev(r)->T_r, pos_t(r-l), xmin);
			if (params.onlybest) {
				//if (xmin > best.xmin) {
				if (m.sim > best.sim) {
					best = m;
				}
			} else {
				if (m.sim > params.tThres) {
					mappings.push_back(m);
				}
			}

			// Prepare for the next step by moving `l` to the right.
			if (++diff_hist[l->kmer_ord] > 0)
				--xmin;

			assert(xmin >= 0);
		}
		assert (xmin == 0);

		if (params.onlybest) { // && best.J > params.tThres)
			mappings.push_back(best);
		}

		return mappings;
	}

	// Return only reasonable matches (i.e. those that are not J-dominated by
	// another overlapping match). Runs in O(|all|).
	vector<Mapping> filter_reasonable(const vector<Mapping> &all, const pos_t P_len) {
		vector<Mapping> reasonable;
		deque<Mapping> recent;

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
				if (recent.back().matches != -1) {
					// Take the leftmost mapping.
					reasonable.push_back(recent.back());
					// Mark the next closeby mappings as not reasonable
					for (auto it=recent.rbegin(); it!=recent.rend() && it->T_l - recent.back().T_l < sep; ++it)
						it->matches = -1;
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
				recent.front().matches = -1;
		}

		// 5. Add the last mapping if it is reasonable
		if (!recent.empty() && recent.back().matches != -1)
			reasonable.push_back(recent.back());

		return reasonable;
	}

  public:
	SweepMap(const SketchIndex &tidx, const params_t &params, Timers *T, Counters *C)
		: tidx(tidx), params(params), T(T), C(C) {
			if (params.tThres < 0.0 || params.tThres > 1.0) {
				cerr << "tThres = " << params.tThres << " outside of [0,1]." << endl;
				exit(1);
			}
		}

	void map(const string &pFile, const blmers_t &bLstmers) {
		T->start("mapping");
		// Load pattern sequences in batches
		T->start("query_reading");
		read_fasta_klib(pFile, [this, &bLstmers](kseq_t *seq) {
			T->stop("query_reading");
			T->start("sketching");
			Sketch p = buildFMHSketch(seq->seq.s, params.k, params.hFrac, bLstmers);
			T->stop("sketching");

			string seqID = seq->name.s;
			pos_t P_sz = (pos_t)seq->seq.l;

			C->inc("read_len", P_sz);
			vector<int> p_hist;

			T->start("seeding");
			auto seeds = choose_seeds(p, &p_hist);
			T->stop("seeding");

			T->start("matching");
			auto L = match_seeds(p, seeds, &p_hist);
			T->stop("matching");
			//print_matches(L);

			T->start("sweep");
			//print_sketches(seqID, p);
			auto mappings = sweep(p_hist, L, (pos_t)p.size(), P_sz);
			T->stop("sweep");

			T->start("postproc");
			auto reasonable_mappings = params.overlaps ? mappings : filter_reasonable(mappings, P_sz);

			for (auto &m: reasonable_mappings) {
				m.print_paf(params, seqID, tidx.T_sz, tidx.name);
				m.map_time = T->secs("sweep") / (double)mappings.size();
				C->inc("similarity", int(10000.0*m.sim));
				C->inc("mappings");
				C->inc("sketched_kmers", m.p_sz);
				C->inc("matches", m.matches);
			}
			C->inc("reads");
			if (reasonable_mappings.empty())
				C->inc("unmapped_reads");

			T->stop("postproc");
			T->start("query_reading");
		});
		T->stop("query_reading");
		T->stop("mapping");
	}

	void print_report(const Counters &C, const Timers &T) {
		cerr << fixed << setprecision(1);
		// TODO: report average Jaccard similarity
		cerr << "Params:" << endl;
		cerr << " | reference:             " << params.tFile << " (" << C.count("T_sz") << " nb)" << endl;
		cerr << " | queries:               " << params.pFile << endl;
		cerr << " | k:                     " << params.k << endl;
		cerr << " | w:                     " << params.w << endl;
		cerr << " | blacklist file:        " << (params.bLstFl.size() ? params.bLstFl : "-") << endl;
		//cerr << " | hFrac:                 " << params.hFrac << endl;
		cerr << " | max_seeds (S):         " << params.max_seeds << endl;
		cerr << " | max_matches (M):       " << params.max_matches << endl;
		cerr << " | onlybest:              " << params.onlybest << endl;
		cerr << " | tThres:                " << params.tThres << endl;
		cerr << "Stats:" << endl;
		cerr << " | Total reads:           " << C.count("reads") << " (~" << 1.0*C.count("read_len") / C.count("reads") << " nb per read)" << endl;
		cerr << " | Sketched read kmers:   " << C.count("sketched_kmers") << " (" << C.frac("sketched_kmers", "reads") << " per read)" << endl;
		cerr << " | Kmer matches:          " << C.count("matches") << " (" << C.frac("matches", "reads") << " per read)" << endl;
		cerr << " | Seed limit reached:    " << C.count("seeds_limit_reached") << " (" << C.perc("seeds_limit_reached", "reads") << "%)" << endl;
		cerr << " | Matches limit reached: " << C.count("matches_limit_reached") << " (" << C.perc("matches_limit_reached", "reads") << "%)" << endl;
		cerr << " | Unmapped reads:        " << C.count("unmapped_reads") << " (" << C.perc("unmapped_reads", "reads") << "%)" << endl;
		cerr << " | Average similarity:    " << C.frac("similarity", "mappings") / 10000.0 << endl;
		cerr << "Total time [sec]:         " << setw(5) << right << T.secs("total")     << " (" << setw(4) << right << T.secs("total") / C.count("reads") << " per read)" << endl;
		cerr << " | Index:                 " << setw(5) << right << T.secs("indexing")  << " (" << setw(4) << right << T.perc("indexing", "total") << "\% of total)" << endl;
		cerr << " |  | read:                  " << setw(5) << right << T.secs("index_reading")  << " (" << setw(4) << right << T.perc("index_reading", "indexing") << "\% of indexing)" << endl;
		cerr << " |  | sketch:                " << setw(5) << right << T.secs("index_sketching")<< " (" << setw(4) << right << T.perc("index_sketching", "indexing") << "\%)" << endl;
		cerr << " |  | initialize:            " << setw(5) << right << T.secs("index_initializing")<< " (" << setw(4) << right << T.perc("index_initializing", "indexing") << "\%)" << endl;
		cerr << " | Map:                   " << setw(5) << right << T.secs("mapping")   << " (" << setw(4) << right << T.perc("mapping", "total") << "\% of total)" << endl;
		cerr << " |  | read queries:          " << setw(5) << right << T.secs("query_reading")  << " (" << setw(4) << right << T.perc("query_reading", "mapping") << "\% of mapping)" << endl;
		cerr << " |  | sketch reads:          " << setw(5) << right << T.secs("sketching") << " (" << setw(4) << right << T.perc("sketching", "mapping") << "\%)" << endl;
		cerr << " |  | seeding the read:      " << setw(5) << right << T.secs("seeding")  << " (" << setw(4) << right << T.perc("seeding", "mapping") << "\%)" << endl;
		cerr << " |  |  | collect seed info:     " << setw(5) << right << T.secs("collect_seed_info")   << " (" << setw(4) << right << T.perc("collect_seed_info", "seeding") << "\% of seeding)" << endl;
		cerr << " |  |  | sort seeds:            " << setw(5) << right << T.secs("sort_seeds")   << " (" << setw(4) << right << T.perc("sort_seeds", "seeding") << "\%)" << endl;
		cerr << " |  | matching seeds:        " << setw(5) << right << T.secs("matching")  << " (" << setw(4) << right << T.perc("matching", "mapping") << "\%)" << endl;
		cerr << " |  |  | collect matches:       " << setw(5) << right << T.secs("collect_matches")   << " (" << setw(4) << right << T.perc("collect_matches", "matching") << "\% of matching)" << endl;
		cerr << " |  |  | sort matches:          " << setw(5) << right << T.secs("sort_matches")   << " (" << setw(4) << right << T.perc("sort_matches", "matching") << "\%)" << endl;
		cerr << " |  | sweep:                 " << setw(5) << right << T.secs("sweep")     << " (" << setw(4) << right << T.perc("sweep", "mapping") << "\%)" << endl;
		cerr << " |  | post proc:             " << setw(5) << right << T.secs("postproc")  << " (" << setw(4) << right << T.perc("postproc", "mapping") << "\%)" << endl;
//		cerr << "Total sketching time:     " << setw(5) << right << FMH_time.secs()     << " (" << setw(4) << right << FMH_time.secs() / T.secs("total") << "\%)" << endl; // << " (" << setw(4) << right << T.perc("postproc", "mapping") << "\%)" << endl;
	}
};

#endif

	// Returns the Jaccard score of the window [l,r)
	//inline double J(size_t p_sz, const vector<Match>::const_iterator l, const vector<Match>::const_iterator r, int xmin) {
	//	int s_sz = prev(r)->t_pos - l->t_pos + 1;
	//	assert(s_sz >= 0);
	//	//if (s_sz < 0) return 0;
	//	assert(p_sz + s_sz - xmin > 0);
	//	double scj = double(xmin) / double(p_sz + s_sz - xmin);
	//	if (!(0.0 <= scj && scj <= 1.0))
	//		cerr << "ERROR: scj=" << scj << ", xmin=" << xmin << ", p_sz=" << p_sz << ", s_sz=" << s_sz << ", l=" << l->t_pos << ", r=" << r->t_pos << endl;
	//	//assert (0.0 <= scj && scj <= 1.0);
	//	return scj;
	//}