#include <algorithm>
#include <iostream>
#include <iomanip>
#include <deque>
#include <set>

#include "Sketch.h"
#include "IO.h"

const float EPS = 1e-7;
unordered_map<uint64_t, char> bLstmers;

clock_t sorting_time = 0;
clock_t total_match_kmers_time = 0;

struct Match {
	uint32_t kmer_ord;
	uint32_t P_l;
	uint32_t T_r;
	uint32_t t_pos;
	bool strand;  // 0 - forward, 1 - reverse

	Match() {}
//	Match(uint32_t _kmer_ord, uint32_t _T_pos, uint32_t _t_pos)
//		: kmer_ord(_kmer_ord), T_r(_T_pos), t_pos(_t_pos) {}

	bool operator<(const Match &other) {
		return T_r < other.T_r;
	}

	char get_strand() const {
		return strand ? '-' : '+';
	}
};

struct Mapping {
	uint32_t k; 	   // kmer size
	uint32_t P_sz;     // pattern size |P| bp 
	uint32_t p_sz;     // pattern sketch size |p| kmers
	uint32_t matches;  // L.size() -- total number of matches in `t' 
	uint32_t T_l;      // the position of the leftmost nucleotide of the mapping
	uint32_t T_r;      // the position of the rightmost nucleotide of the mapping
	uint32_t xmin;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	uint32_t dT_l;      // delta to be applied before output
	uint32_t dT_r;      // -- || --
	double J;          // Jaccard score/similarity [0;1]
	double map_time;

	Mapping(uint32_t k=0, uint32_t P_sz=0, uint32_t p_sz=0, uint32_t matches=0, uint32_t T_l=0, uint32_t T_r=0, uint32_t xmin=0, uint32_t dT_l=0, uint32_t dT_r=0, double J=0.0)
		: k(k), P_sz(P_sz), p_sz(p_sz), matches(matches), T_l(T_l), T_r(T_r), xmin(xmin), dT_l(dT_l), dT_r(dT_r), J(J) {}

	void print() {
		cerr << "k=" << k 
			<< " |P|=" << P_sz
			<< " |p|=" << p_sz
			<< " matches=" << matches
			<< " T_l=" << T_l
			<< " T_r=" << T_r
			<< " xmin=" << xmin
			<< " dT_l=" << dT_l
			<< " dT_r=" << dT_r
			<< " J=" << J
			<< " map_time=" << map_time
			<< endl;
	}
};

template<typename TT> auto prev(const typename TT::iterator &it) {
    auto pr = it; return --pr;
}

//template<typename TT> auto next(const typename TT::iterator &it) {
//    auto pr = it; return ++pr;
//}

class Sweep {
	const mm_idx_t *tidx;
	const params_t &params;

	void unpack(uint64_t idx, Match *m) {
		m->T_r    = (uint32_t)(idx >> 32);
		m->t_pos  = (uint32_t)(((idx << 32) >> 32) >> 1);
		m->strand = (uint32_t(idx) << 31) >> 31;
		//cerr << "strand: " << m->strand << endl;
	}

	// return if the kmer has not been seen before
	bool add2hist(const uint64_t kmer_hash,
			vector<int32_t> *hist, unordered_map<uint64_t, uint32_t> *hash2ord, vector<uint64_t> *ord2hash, uint32_t *kmer_ord) {
		auto ord_it = hash2ord->find(kmer_hash);
		if (ord_it != hash2ord->end()) {
			*kmer_ord = ord_it->second;
			++(*hist)[*kmer_ord];
			return false;
		} else {
			*kmer_ord = hist->size();
			hash2ord->insert({kmer_hash, *kmer_ord});  //(*hash2ord)[kmer_hash] = kmer_ord;
			ord2hash->push_back(kmer_hash);
			hist->push_back(1);
			assert(hist->size() == ord2hash->size());
			return true;
		}
	}

	struct idx_result_t {
		int32_t kmer_pos;
		uint32_t kmer_ord;
		const uint64_t *idx_p;
		int nHits;

		idx_result_t(int32_t kmer_pos, uint32_t kmer_ord, const uint64_t *idx_p, int nHits) :
			kmer_pos(kmer_pos), kmer_ord(kmer_ord), idx_p(idx_p), nHits(nHits) {}

		bool operator<(const idx_result_t &other) {
			return nHits < other.nHits;
		}
	};

	// Initializes the histogram of the pattern and the list of matches
	void match_kmers(
			// input
			const Sketch& p,
			// output
			vector<int32_t> *p_hist,
			vector<Match> *L,
			vector<uint64_t> *L2) {
		unordered_map<uint64_t, uint32_t> hash2ord;
		vector<uint64_t> ord2hash;
		L->reserve(P_MULTIPLICITY * p.size());

		vector<idx_result_t> match_lists;
		match_lists.reserve(p.size());

		uint64_t kmer_hash, prev_kmer_hash = 0;
		for(auto p_it = p.begin(); p_it != p.end(); ++p_it, prev_kmer_hash = kmer_hash) {
			auto kmer_pos = p_it->first;
			kmer_hash = p_it->second;
			uint32_t kmer_ord;

			// Add xor'ed consecuent pattern kmers to p_hist 
			if (p_it != p.begin()) {
				auto xor_hash = prev_kmer_hash ^ kmer_hash;
				(void)add2hist(xor_hash, p_hist, &hash2ord, &ord2hash, &kmer_ord);
			}

			// Add pattern kmers from the pattern to p_hist 
			if (add2hist(kmer_hash, p_hist, &hash2ord, &ord2hash, &kmer_ord)) {
				// Add kmer matches to L
				int nHits;
				auto *idx_p = mm_idx_get(tidx, kmer_hash, &nHits);
				if (nHits > 0)
					match_lists.push_back(idx_result_t(kmer_pos, kmer_ord, idx_p, nHits));
			}
		}

		//Sort by number of hits and get MAX_SEEDS of kmers with the lowest number of hits.
		sort(match_lists.begin(), match_lists.end());
		int total_hits = 0;
		//const int MAX_SEEDS = 1000, MAX_TOTAL_HITS = 100000;
		for (int seed=0; seed<(int)match_lists.size(); seed++) {
			if (seed > params.max_seeds) {
				++seeds_limit_reached;
				break;
			}
			auto &res = match_lists[seed];
			//cerr << "seed=" << seed << ", nHits=" << res.nHits << endl;
			for (int i = 0; i < res.nHits; ++i, ++res.idx_p) {					// Iterate over all occurrences
				Match m;
				m.P_l = res.kmer_pos;
				m.kmer_ord = res.kmer_ord;
				unpack(*res.idx_p, &m);
				L->push_back(m); // Push (k-mer ord in p, k-mer position in reference, k-mer position in sketch) pair
			}
			if ((total_hits += res.nHits) > params.max_matches) {
				++matches_limit_reached;
				break;
			}
		}

		clock_t _ = clock();
		//Sort L by ascending positions in reference
		sort(L->begin(), L->end());
		sorting_time += clock() - _;

		// Add elastic kmers
		if (params.elastic == elastic_t::consecutive) {
			if (L->size() == 0)
				return;

			// Add xor'ed consecuent text kmers to L2
			L2->reserve(L->size());
			L2->push_back(0);
			for (int i=1; i<L->size(); i++) {
				auto xor_h = ord2hash[ (*L)[i-1].kmer_ord ] ^ ord2hash[ (*L)[i].kmer_ord ];
				auto it = hash2ord.find(xor_h); 
				L2->push_back(it != hash2ord.end() ? it->second : p_hist->size());
			}

			if (L->size() != L2->size())
				cerr << L->size() << " != " << L2->size() << endl;
			assert(L->size() == L2->size());
		}
	}

	// Returns the Jaccard score of the window [l,r)
	double J(size_t p_sz, const vector<Match>::iterator l, const vector<Match>::iterator r, int xmin) {
		int s_sz = prev(r)->t_pos - l->t_pos + 1;
		if (s_sz < 0)
			return 0;
		if (params.elastic == elastic_t::consecutive) {
			p_sz *= 2, s_sz *= 2;  // to account for L2
		}
		assert(p_sz + s_sz - xmin > 0);
		double scj = 1.0 * xmin / (p_sz + s_sz - xmin);
		//if (scj >= 0.9)
		//	cerr << "J: p_sz=" << p_sz << ", s_sz=" << s_sz << ", xmin=" << xmin << ", l=" << l->t_pos << ", r=" << r->t_pos << ", scj=" << scj << endl;
		if (!(0.0 <= scj && scj <= 1.0))
			cerr << "ERROR: scj=" << scj << ", xmin=" << xmin << ", p_sz=" << p_sz << ", s_sz=" << s_sz << ", l=" << l->t_pos << ", r=" << r->t_pos << endl;
		//assert (0.0 <= scj && scj <= 1.0);
		return scj;
	}

	// Returns true if the window [l,r) should be extended to the right
	inline bool should_extend_right(
			const Sketch &p,
			const vector<int32_t> &hist,
			const vector<Match> &L,
			const vector<Match>::iterator l, const vector<Match>::iterator r,
			const int xmin, const int Plen_nucl, const int k) {
		if (r == L.end())
			return false;
		bool extension_stays_within_window = r->T_r + k <= l->T_r + Plen_nucl;
		return extension_stays_within_window;
		//if (extension_stays_within_window)
		//	return true;
		//bool extension_improves_jaccard = J(p, l, r, xmin) < J(p, l, next(r), xmin + (hist[r->kmer_ord] > 0));
		//assert(!extension_improves_jaccard);
		//return extension_improves_jaccard;
	}

  public:
	// stats
	int seeds_limit_reached = 0;
	int matches_limit_reached = 0;

	Sweep(const mm_idx_t *tidx, const params_t &params)
		: tidx(tidx), params(params) {
			if (params.tThres < 0.0 || params.tThres > 1.0) {
				cerr << "tThres = " << params.tThres << " outside of [0,1]." << endl;
				exit(1);
			}
		}

	const vector<Mapping> map(const Sketch& p, const uint32_t Plen_nucl) {
		vector<int32_t> diff_hist;  // rem[kmer_hash] = #occurences in `p` - #occurences in `s`
		vector<Match> L;    		// for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
		vector<uint64_t> L2;		// for all consecutive matches in L
		vector<Mapping> mappings;	// List of tripples <i, j, score> of matches

		multiset<int32_t> P_l_set;

		int xmin = 0;
		Mapping best;
		best.k = params.k;
		best.P_sz = Plen_nucl;
		best.p_sz = p.size();

		//cerr << "Mapping " << p.size() << " kmers of length " << Plen_nucl << endl;

		clock_t start_match_kmers = clock();
		match_kmers(p, &diff_hist, &L, &L2);
		total_match_kmers_time += clock() - start_match_kmers;


		// Increase the left point end of the window [l,r) one by one.
		// O(matches)
		int i = 0, j = 0;
		for(auto l = L.begin(), r = L.begin(); l != L.end(); ++l, ++i) {
			// Increase the right end of the window [l,r) until it gets out.
			for(; should_extend_right(p, diff_hist, L, l, r, xmin, Plen_nucl, params.k); ++r, ++j) {
				P_l_set.insert(r->P_l);
				// If taking this kmer from T increases the intersection with P.
				if (--diff_hist[r->kmer_ord] >= 0)
					++xmin;
				assert (l->T_r <= r->T_r);

				if (params.elastic == elastic_t::consecutive) {
					auto pair_ord = L2[r-L.begin()];
					if (pair_ord<diff_hist.size() && --diff_hist[pair_ord] >= 0)
						++xmin;
				}
			}

			auto curr_J = J(p.size(), l, r, xmin);
			auto m = Mapping(params.k, Plen_nucl, p.size(), L.size(), l->T_r, prev(r)->T_r, xmin, 0, 0, curr_J);  // TODO: create only if curr_J is high enough

			if (params.alignment_edges == alignment_edges_t::fine) {
				if (P_l_set.size() > 0) {
					m.dT_l = - *P_l_set.begin();
					m.dT_r = Plen_nucl - (*P_l_set.rbegin()+params.k);
				}
			}

			if (params.onlybest) {
				if (curr_J > best.J) {
					best = m;
				}
			} else {
				if (curr_J > params.tThres) {
					mappings.push_back(m);
				}
			}

			P_l_set.erase(P_l_set.find(l->P_l));
			// Prepare for the next step by moving `l` to the right.
			if (++diff_hist[l->kmer_ord] > 0)
				--xmin;

			if (params.elastic == elastic_t::consecutive) {
				auto pair_ord = L2[l-L.begin()];
				if (pair_ord<diff_hist.size() && ++diff_hist[pair_ord] > 0)
					--xmin;
			}

			assert(xmin >= 0);
		}
		assert (xmin == 0);

		if (params.onlybest) { // && best.J > params.tThres)
			mappings.push_back(best);
		}

		return mappings;
	}
};

// Return only reasonable matches (i.e. those that are not J-dominated by
// another overlapping match). Runs in O(|all|).
vector<Mapping> filter_reasonable(const params_t &params, const vector<Mapping> &all, const uint32_t Plen_nucl) {
	vector<Mapping> reasonable;
	deque<Mapping> recent;

	// Minimal separation between mappings to be considered reasonable
	uint32_t sep = (1.0 - params.tThres) * Plen_nucl;

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
		while(!recent.empty() && recent.front().J < next.J - EPS)
			recent.pop_front();
		assert(recent.empty() || (recent.back().J >= recent.front().J - EPS && recent.front().J >= next.J - EPS));

		// 3. Add the next mapping to the front
		recent.push_front(next);	 

		// 4. If there is another mapping in the deque, it is near and better.
		// Mark the 
		if (recent.size() > 1)
			recent.front().matches = -1;
	}

	// 5. Add the last mapping if it is reasonable
	if (!recent.empty() && recent.back().matches != -1)
		reasonable.push_back(recent.back());

	return reasonable;
}

// Extends the mappings to the left and right
void normalize_mappings(const params_t &params, vector<Mapping> *mappings, const uint32_t pLen, const string &seqID, const uint32_t T_sz, const string &text) {
	for(auto &m: *mappings) {
		if (params.alignment_edges == alignment_edges_t::extend_equally) {
			int span   = m.T_r - m.T_l + 1;
			//assert(span <= m.P_sz);
			int shift  = (m.P_sz - span) / 2;       assert(shift >= 0);
			m.T_l -= shift; m.T_l = max((uint32_t)0, m.T_l);
			m.T_r += shift; m.T_r = min(T_sz-1, m.T_r);
		} else if (params.alignment_edges == alignment_edges_t::fine) {
			m.T_l += m.dT_l;
			m.T_r += m.dT_r;
		}
	}
}

// Outputs all given mappings
void mappings2paf(const params_t &params, const vector<Mapping>& res, const uint32_t P_sz, const uint32_t pLen, const string &seqID, const uint32_t T_sz, const char *T_name, const string &text) {
	for(auto m: res) {
		//m.print();
		// --- https://github.com/lh3/miniasm/blob/master/PAF.md ---
		cout << seqID << "\t"  // Query sequence name
			<< P_sz << "\t"   // query sequence length
			<< 0 << "\t"   // query start (0-based; closed)
			<< P_sz << "\t"   // query end (0-based; open)
			<< "+" << "\t"   // m.get_strand() //strand; TODO
			<< T_name << "\t"   // reference name
			<< T_sz << "\t"  // target sequence length
			<< m.T_l << "\t"  // target start on original strand (0-based)
			<< m.T_r << "\t"  // target start on original strand (0-based)
			<< -1 << "\t"  // Number of residue matches
			<< -1 << "\t"  // Alignment block length
			<< 255 << "\t"  // Mapping quality (0-255; 255 for missing)
		// ----- end of required PAF fields -----
			<< "k:i:" << m.k << "\t"
			//<< "P:i:" << m.P_sz << "\t"  // redundant
			<< "p:i:" << m.p_sz << "\t"
			<< "M:i:" << m.matches << "\t" // matches of `p` in `s` [kmers]
			<< "I:i:" << m.xmin << "\t"  // intersection of `p` and `s` [kmers]
			<< "J:f:" << m.J << "\t"  // Jaccard similarity [0; 1]
			<< "t:f:" << m.map_time << endl;
        if (!text.empty())
            cerr << "   text: " << text.substr(m.T_l, m.T_r-m.T_l+1) << endl;
	}
}

int main(int argc, char **argv) {
	double total_sketching_time=0, total_mapping_time=0, total_time=0;
	clock_t total_start = clock();
	int total_reads = 0, unmapped_reads = 0;
	double total_J = 0.0;
	int total_mappings = 0;

	reader_t reader;
	auto res = reader.init(argc, argv);
	const params_t &params = reader.params;
	assert(res == 0);

	Sweep sweep(reader.tidx, params);

	if (!params.paramsFile.empty()) {
		cerr << "Writing parameters to " << params.paramsFile << "..." << endl;
		auto fout = ofstream(params.paramsFile);
		reader.params.print(fout, false);
	} else {
		reader.params.print(cerr, true);
	}

	//cerr << "aligning reads from " << params.pFile << "..." << endl;

	// Load pattern sequences in batches
	while(true) {
		cerr << "Sketching a batch of reads..." << endl;
		clock_t start_sketching = clock();
		if (!lMiniPttnSks(reader.fStr, params.k, reader.tidx->w, bLstmers, reader.pSks) && reader.pSks.empty())
			break;
		total_sketching_time += clock() - start_sketching;

		cerr << "Aligning a batch of reads..." << endl;
		clock_t start_mapping = clock();
		// Iterate over pattern sketches
		for(auto p = reader.pSks.begin(); p != reader.pSks.end(); ++p) {
			auto seqID = get<0>(*p);
            auto P_sz = get<1>(*p);
            auto sks = get<2>(*p);
			clock_t begin_map_time = clock();

            auto mappings = sweep.map(sks, P_sz);
			auto reasonable_mappings = params.overlaps ? mappings : filter_reasonable(params, mappings, P_sz);
			normalize_mappings(params, &reasonable_mappings, sks.size(), seqID, reader.T_sz, reader.text);

			double all_map_time = 1.0 * (clock() - begin_map_time) / CLOCKS_PER_SEC;
			for (auto &m: reasonable_mappings) {
				m.map_time = all_map_time / mappings.size();
				total_J += m.J;
				total_mappings ++;
			}
			mappings2paf(params, reasonable_mappings, P_sz, sks.size(), seqID, reader.T_sz, reader.tidx->seq->name, reader.text);

			// stats
			++total_reads;
			if (mappings.empty()) ++unmapped_reads;
		}
		reader.pSks.clear();
		total_mapping_time += clock() - start_mapping;
	}

	total_time = double(clock() - total_start) / CLOCKS_PER_SEC;
	total_sketching_time /= CLOCKS_PER_SEC;
	total_mapping_time /= CLOCKS_PER_SEC;
	sorting_time /= CLOCKS_PER_SEC;
	total_match_kmers_time /= CLOCKS_PER_SEC;

	// TODO: report average Jaccard similarity
	cerr << "Total reads:          " << total_reads << endl;
	cerr << "Seed limit reached:   " << sweep.seeds_limit_reached << " (" << 100.0 * sweep.seeds_limit_reached / total_reads << ")" << endl;
	cerr << "Matches limit reached:" << sweep.matches_limit_reached << " (" << 100.0 * sweep.matches_limit_reached / total_reads << ")" << endl;
	cerr << "Unmapped reads:       " << unmapped_reads << " (" << 100.0 * unmapped_reads / total_reads << "%)" << endl;
	cerr << "Average J:            " << total_J / total_mappings << endl;
	cerr << endl;
	cerr << "Total time:           " << total_time << endl;
	cerr << "Average time per read:" << total_time / total_reads << endl;
	cerr << "Total read sketching: " << total_sketching_time << " (" << 100*total_sketching_time / total_time << "\% of total)" << endl;
	cerr << "Total mapping:        " << total_mapping_time << " (" << 100*total_mapping_time / total_time << "\% of total)" << endl;
	cerr << "        Sorting time: " << sorting_time << " (" << 100*sorting_time / total_time << "\% of total)" << endl;
	cerr << "      Matching kmers: " << total_match_kmers_time << " (" << 100*total_match_kmers_time / total_time << "\% of total)" << endl;

	return 0;
}