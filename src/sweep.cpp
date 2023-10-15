#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Sketch.h"
#include "IO.h"
#include "Index.h"

unordered_map<uint64_t, char> bLstmers;

struct Match {
	uint32_t kmer_ord;
	uint32_t T_pos;
	uint32_t t_pos;

	Match() {}
	Match(uint32_t _kmer_ord, uint32_t _T_pos, uint32_t _t_pos)
		: kmer_ord(_kmer_ord), T_pos(_T_pos), t_pos(_t_pos) {}

	bool operator<(const Match &other) {
		return T_pos < other.T_pos;
	}
};

struct Mapping {
	uint32_t k; 		// kmer size
	uint32_t P_sz;     // pattern size |P| bp 
	uint32_t p_sz;     // pattern sketch size |p| kmers
	uint32_t matches;  // L.size()
	uint32_t l_T_pos;
	uint32_t r_T_pos;
	uint32_t xmin;
	double J;  			// Jaccard score

	Mapping(uint32_t k=0, uint32_t P_sz=0, uint32_t p_sz=0, uint32_t matches=0, uint32_t l_T_pos=0, uint32_t r_T_pos=0, uint32_t xmin=0, double J=0.0)
		: k(k), P_sz(P_sz), p_sz(p_sz), matches(matches), l_T_pos(l_T_pos), r_T_pos(r_T_pos), xmin(xmin), J(J) {}
};

//This function outputs all given t-homologies
inline void outputMappings(const params_t &params, const vector<Mapping>& res, const uint32_t pLen, const string &seqID, const uint32_t T_sz, const string &text){
	//Iterate over t-homologies
	for(auto m: res) {
		if (params.alignment_edges == alignment_edges_t::extend_equally) {
			int span   = m.r_T_pos - m.l_T_pos + 1;  assert(span <= m.P_sz);
			int shift  = (m.P_sz - span) / 2;       assert(shift >= 0);
			m.l_T_pos -= shift; m.l_T_pos = max((uint32_t)0, m.l_T_pos);
			m.r_T_pos += shift; m.r_T_pos = min(T_sz-1, m.r_T_pos);
		}

		cout << seqID << "\t" << m.k << "\t" << m.P_sz << "\t" << m.p_sz << "\t"
			<< m.matches << "\t" << m.l_T_pos << "\t" << m.r_T_pos << "\t" << m.xmin << "\t" << m.J << endl;
        if (!text.empty())
            cout << "   text: " << text.substr(m.l_T_pos, m.r_T_pos-m.l_T_pos+1) << endl;
	}
}

template<typename TT> auto prev(const typename TT::iterator &it) {
    auto pr = it; return --pr;
}

template<typename TT> auto next(const typename TT::iterator &it) {
    auto pr = it; return ++pr;
}

class Sweep {
	const mm_idx_t *tidx;
	const params_t &params;

	void unpack(uint64_t idx, uint32_t *T_pos, uint32_t *t_pos) {
		*T_pos = (uint32_t)(idx >> 32);
		*t_pos = (uint32_t)(((idx << 32) >> 32) >> 1);
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

	void init(
			// input
			const Sketch& p,
			// output
			vector<int32_t> *p_hist,
			vector<Match> *L,
			vector<uint64_t> *L2) {
		unordered_map<uint64_t, uint32_t> hash2ord;
		vector<uint64_t> ord2hash;
		L->reserve(P_MULTIPLICITY * p.size());

		uint64_t kmer_hash, prev_kmer_hash = 0;
		for(auto p_it = p.begin(); p_it != p.end(); ++p_it, prev_kmer_hash = kmer_hash) {
			kmer_hash = *p_it;
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
				for (int i = 0; i < nHits; ++i, ++idx_p) {					// Iterate over all occurrences
					Match m;
					m.kmer_ord = kmer_ord;
					unpack(*idx_p, &m.T_pos, &m.t_pos);
					L->push_back(m); // Push (k-mer ord in p, k-mer position in reference, k-mer position in sketch) pair
				}
			}
		}

		//Sort L by ascending positions in reference
		sort(L->begin(), L->end());

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

	//	for (auto it: hash2ord)
	//		cout << "hash2ord: " << it.first << " " << it.second << endl;
	//	for (auto it: *L)
	//		cout << "L: " << setw(5) << right << it.kmer_ord
	//			 << " " << setw(10) << right << it.T_pos
	//			 << " " << setw(10) << right << it.t_pos << endl;
	//	for (int i=0; i<hist->size(); i++)
	//		cout << "hist[" << i << "]: " << (*hist)[i] << endl;
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
		bool extension_stays_within_window = r->T_pos + k <= l->T_pos + Plen_nucl;
		return extension_stays_within_window;
		//if (extension_stays_within_window)
		//	return true;
		//bool extension_improves_jaccard = J(p, l, r, xmin) < J(p, l, next(r), xmin + (hist[r->kmer_ord] > 0));
		//assert(!extension_improves_jaccard);
		//return extension_improves_jaccard;
	}

  public:
	Sweep(const mm_idx_t *tidx, const params_t &params)
		: tidx(tidx), params(params) {}

	const vector<Mapping> map(const Sketch& p, const uint32_t Plen_nucl) {
		vector<int32_t> diff_hist;  // rem[kmer_hash] = #occurences in `p` - #occurences in `s`
		vector<Match> L;    		// for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
		vector<uint64_t> L2;		// for all consecutive matches in L
		vector<Mapping> res;		// List of tripples <i, j, score> of matches

		int xmin = 0;
		Mapping best;
		best.k = params.k;
		best.P_sz = Plen_nucl;
		best.p_sz = p.size();

		init(p, &diff_hist, &L, &L2);

		// Increase the left point end of the window [l,r) one by one.
		int i = 0, j = 0;
		for(auto l = L.begin(), r = L.begin(); l != L.end(); ++l, ++i) {
			// Increase the right end of the window [l,r) until it gets out.
			for(; should_extend_right(p, diff_hist, L, l, r, xmin, Plen_nucl, params.k); ++r, ++j) {
				// If taking this kmer from T increases the intersection with P.
				if (--diff_hist[r->kmer_ord] >= 0)
					++xmin;
				assert (l->T_pos <= r->T_pos);

				if (params.elastic == elastic_t::consecutive) {
					auto pair_ord = L2[r-L.begin()];
					if (pair_ord<diff_hist.size() && --diff_hist[pair_ord] >= 0)
						++xmin;
				}
			}
			
			// Update best.
			auto curr_J = J(p.size(), l, r, xmin);
			if (curr_J > best.J)
				best = Mapping(params.k, Plen_nucl, p.size(), L.size(), l->T_pos, prev(r)->T_pos, xmin, curr_J);

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

		res.push_back(best);
		return res;
	}
};

int main(int argc, char **argv){
	reader_t reader;
	auto res = reader.init(argc, argv);
	const params_t &params = reader.params;
	assert(res == 0);

	Sweep sweep(reader.tidx, params);
	reader.params.print(cerr);

	//Load pattern sequences in batches
	while(lMiniPttnSks(reader.fStr, params.k, reader.tidx->w, bLstmers, reader.pSks) || !reader.pSks.empty()){ //TODO: This function still needs to be tested!
		//Iterate over pattern sketches
		for(auto p = reader.pSks.begin(); p != reader.pSks.end(); ++p){
			auto seqID = get<0>(*p);
            auto plen_nucl = get<1>(*p);
            auto sks = get<2>(*p);

			//Calculate an adapted threshold if we have the necessary informations
			//float curr_tThres = (reader.dec != 0 && reader.inter != 0) ? reader.dec * plen_nucl + reader.inter : reader.tThres;

			//Find t-homologies and output them
            auto res = sweep.map(sks, plen_nucl);
			outputMappings(params, res, sks.size(), seqID, reader.T_sz, reader.text);//TODO: Tests for this function need to be adaptated!
		}

		reader.pSks.clear();
	}

	return 0;
}
