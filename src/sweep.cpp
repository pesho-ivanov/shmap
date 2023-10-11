#include <algorithm>
#include <iostream>
#include <iomanip>

#include "Sketch.h"
#include "IO.h"
#include "Index.h"

//The FracMinHash ratio
double hFrac = HASH_RATIO;
unordered_map<uint64_t, char> bLstmers;

//unordered_map<uint65_t, uint32_t> zero_occurences(const Sketch& skP) {
//    unordered_map<uint64_t, uint32_t> occs(skP.size());
//
//    //Fill occp
//    for(Sketch::const_iterator fSkIt = skP.begin(); fSkIt != skP.end(); ++fSkIt)
//        occs[*fSkIt] = 0;
//    
//    return occs;
//}

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

template<typename TT> auto prev(const typename TT::iterator &it) {
    auto pr = it; return --pr;
}

template<typename TT> auto next(const typename TT::iterator &it) {
    auto pr = it; return ++pr;
}

//void print_L(const vector<Match> &L) {
//    int i=0;
//    for (const auto &s: L) {
//        cerr << "L[" << i << "]=(" << s.kmer_ord << ", " << s.T_pos << ", " << s.t_pos << ")" << endl;
//        ++i;
//    }
//}
//
//void print_skP(const Sketch& skP) {
//    int i=0;
//    for (const auto &p: skP) {
//        cerr << "skP[" << i << "]=" << p << endl;
//        ++i;
//    }
//}
//
//const int window2score(const vector<Match> &L, const Sketch& skP, int from_nucl, int to_nucl, int plen_nucl) {
//    int xmin = 0;
//    auto l = L.end(), r = L.begin();
//	for(auto s = L.begin(); s != L.end(); ++s) {
//        if (from_nucl <= s->second && s->second + plen_nucl < to_nucl) {
//            cout << "[" << from_nucl << ", " << to_nucl << "] includes sketch (" << s->first << ", " << s->second << ")" << endl; 
//            if (s->second < l->second) l=s;
//            if (s->second >= r->second) r=s;
//            ++xmin;
//        }
//    }
//
//    if (l == L.end() && r == L.begin())
//        l = r;
//
//    for (auto s=l; s<r; ++s)
//        assert(from_nucl <= s->second && to_nucl <= s->second + plen_nucl);
//
//    auto scj = scoreX1000(skP, xmin, r-l);
//    //auto scj = 0;
//    cout << "[" << from_nucl << ", " << to_nucl << "] includes " << xmin << " sketches,"  << ", r-l=" << (r-l) <<  ", scj=" << scj << endl;
//    return xmin;
//}

void unpack(uint64_t idx, uint32_t *T_pos, uint32_t *t_pos) {
	*T_pos = (uint32_t)(idx >> 32);
	*t_pos = (uint32_t)(((idx << 32) >> 32) >> 1);
	//uint32_t T_pos = (uint32_t)((*idx_p) >> 32); 					// Position in reference
	//uint32_t t_pos = (uint32_t)((( (*idx_p) << 32) >> 32) >> 1); 	// Position in sketch
	//uint32_t t_pos = ( (*idx_p) ^ ((uint64_t)T_pos << 32) ) >> 1; 	// Position in sketch
}

void init(
		// input
		const Sketch& p,
		const mm_idx_t *tidx,
		// output
		vector<int32_t> *hist,
		vector<Match> *L) {

	unordered_map<uint64_t, uint32_t> hash2ord;
	uint32_t kmers = 0;
	L->reserve(P_MULTIPLICITY * p.size());

    for(auto p_it = p.begin(); p_it != p.end(); ++p_it) {
		uint32_t kmer_ord;
		auto kmer_hash = *p_it;
		auto ord_it = hash2ord.find(kmer_hash);
		if (ord_it != hash2ord.end()) {
			++(*hist)[ord_it->second];
		} else {
			kmer_ord = hist->size();
			hash2ord[kmer_hash] = kmer_ord;
			hist->push_back(1);

			int nHits;
			auto *idx_p = mm_idx_get(tidx, kmer_hash, &nHits);
			for (int i = 0; i < nHits; ++i, ++idx_p) {					// Iterate over all occurrences
				Match m;
				m.kmer_ord = kmer_ord;
				unpack(*idx_p, &m.T_pos, &m.t_pos);
				//cout << "kmer_hash=" << m.kmer_ord << ", T_pos=" << m.T_pos << ", t_pos=" << m.t_pos << endl;
				L->push_back(m); // Push (k-mer ord in p, k-mer position in reference, k-mer position in sketch) pair
			}
		}
	}

	//Sort L by ascending positions in reference
	sort(L->begin(), L->end());

//	for (auto it: hash2ord)
//		cout << "hash2ord: " << it.first << " " << it.second << endl;
//	for (auto it: *L)
//		cout << "L: " << setw(5) << right << it.kmer_ord
//			 << " " << setw(10) << right << it.T_pos
//			 << " " << setw(10) << right << it.t_pos << endl;
//	for (int i=0; i<hist->size(); i++)
//		cout << "hist[" << i << "]: " << (*hist)[i] << endl;
}

double J(const Sketch& p, const vector<Match>::iterator l, const vector<Match>::iterator r, int xmin) {
	int s_sz = prev(r)->t_pos - l->t_pos + 1;
	if (s_sz < 0)
		return 0;
    assert (s_sz >= 0);
    assert(p.size() + s_sz - xmin > 0);
    double scj = 1.0 * xmin / (p.size() + s_sz - xmin);
    assert (0 <= scj && scj <= 1.0);
    return scj;
}

inline bool better_extend_right(const Sketch &p, const vector<int32_t> &hist, const vector<Match> &L,
		const vector<Match>::iterator l, const vector<Match>::iterator r,
		int xmin, const int Plen_nucl, const int k) {
	if (r == L.end())
		return false;
	bool extension_stays_within_window = r->T_pos + k <= l->T_pos + Plen_nucl;
	auto r_next = next(r);
	bool extension_improves_jaccard = r_next != L.end() && J(p, l, r, xmin) < J(p, l, r_next, xmin + (hist[r_next->kmer_ord] > 0));
	return extension_stays_within_window || extension_improves_jaccard;
}

struct Mapping {
	uint32_t l_T_pos;
	uint32_t r_T_pos;
	double J;  			// Jaccard score
	Mapping() {}
	Mapping(uint32_t _l_T_pos, uint32_t _r_T_pos, double _J)
		: l_T_pos(_l_T_pos), r_T_pos(_r_T_pos), J(_J) {}
};

//This function outputs all given t-homologies
inline void outputMappings(const vector<Mapping>& res, const uint32_t& pLen, const string &seqID, const string &text){
	//Iterate over t-homologies
	for(const auto &m: res) {
        cout << seqID << "\t" << m.l_T_pos << "\t" << m.r_T_pos << "\t" << m.J << endl;
        if (!text.empty())
            cout << "   text: " << text.substr(m.l_T_pos, m.r_T_pos-m.l_T_pos+1) << endl;
	}
}


const vector<Mapping> sweep(const Sketch& p, const mm_idx_t *tidx, const int Plen_nucl, const int k) {
    vector<int32_t> hist;  		// rem[kmer_hash] = #occurences in `p` - #occurences in `s`
    vector<Match> L;    		// for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
	vector<Mapping> res;		// List of tripples <i, j, score> of matches

    int xmin = 0;
    Mapping best(0, 0, 0);	// <i, j, score>

    init(p, tidx, &hist, &L);
 
    //cerr << "|P| = " << Plen_nucl << ", |p| = " << skP.size() << endl;
    //cerr << "|L| = " << L.size() << endl;
    //cerr << "k = " << k <<endl;

    // Increase the left point end of the window [l,r) one by one.
    int i = 0, j = 0;
	for(auto l = L.begin(), r = L.begin(); l != L.end(); ++l, ++i) {
        // Increase the right end of the window [l,r) until it gets out.
        for(; better_extend_right(p, hist, L, l, r, xmin, Plen_nucl, k); ++r, ++j) {
			// If taking this kmer from T increases the intersection with P 
			if (--hist[r->kmer_ord] >= 0)
				++xmin;
			assert (l->T_pos <= r->T_pos);
        }
		
        auto curr_J = J(p, l, r, xmin);
        //cout << "i=" << i << ", j=" << j << ", l=" << l->second << ", r=" << r->second << ", |s|=" << s_sz << ", xmin=" << xmin << ", |s|=" << (r-l) <<  ", scj=" << scj << endl;

        if (curr_J > best.J) {
            //cout << "s=r-l=" << r-l << endl;
            //cout << "denom: " << skP.size() + (r-l) - xmin << endl;
            best = Mapping(l->T_pos, prev(r)->T_pos, curr_J);
        }

        // Prepare for the next step by moving `l` to the right
		if (++hist[l->kmer_ord] > 0)
			--xmin;

        assert(xmin >= 0);
    }
	assert (xmin == 0);

    //auto scj = window2score(L, skP, 4681377, 4683857, plen_nucl);
    //cout << "l=" << 4681377 << ", r=" << 4683857 << ", xmin=" << xmin << ", r-l=" << 4683857-4681377 <<  ", scj=" << scj << endl;

    res.push_back(best);
    return res;
}

int main(int argc, char **argv){
	bool normalize = NORM_FLAG_DEFAULT; 		//Flag to save that scores are to be normalized
	bool noNesting = NESTING_FLAG_DEFAULT; 		//Flag to state if we are interested in nested results
	uint32_t kmerLen = K; 						//The k-mer length
	uint32_t w = W; 							//The window size
	uint32_t comWght = DEFAULT_WEIGHT; 			//Scoring weights
	float uniWght = DEFAULT_WEIGHT;
	float tThres = T; 							//The t-homology threshold
	float dec = 0; 								//Intercept and decent to interpolate thresholds
	float inter = 0;
	string pFile, tFile, bLstFl = "highAbundKmersMiniK15w10Lrgr100BtStrnds.txt";  //Input file names
	string seq; 								//An input sequence
	ifstream fStr; 								//A file stream

	//A hash table to store black listed k-mers
	// unordered_map<uint64_t, char> bLstmers;

	mm_idxopt_t iopt; 							//An index option struct
	mm_mapopt_t mopt; 							//A mapping options struct
	mm_idx_reader_t *r; 						//An index reader
	const mm_idx_t *tidx; 						//A pointer to the text index
	const mm_idx_t *pidx; 						//A pointer to the pattern index
	vector<tuple<string, uint32_t, Sketch>> pSks; //A vector of pattern sketches
	vector<tuple<string, uint32_t, Sketch>>::const_iterator p; //An iterator to iterate over pattern sketches

	//Parse arguments
	if(!prsArgs(argc, argv, pFile, tFile, kmerLen, w, hFrac, bLstFl, comWght, uniWght, tThres, normalize, dec, inter, noNesting)){//TODO: Tests for this function need to be adapted!
		dsHlp(); //Display help message
		return 1;
	}

	//Set index options to default
	mm_set_opt(0, &iopt, &mopt);
	iopt.k = kmerLen; //Adjust k if necessary
	iopt.w = w;
	r = mm_idx_reader_open(tFile.c_str(), &iopt, INDEX_DEFAULT_DUMP_FILE);  //Open an index reader //TODO: We do not allow yet to use a prebuilt index

	//Check if index could be opened successfully
	if(r == NULL){
		cerr << "ERROR: Text sequence file could not be read" << endl;
		return -1;
	}

	bLstmers = readBlstKmers(bLstFl);

	//Construct index of reference
	if((tidx = mm_idx_reader_read(r, 1)) == 0){//TODO: Make use of multithreading here!
		cerr << "ERROR: Text index cannot be read" << endl;
		return -1;
	}

	//For simplicity we assume that an index always consists of only one part
	if(mm_idx_reader_read(r, 1) != 0){
		cerr << "ERROR: Text index consists of several parts! We cannot handle this yet" << endl;
		return -1; 
	}

	//Open index reader to read pattern
	r = mm_idx_reader_open(pFile.c_str(), &iopt, INDEX_DEFAULT_DUMP_FILE);

	//Check if index could be opened successfully
	if(r == NULL){
		cerr << "ERROR: Pattern sequence file could not be read" << endl;
		return -1;
	}

	//Construct index of reference
	if((pidx = mm_idx_reader_read(r, 1)) == 0){//TODO: Make use of multithreading here!
		cerr << "ERROR: Pattern index cannot be read" << endl;
		return -1;
	}

	//For simplicity we assume that an index always consists of only one part
	if(mm_idx_reader_read(r, 1) != 0){
		cerr << "ERROR: Pattern index consists of several parts! We cannot handle this yet" << endl;
		return -1; 
	}

	fStr.open(pFile);  //Open stream to read in patterns

    string text;
    //readFASTA(tFile, text);
	//Sketch tsk = buildMiniSketch(text, kmerLen, tidx->w, bLstmers);

    //cerr << "|T| = " << text.size() << ", |t| = " << tsk.size() << endl;

	//Load pattern sequences in batches
	// while(lPttnSks(fStr, kmerLen, hFrac, bLstmers, pSks) || !pSks.empty()){//TODO: Test for this function need to be adaptated!
	while(lMiniPttnSks(fStr, kmerLen, tidx->w, bLstmers, pSks) || !pSks.empty()){//TODO: This function still needs to be tested!
		//Iterate over pattern sketches
		for(p = pSks.begin(); p != pSks.end(); ++p){
			auto seqID = get<0>(*p);
            auto plen_nucl = get<1>(*p);
            auto sks = get<2>(*p);

			//Calculate an adapted threshold if we have the necessary informations
			if(dec != 0 && inter != 0) tThres = dec * plen_nucl + inter;

			//Find t-homologies and output them
            auto res = sweep(sks, tidx, plen_nucl, iopt.k);
			outputMappings(res, sks.size(), seqID, text);//TODO: Tests for this function need to be adaptated!
		}

		//Remove processed pattern sketches
		pSks.clear();
	}

	return 0;
}
