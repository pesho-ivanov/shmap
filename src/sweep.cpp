#include <algorithm>

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

template<typename TT>
const TT::iterator prev(const typename TT::iterator &it) {
    auto pr = it;
    return --pr;
}

typedef vector<pair<uint32_t, uint32_t>> L_t;

void print_L(const L_t &L) {
    int i=0;
    for (const auto &s: L) {
        cerr << "L[" << i << "]=(" << s.first << ", " << s.second << ")" << endl;
        ++i;
    }
}

void print_skP(const Sketch& skP) {
    int i=0;
    for (const auto &p: skP) {
        cerr << "skP[" << i << "]=" << p << endl;
        ++i;
    }
}

//const int window2score(const L_t &L, const Sketch& skP, int from_nucl, int to_nucl, int plen_nucl) {
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

void init(
		// input
		const Sketch& p,
		const mm_idx_t *tidx,
		// output
		unordered_map<uint32_t, int32_t> *hist,
		L_t *L) {

	unordered_map<uint64_t, uint32_t> hash2ord;
	uint32_t kmers = 0;
	L->reserve(P_MULTIPLICITY * p.size());

    for(auto p_it = p.begin(); p_it != p.end(); ++p_it) {
		uint32_t kmer_ord;
		auto kmer_hash = *p_it;
		auto ord_it = hash2ord.find(kmer_hash);
		if (ord_it != hash2ord.end()) {
			kmer_ord = ord_it->second;
			++(hist->at(kmer_ord));
		} else {
			kmer_ord = kmers;
			assert(hist->find(kmer_ord) == hist->end());
			(*hist)[kmer_ord] = 1;
			hash2ord[kmer_hash] = kmers;
			++kmers;

			int nHits;
			auto *idx_p = mm_idx_get(tidx, kmer_hash, &nHits);
			for (int i = 0; i < nHits; ++i, ++idx_p)// Iterate over all occurrences
				L->push_back(make_pair(kmer_ord, ((uint32_t)(*idx_p))>>1));  // Push (hash value, k-mer position in sketch) pair
		}
	}

	//Sort L by ascending positions in reference
	sort(L->begin(), L->end(), [](const pair<uint32_t, uint32_t>& hpa, const pair<uint32_t, uint32_t>& hpb) {
		return hpa.second < hpb.second;
	});

	//for (auto it: *hash2ord)
	//	cout << "hash2ord: " << it.first << " " << it.second << endl;
	//for (auto it: *L)
	//	cout << "L: " << it.first << " " << it.second << endl;
	//for (auto it: *hist)
	//	cout << "hist: " << it.first << " " << it.second << endl;
}

int32_t scoreX1000(const Sketch& p, int s_sz, int xmin) {
    assert (s_sz >= 0);
    assert(p.size() + s_sz - xmin > 0);
    auto scj = 1000 * xmin / (p.size() + s_sz - xmin);
    //cout << scj << endl;
    assert (0 <= scj && scj <= 1000);
    return scj;
}

const vector<Thomology> sweep(const Sketch& p, const mm_idx_t *tidx, const int Plen_nucl, const int k) {
    unordered_map<uint32_t, int32_t> hist;  // rem[kmer_hash] = #occurences in `p` - #occurences in `s`
    vector<pair<uint32_t, uint32_t>> s;    // for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
	vector<Thomology> res;                 // List of tripples <i, j, score> of matches

    init(p, tidx, &hist, &s);

    int xmin = 0;
    Thomology best(0, 0, 0);            // <i, j, score>
 
    //cerr << "|P| = " << Plen_nucl << ", |p| = " << skP.size() << endl;
    //cerr << "|L| = " << L.size() << endl;
    //cerr << "k = " << k <<endl;

    // Increase the left point end of the window [l,r) one by one.
    int i = 0, j = 0;
	for(auto l = s.begin(), r = s.begin(); l != s.end(); ++l, ++i) {
        // Increase the right end of the window [l,r) until it gets out.
        for(; r != s.end() && r->second + k <= l->second + Plen_nucl; ++r, ++j) {
			// If taking this kmer from T increases the intersection with P 
			if (--hist.at(r->first) >= 0)
				++xmin;
			assert (l->second <= r->second);
        }

		int s_sz = r-l;  // TODO: fix! S is from T, not L
        auto scj = scoreX1000(p, s_sz, xmin);
        //cout << "i=" << i << ", j=" << j << ", l=" << l->second << ", r=" << r->second << ", |s|=" << s_sz << ", xmin=" << xmin << ", |s|=" << (r-l) <<  ", scj=" << scj << endl;

        // If better than best
        if (scj > get<2>(best)) {
            //cout << "s=r-l=" << r-l << endl;
            //cout << "denom: " << skP.size() + (r-l) - xmin << endl;
            best = Thomology(l->second, prev(r)->second, scj);
        }

        // Prepare for the next step by moving `l` to the right
		if (++hist.at(l->first) > 0)
			--xmin;

        assert(xmin >= 0);
    }
	assert (xmin == 0);

    //auto scj = window2score(L, skP, 4681377, 4683857, plen_nucl);
    //cout << "l=" << 4681377 << ", r=" << 4683857 << ", xmin=" << xmin << ", r-l=" << 4683857-4681377 <<  ", scj=" << scj << endl;

    res.push_back(best);
    return res;
}

// Score for the current window
//double scj = 1.0*xmin / (2*len - xmin);
//cout << "Curr: " << l->second << " " << r->second << " " << xmin << endl;
//cout << "Best: " << get<0>(best) << " " << get<1>(best) << " " << get<2>(best) << endl;

int main(int argc, char **argv){
	//Flag to save that scores are to be normalized
	bool normalize = NORM_FLAG_DEFAULT;
	//Flag to state if we are interested in nested results
	bool noNesting = NESTING_FLAG_DEFAULT;
	//The k-mer length
	uint32_t kmerLen = K;
	//The window size
	uint32_t w = W;
	//Scoring weights
	uint32_t comWght = DEFAULT_WEIGHT;
	float uniWght = DEFAULT_WEIGHT;
	//The t-homology threshold
	float tThres = T;
	//Intercept and decent to interpolate thresholds
	float dec = 0;
	float inter = 0;
	//Input file names
	string pFile, tFile, bLstFl = "highAbundKmersMiniK15w10Lrgr100BtStrnds.txt";
	//An input sequence
	string seq;
	//A file stream
	ifstream fStr;

	//A hash table to store black listed k-mers
	// unordered_map<uint64_t, char> bLstmers;

	//An index option struct
	mm_idxopt_t iopt;
	//A mapping options struct
	mm_mapopt_t mopt;
	//An index reader
	mm_idx_reader_t *r;
	//A pointer to the text index
	const mm_idx_t *tidx;
	//A pointer to the pattern index
	const mm_idx_t *pidx;
	//A vector of pattern sketches
	vector<tuple<string, uint32_t, Sketch>> pSks;
	//An iterator to iterate over pattern sketches
	vector<tuple<string, uint32_t, Sketch>>::const_iterator p;

	//Parse arguments
	if(!prsArgs(argc, argv, pFile, tFile, kmerLen, w, hFrac, bLstFl, comWght, uniWght, tThres, normalize, dec, inter, noNesting)){//TODO: Tests for this function need to be adapted!
		//Display help message
		dsHlp();
		return 1;
	}

	//Set index options to default
	mm_set_opt(0, &iopt, &mopt);
	//Adjust k if necessary
	iopt.k = kmerLen;
	iopt.w = w;
	//Open an index reader //TODO: We do not allow yet to use a prebuilt index
	r = mm_idx_reader_open(tFile.c_str(), &iopt, INDEX_DEFAULT_DUMP_FILE);

	//Check if index could be opened successfully
	if(r == NULL){
		cerr << "ERROR: Text sequence file could not be read" << endl;
		return -1;
	}

	bLstmers = readBlstKmers(bLstFl);

	// TODO: ignore blacklisted kmers in the index

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

	//Open stream to read in patterns
	fStr.open(pFile);

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
			outputHoms(res, normalize, sks.size(), seqID, text);//TODO: Tests for this function need to be adaptated!
			//outputHomsSeqs(res, normalize, sks.size());//TODO: Tests for this function need to be adaptated!
		}

		//Remove processed pattern sketches
		pSks.clear();
	}

	return 0;
}
