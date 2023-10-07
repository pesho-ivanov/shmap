#include "Sketch.cpp"
#include "IO.cpp"
#include "Thomology.cpp"
#include "Index.cpp"

//The FracMinHash ratio
double hFrac = HASH_RATIO;
unordered_map<uint64_t, char> bLstmers;

unordered_map<uint64_t, uint32_t> occurences(const Sketch& skP) {
    unordered_map<uint64_t, uint32_t> occp(skP.size());

    //Fill occp
    for(Sketch::const_iterator fSkIt = skP.begin(); fSkIt != skP.end(); ++fSkIt){
        if(occp.contains(*fSkIt)){
            ++occp[*fSkIt];
        } else{
            occp[*fSkIt] = 1;
        }
    }
    
    return occp;
}

unordered_map<uint64_t, uint32_t> zero_occurences(const Sketch& skP) {
    unordered_map<uint64_t, uint32_t> occs(skP.size());

    //Fill occp
    for(Sketch::const_iterator fSkIt = skP.begin(); fSkIt != skP.end(); ++fSkIt)
        occs[*fSkIt] = 0;
    
    return occs;
}

template<typename TT>
const TT::iterator prev(const typename TT::iterator &it) {
    auto pr = it;
    return --pr;
}

typedef vector<pair<uint64_t, uint32_t>> L_t;

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

int32_t scoreX1000(const Sketch& skP, int xmin, int S_sz) {
    assert (S_sz >= 0);
    assert(skP.size() + S_sz - xmin > 0);
    auto scj = 1000 * xmin / (skP.size() + S_sz - xmin);
    //cout << scj << endl;
    assert (0 <= scj && scj <= 1000);
    return scj;
}

const int window2score(const L_t &L, const Sketch& skP, int from_nucl, int to_nucl, int plen_nucl) {
    int xmin = 0;
    auto l = L.end(), r = L.begin();
	for(auto s = L.begin(); s != L.end(); ++s) {
        if (from_nucl <= s->second && s->second + plen_nucl < to_nucl) {
            cout << "[" << from_nucl << ", " << to_nucl << "] includes sketch (" << s->first << ", " << s->second << ")" << endl; 
            if (s->second < l->second) l=s;
            if (s->second >= r->second) r=s;
            ++xmin;
        }
    }

    if (l == L.end() && r == L.begin())
        l = r;

    for (auto s=l; s<r; ++s)
        assert(from_nucl <= s->second && to_nucl <= s->second + plen_nucl);

    auto scj = scoreX1000(skP, xmin, r-l);
    //auto scj = 0;
    cout << "[" << from_nucl << ", " << to_nucl << "] includes " << xmin << " sketches,"  << ", r-l=" << (r-l) <<  ", scj=" << scj << endl;
    return xmin;
}

const vector<Thomology> sweep(const Sketch& skP, const mm_idx_t *tidx, const int Plen_nucl, const int k) {
    unordered_map<uint64_t, uint32_t> occp;  // occp[kmer_hash] = #occurences in P
    unordered_map<uint64_t, uint32_t> occs;  // occp[kmer_hash] = #occurences in s = T[l,r]
    vector<pair<uint64_t, uint32_t>> L;      // for all kmers from P in T: <kmer_hash, pos_in_T> * |P| sorted by pos_in_T
	vector<Thomology> res;                   // List of tripples <i, j, score> of matches

    occp = occurences(skP);
    L    = genL(occp, tidx);
    occs = zero_occurences(skP);

    //print_skP(skP);
    //print_L(L);

    int xmin = 0;
    Thomology best(-1, -1, 0);            // <i, j, score>
 
    //cerr << "|P| = " << Plen_nucl << ", |p| = " << skP.size() << endl;
    //cerr << "|L| = " << L.size() << endl;
    //cerr << "k = " << k <<endl;

    // Increase the left point end of the window [l,r) one by one.
    int i = 0, j = 0;
	for(auto l = L.begin(), r = L.begin(); l != L.end(); ++l, ++i) {
        // Increase the right end of the window [l,r) until it gets out.
        for(; r != L.end() && l->second + Plen_nucl >= r->second + k; ++r, ++j) {
            // If taking this kmer from T increases the intersection with P 
            if (++occs[r->first] <= occp[r->first])
                ++xmin;
            assert (l->second <= r->second);
        }

        auto scj = scoreX1000(skP, xmin, r-l);
        //cerr << "i=" << i << ", j=" << j << ", l=" << l->second << ", r=" << r->second << ", r-l=" << (prev(r)->second - l->second) << ", xmin=" << xmin << ", |s|=" << (r-l) <<  ", scj=" << scj << endl;

        // If better than best
        if (scj > get<2>(best)) {
            //cout << "s=r-l=" << r-l << endl;
            //cout << "denom: " << skP.size() + (r-l) - xmin << endl;
            best = Thomology(l->second, prev(r)->second, scj);
        }

        // Prepare for the next step 
        if (--occs[l->first] < occp[l->first])
            --xmin;

        assert(xmin >= 0);
    }

    for (auto it: occs)
        assert(it.second == 0);

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
