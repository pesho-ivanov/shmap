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
    unordered_map<uint64_t, uint32_t> occt(skP.size());

    //Fill occp
    for(Sketch::const_iterator fSkIt = skP.begin(); fSkIt != skP.end(); ++fSkIt)
        occt[*fSkIt] = 0;
    
    return occt;
}

const vector<Thomology> sweep(const Sketch& skP, const mm_idx_t *tidx, const int len) {
    unordered_map<uint64_t, uint32_t> occp;  // occp[kmer_hash] = #occurences in P
    unordered_map<uint64_t, uint32_t> occt;  // ---||--- in T
    vector<pair<uint64_t, uint32_t>> L;      // for all kmers from P in T: <kmer_hash, pos_in_T> * |P|
    unordered_map<uint64_t, char> cnt;       // cnt[kmer_hash] -> #kmer_hash in T[curr, curr+len)
	vector<Thomology> res;                   // List of tripples <i, j, score> of matches

    occp = occurences(skP);
    occt = zero_occurences(skP);
    L    = genL(occp, tidx);
    int xmin = 0;
    Thomology best(-1, -1, -1.0);            // <i, j, score>
 
	for(auto l = L.begin(), r = L.begin(); l != L.end(); ++l) {
        // Move the right end to the right.
        for(; r != L.end() && r->second <= l->second + len; ++r) 
            if (++occt[r->first] <= occp[r->first])
                ++xmin;

        // If better than best
        if (xmin > get<2>(best))
            best = Thomology(l->second, r->second, xmin);

        // Prepare for the next step 
        if (--occt[l->first] < occp[l->first])
            --xmin;

        assert(xmin >= 0);
    }

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

	//Load pattern sequences in batches
	// while(lPttnSks(fStr, kmerLen, hFrac, bLstmers, pSks) || !pSks.empty()){//TODO: Test for this function need to be adaptated!
	while(lMiniPttnSks(fStr, kmerLen, tidx->w, bLstmers, pSks) || !pSks.empty()){//TODO: This function still needs to be tested!
		//Iterate over pattern sketches
		for(p = pSks.begin(); p != pSks.end(); ++p){
			auto seqID = get<0>(*p);
            auto len = get<1>(*p);
            auto sks = get<2>(*p);
			//Only output pattern sequence name if there is more than one sequence
			if(pSks.size() > 1) cout << seqID << endl;

			//Calculate an adapted threshold if we have the necessary informations
			if(dec != 0 && inter != 0) tThres = dec * len + inter;

			//Find t-homologies and output them
			outputHoms(sweep(sks, tidx, len), normalize, sks.size());//TODO: Tests for this function need to be adaptated!
		}

		//Remove processed pattern sketches
		pSks.clear();
	}

	return 0;
}
