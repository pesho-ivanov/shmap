#include <queue>

#include "Sketch.cpp"
#include "IO.cpp"
#include "Thomology.cpp"
#include "Index.cpp"

//The FracMinHash ratio
double hFrac = HASH_RATIO;

//A hash table to store black listed k-mers
unordered_map<uint64_t, char> bLstmers;

using Cost = float;
using State = pair<uint32_t, uint32_t>;           // <text pos, pattern pos>
using StateWithCost = tuple<State, Cost, Cost>;   // <State,  g, f>

Cost jaccard_h(State s) {
    return 0.0;
}

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

const vector<Thomology> astar(const Sketch& skP, const mm_idx_t *tidx) {
    //Initialize occp. occp[kmer_hash] = #occurences in P
    unordered_map<uint64_t, uint32_t> occp = occurences(skP);

    //Hash position array of all hashes and their positions in the text sketch which also appear inside the pattern
    // <kmer_hash, pos_in_T> x |P|
    vector<pair<uint64_t, uint32_t>> L;
    L = genL(occp, tidx);

    // List of tripples <i, j, score> of matches
	vector<Thomology> res;

    std::priority_queue<StateWithCost> Q;

    for (int i=0; i<L.size(); i++) {
        State u(i, 0);
        Cost h = jaccard_h(u);
        Q.push(StateWithCost(u, 0.0, h));
    }

    while(!Q.empty()) {
        StateWithCost u = Q.top();
        Q.pop(); 

        st = get<0>(u);
        g = get<1>(u);
        //f = get<2>(u);

        {
            // move through text only (+1, +0)
            Cost edge_cost = 1.0; //skP[st.second];
            Cost next_g = g + edge_cost;
            Cost h = jaccard_h(next_st);
            State next_st(st.first+1, st.second);
            StateWithCost v(next_st, next_g, next_g + h);
            Q.push(v);
        }
        
        //{
        //    // move through pattern only (+0, +1)
        //    State next_st(get<0>(st)+1, get<1>(st));
        //}
        
        //{
        //    // match pattern and text (+1, +1)
        //}

    }

    //res.push_back(make_tuple(L[i].second, L[j].second, *rowIt));

    return res;
}

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

	//Load high abundance k-mers
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
	while(lMiniPttnSks(fStr, kmerLen, tidx->w, bLstmers, pSks) || !pSks.empty()){//TODO: This function still needs to be tested!
		//Iterate over pattern sketches
        
		for(p = pSks.begin(); p != pSks.end(); ++p){
			//Only output pattern sequence name if there is more than one sequence
			if(pSks.size() > 1) cout << get<0>(*p) << endl;
            
			//Calculate an adapted threshold if we have the necessary informations
			if(dec != 0 && inter != 0) tThres = dec * get<1>(*p) + inter;

            outputHoms(astar(get<2>(*p)));

			//Find t-homologies and output them
			//outputHoms(findThoms(get<2>(*p), tidx, comWght, uniWght, tThres, noNesting), normalize, get<2>(*p).size());//TODO: Tests for this function need to be adaptated!
		}

		//Remove processed pattern sketches
		pSks.clear();
	}

	return 0;
}
