#ifndef IO_HPP
#define IO_HPP

#include <getopt.h>
#include <fstream>
#include <map>

#define INDEX_DEFAULT_DUMP_FILE "indexDump.idx"
//#define T 0.9
#define DEFAULT_WEIGHT 1
using Thomology = tuple<uint32_t, uint32_t, int32_t>;

#define T_HOM_OPTIONS "p:s:k:w:e:a:b:S:M:t:d:i:z:nxoh"

//#define MIN_PARAM_NB 6
#define MAX_RATIO 1.0
//#define NORM_FLAG_DEFAULT false
#define PATTERN_BATCH_SIZE 250000
#define STRING_BUFFER_SIZE_DEFAULT 50

enum class elastic_t {
    wrong = -1,
    off = 0,
    consecutive = 2,
    random = 9
}; 

elastic_t str2elastic(string t) {
    if (t == "off") return elastic_t::off;
    if (t == "consecutive") return elastic_t::consecutive;
    if (t == "random") return elastic_t::random;
	return elastic_t::wrong;
}

std::string getElasticDescription(elastic_t e) {
    switch (e) {
        case elastic_t::off:
            return "off";
        case elastic_t::consecutive:
            return "consecutive";
        case elastic_t::random:
            return "random";
        default:	
            return "wrong";
    }
}

enum class alignment_edges_t {
    wrong = -1,
    sketch_edges = 0,
    extend_equally = 1,
    fine = 2
}; 

alignment_edges_t str2alignment_edges(string t) {
    if (t == "sketch_edges") return alignment_edges_t::sketch_edges;
    if (t == "extend_equally") return alignment_edges_t::extend_equally;
    if (t == "fine") return alignment_edges_t::fine;
	return alignment_edges_t::wrong;
}

std::string getAlignmentEdgesDescription(alignment_edges_t ae) {
    switch (ae) {
        case alignment_edges_t::sketch_edges:
            //return "Use the beginning of the first matched sketch and the end of the last matched sketch.";
            return "sketch_edges";
        case alignment_edges_t::extend_equally:
            //return "Make the length equal to |P| by extending the sketch edges equally.";
            return "extend_equally";
        case alignment_edges_t::fine:
            return "fine";
        default:	
            return "wrong";
    }
}

struct params_t {
	bool normalize; 		//Flag to save that scores are to be normalized
	bool onlybest;          // Output up to one (best) mapping (if above the threshold)
	bool overlaps;          // Permit overlapping mappings 
	uint32_t k;  						//The k-mer length
	uint32_t w;  							//The window size
	elastic_t elastic;
	alignment_edges_t alignment_edges;
	//double hFrac;  					//The FracMinHash ratio
	uint32_t max_seeds;  							//Maximum seeds in a sketch
	uint32_t max_matches;  							//Maximum seed matches in a sketch
	float tThres;  							//The t-homology threshold
	float dec; 								//Intercept and decent to interpolate thresholds
	float inter; 
	string pFile, tFile;
	string paramsFile;
	string bLstFl; // = "highAbundKmersMiniK15w10Lrgr100BtStrnds.txt";  //Input file names

	params_t() {
		normalize = false; 		//Flag to save that scores are to be normalized
		onlybest = false;
		overlaps = false;          // 
		k = K; 						//The k-mer length
		w = W; 							//The window size
		elastic = elastic_t::off;
		alignment_edges = alignment_edges_t::fine;
		//hFrac = HASH_RATIO;
		max_seeds = 10000;
		max_matches = 1000000;
		tThres = 0.9; 							//The t-homology threshold
		dec = 0; 								//Intercept and decent to interpolate thresholds
		inter = 0;
		bLstFl = ""; //"highAbundKmersMiniK15w10Lrgr100BtStrnds.txt";  //Input file names
	}

	void print(std::ostream& out = std::cout, bool human = true) {
		vector<pair<string, string>> m;
		m.push_back({"normalize", std::to_string(normalize)});
		m.push_back({"onlybest", std::to_string(onlybest)});
		m.push_back({"overlaps", std::to_string(overlaps)});
		m.push_back({"k", std::to_string(k)});
		m.push_back({"w", std::to_string(w)});
		m.push_back({"elastic", getElasticDescription(elastic)});
		m.push_back({"alignment_edges", getAlignmentEdgesDescription(alignment_edges)});
		//m.push_back({"hFrac", std::to_string(hFrac)});
		m.push_back({"tThres", std::to_string(tThres)});
		m.push_back({"dec", std::to_string(dec)});
		m.push_back({"inter", std::to_string(inter)});
		m.push_back({"pFile", pFile});
		m.push_back({"tFile", tFile});
		m.push_back({"paramsFile", paramsFile});
		m.push_back({"bLstFl", bLstFl});

		if (human) {
			out << "Parameters:" << endl;
			for (auto& p : m)
				out << std::setw(20) << std::right << p.first << ": " << p.second << endl;
		} else {
			for (auto& p : m)
				out << p.first << "\t";
			out << endl;
			for (auto& p : m)
				out << p.second << "\t";
		}
	}
};

struct reader_t {
	params_t params;
	ifstream fStr; 								//A file stream

	mm_idxopt_t iopt; 							//An index option struct
	mm_mapopt_t mopt; 							//A mapping options struct
	const mm_idx_t *tidx; 						//A pointer to the text index
	//const mm_idx_t *pidx; 						//A pointer to the pattern index
	vector<tuple<string, uint32_t, Sketch>> pSks; //A vector of pattern sketches

	uint32_t T_sz;
	string text;

	reader_t() {}
	bool init_and_index(int argc, char **argv);
};

//This function prints usage infos
inline void dsHlp() {
	cerr << "sweep [-hn] [-p PATTERN_FILE] [-s TEXT_FILE] [-k KMER_LEN] [-r HASH_RATIO] [-b BLACKLIST] [-c COM_HASH_WGHT] [-u UNI\
	_HASH_WGHT] [-t HOM_THRES] [-d DECENT] [-i INTERCEPT]" << endl << endl;
	cerr << "Find sketch-based pattern similarity in text." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -p   --pattern  Pattern sequences file (FASTA format)" << endl;
	cerr << "   -s   --text     Text sequence file (FASTA format)" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -k   --ksize             K-mer length to be used for sketches (default " << K << ")" << endl;
	cerr << "   -w   --windowsize        Window size for minimizer sketching approach (default " << W << ")" << endl;
	cerr << "   -e   --elastic           Elastic pairs of kmers {off, consecutive, random} [off]" << endl;
	cerr << "   -a   --alignment_edges   Alignment interval {sketch_edges, extend_equally, fine} [fine]" << endl;
	//cerr << "   -r   --hashratio         FracMin hash ratio to be used for sketches (default " << HASH_RATIO << ")" << endl;
	cerr << "   -b   --blacklist         File containing hashes to ignore for sketch calculation" << endl;
	cerr << "   -S   --max_seeds         Max seeds in a sketch" << endl;
	cerr << "   -M   --max_matches       Max seed matches in a sketch" << endl;
	cerr << "   -t   --hom_thres         Homology threshold (default " << 0.9 << ")" << endl;
	cerr << "   -d   --decent            Decent required for dynamic threshold selection" << endl;
	cerr << "   -i   --intercept         Intercept required for dynamic threshold selection" << endl;
	cerr << "   -z   --params     		 Output file with parameters (tsv)" << endl << endl;
	cerr << "   -o   --overlaps     		 Permit overlapping mappings" << endl << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -n   --normalize  Normalize scores by length" << endl;
	cerr << "   -x   --onlybest   Output the best alignment if above threshold (otherwise none)" << endl;
	cerr << "   -h   --help       Display this help message" << endl;
}

//This function parses the program parameters. Returns false if given arguments are not valid
bool prsArgs(int& nArgs, char** argList, params_t *params) {
	int option_index = 0, a;

	static struct option long_options[] = {
        {"pattern",            required_argument,  0, 'p'},
        {"text",               required_argument,  0, 's'},
        {"ksize",              required_argument,  0, 'k'},
        {"windowsize",         required_argument,  0, 'w'},
        {"elastic",            required_argument,  0, 'e'},
        {"alignment_edges",    required_argument,  0, 'a'},
        {"hashratio",          required_argument,  0, 'r'},
        {"blacklist",          required_argument,  0, 'b'},
        {"max_seeds",          required_argument,  0, 'S'},
        {"max_matches",        required_argument,  0, 'M'},
        //{"commonhashweight",   required_argument,  0, 'c'},
        //{"uniquehashweight",   required_argument,  0, 'u'},
        {"hom_thres",          required_argument,  0, 't'},
        {"decent",             required_argument,  0, 'd'},
        {"intercept",          required_argument,  0, 'i'},
        {"params",             required_argument,  0, 'z'},
        {"normalize",          no_argument,        0, 'n'},
        {"onlybest",           no_argument,        0, 'x'},
        {"overlaps",           no_argument,        0, 'o'},
        {"help",               no_argument,        0, 'h'},
        {0,                    0,                  0,  0 }
    };

    //Parse all parameters given
	while ((a = getopt_long(nArgs, argList, T_HOM_OPTIONS, long_options, &option_index)) != -1) {
		//Assign parameter values
		switch(a) {
			case 'p':
				//Save input sequence
				params->pFile = optarg;
				break;
			case 's':
				//Save input sequence
				params->tFile = optarg;
				break;
			case 'z':
				//Save input sequence
				params->paramsFile = optarg;
				break;
			case 'k':
				//A k-mer length should be positive
				if(atoi(optarg) <= 0) {
					cerr << "ERROR: K-mer length not applicable" << endl;
					return false;
				}

				params->k = atoi(optarg);
				break;
			case 'S':
				if(atoi(optarg) <= 0) {
					cerr << "ERROR: The number of seeds should be positive." << endl;
					return false;
				}
				params->max_seeds = atoi(optarg);
				break;
			case 'M':
				if(atoi(optarg) <= 0) {
					cerr << "ERROR: The number of seed matches should be positive." << endl;
					return false;
				}
				params->max_matches = atoi(optarg);
				break;
			case 'e':
				//Elastic kmers
				params->elastic = str2elastic(optarg);
				if(params->elastic == elastic_t::wrong) {
					cerr << "ERROR: Elastic parameter not from the list." << endl;
					return false;
				}
				break;
			case 'a':
			 	// alignment edges
				params->alignment_edges = str2alignment_edges(optarg);
				if(params->alignment_edges == alignment_edges_t::wrong) {
					cerr << "ERROR: Alignment edges not from the list." << endl;
					return false;
				}
				break;
			case 'w':
				//Window size should be positive
				if(atoi(optarg) <= 0) {
					cerr << "ERROR: Window size not applicable" << endl;
					return false;
				}

				params->w = atoi(optarg);
				break;
//			case 'r':
//				//Check if given value is reasonable to represent a ratio
//				if(atof(optarg) <= 0 || atof(optarg) > MAX_RATIO) {
//					cerr << "ERROR: Given hash ratio not applicable" << endl;
//
//					return false;
//				}
//
//				params->hFrac = atof(optarg);
//				break;
			case 'b':
				//Save blacklist file name
				params->bLstFl = optarg;
				break;
			//case 'c':
			//	//Weights should be positive
			//	if(atoi(optarg) <= 0) {
			//		cerr << "ERROR: Common hash weight not applicable" << endl;

			//		return false;
			//	}

			//	reader->cw = atoi(optarg);
			//	break;
			//case 'u':
			//	//Weights should be positive
			//	if(atof(optarg) <= 0) {
			//		cerr << "ERROR: Unique hash weight not applicable" << endl;

			//		return false;
			//	}

			//	reader->uw = atof(optarg);
			//	break;
			case 't':
				params->tThres = atof(optarg);
				break;
			case 'd':
				params->dec = atof(optarg);
				break;
			case 'i':
				params->inter = atof(optarg);
				break;
			case 'n':
				params->normalize = true;
				break;
			case 'x':
				params->onlybest = true;
				break;
			case 'o':
				params->overlaps = true;
				break;
			case 'h':
				return false;
			default:
				cerr << "Unknown option " << a << " '" << char(a) << "', " << optarg << endl ;
				break;
		}
	}

	return !params->pFile.empty() && !params->tFile.empty();
}

//This function reads a file in FASTA format and returns true on success
//NOTE: Using this function, we first read in a complete sequence before calculating a sketch from it. Calculating a sketch could
//		also be done on the fly while reading the file. This possibility would be more memory saving most likely, because we would
//		not have to store the whole sequence in memory first. However, we had no chance to estimate the sketch size. Thus, the vec-
//		tor to store the sketch would have to be elongated several times which might be time consuming. On the other hand, it might
//		also be time consuming to handle the complete sequence in memory first. Eventually, it will depend on a try to find out
//		what is the better alternative. Could also be that it does not matter at all!
bool readFASTA(const string& filePath, string& seq) {
	bool headerRead = false, lnBrkDiscvd = false;
	char c;

	//Open the file
	ifstream fStr(filePath);

	//Check if the file is open
	if(!fStr.is_open()) return false;

	//Read in file character by character
	while(fStr.get(c)) {
		//We are done if we find a second header (which can only start after at least one line break) in the file
		if(c == '>' && headerRead && lnBrkDiscvd) break;

		//Header lines are skipped
		if(c == '>') {
			headerRead = true;
			continue;
		}

		//The first line break indicates that we have left the header line
		if(c == '\n') {
			lnBrkDiscvd = true;
			continue;
		}

		//There is no sequence to load in the header line
		if(headerRead && !lnBrkDiscvd) continue;

		//We are only interested in unambigous, unmasked nucleotides
		if(c == 'A' || c == 'C' || c == 'G' || c == 'T') seq += c;
	}

	//Close file
	fStr.close();

	return true;
}

//This function reads in batches of FASTA sequence entries from file and transforms them into minimap sketches. Returns false if end
// of file was reached.
bool lMiniPttnSks(ifstream& fStr, const uint32_t& k, const uint32_t& w, const unordered_map<uint64_t, char>& blmers, 
	vector<tuple<string, uint32_t, Sketch>>& pSks) {
	bool headerRead, idRead = false, lnBrkDiscvd = false;
	char c;
	string seq, seqID;

	//Check if the file is open
	if(!fStr.is_open()) return false;

	//If this function is called iteratively, the '>' of the next entry has already been read
	headerRead = fStr.gcount() != 0;

	//Read in file character by character
	while(fStr.get(c)) {
		//An entry's sequence is complete if we find a second header (which can only start after at least one line break) in the 
		//file
		if(c == '>' && headerRead && lnBrkDiscvd) {
			//Add sequence's sketch, length and id to result vector
			//auto sketch = buildMiniSketch(seq, k, w, blmers);
			mm128_v a = {0,0,0};
			mm_sketch(0, seq.c_str(), seq.size(), w, k, 0, 0, &a);
			auto sketch2 = convertToSketch(&a);
			//if (sketch.size() != sketch2.size()) {
			//	cerr << "ERROR: sketch size mismatch" << endl;
			//	cerr << "sketch:  " << sketch.size() << endl;
			//	cerr << "sketch2: " << sketch2.size() << endl;
			//	//exit(1);
			//}
			//for (size_t i = 1; i < sketch.size(); ++i) {
			//	if (sketch[i] != sketch2[i]) {
			//		cerr << "ERROR: sketch mismatch" << endl;
			//		cerr << "sketch:  " << sketch[i].first << ", " << sketch[i].second << endl;
			//		cerr << "sketch2: " << sketch2[i].first << ", " << sketch2[i].second << endl;
			//		//exit(1);
			//	}
			//}	
			pSks.push_back(make_tuple(seqID, seq.length(), sketch2));//TODO: This function still needs to be tested!
            //cout << "pattern: " << seq << endl;
			//Clear sequence id
			seqID.clear();
			//Clear sequence
			seq.clear();

			//Check if enough sequences have been read
			if(pSks.size() == PATTERN_BATCH_SIZE) return true;

			//Reset id-read flag
			idRead = false;
			//Reset line-break-discovered flag
			lnBrkDiscvd = false;
			continue;
		}

		//Note if we have completely read the sequence id
		idRead = idRead || (headerRead && c == ' ' && !lnBrkDiscvd);
		//Note if we have found the first line break after a new header started
		lnBrkDiscvd = lnBrkDiscvd || c == '\n';

		//Update sequence id if we are still reading it
		if(headerRead && !lnBrkDiscvd && !idRead) {
			seqID += c;
			continue;
		}

		//Note if we have found the beginning of a header
		headerRead = headerRead || (c == '>');

		//There is no sequence to load in the header line
		if(headerRead && !lnBrkDiscvd) continue;

		//We are only interested in unambigous, unmasked nucleotides
		if(c == 'A' || c == 'C' || c == 'G' || c == 'T') seq += c;
	}

	//Add last entry's sketch and sequence id to result vector if it is not empty
	if(!seq.empty()) {
        pSks.push_back(make_tuple(seqID, seq.length(), buildMiniSketch(seq, k, w, blmers)));
        //cout << "pattern: " << seq << endl;
    }

	return false;
}

//This function reads 64-bit numbers from file and returns them as a hash table
const unordered_map<uint64_t, char> readBlstKmers(const string& fname) {
	char nb[STRING_BUFFER_SIZE_DEFAULT];
	char* end;
	ifstream iFile(fname);
	unordered_map<uint64_t, char> numbers;

	if(!iFile.is_open()) {
		cerr << "readBlstKmers: WARNING input file could not be opened. Hash table will be empty!" << endl;
		return numbers;
	}

	//Read kmer hashes
	while(iFile.getline(nb, STRING_BUFFER_SIZE_DEFAULT))
		numbers[strtoull(nb, &end, 0)] = 1;

	return numbers;
}
	
bool reader_t::init_and_index(int argc, char **argv) {
	//A hash table to store black listed k-mers
	// unordered_map<uint64_t, char> bLstmers;
	mm_idx_reader_t *r; 						//An index reader

	cerr << "Reading index " << params.tFile << "..." << endl;

	//Parse arguments
	if(!prsArgs(argc, argv, &params)) {//TODO: Tests for this function need to be adapted!
		dsHlp(); //Display help message
		return 0;
	}

	//Set index options to default
	mm_set_opt(0, &iopt, &mopt);
	iopt.k = params.k; //Adjust k if necessary
	iopt.w = params.w;
	r = mm_idx_reader_open(params.tFile.c_str(), &iopt, INDEX_DEFAULT_DUMP_FILE);  //Open an index reader //TODO: We do not allow yet to use a prebuilt index

	//Check if index could be opened successfully
	if(r == NULL) {
		cerr << "ERROR: Text sequence file could not be read" << endl;
		return 1;
	}

	if (!params.bLstFl.empty()) {
		bLstmers = readBlstKmers(params.bLstFl);
	}

	//Construct index of reference
	if((tidx = mm_idx_reader_read(r, 1)) == 0) { //TODO: Make use of multithreading here!
		cerr << "ERROR: Text index cannot be read" << endl;
		return 1;
	}

	//For simplicity we assume that an index always consists of only one part
	if(mm_idx_reader_read(r, 1) != 0) {
		cerr << "ERROR: Text index consists of several parts! We cannot handle this yet" << endl;
		return 1; 
	}

	//Open index reader to read pattern
	//r = mm_idx_reader_open(params.pFile.c_str(), &iopt, INDEX_DEFAULT_DUMP_FILE);

	////Check if index could be opened successfully
	//if(r == NULL) {
	//	cerr << "ERROR: Pattern sequence file could not be read" << endl;
	//	return 1;
	//}

	////Construct index of reference
	//if((pidx = mm_idx_reader_read(r, 1)) == 0) {//TODO: Make use of multithreading here!
	//	cerr << "ERROR: Pattern index cannot be read" << endl;
	//	return 1;
	//}

	////For simplicity we assume that an index always consists of only one part
	//if(mm_idx_reader_read(r, 1) != 0) {
	//	cerr << "ERROR: Pattern index consists of several parts! We cannot handle this yet" << endl;
	//	return 1; 
	//}

	T_sz = tidx->seq->len;  //The length of the text
	fStr.open(params.pFile);  //Open stream to read in patterns

	//readFASTA(tFile, text);
	//Sketch tsk = buildMiniSketch(text, kmerLen, tidx->w, bLstmers);
	//cerr << "|T| = " << text.size() << ", |t| = " << tsk.size() << endl;
	return 0;
}

#endif