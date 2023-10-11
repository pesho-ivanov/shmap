#ifndef IO_HPP
#define IO_HPP

#include <getopt.h>
#include <fstream>

#include "Measures.h"

#define T 0
#define DEFAULT_WEIGHT 1
using Thomology = tuple<uint32_t, uint32_t, int32_t>;

#define OPTIONS "a:b:ilh"
#define T_HOM_OPTIONS "p:s:k:w:r:b:c:u:t:d:i:nNh"
#define MIN_PARAM_NB 6
#define MAX_RATIO 1.0
#define NORM_FLAG_DEFAULT false
#define NESTING_FLAG_DEFAULT true
#define PATTERN_BATCH_SIZE 250000
#define STRING_BUFFER_SIZE_DEFAULT 50

//This function prints usage infos
inline void dspHlp(){
	cerr << "CalcSim [-hil] [-a SEQ1] [-b SEQ2]" << endl << endl;
	cerr << "Calculating sequence similarity based on their sketches." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -a   --seqa  First input sequence" << endl;
	cerr << "   -b   --seqb  Second input sequence" << endl << endl;
	cerr << "Required parameters without argument:" << endl;
	cerr << "   -i   --intersim     Calculate similarity based on intersection measure" << endl;
	cerr << "   -l   --algnHshsSim  Calculate similarity based on aligning hashes measure" << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -h   --help      Display this help message" << endl;
}

//This function prints usage infos
inline void dsHlp(){
	cerr << "eskemap [-hn] [-p PATTERN_FILE] [-s TEXT_FILE] [-k KMER_LEN] [-r HASH_RATIO] [-b BLACKLIST] [-c COM_HASH_WGHT] [-u UNI\
	_HASH_WGHT] [-t HOM_THRES] [-d DECENT] [-i INTERCEPT] [-N NESTING]" << endl << endl;
	cerr << "Find sketch-based pattern similarity in text." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -p   --pattern  Pattern sequences file (FASTA format)" << endl;
	cerr << "   -s   --text     Text sequence file (FASTA format)" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -k   --ksize             K-mer length to be used for sketches (default " << K << ")" << endl;
	cerr << "   -w   --windowsize        Window size for minimizer sketching approach (default " << W << ")" << endl;
	cerr << "   -r   --hashratio         FracMin hash ratio to be used for sketches (default " << HASH_RATIO << ")" << endl;
	cerr << "   -b   --blacklist         File containing hashes to ignore for sketch calculation" << endl;
	cerr << "   -c   --commonhashweight  Weight to reward common hashes (default " << DEFAULT_WEIGHT << ")" << endl;
	cerr << "   -u   --uniquehashweight  Weight to punish unique hashes (default " << DEFAULT_WEIGHT << ")" << endl;
	cerr << "   -t   --hom_thres         Homology threshold (default " << T << ")" << endl;
	cerr << "   -d   --decent            Decent required for dynamic threshold selection" << endl;
	cerr << "   -i   --intercept         Intercept required for dynamic threshold selection" << endl;
	cerr << "Optional parameters without argument:" << endl;
	cerr << "   -n   --normalize  Normalize scores by length" << endl;
	cerr << "   -N   --nesting    Output nested homologies" << endl;
	cerr << "   -h   --help       Display this help message" << endl;
}

//inline void outputHomsWithNucl(const vector<Thomology>& homs, const uint32_t& pLen){
//	//Iterate over t-homologies
//	for(vector<Thomology>::const_iterator h = homs.begin(); h != homs.end(); ++h){
//		//Normalize score before reporting if requested
//        cout << "i: " << get<0>(*h) << " j: " << get<1>(*h) << " score: " << (double) get<2>(*h) / max(pLen, get<1>(*h) - 
//            get<0>(*h) + 1) << endl;
//	}
//}

////This function parses the program parameters. Returns false if given arguments are not valid
//const bool prsArgs(int& nArgs, char** argList, string& seqa, string& seqb, Measure& msr);
//
////This function parses the program parameters. Returns false if given arguments are not valid
//const bool prsArgs(int& nArgs, char** argList, string& pFl, string& tFl, uint32_t& k, uint32_t& w, double& hFrac, string& blFl, 
//	uint32_t& cw, float& uw, float& tThres, bool& norm, float& dec, float& inter, bool& noNesting);
//
////This function reads a file in FASTA format and returns true on success
//const bool readFASTA(const string& filePath, string& seq);
//
////This function reads in batches of FASTA sequence entries from file and transforms them into sketches. Returns false if end of file
////was reached.
//const bool lPttnSks(ifstream& fStr, const uint32_t& k, const double& hFrac, const unordered_map<uint64_t, char>& bLstmers, 
//	vector<tuple<string, uint32_t, Sketch>>& pSks);
//
////This function reads in batches of FASTA sequence entries from file and transforms them into minimap sketches. Returns false if end
//// of file was reached.
//const bool lMiniPttnSks(ifstream& fStr, const uint32_t& k, const uint32_t& w, const unordered_map<uint64_t, char>& blmers, 
//	vector<tuple<string, uint32_t, Sketch>>& pSks);
//
////This function reads 64-bit numbers from file and returns them as a hash table
//const unordered_map<uint64_t, char> readBlstKmers(const string& fname);

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& seqa, string& seqb, Measure& msr){
	bool seqaGvn = false, seqbGvn = false;
	int option_index = 0, a;

	//Check if enough arguments are given at all
	if(nArgs < MIN_PARAM_NB) return false;

	static struct option long_options[] = {
        {"seqa",        required_argument,  0, 'a'},
        {"seqb",        required_argument,  0, 'b'},
        {"intersim",    no_argument,        0, 'i'},
        {"algnHshsSim", no_argument,        0, 'l'},
        {"help",        no_argument,        0, 'h'},
        {0,             0,                  0,  0 }
    };

    //Parse all parameters given
	while ((a = getopt_long(nArgs, argList, OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'a':
				//Save input sequence
				seqa = optarg;
				//Note that we have seen one input sequence
				seqaGvn = true;
				break;
			case 'b':
				//Save input sequence
				seqb = optarg;
				//Note that we have seen one input sequence
				seqbGvn = true;
				break;
			case 'i':
				//Check if another measure was also given as a parameter
				if(msr != none && msr != intersec){
					cerr << "ERROR: Only one measure can be calculated at a time!" << endl;

					return false;
				}

				msr = intersec;
				break;
			case 'l':
				//Check if another measure was also given as a parameter
				if(msr != none && msr != algnWoutOffs){
					cerr << "ERROR: Only one measure can be calculated at a time!" << endl;

					return false;
				}

				msr = algnWoutOffs;
				break;
			case 'h':
				return false;
			default:
				break;
		}
	}

	return seqaGvn && seqbGvn && msr != none;
}

//This function parses the program parameters. Returns false if given arguments are not valid
const bool prsArgs(int& nArgs, char** argList, string& pFl, string& tFl, uint32_t& k, uint32_t& w, double& hFrac, string& blFl, 
	uint32_t& cw, float& uw, float& tThres, bool& norm, float& dec, float& inter, bool& noNesting){
	int option_index = 0, a;

	static struct option long_options[] = {
        {"pattern",            required_argument,  0, 'p'},
        {"text",               required_argument,  0, 's'},
        {"ksize",              required_argument,  0, 'k'},
        {"windowsize",         required_argument,  0, 'w'},
        {"hashratio",          required_argument,  0, 'r'},
        {"blacklist",          required_argument,  0, 'b'},
        {"commonhashweight",   required_argument,  0, 'c'},
        {"uniquehashweight",   required_argument,  0, 'u'},
        {"hom_thres",          required_argument,  0, 't'},
        {"decent",             required_argument,  0, 'd'},
        {"intercept",          required_argument,  0, 'i'},
        {"normalize",          no_argument,        0, 'n'},
        {"nesting",            no_argument,        0, 'N'},
        {"help",               no_argument,        0, 'h'},
        {0,                    0,                  0,  0 }
    };

    //Parse all parameters given
	while ((a = getopt_long(nArgs, argList, T_HOM_OPTIONS, long_options, &option_index)) != -1){
		//Assign parameter values
		switch(a){
			case 'p':
				//Save input sequence
				pFl = optarg;
				break;
			case 's':
				//Save input sequence
				tFl = optarg;
				break;
			case 'k':
				//A k-mer length should be positive
				if(atoi(optarg) <= 0){
					cerr << "ERROR: K-mer length not applicable" << endl;
					return false;
				}

				k = atoi(optarg);
				break;
			case 'w':
				//Window size should be positive
				if(atoi(optarg) <= 0){
					cerr << "ERROR: Window size not applicable" << endl;
					return false;
				}

				w = atoi(optarg);
				break;
			case 'r':
				//Check if given value is reasonable to represent a ratio
				if(atof(optarg) <= 0 || atof(optarg) > MAX_RATIO){
					cerr << "ERROR: Given hash ratio not applicable" << endl;

					return false;
				}

				hFrac = atof(optarg);
				break;
			case 'b':
				//Save blacklist file name
				blFl = optarg;
				break;
			case 'c':
				//Weights should be positive
				if(atoi(optarg) <= 0){
					cerr << "ERROR: Common hash weight not applicable" << endl;

					return false;
				}

				cw = atoi(optarg);
				break;
			case 'u':
				//Weights should be positive
				if(atof(optarg) <= 0){
					cerr << "ERROR: Unique hash weight not applicable" << endl;

					return false;
				}

				uw = atof(optarg);
				break;
			case 't':
				tThres = atof(optarg);
				break;
			case 'd':
				dec = atof(optarg);
				break;
			case 'i':
				inter = atof(optarg);
				break;
			case 'n':
				norm = true;
				break;
			case 'N':
				noNesting = false;
				break;
			case 'h':
				return false;
			default:
				break;
		}
	}

	return !pFl.empty() && !tFl.empty();
}

//This function reads a file in FASTA format and returns true on success
//NOTE: Using this function, we first read in a complete sequence before calculating a sketch from it. Calculating a sketch could
//		also be done on the fly while reading the file. This possibility would be more memory saving most likely, because we would
//		not have to store the whole sequence in memory first. However, we had no chance to estimate the sketch size. Thus, the vec-
//		tor to store the sketch would have to be elongated several times which might be time consuming. On the other hand, it might
//		also be time consuming to handle the complete sequence in memory first. Eventually, it will depend on a try to find out
//		what is the better alternative. Could also be that it does not matter at all!
const bool readFASTA(const string& filePath, string& seq){
	bool headerRead = false, lnBrkDiscvd = false;
	char c;

	//Open the file
	ifstream fStr(filePath);

	//Check if the file is open
	if(!fStr.is_open()) return false;

	//Read in file character by character
	while(fStr.get(c)){
		//We are done if we find a second header (which can only start after at least one line break) in the file
		if(c == '>' && headerRead && lnBrkDiscvd) break;

		//Header lines are skipped
		if(c == '>'){
			headerRead = true;
			continue;
		}

		//The first line break indicates that we have left the header line
		if(c == '\n'){
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

//This function reads in batches of FASTA sequence entries from file and transforms them into sketches. Returns false if end of file
//was reached.
//const bool lPttnSks(ifstream& fStr, const uint32_t& k, const double& hFrac, const unordered_map<uint64_t, char>& bLstmers, 
//	vector<tuple<string, uint32_t, Sketch>>& pSks){
//	bool headerRead, idRead = false, lnBrkDiscvd = false;
//	char c;
//	string seq, seqID;
//
//	//Check if the file is open
//	if(!fStr.is_open()) return false;
//
//	//If this function is called iteratively, the '>' of the next entry has already been read
//	headerRead = fStr.gcount() != 0;
//
//	//Read in file character by character
//	while(fStr.get(c)){
//		//An entry's sequence is complete if we find a second header (which can only start after at least one line break) in the 
//		//file
//		if(c == '>' && headerRead && lnBrkDiscvd){
//			//Add sequence's sketch, length and id to result vector
//			pSks.push_back(make_tuple(seqID, seq.length(), buildSketch(seq, k, hFrac, bLstmers)));
//			//Clear sequence id
//			seqID.clear();
//			//Clear sequence
//			seq.clear();
//
//			//Check if enough sequences have been read
//			if(pSks.size() == PATTERN_BATCH_SIZE) return true;
//
//			//Reset id-read flag
//			idRead = false;
//			//Reset line-break-discovered flag
//			lnBrkDiscvd = false;
//			continue;
//		}
//
//		//Note if we have completely read the sequence id
//		idRead = idRead || (headerRead && c == ' ' && !lnBrkDiscvd);
//		//Note if we have found the first line break after a new header started
//		lnBrkDiscvd = lnBrkDiscvd || c == '\n';
//
//		//Update sequence id if we are still reading it
//		if(headerRead && !lnBrkDiscvd && !idRead){
//			seqID += c;
//			continue;
//		}
//
//		//Note if we have found the beginning of a header
//		headerRead = headerRead || (c == '>');
//
//		//There is no sequence to load in the header line
//		if(headerRead && !lnBrkDiscvd) continue;
//
//		//We are only interested in unambigous, unmasked nucleotides
//		if(c == 'A' || c == 'C' || c == 'G' || c == 'T') seq += c;
//	}
//
//	//Add last entry's sketch and sequence id to result vector if it is not empty
//	if(!seq.empty()) pSks.push_back(make_tuple(seqID, seq.length(), buildSketch(seq, k, hFrac, bLstmers)));
//
//	return false;
//}

//This function reads in batches of FASTA sequence entries from file and transforms them into minimap sketches. Returns false if end
// of file was reached.
const bool lMiniPttnSks(ifstream& fStr, const uint32_t& k, const uint32_t& w, const unordered_map<uint64_t, char>& blmers, 
	vector<tuple<string, uint32_t, Sketch>>& pSks){
	bool headerRead, idRead = false, lnBrkDiscvd = false;
	char c;
	string seq, seqID;

	//Check if the file is open
	if(!fStr.is_open()) return false;

	//If this function is called iteratively, the '>' of the next entry has already been read
	headerRead = fStr.gcount() != 0;

	//Read in file character by character
	while(fStr.get(c)){
		//An entry's sequence is complete if we find a second header (which can only start after at least one line break) in the 
		//file
		if(c == '>' && headerRead && lnBrkDiscvd){
			//Add sequence's sketch, length and id to result vector
			pSks.push_back(make_tuple(seqID, seq.length(), buildMiniSketch(seq, k, w, blmers)));//TODO: This function still needs to be tested!
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
		if(headerRead && !lnBrkDiscvd && !idRead){
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
const unordered_map<uint64_t, char> readBlstKmers(const string& fname){
	char nb[STRING_BUFFER_SIZE_DEFAULT];
	char* end;
	ifstream iFile(fname);
	unordered_map<uint64_t, char> numbers;

	//Check if file could be opened
	if(!iFile.is_open()){
		cerr << "readBlstKmers: WARNING input file could not be opened. Hash table will be empty!" << endl;

		return numbers;
	}

	//Read numbers
	while(iFile.getline(nb, STRING_BUFFER_SIZE_DEFAULT)) numbers[strtoull(nb, &end, 0)] = 1;

	return numbers;
}

#endif