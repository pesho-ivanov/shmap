#pragma once

#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <string>
#include <vector>

#include <zlib.h>  
#include "kseq.h"
//#include "cxxopts.hpp"
#include "utils.h"

namespace sweepmap {

KSEQ_INIT(gzFile, gzread)  

using std::cerr;
using std::string;
using std::pair;
using std::ifstream;
using std::endl;

#define T_HOM_OPTIONS "p:s:k:r:S:M:t:z:aonxh"

struct params_t {
	// required
	string pFile, tFile;

	// with an argument:
	int k;							// The k-mer length
	double hFrac;					// The FracMinHash ratio
	int max_seeds; 					// Maximum seeds in a sketch
	int max_matches; 				// Maximum seed matches in a sketch
	double tThres; 					// The t-homology threshold
	string paramsFile;

	// no arguments
	bool sam; 				// Output in SAM format (PAF by default)
	bool overlaps;			// Permit overlapping mappings 
	bool normalize; 		// Flag to save that scores are to be normalized
	bool onlybest;			// Output up to one (best) mapping (if above the threshold)

	params_t() :
		k(15), hFrac(0.05), max_seeds(10000), max_matches(1000000), tThres(0.9),
		sam(false), overlaps(false), normalize(false), onlybest(false) {}

	void print(std::ostream& out, bool human) const {
		std::vector<pair<string, string>> m;
		m.push_back({"pFile", pFile});
		m.push_back({"tFile", tFile});
		m.push_back({"k", std::to_string(k)});
		m.push_back({"hFrac", std::to_string(hFrac)});
		m.push_back({"max_seeds", std::to_string(max_seeds)});
		m.push_back({"max_matches", std::to_string(max_matches)});
		m.push_back({"tThres", std::to_string(tThres)});
		m.push_back({"paramsFile", paramsFile});

		m.push_back({"sam", std::to_string(sam)});
		m.push_back({"overlaps", std::to_string(overlaps)});
		m.push_back({"normalize", std::to_string(normalize)});
		m.push_back({"onlybest", std::to_string(onlybest)});

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

	void print_display(std::ostream& out) {
		out << "Params:" << endl;
		out << " | reference:             " << tFile << endl;
		out << " | queries:               " << pFile << endl;
		out << " | k:                     " << k << endl;
		out << " | hFrac:                 " << hFrac << endl;
		out << " | max_seeds (S):         " << max_seeds << endl;
		out << " | max_matches (M):       " << max_matches << endl;
		out << " | sam:                   " << sam << endl;
		out << " | overlaps:              " << overlaps << endl;
		out << " | onlybest:              " << onlybest << endl;
		out << " | tThres:                " << tThres << endl;
	}

	void dsHlp() {
		cerr << "sweep [-hn] [-p PATTERN_FILE] [-s TEXT_FILE] [-k KMER_LEN] [-r HASH_RATIO] [-b BLACKLIST] [-c COM_HASH_WGHT] [-u UNI\
		_HASH_WGHT] [-t HOM_THRES] [-d DECENT] [-i INTERCEPT]" << endl;
		cerr << endl;
		cerr << "Find sketch-based pattern similarity in text." << endl;
		cerr << endl;
		cerr << "Required parameters:" << endl;
		cerr << "   -p   --pattern           Pattern sequences file (FASTA format)" << endl;
		cerr << "   -s   --text              Text sequence file (FASTA format)" << endl;
		cerr << endl;
		cerr << "Optional parameters with an argument:" << endl;
		cerr << "   -k   --ksize             K-mer length to be used for sketches" << endl;
		cerr << "   -r   --ratio   			 FracMinHash ratio in [0; 1] [0.1]" << endl;
		cerr << "   -S   --max_seeds         Max seeds in a sketch" << endl;
		cerr << "   -M   --max_matches       Max seed matches in a sketch" << endl;
		cerr << "   -t   --hom_thres         Homology threshold" << endl;
		cerr << "   -z   --params     		 Output file with parameters (tsv)" << endl;
		cerr << endl;
		cerr << "Optional parameters without an argument:" << endl;
		cerr << "   -a                       Output in SAM format (PAF by default)" << endl;
		cerr << "   -o   --overlaps          Permit overlapping mappings" << endl;
		cerr << "   -n   --normalize         Normalize scores by length" << endl;
		cerr << "   -x   --onlybest          Output the best alignment if above threshold (otherwise none)" << endl;
		cerr << "   -h   --help              Display this help message" << endl;
	}

	//This function parses the program parameters. Returns false if given arguments are not valid
	bool prsArgs(int& nArgs, char** argList) {
		static struct option long_options[] = {
			{"pattern",            required_argument,  0, 'p'},
			{"text",               required_argument,  0, 's'},
			{"ksize",              required_argument,  0, 'k'},
			{"hashratio",          required_argument,  0, 'r'},
			{"max_seeds",          required_argument,  0, 'S'},
			{"max_matches",        required_argument,  0, 'M'},
			{"hom_thres",          required_argument,  0, 't'},
			{"params",             required_argument,  0, 'z'},
			{"overlaps",           no_argument,        0, 'o'},
			{"normalize",          no_argument,        0, 'n'},
			{"onlybest",           no_argument,        0, 'x'},
			{"help",               no_argument,        0, 'h'},
			{0,                    0,                  0,  0 }
		};

		int option_index = 0, a;
		while ((a = getopt_long(nArgs, argList, T_HOM_OPTIONS, long_options, &option_index)) != -1) {
			switch(a) {
				case 'p':
					pFile = optarg;
					break;
				case 's':
					tFile = optarg;
					break;
				case 'k':
					if(atoi(optarg) <= 0) {
						cerr << "ERROR: K-mer length not applicable" << endl;
						return false;
					}
					k = atoi(optarg);
					break;
				case 'r':
					if(atof(optarg) <= 0 || atof(optarg) > 1.0) {
						cerr << "ERROR: Given hash ratio " << optarg << " not applicable" << endl;
						return false;
					}
					hFrac = atof(optarg);
					break;
				case 'S':
					if(atoi(optarg) <= 0) {
						cerr << "ERROR: The number of seeds should be positive." << endl;
						return false;
					}
					max_seeds = atoi(optarg);
					break;
				case 'M':
					if(atoi(optarg) <= 0) {
						cerr << "ERROR: The number of seed matches should be positive." << endl;
						return false;
					}
					max_matches = atoi(optarg);
					break;
				case 't':
					tThres = atof(optarg);
					break;
				case 'z':
					paramsFile = optarg;
					break;
				case 'a':
					sam = true;
					break;
				case 'o':
					overlaps = true;
					break;
				case 'n':
					normalize = true;
					break;
				case 'x':
					onlybest = true;
					break;
				case 'h':
					return false;
				default:
					cerr << "Unknown option " << a << " '" << char(a) << "', " << optarg << endl ;
					break;
			}
		}

		return !pFile.empty() && !tFile.empty();
	}
};

// seq->name.s, seq->comment.l, seq->comment.s, seq->seq.s, seq->qual.l
void read_fasta_klib(const std::string& filename, std::function<void(kseq_t*)> callback);

} // namespace sweepmap