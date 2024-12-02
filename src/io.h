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

#define T_HOM_OPTIONS "p:s:k:r:S:M:m:t:z:d:aonxhbB"

struct params_t {
	// required
	string pFile, tFile;

	// with an argument:
	int k;							// The k-mer length
	double hFrac;					// The FracMinHash ratio
	int max_seeds; 					// Maximum seeds in a sketch
	int max_matches; 				// Maximum seed matches in a sketch
	double theta; 					// The t-homology threshold
	string paramsFile;
	string mapper;                  // The name of the mapper
	double best_score_delta;         // Score difference threshold for considering mappings ambiguous

	// no arguments
	bool sam; 				// Output in SAM format (PAF by default)
	bool overlaps;			// Permit overlapping mappings 
	bool normalize; 		// Flag to save that scores are to be normalized
	bool onlybest;			// Output up to one (best) mapping (if above the threshold)

	// for degradation evals
	bool no_bucket_pruning;
	bool one_sweep;

    params_t() :
        k(15), hFrac(0.05), max_seeds(-1), max_matches(-1), theta(0.9), mapper("shmap"), best_score_delta(0.015),
        sam(false), overlaps(false), normalize(false), onlybest(false), no_bucket_pruning(false), 
        one_sweep(false) {}

    void print(std::ostream& out, bool human) const {
        std::vector<pair<string, string>> m;
        m.push_back({"pFile", pFile});
        m.push_back({"tFile", tFile});
        m.push_back({"k", std::to_string(k)});
        m.push_back({"hFrac", std::to_string(hFrac)});
        m.push_back({"max_seeds", std::to_string(max_seeds)});
        m.push_back({"max_matches", std::to_string(max_matches)});
        m.push_back({"tThres", std::to_string(theta)});
        m.push_back({"paramsFile", paramsFile});
        m.push_back({"mapper", mapper});
        m.push_back({"best_score_delta", std::to_string(best_score_delta)});  // Add to output

        m.push_back({"sam", std::to_string(sam)});
        m.push_back({"overlaps", std::to_string(overlaps)});
        m.push_back({"normalize", std::to_string(normalize)});
        m.push_back({"onlybest", std::to_string(onlybest)});
        m.push_back({"no-bucket-pruning", std::to_string(no_bucket_pruning)});
        m.push_back({"one-sweep", std::to_string(one_sweep)});

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
        out << " | reads:                 " << pFile << endl;
        out << " | algorithm:             " << mapper << endl;
        out << " | k:                     " << k << endl;
        out << " | hFrac:                 " << hFrac << endl;
        out << " | max_seeds              " << max_seeds << endl;
        out << " | max_matches:           " << max_matches << endl;
        out << " | sam:                   " << sam << endl;
        out << " | overlaps:              " << overlaps << endl;
        out << " | onlybest:              " << onlybest << endl;
        out << " | no-bucket-pruning:     " << no_bucket_pruning << endl;
        out << " | one-sweep:             " << one_sweep << endl;
        out << " | tThres:                " << theta << endl;
        out << " | best_score_delta:      " << best_score_delta << endl;  // Add to display
    }

    void dsHlp() {
        cerr << "sweep [-hn] [-p PATTERN_FILE] [-s TEXT_FILE] [-k KMER_LEN] [-r HASH_RATIO] [-b BLACKLIST] [-c COM_HASH_WGHT] [-u UNI\
        _HASH_WGHT] [-t HOMOLOGY_THRESHOLD]" << endl;
        cerr << endl;
        cerr << "Find sketch-based pattern similarity in text." << endl;
        cerr << endl;
        cerr << "Required parameters:" << endl;
        cerr << "   -p   --pattern           Pattern sequences file (FASTA format)" << endl;
        cerr << "   -s   --text              Text sequence file (FASTA format)" << endl;
        cerr << endl;
        cerr << "Optional parameters with an argument:" << endl;
        cerr << "   -m   --mapper            Mapper name {sweep, bucket} (default: sweep)" << endl;
        cerr << "   -k   --ksize             K-mer length to be used for sketches" << endl;
        cerr << "   -r   --ratio             FracMinHash ratio in [0; 1] [0.1]" << endl;
        cerr << "   -S   --max_seeds         Max seeds in a sketch" << endl;
        cerr << "   -M   --max_matches       Max seed matches in a sketch" << endl;
        cerr << "   -t   --threshold         Homology percentage threshold [0, 1]" << endl;
        cerr << "   -z   --params            Output file with parameters (tsv)" << endl;
        cerr << "   -d   --delta             Score difference threshold for ambiguous mappings [0.015]" << endl;  // Add help text
        cerr << endl;
        cerr << "Optional parameters without an argument:" << endl;
        cerr << "   -a                       Output in SAM format (PAF by default)" << endl;
        cerr << "   -o   --overlaps          Permit overlapping mappings" << endl;
        cerr << "   -n   --normalize         Normalize scores by length" << endl;
        cerr << "   -x   --onlybest          Output the best alignment if above threshold (otherwise none)" << endl;
        cerr << "   -b   --no-bucket-pruning Disables bucket pruning" << endl;
        cerr << "   -B   --one-sweep         Disregards the seed heuristic and runs one sweepmap on all matches" << endl;
        cerr << "   -h   --help              Display this help message" << endl;
    }

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
            {"mapper",             required_argument,  0, 'm'},
            {"threshold",          required_argument,  0, 't'},
            {"delta",              required_argument,  0, 'd'},  // Add option
            {"overlaps",           no_argument,        0, 'o'},
            {"normalize",          no_argument,        0, 'n'},
            {"onlybest",           no_argument,        0, 'x'},
            {"no-bucket-pruning",  no_argument,        0, 'b'},
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
				case 'm':
					mapper = optarg;
					break;
				case 'k':
					if(atoi(optarg) <= 0) {
						cerr << "ERROR: K-mer length not applicable" << endl;
                        return false;
                    }
					k = atoi(optarg);
                    break;
                case 'd':
                    if(atof(optarg) < 0.0) {
                        cerr << "ERROR: The delta threshold should be non-negative." << endl;
                        return false;
                    }
                    best_score_delta = atof(optarg);
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
					if(atoi(optarg) < 0 || atoi(optarg) > 1.0) {
						cerr << "ERROR: The threshold t should be between 0 and 1." << endl;
						return false;
					}
					theta = atof(optarg);
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
				case 'b':
					no_bucket_pruning = true;
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

struct ParsedQueryId {
	bool valid;
	string segm_id;
	rpos_t start_pos;
	rpos_t end_pos;
	char strand;
	
	static ParsedQueryId parse(const string& query_id) {
		// query_id is in the form "S1_21!NC_060948.1!57693539!57715501!+"
		ParsedQueryId result{false, "", 0, 0, '?'};
		
		// Skip if query_id doesn't contain correct mapping info
		if (query_id.find('!') == string::npos) return result;

		// Split query_id on '!' character
		std::vector<string> parts;
		size_t start = 0;
		size_t end = query_id.find('!');
		while (end != string::npos) {
			parts.push_back(query_id.substr(start, end - start));
			start = end + 1;
			end = query_id.find('!', start);
		}
		parts.push_back(query_id.substr(start));

		// Need 5 parts: read_id, segment, start, end, strand
		if (parts.size() != 5) return result;

		result.segm_id = parts[1];
		result.start_pos = stoll(parts[2]);
		result.end_pos = stoll(parts[3]);
		result.strand = parts[4][0];
		result.valid = true;
		
		return result;
	}
};
	

// seq->name.s, seq->comment.l, seq->comment.s, seq->seq.s, seq->qual.l
void read_fasta_klib(const std::string& filename, std::function<void(const std::string&, const std::string&)> callback);

} // namespace sweepmap