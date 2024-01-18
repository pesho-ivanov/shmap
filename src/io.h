#ifndef IO_HPP
#define IO_HPP

#include <fstream>
#include <functional>
#include <getopt.h>
#include <iomanip>
#include <string>
#include <vector>

#include <zlib.h>  
#include "kseq.h"  
#include "utils.h"

using std::cerr;
using std::cout;
using std::vector;
using std::string;
using std::pair;
using std::ifstream;
using std::endl;

#define INDEX_DEFAULT_DUMP_FILE "indexDump.idx"
#define DEFAULT_WEIGHT 1

#define T_HOM_OPTIONS "p:s:k:w:r:e:a:b:S:M:t:z:nxoh"

//#define MIN_PARAM_NB 6
//#define MAX_RATIO 1.0
//#define NORM_FLAG_DEFAULT false
#define PATTERN_BATCH_SIZE 250000
#define STRING_BUFFER_SIZE_DEFAULT 50

//The hash value threshold
#define MAX_HASH 26214
//The FracMinHash ratio
//#define HASH_RATIO 0.1
//Size of the considered alphabet
#define ALPHABET_SIZE 4
//Number of times we expect p to occur in t
#define P_MULTIPLICITY 2

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

string getElasticDescription(elastic_t e) {
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

string getAlignmentEdgesDescription(alignment_edges_t ae) {
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
	int k;  						//The k-mer length
	int w;  							//The window size
	elastic_t elastic;
	alignment_edges_t alignment_edges;
	double hFrac;  					//The FracMinHash ratio
	int max_seeds;  							//Maximum seeds in a sketch
	int max_matches;  							//Maximum seed matches in a sketch
	double tThres;  							//The t-homology threshold
	string pFile, tFile;
	string paramsFile;
	string bLstFl; // = "highAbundKmersMiniK15w10Lrgr100BtStrnds.txt";  //Input file names

	params_t() {
		normalize = false; 		//Flag to save that scores are to be normalized
		onlybest = false;
		overlaps = false;          // 
		k = 15; 						//The k-mer length
		w = 10; 						//The window size
		elastic = elastic_t::off;
		alignment_edges = alignment_edges_t::extend_equally;
		hFrac = 0.05;
		max_seeds = 10000;
		max_matches = 1000000;
		tThres = 0.9; 							//The t-homology threshold
		bLstFl = ""; //"highAbundKmersMiniK15w10Lrgr100BtStrnds.txt";  //Input file names
	}

	void print(std::ostream& out = std::cout, bool human = true) {
		vector<pair<string, string>> m;
		m.push_back({"normalize", std::to_string(normalize)});
		m.push_back({"onlybest", std::to_string(onlybest)});
		m.push_back({"overlaps", std::to_string(overlaps)});
		m.push_back({"k", std::to_string(k)});
		m.push_back({"w", std::to_string(w)});
		m.push_back({"hFrac", std::to_string(hFrac)});
		m.push_back({"elastic", getElasticDescription(elastic)});
		m.push_back({"alignment_edges", getAlignmentEdgesDescription(alignment_edges)});
		m.push_back({"tThres", std::to_string(tThres)});
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

//This function prints usage infos
inline void dsHlp() {
	cerr << "sweep [-hn] [-p PATTERN_FILE] [-s TEXT_FILE] [-k KMER_LEN] [-r HASH_RATIO] [-b BLACKLIST] [-c COM_HASH_WGHT] [-u UNI\
	_HASH_WGHT] [-t HOM_THRES] [-d DECENT] [-i INTERCEPT]" << endl << endl;
	cerr << "Find sketch-based pattern similarity in text." << endl << endl;
	cerr << "Required parameters:" << endl;
	cerr << "   -p   --pattern  Pattern sequences file (FASTA format)" << endl;
	cerr << "   -s   --text     Text sequence file (FASTA format)" << endl << endl;
	cerr << "Optional parameters with required argument:" << endl;
	cerr << "   -k   --ksize             K-mer length to be used for sketches" << endl;
	cerr << "   -w   --windowsize        Window size for minimizer sketching approach" << endl;
	cerr << "   -e   --elastic           Elastic pairs of kmers {off, consecutive, random} [off]" << endl;
	cerr << "   -a   --alignment_edges   Alignment interval {sketch_edges, extend_equally, fine} [fine]" << endl;
	//cerr << "   -c   --hashratio         FracMin hash ratio to be used for sketches (default " << HASH_RATIO << ")" << endl;
	cerr << "   -b   --blacklist         File containing hashes to ignore for sketch calculation" << endl;
	cerr << "   -S   --max_seeds         Max seeds in a sketch" << endl;
	cerr << "   -M   --max_matches       Max seed matches in a sketch" << endl;
	cerr << "   -t   --hom_thres         Homology threshold" << endl;
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
        {"hashratio",          required_argument,  0, 'r'},
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
			case 'r':
				//Check if given value is reasonable to represent a ratio
				if(atof(optarg) <= 0 || atof(optarg) > 1.0) {
					cerr << "ERROR: Given hash ratio " << optarg << " not applicable" << endl;
					return false;
				}
				params->hFrac = atof(optarg);
				break;
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

KSEQ_INIT(gzFile, gzread)  

// seq->name.s, seq->comment.l, seq->comment.s, seq->seq.s, seq->qual.l
void read_fasta_klib(const std::string& filename, std::function<void(kseq_t*)> callback) {
    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int l;

    while ((l = kseq_read(seq)) >= 0) {
        callback(seq);
    }

    kseq_destroy(seq);
    gzclose(fp);
}

#endif