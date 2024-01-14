#ifndef UTILS_HPP
#define UTILS_HPP

#include <chrono>
#include <iostream>
#include <string>
#include <unordered_map>

using std::string;

using hash_t     = uint64_t;
using pos_t      = int32_t;
using kmer_num_t = int32_t;

const double EPS = 1e-7;

//Character constants used for nucleotide bases
#define NUCL_BASE_A 'A'
#define NUCL_BASE_C 'C'
#define NUCL_BASE_G 'G'
#define NUCL_BASE_T 'T'

//Complements of nucleotide bases
#define CMPL_BASE_A 'T'
#define CMPL_BASE_C 'G'
#define CMPL_BASE_G 'C'
#define CMPL_BASE_T 'A'

class Timer {
private:
    std::chrono::time_point<std::chrono::high_resolution_clock> start_time_point_, end_time_point_;
    double accumulated_time_;
    bool running_;

public:
    Timer() : accumulated_time_(0.0), running_(false) {}

    void start() {
        if (!running_) {
            start_time_point_ = std::chrono::high_resolution_clock::now();
            running_ = true;
        }
    }

    void stop() {
		assert(running_);
        if (running_) {
            end_time_point_ = std::chrono::high_resolution_clock::now();
            accumulated_time_ += std::chrono::duration<double>(end_time_point_ - start_time_point_).count();
            running_ = false;
        }
    }

    double get_secs() const {
		assert(!running_);
		return accumulated_time_;
    }
};

class Timers {
private:
    std::unordered_map<std::string, Timer> timers_;

public:
    void start(const std::string& name) {
        timers_[name].start();
    }

    void stop(const std::string& name) {
        auto it = timers_.find(name);
        if (it != timers_.end()) {
            it->second.stop();
        }
    }

    double get_secs(const std::string& name) const {
        auto it = timers_.find(name);
		assert(it != timers_.end());
        if (it != timers_.end()) {
            return it->second.get_secs();
        }
        return 0.0;
    }

    double get_perc(const std::string& name, const std::string& total) const {
        auto it = timers_.find(name);
		assert(it != timers_.end());
        auto it_total = timers_.find(total);
		assert(it_total != timers_.end());
        if (it != timers_.end()) {
			if (it_total != timers_.end()) {
				assert(it->second.get_secs() <= it_total->second.get_secs());
				return it->second.get_secs() / it_total->second.get_secs() * 100.0;
			}
        }
        return 0.0;
    }
};

struct Match {
	kmer_num_t kmer_ord;
	pos_t P_l;
	pos_t T_r;
	pos_t t_pos;
	//bool strand;  // 0 - forward, 1 - reverse

	//char get_strand() const {
	//	return strand ? '-' : '+';
	//}
};

struct Mapping {
	int k; 	   // kmer size
	pos_t P_sz;     // pattern size |P| bp 
	pos_t p_sz;     // pattern sketch size |p| kmers
	int matches;  // L.size() -- total number of matches in `t' 
	pos_t T_l;      // the position of the leftmost nucleotide of the mapping
	pos_t T_r;      // the position of the rightmost nucleotide of the mapping
	int xmin;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	pos_t dT_l;      // delta to be applied before output
	pos_t dT_r;      // -- || --
	double J;          // Jaccard score/similarity [0;1]
	double map_time;

	Mapping(int k=0, pos_t P_sz=0, pos_t p_sz=0, int matches=0, pos_t T_l=0, pos_t T_r=0, int xmin=0, pos_t dT_l=0, pos_t dT_r=0, double J=0.0)
		: k(k), P_sz(P_sz), p_sz(p_sz), matches(matches), T_l(T_l), T_r(T_r), xmin(xmin), dT_l(dT_l), dT_r(dT_r), J(J) {}

	//void print() {
	//	std::cerr << "k=" << k 
	//		<< " |P|=" << P_sz
	//		<< " |p|=" << p_sz
	//		<< " matches=" << matches
	//		<< " T_l=" << T_l
	//		<< " T_r=" << T_r
	//		<< " xmin=" << xmin
	//		<< " dT_l=" << dT_l
	//		<< " dT_r=" << dT_r
	//		<< " J=" << J
	//		<< " map_time=" << map_time
	//		<< std::endl;
	//}
};

template<typename TT> auto prev(const typename TT::iterator &it) {
    auto pr = it; return --pr;
}

//template<typename TT> auto next(const typename TT::iterator &it) {
//    auto pr = it; return ++pr;
//}

//This function returns a numerical value for each nucleotide base
inline hash_t getBaseNb(const char& c) {
	switch(c) {
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		default:
			std::cerr << "WARNING: Unsupported character in input sequence detected!" << std::endl;
			return 0;
	}
}

inline char to_upper(char c) {
	return (c >= 'a' && c <= 'z') ? c - 'a' + 'A' : c;
}

//This function calculates the reverse complement of a DNA sequence
string revComp(const string &seq) {
	//The result string
	string revSeq;

	//Go through the query from the end to the beginning
	for(int i = (int)seq.length() - 1; i >= 0; --i) {
		//Check which base we are dealing with and append its complement
		switch(seq[i]) {
			case NUCL_BASE_A:
				revSeq += CMPL_BASE_A;
				break;
			case NUCL_BASE_C:
				revSeq += CMPL_BASE_C;
				break;
			case NUCL_BASE_G:
				revSeq += CMPL_BASE_G;
				break;
			case NUCL_BASE_T:
				revSeq += CMPL_BASE_T;
				break;
			default:
				revSeq += "N";
				break;
		}
	}

	return revSeq;
}


//This function calculates a hash for the given numerical representation of a k-mer and a mask of the form (4**k)-1 where k is the k-mer length;
//it originates from the minimap2 source code (the function hash64 in sketch.c)
static inline uint64_t getHash(uint64_t kmer, uint64_t mask) {
	kmer = (~kmer + (kmer << 21)) & mask; // kmer = (kmer << 21) - kmer - 1;
	kmer = kmer ^ kmer >> 24;
	kmer = ((kmer + (kmer << 3)) + (kmer << 8)) & mask; // kmer * 265
	kmer = kmer ^ kmer >> 14;
	kmer = ((kmer + (kmer << 2)) + (kmer << 4)) & mask; // kmer * 21
	kmer = kmer ^ kmer >> 28;
	kmer = (kmer + (kmer << 31)) & mask;
	return kmer;
}

//This function calculates the numerical representation of a k-mer
uint64_t calcKmerNb(const string& kmer) {
	uint64_t kmerNb = 0;
	for(string::const_iterator c = kmer.begin(); c != kmer.end(); ++c)
		kmerNb = (kmerNb << 2) + getBaseNb(*c);
	return kmerNb;
}

#endif