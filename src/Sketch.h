#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <cstdint>
#include <utility>
#include <list>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include <math.h>

#include "../minimap2/minimap.h"
#include "../minimap2/index.c"

//The k-mer size
#define K 9
//The window size
#define W 10
//The hash value threshold
#define MAX_HASH 26214
//The FracMinHash ratio
#define HASH_RATIO 0.1
//Size of the considered alphabet
#define ALPHABET_SIZE 4
//Number of times we expect p to occur in t
#define P_MULTIPLICITY 2
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

using namespace std;
//A PairSketch is a list of offset hash pairs
using PairSketch = list<pair<uint32_t, uint64_t>>;
//A Sketch is a list of hashes
using Sketch = vector<uint64_t>;

//A compare function to sort hashes in a sketch
inline const bool smHshsFrst(const pair<uint32_t, uint64_t>& left, const pair<uint32_t, uint64_t>& right){ return left.second < right.second; }

//This function returns a numerical value for each nucleotide base
inline uint64_t getBaseNb(const char& c){
	switch(c){
		case 'A':
			return 0;
		case 'C':
			return 1;
		case 'G':
			return 2;
		case 'T':
			return 3;
		default:
			cerr << "WARNING: Unsupported character in input sequence detected!" << endl;
			return 0;
	}
}

//This function calculates the reverse complement of a DNA sequence
string revComp(const string &seq){
	//The result string
	string revSeq;

	//Go through the query from the end to the beginning
	for(int32_t i = seq.length() - 1; i >= 0; --i){
		//Check which base we are dealing with and append its complement
		switch(seq[i]){
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
static inline uint64_t getHash(uint64_t kmer, uint64_t mask){
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
uint64_t calcKmerNb(const string& kmer){
	uint64_t kmerNb = 0;

	for(string::const_iterator c = kmer.begin(); c != kmer.end(); ++c) kmerNb = (kmerNb << 2) + getBaseNb(*c);

	return kmerNb;
}

//This function builds a minimap2 sketch of a sequence by querying from a prebuilt minimap index; this function is influenced by the
// code of "The minimizer Jaccard estimator is biased and inconsistent." from Belbasi et al. (function 
//"winnowed_minimizers_linear(perm,windowSize)" in file "winnowed_minimizers.py").
const Sketch buildMiniSketch(const string& seq, const uint32_t& k, const uint32_t& w, const unordered_map<uint64_t, char>& blmers){
	uint32_t lastIdx = UINT32_MAX;
	int32_t windowBorder;
	const uint64_t mask = pow(ALPHABET_SIZE, k) - 1;
	uint64_t kmerHash, kmerBitSeq, revKmerBitSeq;
	//This stores pairs of k-mer starting positions and their hashes within the current window
	vector<pair<int32_t, uint64_t>> windowKmers;
	Sketch sk;

	//If the sequence is smaller than k we are done
	if(seq.length() < k){
		cerr << "WARNING: Length of input sequence " << seq << " is smaller than k (k=" << k << ")" << endl;

		return sk;
	}

	//Reserve as much space as is approximately needed to store the sketch (which hopefully saves some time)
	sk.reserve(seq.length() * 0.03);

	//Iterate over k-mer starting positions in sequence
	for(int32_t i = 0; i < seq.length() - k + 1; ++i){
		kmerBitSeq = calcKmerNb(seq.substr(i, k));
		revKmerBitSeq = calcKmerNb(revComp(seq.substr(i, k)));
		windowBorder = i - (w - 1);

		//We do not consider k-mers which are their own reverse complement
		if(kmerBitSeq == revKmerBitSeq) continue;

		//Depending on which is lexicographically smaller, we consider either the k-mer or its reverse complement
		if(kmerBitSeq < revKmerBitSeq){
			//Calculate k-mer's hash
			kmerHash = getHash(kmerBitSeq, mask);
		} else{
			//Calculate k-mer's hash
			kmerHash = getHash(revKmerBitSeq, mask);
		}
		
		//Remove all pairs with a larger hash from the back
		while(!windowKmers.empty() && windowKmers.back().second > kmerHash) windowKmers.pop_back();

		//Add current k-mer hash to windowKmers
		windowKmers.push_back(make_pair(i, kmerHash));

		//Remove pairs of k-mers which are no longer inside the window
		while(!windowKmers.empty() && windowKmers.front().first < windowBorder) windowKmers.erase(windowKmers.begin());

		//Choose a minimizer as soon as we have the first full window of k-mers and make sure we do the same minimizer a second time
		if(windowBorder >= 0 && !windowKmers.empty()){
			if(lastIdx != windowKmers.front().first && !blmers.contains(windowKmers.front().second)){
				lastIdx = windowKmers.front().first;
				sk.push_back(windowKmers.front().second);
			}

			//If the same k-mer appears several times inside the window and it has the smallest hash we want to save all occurrences
			while(windowKmers.size() > 1 && windowKmers.front().second == windowKmers[1].second){
				windowKmers.erase(windowKmers.begin());
				lastIdx = windowKmers.front().first;

				//Blacklisted k-mers are not added
				if(!blmers.contains(windowKmers.front().second)) sk.push_back(windowKmers.front().second);
			}
		}
	}

	//In case we have never seen a full window of k-mers take the one with the smallest hash seen for the sketch
	if(windowBorder < 0 && !windowKmers.empty()){
		//Blacklisted k-mers are not added
		if(!blmers.contains(windowKmers.front().second)) sk.push_back(windowKmers.front().second);

		//If the same k-mer appears several times inside the window and it has the smallest hash we want to save all occurrences
		while(windowKmers.size() > 1 && windowKmers.front().second == windowKmers[1].second){
			windowKmers.erase(windowKmers.begin());

			//Blacklisted k-mers are not added
			if(!blmers.contains(windowKmers.front().second)) sk.push_back(windowKmers.front().second);
		}
	}

	//Resize sketch (just for case we have allocated way too much memory)
	sk.shrink_to_fit();

	return sk;
}

#endif