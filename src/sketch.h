#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <bit>
//#include <bitset>
#include <cstdint>
#include <utility>
#include <list>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <math.h>

#include "utils.h"
#include "io.h"
#include <climits>

using namespace std;

using abs_hash_t = pair<pos_t, hash_t>;
using abs_ord_t  = pair<pos_t, kmer_num_t>; 

using Sketch     = vector<abs_hash_t>;  // (kmer's left 0-based position, kmer hash)

using std::rotl;

static hash_t LUT_fw[256], LUT_rc[256];
static Timer FMH_time;

void print_sketches(const string &seqID, const Sketch &sks) {
	cout << seqID << endl;
	for (const auto& sk : sks) {
		cout << "  " << sk.first << ", " << sk.second << endl;
	}
}

void initialize_LUT() {
	LUT_fw['a'] = LUT_fw['A'] = 0x3c8b'fbb3'95c6'0474; // Daniel's
	//LUT_fw['a'] = LUT_fw['A'] = 0x3c8bfbb395c60470;  // Ragnar's
	LUT_fw['c'] = LUT_fw['C'] = 0x3193'c185'62a0'2b4c; // Daniel's
	LUT_fw['g'] = LUT_fw['G'] = 0x2032'3ed0'8257'2324; // Daniel's
	LUT_fw['t'] = LUT_fw['T'] = 0x2955'49f5'4be2'4456; // Daniel's
	//LUT_fw['t'] = LUT_fw['T'] = 0x2d2a04e675310c18;  // Ragnar's
	
	//l[b'A' as usize] = 0x3c8bfbb395c60474u64;
    //l[b'C' as usize] = 0x3193c18562a02b4cu64;
    //l[b'G' as usize] = 0x20323ed082572324u64;
    //l[b'T' as usize] = 0x295549f54be24456u64;

	LUT_rc['a'] = LUT_rc['A'] = LUT_fw['T'];
	LUT_rc['c'] = LUT_rc['C'] = LUT_fw['G'];
	LUT_rc['g'] = LUT_rc['G'] = LUT_fw['C'];
	LUT_rc['t'] = LUT_rc['T'] = LUT_fw['A'];
}

const Sketch buildFMHSketch(const string& s, int k, double hFrac, const unordered_map<hash_t, char>& blmers) {
	FMH_time.start();

	Sketch sk;
	sk.reserve((int)(1.1 * (double)s.size() * hFrac));

	if ((int)s.size() < k) return sk;

	hash_t h, h_fw = 0, h_rc = 0;
	hash_t hThres = hash_t(hFrac * double(std::numeric_limits<hash_t>::max()));
	int r;

	for (r=0; r<k; r++) {
		h_fw ^= rotl(LUT_fw[(int)s[r]], k-r-1);
		h_rc ^= rotl(LUT_rc[(int)s[r]], r);
	}

	while(true) {
		h = h_fw ^ h_rc;
//		cerr << s << " " << r << " " << std::hex << std::setfill('0') << std::setw(16) << h
//							<< " = " << std::hex << std::setfill('0') << std::setw(16) << h_fw
//							<< " ^ " << std::hex << std::setfill('0') << std::setw(16) << h_rc << endl;
		if (h < hThres)  // optimize to only look at specific bits
			//sk.emplace_back(r, h);
			sk.emplace_back(r, h_fw);
		
		if (r >= (int)s.size()) break;

		h_fw = rotl(h_fw, 1) ^ rotl(LUT_fw[(int)s[r-k]], k) ^ LUT_fw[(int)s[r]];
		h_rc = rotr(h_rc, 1) ^ rotr(LUT_rc[(int)s[r-k]], 1) ^ rotl(LUT_rc[(int)s[r]], k-1);

		++r;
	}

	FMH_time.stop();
	return sk;
}

struct kmer_hits_t {
	pos_t P_l;
	kmer_num_t kmer_num;
	vector<abs_ord_t> kmers_in_T;

	kmer_hits_t(pos_t P_l, kmer_num_t kmer_num, const vector<abs_ord_t> &kmers_in_T) :
		P_l(P_l), kmer_num(kmer_num), kmers_in_T(kmers_in_T) {}
};

struct SketchIndex {
	//double s;
	//const params_t &params;
	pos_t T_sz;
	string name;
	unordered_map<hash_t, vector<abs_ord_t>> h2pos;

	void populate_h2pos(const Sketch& sketch) {
		//print_sketches(name, sketch);
		for (kmer_num_t tpos = 0; tpos < (int)sketch.size(); ++tpos) {
			const abs_hash_t& abs_hash = sketch[tpos];
			h2pos[abs_hash.second].push_back(abs_ord_t(abs_hash.first, tpos));
		}
	}

	SketchIndex(const Sketch& sketch, pos_t T_sz)
		: T_sz(T_sz) {
		populate_h2pos(sketch);
	}

	SketchIndex(const string &name, const string &ref, const params_t &params, const unordered_map<hash_t, char>& bLstmers)
		: T_sz((pos_t)ref.size()), name(name) {
		Sketch sketch = buildFMHSketch(ref, params.k, params.hFrac, bLstmers);
		populate_h2pos(sketch);
	}

	SketchIndex(const string &tFile, const params_t &params, const unordered_map<hash_t, char>& bLstmers) {
		string ref;
		cerr << "Reading index " << tFile << "..." << endl;
		if (readFASTA(tFile, &ref)) {
			cerr << "ERROR: Reference file could not be read" << endl;
		}
		T_sz = (pos_t)ref.size();

		ifstream fStr(tFile);
		string line;
		getline(fStr, line);
		name = line.substr(1);

		//Sketch sketch = buildMiniSketch(ref, params.k, params.w, bLstmers);
		Sketch sketch = buildFMHSketch(ref, params.k, params.hFrac, bLstmers);
		populate_h2pos(sketch);
	}
};

//This function reads in batches of FASTA sequence entries from file and transforms them into minimap sketches. Returns false if end
// of file was reached.
bool SketchReads(ifstream& fStr, const params_t &params, const unordered_map<hash_t, char>& blmers, 
	vector<std::tuple<string, pos_t, Sketch>> *pSks) {
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
			//auto sketch = buildMiniSketch(seq, params.k, params.w, blmers);
			auto sketch = buildFMHSketch(seq, params.k, params.hFrac, blmers);
			pSks->push_back(make_tuple(seqID, seq.length(), sketch));//TODO: This function still needs to be tested!
            //cout << "pattern: " << seq << endl;
			//Clear sequence id
			seqID.clear();
			//Clear sequence
			seq.clear();

			//Check if enough sequences have been read
			if(pSks->size() == PATTERN_BATCH_SIZE) return true;

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
        //pSks->push_back(make_tuple(seqID, seq.length(), buildMiniSketch(seq, params.k, params.w, blmers)));
        pSks->push_back(make_tuple(seqID, seq.length(), buildFMHSketch(seq, params.k, params.hFrac, blmers)));
    }

	return false;
}

//This function builds a minimap2 sketch of a sequence by querying from a prebuilt minimap index; this function is influenced by the
// code of "The minimizer Jaccard estimator is biased and inconsistent." from Belbasi et al. (function 
//"winnowed_minimizers_linear(perm,windowSize)" in file "winnowed_minimizers.py").
//const Sketch buildMiniSketch(const string& seq, const int& k, const int& w, const unordered_map<hash_t, char>& blmers) {
//	pos_t lastIdx = std::numeric_limits<pos_t>::max();
//	pos_t windowBorder;
//	const hash_t mask = static_cast<hash_t>(pow(ALPHABET_SIZE, k) - 1);
//	hash_t kmerHash, kmerBitSeq, revKmerBitSeq;
//	//This stores pairs of k-mer starting positions and their hashes within the current window
//	vector<pair<pos_t, hash_t>> windowKmers; // TODO: this can be sped up by using a fixed-size cyclic array with two pointers for begin and end
//	Sketch sk;
//
//	//If the sequence is smaller than k we are done
//	if((int)seq.length() < k) {
//		cerr << "WARNING: Length of input sequence " << seq << " is smaller than k (k=" << k << ")" << endl;
//		return sk;
//	}
//
//	//Reserve as much space as is approximately needed to store the sketch (which hopefully saves some time)
//	sk.reserve((int)((double)seq.length() * 0.03));
//
//	//Iterate over k-mer starting positions in sequence
//	for(pos_t i = 0; i < (int)seq.length() - k + 1; ++i) {
//		// TODO: don't take the substr each time, work only with 1 letter at a time
//		kmerBitSeq = calcKmerNb(seq.substr(i, k));
//		revKmerBitSeq = calcKmerNb(revComp(seq.substr(i, k)));
//		windowBorder = i - (w - 1);
//
//		//We do not consider k-mers which are their own reverse complement
//		if(kmerBitSeq == revKmerBitSeq) continue;
//
//		//Depending on which is lexicographically smaller, we consider either the k-mer or its reverse complement
//		// TODO: rewrite shorter
//		if(kmerBitSeq < revKmerBitSeq) {
//			//Calculate k-mer's hash
//			kmerHash = getHash(kmerBitSeq, mask);
//		} else{
//			//Calculate k-mer's hash
//			kmerHash = getHash(revKmerBitSeq, mask);
//		}
//		
//		//Remove all pairs with a larger hash from the back
//		while(!windowKmers.empty() && windowKmers.back().second > kmerHash)
//			windowKmers.pop_back();
//
//		//Add current k-mer hash to windowKmers
//		windowKmers.push_back(make_pair(i+k, kmerHash));
//
//		//Remove pairs of k-mers which are no longer inside the window
//		while(!windowKmers.empty() && windowKmers.front().first < windowBorder)
//			windowKmers.erase(windowKmers.begin());
//
//		//Choose a minimizer as soon as we have the first full window of k-mers and make sure we do the same minimizer a second time
//		if(windowBorder >= 0 && !windowKmers.empty()) {
//			if(lastIdx != windowKmers.front().first && !blmers.contains(windowKmers.front().second)) {
//				lastIdx = windowKmers.front().first;
//				sk.push_back(windowKmers.front());
//			}
//
//			//If the same k-mer appears several times inside the window and it has the smallest hash we want to save all occurrences
//			while(windowKmers.size() > 1 && windowKmers.front().second == windowKmers[1].second) {
//				windowKmers.erase(windowKmers.begin());
//				lastIdx = windowKmers.front().first;
//
//				//Blacklisted k-mers are not added
//				if(!blmers.contains(windowKmers.front().second))
//					sk.push_back(windowKmers.front());
//			}
//		}
//	}
//
//	//In case we have never seen a full window of k-mers take the one with the smallest hash seen for the sketch
//	if(windowBorder < 0 && !windowKmers.empty()) {
//		//Blacklisted k-mers are not added
//		if(!blmers.contains(windowKmers.front().second)) sk.push_back(windowKmers.front());
//
//		//If the same k-mer appears several times inside the window and it has the smallest hash we want to save all occurrences
//		while(windowKmers.size() > 1 && windowKmers.front().second == windowKmers[1].second) {
//			windowKmers.erase(windowKmers.begin());
//
//			//Blacklisted k-mers are not added
//			if(!blmers.contains(windowKmers.front().second)) sk.push_back(windowKmers.front());
//		}
//	}
//
//	//Resize sketch (just for case we have allocated way too much memory)
//	sk.shrink_to_fit();
//
//	return sk;
//}

#endif