#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <bit>
#include <cstdint>
#include <utility>
#include <list>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>
#include <math.h>

#include <ankerl/unordered_dense.h>
//#include "phmap.hpp"

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

	LUT_rc['a'] = LUT_rc['A'] = LUT_fw['T'];
	LUT_rc['c'] = LUT_rc['C'] = LUT_fw['G'];
	LUT_rc['g'] = LUT_rc['G'] = LUT_fw['C'];
	LUT_rc['t'] = LUT_rc['T'] = LUT_fw['A'];
}

//blmers_t blmers;
//string s = "ACGTTAG";
//Sketch sk1 = buildFMHSketch(s, 5, 1.0, blmers);
//Sketch sk2 = buildFMHSketch(revComp(s), 5, 1.0, blmers);
//return 0;
const Sketch buildFMHSketch(const string& s, int k, double hFrac, const blmers_t& blmers) {
	FMH_time.start();

	Sketch sk;
	sk.reserve((int)(1.1 * (double)s.size() * hFrac));

	if ((int)s.size() < k) return sk;

	//hash_t h, h_fw = 0, h_rc = 0;
	hash_t h_fw = 0;
	hash_t hThres = hash_t(hFrac * double(std::numeric_limits<hash_t>::max()));
	int r;

	for (r=0; r<k; r++) {
		h_fw ^= rotl(LUT_fw[(int)s[r]], k-r-1);
		//h_rc ^= rotl(LUT_rc[(int)s[r]], r);
	}

	while(true) {
		//h = h_fw ^ h_rc;
//		cerr << s << " " << r << " " << std::hex << std::setfill('0') << std::setw(16) << h
//							<< " = " << std::hex << std::setfill('0') << std::setw(16) << h_fw
//							<< " ^ " << std::hex << std::setfill('0') << std::setw(16) << h_rc << endl;
		if (h_fw < hThres)  // optimize to only look at specific bits
			//sk.emplace_back(r, h);
			sk.emplace_back(r, h_fw);
		
		if (r >= (int)s.size()) break;

		h_fw = rotl(h_fw, 1) ^ rotl(LUT_fw[(int)s[r-k]], k) ^ LUT_fw[(int)s[r]];
		//h_rc = rotr(h_rc, 1) ^ rotr(LUT_rc[(int)s[r-k]], 1) ^ rotl(LUT_rc[(int)s[r]], k-1);

		++r;
	}

	FMH_time.stop();
	return sk;
}

struct kmer_hits_t {
	pos_t P_l;
	kmer_num_t kmer_num;
	const vector<abs_ord_t> *kmers_in_T;

	kmer_hits_t() {}
	kmer_hits_t(pos_t P_l, kmer_num_t kmer_num, const vector<abs_ord_t> *kmers_in_T) :
		P_l(P_l), kmer_num(kmer_num), kmers_in_T(kmers_in_T) {}
};

//struct EmptyHash {
//    size_t operator()(hash_t key) const {
//        return key; // ^ (key >> 32);
//    }
//};

struct SketchIndex {
	pos_t T_sz;
	string name;
	//unordered_map<hash_t, vector<abs_ord_t>> h2pos;
	ankerl::unordered_dense::map<hash_t, vector<abs_ord_t>> h2pos;
	//gtl::flat_hash_map<hash_t, vector<abs_ord_t>> h2pos;
	//gtl::node_hash_map<hash_t, vector<abs_ord_t>> h2pos;

	void populate_h2pos(const Sketch& sketch) {
		//print_sketches(name, sketch);
		for (kmer_num_t tpos = 0; tpos < (int)sketch.size(); ++tpos) {
			const abs_hash_t& abs_hash = sketch[tpos];
			h2pos[abs_hash.second].push_back(abs_ord_t(abs_hash.first, tpos));
		}
	}

	SketchIndex(const Sketch& sketch, pos_t T_sz, const string &name)
		: T_sz(T_sz), name(name) {
		populate_h2pos(sketch);
	}

	SketchIndex(const string &name, const string &ref, const params_t &params, const blmers_t &bLstmers)
		: T_sz((pos_t)ref.size()), name(name) {
		Sketch sketch = buildFMHSketch(ref, params.k, params.hFrac, bLstmers);
		populate_h2pos(sketch);
	}

	//SketchIndex(const string &tFile, const params_t &params, const blmers_t &bLstmers) {
	//	read_fasta_klib(tFile, [this, &params, &bLstmers](kseq_t *seq) {
	//		assert(name.empty());  // TODO: support multi-sequence files
	//		name = seq->name.s;
	//		T_sz = seq->seq.l;
	//		Sketch sketch = buildFMHSketch(seq->seq.s, params.k, params.hFrac, bLstmers);
	//		populate_h2pos(sketch);
	//	});
	//}
};

#endif