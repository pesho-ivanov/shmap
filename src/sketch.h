#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <climits>
#include <string>
#include <iostream>
#include <vector>
#include <unordered_map>

#include <ankerl/unordered_dense.h>
//#include "phmap.hpp"

#include "utils.h"
#include "io.h"

using namespace std;

using hash_t     = uint64_t;
using pos_t      = int32_t;
using blmers_t   = std::unordered_map<hash_t, char>;

struct kmer_with_pos_t {
	pos_t r;
	hash_t kmer;

	bool operator<(const kmer_with_pos_t& rhs) const {
		return kmer < rhs.kmer;
	}
};

//using abs_hash_t = pair<pos_t, hash_t>;
//using Hit  = pair<pos_t, pos_t>; 

struct Hit {
	pos_t r;
	pos_t pos;
	int segm_id;

	Hit(pos_t r, pos_t pos, int segm_id) : r(r), pos(pos), segm_id(segm_id) {}
};

using Sketch     = vector<kmer_with_pos_t>;  // (kmer hash, kmer's left 0-based position)

using std::rotl;

static hash_t LUT_fw[256], LUT_rc[256];
//static Timer FMH_time;

//void print_sketches(const string &seqID, const Sketch &sks) {
//	cout << seqID << endl;
//	for (const auto& sk : sks) {
//		cout << "  " << sk.r << ", " << sk.kmer << endl;
//	}
//}

void initialize_LUT() {
	// https://gist.github.com/Daniel-Liu-c0deb0t/7078ebca04569068f15507aa856be6e8
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
// TODO: use either only forward or only reverse
const Sketch buildFMHSketch(const string& s, int k, double hFrac, const blmers_t& blmers) {
	// TODO: accept char*
	Sketch sk;
	sk.reserve((int)(1.1 * (double)s.size() * hFrac));

	if ((int)s.size() < k) return sk;

	hash_t h, h_fw = 0, h_rc = 0;
	//hash_t h_fw = 0;
	hash_t hThres = hash_t(hFrac * double(std::numeric_limits<hash_t>::max()));
	int r;

	for (r=0; r<k; r++) {
		h_fw ^= rotl(LUT_fw[(int)s[r]], k-r-1);
		h_rc ^= rotl(LUT_rc[(int)s[r]], r);
	}

	while(true) {
		// HACK! the least significant bit is quite random and does not correlate with (h < hThres)
		h = std::countr_zero(h_fw) < std::countr_zero(h_rc) ? h_fw : h_rc;

		//cerr << s << " " << r << " " << std::hex << std::setfill('0') << std::setw(16) << h
		//					<< " = " << std::hex << std::setfill('0') << std::setw(16) << h_fw
		//					<< " ^ " << std::hex << std::setfill('0') << std::setw(16) << h_rc << endl;
		if (h < hThres) { // optimize to only look at specific bits
			sk.emplace_back(r, h);
		}
					
		if (r >= (int)s.size()) break;

		h_fw = rotl(h_fw, 1) ^ rotl(LUT_fw[(int)s[r-k]], k) ^ LUT_fw[(int)s[r]];
		h_rc = rotr(h_rc, 1) ^ rotr(LUT_rc[(int)s[r-k]], 1) ^ rotl(LUT_rc[(int)s[r]], k-1);

		++r;
	}

	return sk;
}

const Sketch buildFMHSketch_onlyfw(const string& s, int k, double hFrac, int max_sketch_size=-1) {
	Sketch sk;
	sk.reserve((int)(1.1 * (double)s.size() * hFrac));

	if ((int)s.size() < k) return sk;

	hash_t h_fw = 0, h_rc = 0;
	hash_t hThres = hash_t(hFrac * double(std::numeric_limits<hash_t>::max()));
	int r;

	for (r=0; r<k; r++) {
		h_fw ^= rotl(LUT_fw[(int)s[r]], k-r-1);
		h_rc ^= rotl(LUT_rc[(int)s[r]], r);
	}

	while(true) {
		if (std::countr_zero(h_fw) < std::countr_zero(h_rc))
			if (h_fw < hThres) {
				sk.emplace_back(r, h_fw);
				if (max_sketch_size != -1 && (int)sk.size() >= max_sketch_size) break;
			}
					
		if (r >= (int)s.size()) break;

		h_fw = rotl(h_fw, 1) ^ rotl(LUT_fw[(int)s[r-k]], k) ^ LUT_fw[(int)s[r]];
		h_rc = rotr(h_rc, 1) ^ rotr(LUT_rc[(int)s[r-k]], 1) ^ rotl(LUT_rc[(int)s[r]], k-1);

		++r;
	}

	return sk;
}

struct RefSegment {
	string name;
	string seq;  // TODO: change type to char* and store compressed (2 bits per base)
	RefSegment(const string &name, const string &seq) : name(name), seq(seq) {}
};

struct SketchIndex {
	vector<RefSegment> T;
	pos_t total_size;  // total size of all segments // TODO: change type to size_t
	const params_t &params;
	//unordered_map<hash_t, vector<Hit>> h2pos;
	ankerl::unordered_dense::map<hash_t, vector<Hit>> h2pos;
	//gtl::flat_hash_map<hash_t, vector<Hit>> h2pos;
	//gtl::node_hash_map<hash_t, vector<Hit>> h2pos;
	//ankerl::unordered_dense::map<hash_t, Hit> h2singlepos;

	void print_hist() {
		vector<int> hist(10, 0);
		int kmers = 0, different_kmers = 0, max_occ = 0;
		for (const auto& h2p : h2pos) {
			int occ = h2p.second.size();
			kmers += occ;
			++different_kmers;
			if (occ >= (int)hist.size()-1) {
				hist.back() += occ;
				if (occ > max_occ)
					max_occ = occ;
			} else 
				hist[occ] += occ;
		}

		cerr << fixed << setprecision(2);
		cerr << "Histogram of " << kmers << " kmers ("
			<< 100.0*double(different_kmers)/double(kmers) << "\% different) covering "
			<< 100.0*double(params.k)*double(kmers)/double(total_size) << "\% of the " << total_size << "nb index" << endl;
		for (size_t i=0; i<hist.size(); ++i)
			if (hist[i] > 0)
				cerr << setw(5) << right << i << (i<hist.size()-1?" ":"+") << "occ: " << setw(9) << right << hist[i] << " kmers (" << 100.0*double(hist[i])/double(kmers) << "\%)" << endl;
		cerr << "The most frequent kmer occurs " << max_occ << " times." << endl;
	}

	void populate_h2pos(const Sketch& sketch, int segm_id) {
		//print_sketches(name, sketch);
		for (size_t tpos = 0; tpos < sketch.size(); ++tpos) {
			const kmer_with_pos_t& abs_hash = sketch[tpos];
			h2pos[abs_hash.kmer].push_back(Hit(abs_hash.r, tpos, segm_id));
		}
	}

	void add_segment(const Sketch& sketch, const string &name, char *T_segm) {
		T.push_back(RefSegment(name, T_segm));
		total_size += this->T.back().seq.size();
		populate_h2pos(sketch, T.size()-1);
	}

	SketchIndex(const params_t &params, const blmers_t &bLstmers)
		: total_size(0), params(params) {
	}

	//SketchIndex(const string &name, const string &ref, const params_t &params, const blmers_t &bLstmers)
	//	: T_sz((pos_t)ref.size()), name(name), params(params) {
	//	Sketch sketch = buildFMHSketch(ref, params.k, params.hFrac, bLstmers);
	//	populate_h2pos(sketch);
	//}

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

string revCompl(string s) {
    for (auto& c : s) {
		switch (c) {
			case 'A': c = 'T'; break;
			case 'C': c = 'G'; break;
			case 'G': c = 'C'; break;
			case 'T': c = 'A'; break;
			default: break;
		}
	}
    std::reverse(s.begin(), s.end());
    return s;
}

#endif