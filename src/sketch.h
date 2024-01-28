#ifndef SKETCH_HPP
#define SKETCH_HPP

#include <climits>
#include <string>
#include <iostream>
#include <vector>

#include <ankerl/unordered_dense.h>

#include "utils.h"
#include "io.h"

namespace sweepmap {

using hash_t     = uint64_t;
using pos_t      = int32_t;

struct Kmer {
	pos_t r;
	hash_t h;
	bool strand;  // false: forward, true: reverse
	Kmer(pos_t r, hash_t h, bool strand) : r(r), h(h), strand(strand) {}
};

struct Hit {
	pos_t r;
	bool strand;
	pos_t pos;
	int segm_id;
	Hit(const Kmer &kmer, pos_t pos, int segm_id)
		: r(kmer.r), strand(kmer.strand), pos(pos), segm_id(segm_id) {}
};

using Sketch = std::vector<Kmer>;  // (kmer hash, kmer's left 0-based position)

static hash_t LUT_fw[256], LUT_rc[256];

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

//string s = "ACGTTAG";
//Sketch sk1 = buildFMHSketch(s, 5, 1.0, blmers);
//Sketch sk2 = buildFMHSketch(revComp(s), 5, 1.0, blmers);
//return 0;
// TODO: use either only forward or only reverse
// TODO: accept char*
const Sketch buildFMHSketch(const string& s, int k, double hFrac) {
	Sketch sk;
	sk.reserve((int)(1.1 * (double)s.size() * hFrac));

	if ((int)s.size() < k) return sk;

	hash_t h, h_fw = 0, h_rc = 0;
	hash_t hThres = hash_t(hFrac * double(std::numeric_limits<hash_t>::max()));
	int r;

	for (r=0; r<k; r++) {
		h_fw ^= std::rotl(LUT_fw[(int)s[r]], k-r-1);
		h_rc ^= std::rotl(LUT_rc[(int)s[r]], r);
	}

	while(true) {
		// HACK! the least significant bit is quite random and does not correlate with (h < hThres)
		// TODO: tie break
		const auto first_diff_bit = 1 << std::countr_zero(h_fw ^ h_rc);
		const bool strand         = h_fw & first_diff_bit;
		h = strand ? h_rc : h_fw;

		//cerr << s << " " << r << " " << std::hex << std::setfill('0') << std::setw(16) << h
		//					<< " = " << std::hex << std::setfill('0') << std::setw(16) << h_fw
		//					<< " ^ " << std::hex << std::setfill('0') << std::setw(16) << h_rc << endl;
		if (h < hThres) // optimize to only look at specific bits
			sk.push_back(Kmer(r, h, strand));
					
		if (r >= (int)s.size()) break;

		h_fw = std::rotl(h_fw, 1) ^ std::rotl(LUT_fw[(int)s[r-k]], k) ^ LUT_fw[(int)s[r]];
		h_rc = std::rotr(h_rc, 1) ^ std::rotr(LUT_rc[(int)s[r-k]], 1) ^ std::rotl(LUT_rc[(int)s[r]], k-1);

		++r;
	}

	return sk;
}

struct RefSegment {
	string name;
	int sz;
	RefSegment(const string &name, const int sz) : name(name), sz(sz) {}
};

struct SketchIndex {
	std::vector<RefSegment> T;
	pos_t total_size;  // total size of all segments // TODO: change type to size_t
	const params_t &params;
	ankerl::unordered_dense::map<hash_t, std::vector<Hit>> h2pos;

	void print_kmer_hist() {
		std::vector<int> hist(10, 0);
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

		cerr << std::fixed << std::setprecision(2);
		cerr << "Histogram of " << kmers << " kmers ("
			<< 100.0*double(different_kmers)/double(kmers) << "\% different) covering "
			<< 100.0*double(params.k)*double(kmers)/double(total_size) << "\% of the " << total_size << "nb index" << endl;
		for (size_t i=0; i<hist.size(); ++i)
			if (hist[i] > 0)
				cerr << std::setw(5) << std::right << i << (i<hist.size()-1?" ":"+") << "occ: " << std::setw(9) << std::right << hist[i] << " kmers (" << 100.0*double(hist[i])/double(kmers) << "\%)" << endl;
		cerr << "The most frequent kmer occurs " << max_occ << " times." << endl;
	}

	void apply_blacklist(int blacklist_threshold) {
		for (auto hits: h2pos)
			if (hits.second.size() > (size_t)blacklist_threshold) {
				h2pos.erase(hits.first);  // TODO: use the iterator instead
			}
	}

	void populate_h2pos(const Sketch& sketch, int segm_id) {
		h2pos.reserve(sketch.size());
		for (size_t tpos = 0; tpos < sketch.size(); ++tpos) {
			const Kmer& kmer = sketch[tpos];
			h2pos[kmer.h].push_back(Hit(kmer, tpos, segm_id));
		}
	}

	void add_segment(const Sketch& sketch, const string &name, int T_sz) {
		T.push_back(RefSegment(name, T_sz));
		total_size += T_sz;
		populate_h2pos(sketch, T.size()-1);
	}

	SketchIndex(const params_t &params) : total_size(0), params(params) {}
};

} // namespace sweepmap

#endif