#ifndef SWEEPMAP_SKETCH_HPP
#define SWEEPMAP_SKETCH_HPP

#include <climits>
#include <iostream>
#include <string>
#include <vector>

#include "utils.h"

namespace sweepmap {

struct Kmer {
	pos_t r;
	hash_t h;
	bool strand;  // false: forward, true: reverse
	Kmer(pos_t r, hash_t h, bool strand) : r(r), h(h), strand(strand) {}
};

struct Hit {  // TODO: compress all in 32bit
	pos_t r;
	bool strand;
	//pos_t pos;
		 //pos(pos),
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
const Sketch buildFMHSketch(const std::string& s, int k, double hFrac) {
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
		// HACK! the lowest differing bit is not expected to correlate much with (h < hThres)
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

} // namespace sweepmap

#endif