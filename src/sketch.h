#pragma once

#include <climits>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "cmath"

#include "utils.h"
#include "types.h"

namespace sweepmap {

class FracMinHash {
public:
	hash_t LUT_fw[256], LUT_rc[256];
	int k;
	double hFrac;
	Counters *C;
	Timers *T;

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

	// TODO: use either only forward or only reverse
	// TODO: accept char*
	const sketch_t sketch(const std::string& s) const {
		sketch_t kmers;
		kmers.reserve((rpos_t)(1.1 * (double)s.size() * hFrac));

		if (rpos_t(s.size()) < k) return kmers;

		hash_t h, h_fw = 0, h_rc = 0;
		hash_t hThres = hFrac < 1.0 ? hash_t(hFrac * double(std::numeric_limits<hash_t>::max())) : std::numeric_limits<hash_t>::max();
		rpos_t r;

		for (r=0; r<k; r++) {
			h_fw ^= std::rotl(LUT_fw[ size_t(s[r]) ], k-r-1);
			h_rc ^= std::rotl(LUT_rc[ size_t(s[r]) ], r);
		}

		while(true) {
			// HACK! the lowest differing bit is not expected to correlate much with (h < hThres)
			// TODO: tie break
			const auto first_diff_bit = 1 << std::countr_zero(h_fw ^ h_rc);
			const bool strand         = h_fw & first_diff_bit;
			h = strand ? h_rc : h_fw;

			if (h <= hThres) // optimize to only look at specific bits
				kmers.push_back(Kmer(r-1, h, strand));
						
			if (r >= rpos_t(s.size())) break;

			h_fw = std::rotl(h_fw, 1) ^ std::rotl(LUT_fw[ size_t(s[r-k]) ], k) ^ LUT_fw[ size_t(s[r]) ];
			h_rc = std::rotr(h_rc, 1) ^ std::rotr(LUT_rc[ size_t(s[r-k]) ], 1) ^ std::rotl(LUT_rc[ size_t(s[r]) ], k-1);

			++r;
		}

		C->inc("sketched_seqs");
		C->inc("sketched_len", s.size());
		C->inc("original_kmers", kmers.size());
		C->inc("sketched_kmers", kmers.size());

		return kmers;
	}

	FracMinHash(int k, double hFrac, Counters *C, Timers *T)
			: k(k), hFrac(hFrac), C(C), T(T) {
        initialize_LUT();
	}
};

} // namespace sweepmap
