#pragma once

#include <climits>
#include <iostream>
#include <string>
#include <vector>

#include "io.h"
#include "utils.h"

namespace sweepmap {

struct Kmer {
	pos_t r;      // kmer resides [l, r), where l+k=r
	hash_t h;
	bool strand;  // false: forward, true: reverse
	bool black;   // if true, the kmer is blacklisted
	Kmer(pos_t r, hash_t h, bool strand) : r(r), h(h), strand(strand) {}
};

class Sketch {
public:
	using sketch_t = std::vector<Kmer>;

	inline static hash_t LUT_fw[256], LUT_rc[256];
	inline static params_t *params;
	inline static Timers *T;
	inline static Counters *C;

	static void initialize_LUT() {
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

private:
	sketch_t addPairKmers(const sketch_t &kmers, double hFrac) {
		hash_t hThres = hash_t(hFrac * double(std::numeric_limits<hash_t>::max()));
		sketch_t kmers2;
		kmers2.reserve(2*kmers.size());
		for (int i=0; i<(int)kmers.size(); ++i) {
			kmers2.push_back(kmers[i]);
			for (int j=i+1; j<std::min(i+20, (int)kmers.size()); ++j) {
				if (kmers[i].strand == kmers[j].strand) {
					//if (sk[j].r - sk[i].r > 0 && sk[j].r - sk[i].r < 5000) {
					if ((kmers[i].h^kmers[j].h) < hash_t(0.15*hThres)) {
						kmers2.push_back(Kmer(kmers[i].r, kmers[i].h^kmers[j].h, kmers[i].strand));
						//break;
					}
					//}
				}
			}
		}
		C->inc("paired_kmers", kmers2.size() - kmers.size());
		return kmers2;
	}
	
//	const sketch_t buildMinimizerSketch(const std::string &s, int k, int w) {
//		sketch_t kmers;
//		int r, min_r;
//		hash_t min_h, h, h_fw = 0, h_rc = 0;
//		std::deque<Kmer> W;  // kmers in [r-w; r) in increasing order of r
//		for (r=0; r<k; r++) {
//			h_fw ^= std::rotl(LUT_fw[(int)s[r]], k-r-1);
//			h_rc ^= std::rotl(LUT_rc[(int)s[r]], r);
//		}
//		while (true) {
//			const auto first_diff_bit = 1 << std::countr_zero(h_fw ^ h_rc);
//			const bool strand         = h_fw & first_diff_bit;
//			h = strand ? h_rc : h_fw;
//
//			W.push_back(Kmer(r, h, strand));
//
//			if (min_r + w < r) {
//				kmers.push_back(Kmer(r, h, strand));
//				
//				min_r = r;
//			}
//
//			if (r >= (int)s.size()) break;
//
//			h_fw = std::rotl(h_fw, 1) ^ std::rotl(LUT_fw[(int)s[r-k]], k) ^ LUT_fw[(int)s[r]];
//			h_rc = std::rotr(h_rc, 1) ^ std::rotr(LUT_rc[(int)s[r-k]], 1) ^ std::rotl(LUT_rc[(int)s[r]], k-1);
//
//			++r;
//		}
//		return kmers;
//	}
//
	//string s = "ACGTTAG";
	//Sketch sk1 = buildFMHSketch(s, 5, 1.0, blmers);
	//Sketch sk2 = buildFMHSketch(revComp(s), 5, 1.0, blmers);
	//return 0;
	// TODO: use either only forward or only reverse
	// TODO: accept char*
	const sketch_t buildFMHSketch(const std::string& s, int k, double hFrac) {
		sketch_t kmers;
		kmers.reserve((int)(1.1 * (double)s.size() * hFrac));

		if ((int)s.size() < k) return kmers;

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
				kmers.push_back(Kmer(r, h, strand));
						
			if (r >= (int)s.size()) break;

			h_fw = std::rotl(h_fw, 1) ^ std::rotl(LUT_fw[(int)s[r-k]], k) ^ LUT_fw[(int)s[r]];
			h_rc = std::rotr(h_rc, 1) ^ std::rotr(LUT_rc[(int)s[r-k]], 1) ^ std::rotl(LUT_rc[(int)s[r]], k-1);

			++r;
		}

		return kmers;
	}

public:
	sketch_t kmers;   // (kmer hash, kmer's left 0-based position)

	Sketch(const std::string& s) {
		kmers = buildFMHSketch(s, params->k, params->hFrac);
		C->inc("sketched_seqs");
		C->inc("sketched_len", s.size());
		C->inc("original_kmers", kmers.size());
		C->inc("paired_kmers", 0);
		if (params->paired_kmers)
			kmers = addPairKmers(kmers, params->hFrac);
		C->inc("sketched_kmers", kmers.size());
	}

	static void print_stats() {
		cerr << std::fixed << std::setprecision(1);
		cerr << "Sketching:" << endl;
		cerr << " | Sketched sequences:    " << C->count("sketched_seqs") << " (" << C->count("sketched_len") << " nb)" << endl;
		cerr << " | Kmers:                 " << C->count("sketched_kmers") << endl;
		cerr << " |  | original:               " << C->count("original_kmers") << " (" << C->perc("original_kmers", "sketched_kmers") << "%)" << endl;
		cerr << " |  | paired:                 " << C->count("paired_kmers") << " (" << C->perc("paired_kmers", "sketched_kmers") << "%)" << endl;
	}
};

} // namespace sweepmap