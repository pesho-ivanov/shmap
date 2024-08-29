#pragma once

#include <climits>
#include <iostream>
#include <string>
#include <vector>

#include "../ext/edlib.h"
#include "io.h"
#include "utils.h"

namespace sweepmap {

// Kmer -- a kmer with metadata (a position in the sequence, hash, strand)
struct Kmer {
	pos_t r;      // kmer resides [l, r), where l+k=r
	hash_t h;
	bool strand;  // false: forward, true: reverse
	Kmer(pos_t r, hash_t h, bool strand) : r(r), h(h), strand(strand) {}
	bool operator==(const Kmer &other) const { return h == other.h; }
};

using sketch_t = std::vector<Kmer>;

// Hit -- a kmer hit in the reference T
struct Hit {  // TODO: compress into a 32bit field
	pos_t r;            // right end of the kmer [l, r), where l+k=r
	pos_t tpos;		    // position in the reference sketch
	bool strand;
	segm_t segm_id;
	Hit() {}
	Hit(const Kmer &kmer, pos_t tpos, segm_t segm_id)
		: r(kmer.r), tpos(tpos), strand(kmer.strand), segm_id(segm_id) {}
};

// Seed -- a kmer with metadata (a position in the queyr P and number of hits in the reference T)
struct Seed {
	Kmer kmer;
	int r_first, r_last;
	int hits_in_T;
	Seed(const Kmer &kmer, pos_t r_first, pos_t r_last, int hits_in_T) :
		kmer(kmer), r_first(r_first), r_last(r_last), hits_in_T(hits_in_T) {}	
};

// Match -- a pair of a seed and a hit
struct Match {
	Seed seed;
	Hit hit;
	int seed_num; // seed number among the chosen seeds, used for indexing the histogram
	Match(const Seed &seed, const Hit &hit, int seed_num)
		: seed(seed), hit(hit), seed_num(seed_num) {}
	
	inline bool is_same_strand() const {
		return seed.kmer.strand == hit.strand;
	}
};

struct RefSegment {
	sketch_t kmers;
	std::string name;
	std::string seq;   // empty if only mapping and no alignment
	int sz;
	RefSegment(const sketch_t &sk, const std::string &name, const std::string &seq, const int sz)
		: kmers(sk), name(name), seq(seq), sz(sz) {}
};

struct Mapping {
	int k; 	   // kmer size
	pos_t P_sz;     // pattern size |P| bp 
	pos_t seeds;     // number of seeds (subset of the sketch kmers)
	pos_t T_l;      // the position of the leftmost nucleotide of the mapping
	pos_t T_r;      // the position of the rightmost nucleotide of the mapping
	segm_t segm_id;
	pos_t s_sz;      // the position of the rightmost nucleotide of the mapping
	int xmin;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	double J, J2;     // Jaccard similarity in [0;1] for the best and for the second best mapping
	double map_time;
	int mapq;
	char strand;    // '+' or '-'
	bool unreasonable;  // reserved for filtering matches
	std::vector<Match>::const_iterator l, r;

    Mapping() {}
	Mapping(int k, pos_t P_sz, int seeds, pos_t T_l, pos_t T_r, segm_t segm_id, pos_t s_sz, int xmin, int same_strand_seeds, std::vector<Match>::const_iterator l, std::vector<Match>::const_iterator r)
		: k(k), P_sz(P_sz), seeds(seeds), T_l(T_l), T_r(T_r), segm_id(segm_id), s_sz(s_sz), xmin(xmin), J(double(xmin) / std::max(seeds, s_sz)), mapq(255), strand(same_strand_seeds > 0 ? '+' : '-'), unreasonable(false), l(l), r(r) {}

	// --- https://github.com/lh3/miniasm/blob/master/PAF.md ---
    void print_paf(const string &query_id, const RefSegment &segm, std::vector<Match> matches) const {
		int P_start = P_sz, P_end = -1;
		for (auto m = l; m != r; ++m) {
			P_start = std::min(P_start, m->seed.r_first);
			P_end = std::max(P_end, m->seed.r_last);
		}
		if (!(0 <= P_start && P_start <= P_end && P_end <= P_sz))
			std::cerr << "P_start=" << P_start << " P_end=" << P_end << " P_sz=" << P_sz << std::endl;
		assert(0 <= P_start && P_start <= P_end && P_end <= P_sz);

		auto T_l_predicted = std::max(T_l-P_start, 0);  // -P_start -- P_start too big
		auto T_r_predicted = std::min(T_r+(P_sz-P_end), segm.sz);  // +(P_sz-P_end) too big

		std::cout << query_id  			// Query sequence name
			<< "\t" << P_sz     // query sequence length
			<< "\t" << P_start   // query start (0-based; closed)
			<< "\t" << P_end  // query end (0-based; open)
			<< "\t" << strand   // '+' or '-'
			<< "\t" << segm.name // reference segment name
			<< "\t" << segm.sz // T_sz -- target sequence length
			<< "\t" << T_l_predicted  // target start on original strand (0-based)
			<< "\t" << T_r_predicted  // target start on original strand (0-based)
			<< "\t" << P_sz  // TODO: fix; Number of residue matches (number of nucleotide matches)
			<< "\t" << P_sz  // TODO: fix; Alignment block length: total number of sequence matches, mismatches, insertions and deletions in the alignment
			<< "\t" << mapq  // Mapping quality (0-255; 255 for missing)
// ----- end of required PAF fields -----
			<< "\t" << "k:i:" << k
			<< "\t" << "p:i:" << seeds  // sketches
			<< "\t" << "M:i:" << matches.size() // kmer matches in T
			<< "\t" << "s:i:" << s_sz
			<< "\t" << "I:i:" << xmin  // intersection of `p` and `s` [kmers]
			<< "\t" << "J:f:" << J   // Jaccard similarity [0; 1]
			<< "\t" << "J2:f:" << J2   // second best mapping Jaccard similarity [0; 1]
			<< "\t" << "t:f:" << map_time
			<< endl;
	}

	std::string reverseComplement(const std::string& seq) const {
		std::string revComp;
		revComp.reserve(seq.size());
		for (int i = seq.size() - 1; i >= 0; --i) {
			switch (seq[i]) {
				case 'A':
					revComp.push_back('T');
					break;
				case 'C':
					revComp.push_back('G');
					break;
				case 'G':
					revComp.push_back('C');
					break;
				case 'T':
					revComp.push_back('A');
					break;
				default:
					revComp.push_back('N');
					break;
			}
		}
		return revComp;
	}

    int print_sam(const string &query_id, const RefSegment &segm, const int matches, const char *query, const size_t query_size) const {
		int T_start = std::max(T_l-k, 0);
		int T_end = std::max(T_r, T_l-k+P_sz);
		int T_d = T_end - T_start;
		string s = segm.seq.substr(T_start, T_d);
		if (strand == '-')
		   s = reverseComplement(s);
		assert(T_d >= 0);
		auto max_edit_dist = -1; //10000;
		auto cfg = edlibNewAlignConfig(max_edit_dist, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
		EdlibAlignResult result = edlibAlign(query, query_size, s.c_str(), T_d, cfg);
		assert(result.status == EDLIB_STATUS_OK);
		char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
		//printf("query=%s, s=%s, ", query, s.c_str());
		//printf("ed=%d, cigar=%s\n", result.editDistance, cigar);

		int ed = result.editDistance;
		int flag = 0;
		if (strand == '-') flag |= 0x10;
		std::cout << query_id 				// 1 QNAME String [!-?A-~]{1,254} Query template NAME
			<< "\t" << flag 				// 2 FLAG Int [0, 216 − 1] bitwise FLAG
			<< "\t" << segm.name  			// 3 RNAME String \*|[:rname:∧*=][:rname:]* Reference sequence NAME11
			<< "\t" << T_start+1  			// 4 POS Int [0, 231 − 1] 1-based leftmost mapping POSition
			<< "\t" << mapq  				// 5 MAPQ Int [0, 28 − 1] MAPping Quality
			<< "\t" << cigar  				// 6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
			<< "\t" << "="  				// 7 RNEXT String \*|=|[:rname:∧*=][:rname:]* Reference name of the mate/next read
			<< "\t" << 0  					// 8 PNEXT Int [0, 231 − 1] Position of the mate/next read
			<< "\t" << 0  					// 9 TLEN Int [−231 + 1, 231 − 1] observed Template LENgth
			<< "\t" << query  				// 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
			<< "\t" << "!"  				// 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
			<< "\t" << "NM:i:" << ed  		// edit distance
//			<< "\t" << "q:s:" << query		// query string
//			<< "\t" << "s:s:" << s			// window string in T
			<< endl;

		free(cigar);
		edlibFreeAlignResult(result);
		return ed;
	}
};

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

			if (h < hThres) // optimize to only look at specific bits
				kmers.push_back(Kmer(r, h, strand));
						
			if (r >= (int)s.size()) break;

			h_fw = std::rotl(h_fw, 1) ^ std::rotl(LUT_fw[(int)s[r-k]], k) ^ LUT_fw[(int)s[r]];
			h_rc = std::rotr(h_rc, 1) ^ std::rotr(LUT_rc[(int)s[r-k]], 1) ^ std::rotl(LUT_rc[(int)s[r]], k-1);

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