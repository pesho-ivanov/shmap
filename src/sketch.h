#pragma once

#include <climits>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>

#include "../ext/edlib.h"
#include "io.h"
#include "utils.h"

namespace sweepmap {

// Kmer -- a kmer with metadata (a position in the sequence, hash, strand)
struct Kmer {
	rpos_t r;      // kmer resides [l, r), where l+k=r
	hash_t h;
	bool strand;  // false: forward, true: reverse
	Kmer(rpos_t r, hash_t h, bool strand) : r(r), h(h), strand(strand) {}
	bool operator==(const Kmer &other) const { return h == other.h; }
	friend std::ostream& operator<<(std::ostream& os, const Kmer& kmer) {
		os << "Kmer(r=" << kmer.r << ", h=" << kmer.h << ", strand=" << kmer.strand << ")";
		return os;
	}
};

using sketch_t = std::vector<Kmer>;

// Hit -- a kmer hit in the reference T
struct Hit {  // TODO: compress into a 32bit field
	rpos_t r;            // right end of the kmer [l, r), where l+k=r
	rpos_t tpos;		    // position in the reference sketch
	bool strand;
	segm_t segm_id;
	Hit() {}
	Hit(const Kmer &kmer, rpos_t tpos, segm_t segm_id)
		: r(kmer.r), tpos(tpos), strand(kmer.strand), segm_id(segm_id) {}
	friend std::ostream& operator<<(std::ostream& os, const Hit& hit) {
		os << "Hit(r=" << hit.r << ", tpos=" << hit.tpos << ", strand=" << hit.strand << ", segm_id=" << hit.segm_id << ")";
		return os;
	}
};

// Seed -- a kmer with metadata (a position in the queyr P and number of hits in the reference T)
struct Seed {
	Kmer kmer;
	rpos_t hits_in_T;
	qpos_t occs_in_p;
	qpos_t seed_num;
	Seed(const Kmer &kmer, rpos_t hits_in_T, qpos_t seed_num) :
		kmer(kmer), hits_in_T(hits_in_T), seed_num(seed_num) {}	
	friend std::ostream& operator<<(std::ostream& os, const Seed& seed) {
		os << "Seed(" << seed.kmer << ", hits_in_T=" << seed.hits_in_T << ", seed_num=" << seed.seed_num << ")";
		return os;
	}
};

// Match -- a pair of a seed and a hit
struct Match {
	Seed seed;
	Hit hit;
	Match(const Seed &seed, const Hit &hit)
		: seed(seed), hit(hit) {}
	
	inline bool is_same_strand() const {
		return seed.kmer.strand == hit.strand;
	}
	friend std::ostream& operator<<(std::ostream& os, const Match& match) {
		os << "Match(" << match.seed << ", " << match.hit << ")";
		return os;
	}
};

struct RefSegment {
	sketch_t kmers;
	std::string name;
	std::string seq;   // empty if only mapping and no alignment
	rpos_t sz;
	RefSegment(const sketch_t &sk, const std::string &name, const std::string &seq, const rpos_t sz)
		: kmers(sk), name(name), seq(seq), sz(sz) {}
};

struct bucket_t {
	segm_t segm_id;		// segm_id refers to tidx.T[segm_id]
	rpos_t b;   		// b refers to interval [lmax*b, lmax*(b+2)) in tidx.T[segm_id].kmers
	bucket_t() : segm_id(-1), b(-1) {}
	bucket_t(segm_t segm_id, rpos_t b) : segm_id(segm_id), b(b) {}
    bool operator==(const bucket_t &other) const {
        return segm_id == other.segm_id && b == other.b;
    }
	friend std::ostream& operator<<(std::ostream& os, const bucket_t& b) {
		os << "(" << b.segm_id << "," << b.b << ")";
		return os;
	}
};

} // namespace sweepmap

template <>
struct std::hash<sweepmap::bucket_t> {
	std::size_t operator()(const sweepmap::bucket_t& b) const {
		return std::hash<sweepmap::segm_t>()(b.segm_id) ^ (std::hash<sweepmap::rpos_t>()(b.b) << 1);
	}
};

namespace sweepmap {

struct Mapping {
	int k; 	   // kmer size
	const char* query_id;
	qpos_t P_start;
	qpos_t P_end;
	qpos_t P_sz;     // pattern size |P| bp 
	qpos_t p_sz;     // number of seeds (subset of the sketch kmers)
	rpos_t T_l;      // the position of the leftmost nucleotide of the mapping
	rpos_t T_r;      // the position of the rightmost nucleotide of the mapping
	segm_t segm_id;
	const char* segm_name;
	segm_t segm_sz;
	rpos_t s_sz;      // the position of the rightmost nucleotide of the mapping
	qpos_t intersection;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	double J, J2;     // Similarity in [0;1] for the best and for the second best mapping
	double map_time;
	int mapq;
	int same_strand_seeds;  // diff between positive or negative matches
	char strand;    // '+' or '-'
	bool unreasonable;  // reserved for filtering matches
	std::vector<Match>::const_iterator l, r;
//	qpos_t ed;         // edit distance
//	qpos_t ed2;         // second best edit distance

	// internal stats
	qpos_t seeds;               // number of seeds before pruning
	rpos_t max_seed_matches;// number of matches of the most frequent seed
	rpos_t seed_matches;    // number of matches of seeds
	rpos_t total_matches;   // number of matches of all kmers
	double match_inefficiency;   // intersection / total_matches
	rpos_t seeded_buckets;     // the initial number of buckets
	rpos_t final_buckets;   // the number of buckets after pruning with all kmers
	double FPTP;			// false discovery rate: FP / PP
	bucket_t bucket;			 // the bucket where the mapping is found
	bucket_t bucket2;  // the bucket of the second best mapping
	qpos_t intersection2; // number of matches in the second best mapping
	double sigmas_diff;  // how many sigmas is the diff between intersection1 and intersection2
	double gt_J, gt_C, gt_C_bucket;  // Jaccard and Containment index of the ground truth mapping

    Mapping() : J(-1.0) {}
	Mapping(int k, qpos_t P_sz, qpos_t p_sz, rpos_t T_l, rpos_t T_r, segm_t segm_id, qpos_t intersection, double J, int same_strand_seeds, std::vector<Match>::const_iterator l, std::vector<Match>::const_iterator r, bucket_t bucket)
		: k(k), P_sz(P_sz), p_sz(p_sz), T_l(T_l), T_r(T_r), segm_id(segm_id), intersection(intersection), J(J), same_strand_seeds(same_strand_seeds), unreasonable(false), l(l), r(r), bucket(bucket) {
		//if (J > 1.0)
		//	cerr << "s_sz = " << s_sz << ", intersection = " << intersection << ", p_sz = " << p_sz << ", J = " << J << ", r-l=" << r-l << endl;
		assert(J >= -0.0);
		assert(J <= 1.0);
		//query_id = "";
		//segm_name = "";
		segm_sz = -1;
		s_sz = r->hit.tpos - l->hit.tpos + 1;
		mapq = 255;
//		ed = -1;
//		ed2 = -1;
		strand = same_strand_seeds > 0 ? '+' : '-';
		P_start = 0;
		P_end = P_sz-1;

		// stats
		seeds = -1;
		total_matches = -1;
		max_seed_matches = -1;
		seed_matches = -1;
		seeded_buckets = -1;
		final_buckets = -1;
		FPTP = -1.0;
		match_inefficiency = -1.0;
		J2 = -1.0;
		intersection2 = -1;
		sigmas_diff = -1.0;
	}

	// --- https://github.com/lh3/miniasm/blob/master/PAF.md ---
    void print_paf() const {
		std::cout << *this;
	}

	static std::string reverseComplement(const std::string& seq) {
		std::string revComp;
		revComp.reserve(seq.size());
		for (auto i = rpos_t(seq.size()) - 1; i >= 0; --i) {
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

//    int print_sam(const string &query_id, const RefSegment &segm, const rpos_t matches, const char *query, const size_t query_size) const {
//		rpos_t T_start = std::max(T_l-k, rpos_t(0));
//		rpos_t T_end = std::max(T_r, T_l-k+P_sz);
//		rpos_t T_d = T_end - T_start;
//		string s = segm.seq.substr(T_start, T_d);
//		if (strand == '-')
//		   s = reverseComplement(s);
//		assert(T_d >= 0);
//		auto max_edit_dist = -1; //10000;
//		auto cfg = edlibNewAlignConfig(max_edit_dist, EDLIB_MODE_NW, EDLIB_TASK_PATH, NULL, 0);
//		EdlibAlignResult result = edlibAlign(query, query_size, s.c_str(), T_d, cfg);
//		assert(result.status == EDLIB_STATUS_OK);
//		char* cigar = edlibAlignmentToCigar(result.alignment, result.alignmentLength, EDLIB_CIGAR_STANDARD);
//		//printf("query=%s, s=%s, ", query, s.c_str());
//		//printf("ed=%d, cigar=%s\n", result.editDistance, cigar);
//
//		rpos_t ed = result.editDistance;
//		int flag = 0;
//		if (strand == '-') flag |= 0x10;
//		std::cout << query_id 				// 1 QNAME String [!-?A-~]{1,254} Query template NAME
//			<< "\t" << flag 				// 2 FLAG Int [0, 216 − 1] bitwise FLAG
//			<< "\t" << segm.name  			// 3 RNAME String \*|[:rname:∧*=][:rname:]* Reference sequence NAME11
//			<< "\t" << T_start+1  			// 4 POS Int [0, 231 − 1] 1-based leftmost mapping POSition
//			<< "\t" << mapq  				// 5 MAPQ Int [0, 28 − 1] MAPping Quality
//			<< "\t" << cigar  				// 6 CIGAR String \*|([0-9]+[MIDNSHPX=])+ CIGAR string
//			<< "\t" << "="  				// 7 RNEXT String \*|=|[:rname:∧*=][:rname:]* Reference name of the mate/next read
//			<< "\t" << 0  					// 8 PNEXT Int [0, 231 − 1] Position of the mate/next read
//			<< "\t" << 0  					// 9 TLEN Int [−231 + 1, 231 − 1] observed Template LENgth
//			<< "\t" << query  				// 10 SEQ String \*|[A-Za-z=.]+ segment SEQuence
//			<< "\t" << "!"  				// 11 QUAL String [!-~]+ ASCII of Phred-scaled base QUALity+33
//			<< "\t" << "NM:i:" << ed  		// edit distance
////			<< "\t" << "q:s:" << query		// query string
////			<< "\t" << "s:s:" << s			// window string in T
//			<< endl;
//
//		free(cigar);
//		edlibFreeAlignResult(result);
//		return ed;
//	}

	friend std::ostream& operator<<(std::ostream& os, const Mapping& mapping) {
		os << std::fixed << std::setprecision(3);
		os 
			<< mapping.query_id  	// Query sequence name
			<< "\t" << mapping.P_sz     // query sequence length
			<< "\t" << mapping.P_start   // query start (0-based; closed)
			<< "\t" << mapping.P_end  // query end (0-based; open)
			<< "\t" << mapping.strand   // '+' or '-'
			<< "\t" << mapping.segm_name // reference segment name
			<< "\t" << mapping.segm_sz // T_sz -- target sequence length
			<< "\t" << mapping.T_l  // target start on original strand (0-based)
			<< "\t" << mapping.T_r  // target end on original strand (0-based)
			<< "\t" << mapping.P_sz  // TODO: fix; Number of residue matches (number of nucleotide matches)
			<< "\t" << mapping.P_sz  // TODO: fix; Alignment block length: total number of sequence matches, mismatches, insertions and deletions in the alignment
			<< "\t" << mapping.mapq  // Mapping quality (0-255; 255 for missing)
// ----- end of required PAF fields -----
			<< "\t" << "k:i:" 			<< mapping.k
			<< "\t" << "p:i:" 			<< mapping.p_sz  // sketches
			<< "\t" << "s:i:" 			<< mapping.s_sz
			<< "\t" << "I:i:" 			<< mapping.intersection  // intersection of `p` and `s` [kmers]
			<< "\t" << "I2:i:"			<< mapping.intersection2
			<< "\t" << "Idiff:i:"		<< mapping.intersection - mapping.intersection2
			<< "\t" << "Isdiff:f:"		<< mapping.sigmas_diff
			<< "\t" << "J:f:"			<< mapping.J   // Similarity [0; 1]
			<< "\t" << "J2:f:"			<< mapping.J2   // second best mapping similarity [0; 1]
			<< "\t" << "Jdiff:f:"		<< mapping.J - mapping.J2   // second best mapping similarity [0; 1]
			<< "\t" << "gt_J:f:"		<< mapping.gt_J
			<< "\t" << "gt_C:f:"		<< mapping.gt_C
			<< "\t" << "gt_C_bucket:f:"	<< mapping.gt_C_bucket
//			<< "\t" << "ed:i:"			<< mapping.ed
//			<< "\t" << "ed2:i:"			<< mapping.ed2
			<< "\t" << "Seeds:i:"		<< mapping.seeds
			<< "\t" << "MSeedmax:i:"	<< mapping.max_seed_matches
			<< "\t" << "MSeed:i:"		<< mapping.seed_matches
			<< "\t" << "M:i:"			<< mapping.total_matches
			<< "\t" << "Mineff:f:"		<< mapping.match_inefficiency
			<< "\t" << "Bmax:i:"		<< mapping.seeded_buckets
			<< "\t" << "Bfinal:i:"		<< mapping.final_buckets
			<< "\t" << "FPTP:f:"		<< mapping.FPTP
			<< "\t" << "b:s:"			<< mapping.bucket.segm_id << "," << mapping.bucket.b
			<< "\t" << "b2:s:"			<< mapping.bucket2.segm_id << "," << mapping.bucket2.b
			<< "\t" << "strand:i:"		<< mapping.same_strand_seeds
			<< "\t" << "t:f:"			<< mapping.map_time
			<< endl;
		return os;
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