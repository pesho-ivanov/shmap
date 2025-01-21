#pragma once

#include <climits>
#include <iostream>
#include <string>
#include <vector>
#include <iomanip>
#include "cmath"

#include "../ext/edlib.h"
#include "io.h"
#include "utils.h"
#include "buckets.h"

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
using Seeds = std::vector<Seed>;

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
using Matches = std::vector<Match>;

struct RefSegment {
	sketch_t kmers;
	std::string name;
	std::string seq;   // empty if only mapping and no alignment
	rpos_t sz;
	int id;
	RefSegment(const sketch_t &sk, const std::string &name, const std::string &seq, rpos_t sz, int id)
		: kmers(sk), name(name), seq(seq), sz(sz), id(id) {}
};

// the minimal mapping info for PAF output
struct MappingPAF {
	// local information: depending on the specific mapping
	qpos_t P_start;
	qpos_t P_end;
	char strand;     		// '+' or '-'
	const char *segm_name;	// name of the reference segment
	rpos_t segm_sz;  		// number of nucleotides in the segment
	int segm_id;
	rpos_t T_l;      		// the position of the leftmost nucleotide of the mapping
	rpos_t T_r;      		// the position of the rightmost nucleotide of the mapping
	int mapq;        		// from 0 to 60

	// global info: depending on the read
	const char* query_id;   // the read name
	qpos_t P_sz;     		// pattern size |P| bp 

	MappingPAF() : P_start(-1), P_end(-1), strand('?'), segm_name(""), segm_sz(-1), segm_id(-1), T_l(-1), T_r(-1), mapq(255), query_id(""), P_sz(-1) {}

	MappingPAF(qpos_t P_start, qpos_t P_end, char strand, const char *segm_name, rpos_t segm_sz, int segm_id, rpos_t T_l, rpos_t T_r)
		: P_start(P_start), P_end(P_end), strand(strand), segm_name(segm_name), segm_sz(segm_sz), segm_id(segm_id), T_l(T_l), T_r(T_r), mapq(255), query_id(""), P_sz(-1) {}

	friend std::ostream& operator<<(std::ostream& os, const MappingPAF& paf) {
		if (paf.P_sz == 0) throw std::runtime_error("MappingPAF is not set at time of output.");
		os << std::fixed << std::setprecision(3);
		os  << paf.query_id				// query sequence name
			<< "\t" << paf.P_sz			// query sequence length
			<< "\t" << paf.P_start		// query start (0-based; closed)
			<< "\t" << paf.P_end		// query end (0-based; open)
			<< "\t" << paf.strand		// '+' or '-'
			<< "\t" << paf.segm_name	// reference segment name
			<< "\t" << paf.segm_sz 		// T_sz -- target sequence length
			<< "\t" << paf.T_l  		// target start on original strand (0-based)
			<< "\t" << paf.T_r  		// target end on original strand (0-based)
			<< "\t" << paf.P_sz  		// TODO: fix; Number of residue matches (number of nucleotide matches)
			<< "\t" << paf.P_sz  		// TODO: fix; Alignment block length: total number of sequence matches, mismatches, insertions and deletions in the alignment
			<< "\t" << paf.mapq;  		// mapping quality (0-255; 255 for missing)
		return os;
	}
};

struct LocalMappingStats {
	rpos_t s_sz;      // the position of the rightmost nucleotide of the mapping
	qpos_t intersection;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	double map_time;
	int same_strand_seeds;  // diff between positive or negative matches
	bool unreasonable;  // reserved for filtering matches
	std::vector<Match>::const_iterator l, r;
	double J, J2;     // Similarity in [0;1] for the best and for the second best mapping
	qpos_t p_sz;     // number of seeds (subset of the sketch kmers)
	Buckets::Bucket bucket;			 // the bucket where the mapping is found
	Buckets::Bucket bucket2;  // the bucket of the second best mapping
	qpos_t intersection2; // number of matches in the second best mapping
	double sigmas_diff;  // how many sigmas is the diff between intersection1 and intersection2

	LocalMappingStats() {
		s_sz = -1;
		intersection = -1;
		map_time = -1.0;
		same_strand_seeds = -1;
		unreasonable = false;
		J = -1.0;
		J2 = -1.0;
		p_sz = -1;
		intersection2 = -1;
		sigmas_diff = -1.0;
	}

	friend std::ostream& operator<<(std::ostream& os, const LocalMappingStats& stats) {
		// per mapping stats
		os	<< "\t" << "p:i:" 			<< stats.p_sz  // sketches
			<< "\t" << "s:i:" 			<< stats.s_sz
			<< "\t" << "I:i:" 			<< stats.intersection  // intersection of `p` and `s` [kmers]
			<< "\t" << "I2:i:"			<< stats.intersection2
			<< "\t" << "Idiff:i:"		<< stats.intersection - stats.intersection2
			<< "\t" << "Isdiff:f:"		<< stats.sigmas_diff
			<< "\t" << "J:f:"			<< stats.J   // Similarity [0; 1]
			<< "\t" << "J2:f:"			<< stats.J2   // second best mapping similarity [0; 1]
			<< "\t" << "Jdiff:f:"		<< stats.J - stats.J2   // second best mapping similarity [0; 1]
			<< "\t" << "strand:i:"		<< stats.same_strand_seeds
			<< "\t" << "t:f:"			<< stats.map_time
			<< "\t" << "b:s:"			<< stats.bucket.segm_id << "," << stats.bucket.b
			<< "\t" << "b2:s:"			<< stats.bucket2.segm_id << "," << stats.bucket2.b;
		return os;
	}
};

struct GlobalMappingStats {
	// internal stats
	int k; 	   // kmer size
	qpos_t seeds;               // number of seeds before pruning
	rpos_t max_seed_matches;// number of matches of the most frequent seed
	rpos_t seed_matches;    // number of matches of seeds
	rpos_t total_matches;   // number of matches of all kmers
	double match_inefficiency;   // intersection / total_matches
	rpos_t seeded_buckets;     // the initial number of buckets
	rpos_t final_buckets;   // the number of buckets after pruning with all kmers
	double FPTP;			// false discovery rate: FP / PP
	double gt_J, gt_C, gt_C_bucket;  // Jaccard and Containment index of the ground truth mapping
//	qpos_t ed;         // edit distance
//	qpos_t ed2;         // second best edit distance

	GlobalMappingStats() {
		k = -1;
		seeds = -1;
		total_matches = -1;
		max_seed_matches = -1;
		seed_matches = -1;
		seeded_buckets = -1;
		final_buckets = -1;
		FPTP = -1.0;
		match_inefficiency = -1.0;
		gt_J = -1.0;
		gt_C = -1.0;
		gt_C_bucket = -1.0;
//		ed = -1;
//		ed2 = -1;
	}

	friend std::ostream& operator<<(std::ostream& os, const GlobalMappingStats& stats) {
		os << std::fixed << std::setprecision(3);
		os  << "\t" << "k:i:" 					<< stats.k
			<< "\t" << "gt_J:f:"				<< stats.gt_J
			<< "\t" << "gt_C:f:"				<< stats.gt_C
			<< "\t" << "gt_C_bucket:f:"			<< stats.gt_C_bucket
			<< "\t" << "seeds:i:"				<< stats.seeds
			<< "\t" << "max_seed_matches:i:"	<< stats.max_seed_matches
			<< "\t" << "seed_matches:i:"		<< stats.seed_matches
			<< "\t" << "total_matches:i:"		<< stats.total_matches
			<< "\t" << "match_inefficiency:f:"	<< stats.match_inefficiency
			<< "\t" << "seeded_buckets:i:"		<< stats.seeded_buckets
			<< "\t" << "final_buckets:i:"		<< stats.final_buckets
			<< "\t" << "FPTP:f:"				<< stats.FPTP;
		return os;
	}
};

class Mapping {
public:
	MappingPAF paf;
	LocalMappingStats local_stats;
	std::unique_ptr<GlobalMappingStats> global_stats;

    Mapping() {}
	Mapping(const Mapping &m) : paf(m.paf), local_stats(m.local_stats), global_stats(m.global_stats ? std::make_unique<GlobalMappingStats>(*m.global_stats) : nullptr) {}

	void update(qpos_t P_start, qpos_t P_end, rpos_t T_l, rpos_t T_r, const RefSegment &segm, qpos_t intersection, double new_score,
				int same_strand_seeds, std::vector<Match>::const_iterator l, std::vector<Match>::const_iterator r) {
		if (new_score > score()) {
			paf = MappingPAF(P_start, P_end, same_strand_seeds > 0 ? '+' : '-', segm.name.c_str(), segm.sz, segm.id, T_l, T_r);
			local_stats.s_sz = r->hit.tpos - l->hit.tpos + 1;
			local_stats.same_strand_seeds = same_strand_seeds;
			local_stats.intersection = intersection;
			local_stats.J = new_score;
			local_stats.l = l;
			local_stats.r = r;
		}
	}
	//Mapping(qpos_t P_start, qpos_t P_end, rpos_t T_l, rpos_t T_r, const RefSegment &segm, qpos_t intersection, double J, int same_strand_seeds, std::vector<Match>::const_iterator l, std::vector<Match>::const_iterator r)
	//	: 	paf(MappingPAF(P_start, P_end, same_strand_seeds > 0 ? '+' : '-', segm, T_l, T_r)),
	//		global_stats(nullptr) {
	//	local_stats.s_sz = r->hit.tpos - l->hit.tpos + 1;
	//	local_stats.same_strand_seeds = same_strand_seeds;
	//	local_stats.intersection = intersection;
	//};

	//Mapping& operator=(const Mapping&) = default;

	Mapping& operator=(const Mapping& other) {
		if (this != &other) {
			paf = other.paf;
			local_stats = other.local_stats;
			global_stats = other.global_stats ? std::make_unique<GlobalMappingStats>(*other.global_stats) : nullptr;
		}
		return *this;
	}

	double score() const { return local_stats.J; }
	double score2() const { return local_stats.J2; }
	int segm_id() const { return paf.segm_id; }

	void set_score2(double score2) {
		local_stats.J2 = score2;
	}

	qpos_t intersection() const {
		return local_stats.intersection;
	}

	int calc_mapq(double theta, double min_diff) {
		// minimap2: mapQ = 40 (1-f2/f1) min(1, m/10) log f1, where m is #anchors on primary chain
		assert(score() >= 0.0);
		if (score2() < theta)
			set_score2(theta);
		double diff = score() - score2();
		int mapq_strand = (abs(local_stats.same_strand_seeds) < local_stats.intersection/2) ? 5 : 60;
		if (diff > min_diff) {
			return std::min(60, mapq_strand);
		} else if (diff < min_diff/2) {
			return 0;
		} else {
			diff -= min_diff/2;
			return std::min(int(60*diff/(min_diff/2)), mapq_strand);
		}
	}

	double sigmas_diff(qpos_t X, qpos_t Y) {
		if (X==-1) X=0;
		if (Y==-1) Y=0;
		return std::abs(X - Y) / std::sqrt(X + Y);
	}

	void set_local_stats(const RefSegment &segm) {
		//local_stats.s_sz = segm.sz;
		//local_stats.segm_name = segm.name;
	}

	double mapq() const {
		return paf.mapq;
	}

	Buckets::Bucket bucket() const {
		return local_stats.bucket;
	}

	void set_bucket(const Buckets::Bucket &bucket) {
		local_stats.bucket = bucket;
	}

	void set_second_best(const Mapping &m2) {
		local_stats.bucket2 = m2.local_stats.bucket;
		local_stats.intersection2 = m2.local_stats.intersection;

		//cerr << "best_idx: " << best_idx << ", best2_idx: " << best2_idx << endl;
		set_score2( m2.score() );
		local_stats.bucket2 = m2.local_stats.bucket;
		local_stats.intersection2 = m2.local_stats.intersection;
		local_stats.sigmas_diff = sigmas_diff(local_stats.intersection2, local_stats.intersection);
		//assert(m.bucket.segm_id != m.bucket2.segm_id || abs(m.bucket.b - m.bucket2.b) > 1);
		//m.ed2 = edit_distance(m.bucket, P, P_sz, m, kmers);
	}

	void set_global_stats(
			double theta, double min_diff, qpos_t p_sz,
			const char* query_id, qpos_t P_sz, qpos_t k, qpos_t seeds, rpos_t total_matches, rpos_t max_seed_matches, rpos_t seed_matches,
			rpos_t seeded_buckets, rpos_t final_buckets, double FPTP, double map_time) {
		paf.query_id = query_id;
		paf.P_sz = P_sz;
		paf.mapq = calc_mapq(theta, min_diff);

		local_stats.p_sz = p_sz;
		local_stats.map_time = map_time;

		global_stats = std::make_unique<GlobalMappingStats>();
		global_stats->k = k;
		global_stats->seeds = seeds;
		global_stats->total_matches = total_matches;
		global_stats->max_seed_matches = max_seed_matches;
		global_stats->seed_matches = seed_matches;
		global_stats->seeded_buckets = seeded_buckets;
		global_stats->final_buckets = final_buckets;
		global_stats->FPTP = FPTP;
		global_stats->gt_J = -1.0;
		global_stats->gt_C = -1.0;
		global_stats->gt_C_bucket = -1.0;
		global_stats->match_inefficiency = 1.0 * total_matches / local_stats.intersection;
	}

	static bool do_overlap(const Mapping &a, const Mapping &b) {
		if (a.paf.segm_id != b.paf.segm_id)
			return false;
		return a.paf.T_r >= b.paf.T_l && a.paf.T_l <= b.paf.T_r;
	}

	static double overlap(const Mapping &a, const Mapping &b) {
		if (a.paf.segm_id != b.paf.segm_id)
			return -0.0;
		int cap = std::max(0, std::min(a.paf.T_r, b.paf.T_r) - std::max(a.paf.T_l, b.paf.T_l));
		int cup = std::max(a.paf.T_r, b.paf.T_r) - std::min(a.paf.T_l, b.paf.T_l);
		assert(cup >= 0 && cap >=0 && cup >= cap);
		return 1.0 * cap / cup;
	}
		//if (J > 1.0)
		//	cerr << "s_sz = " << s_sz << ", intersection = " << intersection << ", p_sz = " << p_sz << ", J = " << J << ", r-l=" << r-l << endl;
		//assert(J >= -0.0);
		//assert(J <= 1.0);
		//query_id = "";
		//segm_name = "";
		//segm_sz = -1;
		//s_sz = r->hit.tpos - l->hit.tpos + 1;
		//J2 = -1.0;
		//mapq = 255;
		//strand = same_strand_seeds > 0 ? '+' : '-';
		//P_start = 0;
		//P_end = P_sz-1;
	//}

	// --- https://github.com/lh3/miniasm/blob/master/PAF.md ---
    void print_paf(std::ostream& os) const {
		os << *this;
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

	friend std::ostream& operator<<(std::ostream& os, const Mapping& mapping) {
		os << mapping.paf;
		os << mapping.local_stats;
		if (mapping.global_stats != nullptr) os << *mapping.global_stats;
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

		hash_t h_fw = 0, h_rc = 0;
		hash_t hThres = hFrac < 1.0 ? hash_t(hFrac * double(std::numeric_limits<hash_t>::max())) : std::numeric_limits<hash_t>::max();
		rpos_t r;

		for (r=0; r<k; r++) {
			h_fw ^= std::rotl(LUT_fw[ size_t(s[r]) ], k-r-1);
			h_rc ^= std::rotl(LUT_rc[ size_t(s[r]) ], r);
		}

		while(true) {
			// HACK! the lowest differing bit is not expected to correlate much with (h < hThres)
			// TODO: tie break
			auto h = h_rc ^ h_fw;
			if (h <= hThres) {
				//if (h_fw == h_rc) cerr << "h_fw == h_rc == " << h_fw << endl;
				const bool strand = h_fw > h_rc;
				kmers.push_back(Kmer(r-1, h, strand));
			}
						
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