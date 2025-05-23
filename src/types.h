#pragma once

#include <fstream>
#include "gtl/vector.hpp"
#include "unordered_dense/include/ankerl/unordered_dense.h"

namespace sweepmap {

using hash_t     = uint64_t;
using rpos_t     = int32_t;  // reference position; klib anyway doesn't support 64-bit
using qpos_t     = int32_t;  // query position
using segm_t     = int32_t;

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
	std::vector<qpos_t> pmatches;   // positions in `p' of all occurences of `kmer'; sorted in decreasing order
	Seed(const Kmer &kmer, rpos_t hits_in_T, qpos_t occs_in_p, qpos_t seed_num, const std::vector<qpos_t> &pmatches) :
		kmer(kmer), hits_in_T(hits_in_T), occs_in_p(occs_in_p), seed_num(seed_num), pmatches(pmatches) {}	
	friend std::ostream& operator<<(std::ostream& os, const Seed& seed) {
		os << "Seed(" << seed.kmer << ", hits_in_T=" << seed.hits_in_T << ", occs_in_p=" << seed.occs_in_p << ", seed_num=" << seed.seed_num << ")";
		return os;
	}
};

// Match -- a pair of a seed and a hit
struct Match {
	Seed seed;
	Hit hit;
	Match(const Seed &seed, const Hit &hit)
		: seed(seed), hit(hit) {}
	
	inline int codirection() const {
		return seed.kmer.strand == hit.strand ? +1 : -1;
	}
	friend std::ostream& operator<<(std::ostream& os, const Match& match) {
		os << "Match(" << match.seed << ", " << match.hit << ")";
		return os;
	}
};

struct BucketLoc {
	segm_t segm_id;		// segm_id refers to tidx.T[segm_id]
	rpos_t b;   		// b refers to interval [lmax*b, lmax*(b+2)) in tidx.T[segm_id].kmers

	BucketLoc() : segm_id(-1), b(-1) {}
	BucketLoc(segm_t segm_id, rpos_t b) : segm_id(segm_id), b(b) {}

	bool operator==(const BucketLoc &other) const {
		return segm_id == other.segm_id && b == other.b;
	}
	friend std::ostream& operator<<(std::ostream& os, const BucketLoc& b) {
		os << "(" << b.segm_id << "," << b.b << ")";
		return os;
	}
};
struct BucketHash {
	std::size_t operator()(const BucketLoc& b) const {
		return std::hash<segm_t>()(b.segm_id) ^ (std::hash<rpos_t>()(b.b) << 1);
	}
};

//struct LightMatch {
//	qpos_t r;  // TODO: shrink to [0, 2*bucket_halflen)
//	qpos_t seed_num;	// [0, |p_unique|)
//};

struct BucketContent {
	// propagated from Buckets, then updated individually
	int i;   // index of the next kmer for this bucket
	qpos_t seeds;

	// not propagated from Buckets
	qpos_t matches;
	int codirection;
	rpos_t r_min, r_max;

	BucketContent()
		: i(-1), seeds(0), matches(0), codirection(0), r_min(std::numeric_limits<rpos_t>::max()), r_max(-1) {}
	BucketContent(qpos_t matches, qpos_t seeds, int codirection, rpos_t r_min, rpos_t r_max)
		: i(-1), seeds(seeds), matches(matches), codirection(codirection), r_min(r_min), r_max(r_max) {}

	friend std::ostream& operator<<(std::ostream& os, const BucketContent& b) {
		os << "(i=" << b.i << ", seeds=" << b.seeds << ", matches=" << b.matches << ", codirection=" << b.codirection << ", r=[" << b.r_min << "," << b.r_max << "])";
		return os;
	}

	BucketContent operator+=(const BucketContent &other) {
		//ZoneScoped;
		matches += other.matches;
		seeds += other.seeds;
		codirection += other.codirection;
		r_min = std::min(r_min, other.r_min);
		r_max = std::max(r_max, other.r_max);
		return *this;
	}
};


using Seeds 	 = gtl::vector<Seed>;
using Matches 	 = gtl::vector<Match>;
using h2cnt 	 = ankerl::unordered_dense::map<hash_t, qpos_t>;
using h2seed_t   = ankerl::unordered_dense::map<hash_t, Seed>;

inline int codirection(const Kmer &kmer, const Hit &hit) {
	return kmer.strand == hit.strand ? +1 : -1;
}

inline int codirection(const Kmer &kmer1, const Kmer &kmer2) {
	return kmer1.strand == kmer2.strand ? +1 : -1;
}

enum Metric {
	Containment,
	Jaccard,
	bucket_SH,
	bucket_LCS,
};

inline std::string mapping_metric_str(Metric metric) {
	switch (metric) {
		case Metric::Containment: return "Containment";
		case Metric::Jaccard: return "Jaccard";
		case Metric::bucket_SH: return "bucket_SH";
		case Metric::bucket_LCS: return "bucket_LCS";
	}
	throw std::runtime_error("Invalid mapping metric: " + std::to_string(static_cast<int>(metric)));
}

inline Metric mapping_metric_from_str(const std::string& str) {
	if (str == "Containment") return Metric::Containment;
	if (str == "Jaccard") return Metric::Jaccard;
	if (str == "bucket_SH") return Metric::bucket_SH;
	if (str == "bucket_LCS") return Metric::bucket_LCS;
	throw std::runtime_error("Invalid mapping metric: " + str);
}

}  // namespace sweepmap