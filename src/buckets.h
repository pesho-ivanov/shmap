#pragma once

#include <functional>  // for std::hash
#include <vector>
#include "types.h"
#include "tracy/public/tracy/Tracy.hpp"
#include "unordered_dense/include/ankerl/unordered_dense.h"

namespace sweepmap {

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

template<bool abs_pos>
class Buckets {
	using bucket_map_t = ankerl::unordered_dense::map<BucketLoc, BucketContent, BucketHash, std::equal_to<BucketLoc> >;

public:
	qpos_t halflen;  // half-length; the bucket spans [r*len, (r+2)*len)

	int i;  				// index of the next kmer for all buckets
	int seeds; 				// number of seeds in all buckets

	Buckets() : halflen(-1), i(0), seeds(0) {}
	Buckets(qpos_t halflen) : halflen(halflen), i(0), seeds(0) {}

	rpos_t begin(const BucketLoc& b) const {
		return b.b * halflen;
	}

	rpos_t end(const BucketLoc& b) const {
		return (b.b + 2) * halflen;
	}

	void propagate_seeds_to_buckets() {
		//ZoneScoped;
		//cerr << "propagate: i=" << i << " seeds=" << seeds << endl;
		for (auto &bucket : buckets) {
			assert(bucket.second.i == -1);
			bucket.second.i = i;
			bucket.second.seeds = seeds;
		}
	}

    void add_to_pos(const Hit& hit, const BucketContent &content) {
		//ZoneScoped;
		qpos_t b = (abs_pos ? hit.r : hit.tpos) / halflen;
        buckets[ BucketLoc(hit.segm_id, b) ] += content;
        if (b>0) buckets[ BucketLoc(hit.segm_id, b-1) ] += content;
    }

    void add_to_bucket(BucketLoc b, const BucketContent &content) {
		//ZoneScoped;
        buckets[b] += content;
    }

	int size() const {
		return buckets.size();
	}

	std::vector<typename bucket_map_t::value_type> get_sorted_buckets() {
		ZoneScoped;
		std::vector<typename bucket_map_t::value_type> sorted_buckets(buckets.begin(), buckets.end());
		std::sort(sorted_buckets.begin(), sorted_buckets.end(), [](const auto &a, const auto &b) {
			return a.second.matches > b.second.matches;
		});
		return sorted_buckets;
	}

	bucket_map_t buckets;
};

} // namespace sweepmap
