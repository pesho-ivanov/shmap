#pragma once

#include <functional>  // for std::hash
#include <vector>
#include "types.h"
#include "tracy/public/tracy/Tracy.hpp"
#include "unordered_dense/include/ankerl/unordered_dense.h"

#include "index.h"

namespace sweepmap {
    class SketchIndex;
}

namespace sweepmap {

template<bool abs_pos>
class BucketsHash {
	using bucket_map_t = ankerl::unordered_dense::map<BucketLoc, BucketContent, BucketHash, std::equal_to<BucketLoc> >;

public:
	qpos_t halflen;  // half-length; the bucket spans [r*len, (r+2)*len)

	int i;  				// index of the next kmer for all buckets
	int seeds; 				// number of seeds in all buckets

	bucket_map_t buckets;

	BucketsHash() : halflen(-1), i(0), seeds(0) {}
	BucketsHash(qpos_t halflen) : halflen(halflen), i(0), seeds(0) {}

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
};

template<bool abs_pos>
class Buckets {
	static int const MAX_SEGMENTS = 100;
public:
	const SketchIndex &tidx;
	qpos_t min_halflen;
	qpos_t halflen;  // half-length; the bucket spans [r*len, (r+2)*len)

	int i;  				// index of the next kmer for all buckets
	int seeds; 				// number of seeds in all buckets
							//
    std::vector<BucketContent> buckets[MAX_SEGMENTS];   // buckets[segm_id][b]
	std::vector<BucketLoc> non_empty_buckets_with_repeats;

	Buckets(const SketchIndex &tidx, qpos_t min_halflen)
	: tidx(tidx), min_halflen(min_halflen), halflen(-1), i(0), seeds(0) {
		if (tidx.segments() > MAX_SEGMENTS)
			throw std::runtime_error("Number of segments exceeds MAX_SEGMENTS (" + std::to_string(MAX_SEGMENTS) + ")");
		for (int b = 0; b < (int)tidx.segments(); ++b)
			buckets[b].resize(tidx.get_segment(b).sz / min_halflen + 2);
	}
	
	void clear() {
		i = 0;
		seeds = 0;
		for (auto &[segm_id, b] : non_empty_buckets_with_repeats)
			buckets[segm_id][b] = BucketContent();
		non_empty_buckets_with_repeats.clear();
	}
	
	bool set_halflen(qpos_t new_halflen) {
		//cerr << "min_halflen: " << min_halflen << ", new_halflen: " << new_halflen << endl;
		halflen = new_halflen;
		return halflen >= min_halflen;
	}

	rpos_t begin(const BucketLoc& b) const {
		return b.b * halflen;
	}

	rpos_t end(const BucketLoc& b) const {
		return (b.b + 2) * halflen;
	}

	void propagate_seeds_to_buckets() {
		for (auto &[segm_id, b] : non_empty_buckets_with_repeats) {
			//assert(buckets[segm_id][b].i == -1);
			// works even with repeated buckets
			buckets[segm_id][b].i = i;
			buckets[segm_id][b].seeds = seeds;
		}
	}

    void add_to_pos(const Hit& hit, const BucketContent &content) {
		qpos_t b = (abs_pos ? hit.r : hit.tpos) / halflen;
		//cerr << "halflen: " << halflen << ", b: " << b << ", hit.segm_id: " << hit.segm_id << ", buckets[segm_id].size = " << buckets[hit.segm_id].size() << endl;
		assert(hit.segm_id < (int)tidx.segments());
		assert(b < (int)buckets[hit.segm_id].size());
        buckets[hit.segm_id][b] += content;
		non_empty_buckets_with_repeats.push_back(BucketLoc(hit.segm_id, b));
        if (b > 0) {
			buckets[hit.segm_id][b-1] += content;
			non_empty_buckets_with_repeats.push_back(BucketLoc(hit.segm_id, b-1));
		}
    }

    void add_to_bucket(BucketLoc b, const BucketContent &content) {
        buckets[b.segm_id][b.b] += content;
		non_empty_buckets_with_repeats.push_back(b);
    }

	int size() const {
		return -1;  // TODO: implement
	}

	std::vector< std::pair<BucketLoc, BucketContent> > get_sorted_buckets() {
		ZoneScoped;
		
		// Sort non_empty_buckets_with_repeats to prepare for unique operation
		std::sort(non_empty_buckets_with_repeats.begin(), non_empty_buckets_with_repeats.end(), 
			[](const BucketLoc &a, const BucketLoc &b) {
				if (a.segm_id != b.segm_id)
					return a.segm_id < b.segm_id;
				return a.b < b.b;
			});
		
		auto unique_end = std::unique(non_empty_buckets_with_repeats.begin(), non_empty_buckets_with_repeats.end());
		std::vector< std::pair<BucketLoc, BucketContent> > sorted_buckets;
		
		for (auto it = non_empty_buckets_with_repeats.begin(); it != unique_end; ++it) {
			const BucketLoc& loc = *it;
			sorted_buckets.push_back(std::make_pair(loc, buckets[loc.segm_id][loc.b]));
		}
		
		std::sort(sorted_buckets.begin(), sorted_buckets.end(), [](const auto &a, const auto &b) {
			return a.second.matches > b.second.matches;
		});
		return sorted_buckets;
	}
};


} // namespace sweepmap
