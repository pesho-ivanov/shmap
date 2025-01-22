#pragma once

#include <ankerl/unordered_dense.h>
#include <functional>  // for std::hash

#include "utils.h"

namespace sweepmap {

class Buckets {
public:
	Buckets(qpos_t len) : len(len) {}

    struct BucketLoc {
        segm_t segm_id;		// segm_id refers to tidx.T[segm_id]
        rpos_t b;   		// b refers to interval [lmax*b, lmax*(b+2)) in tidx.T[segm_id].kmers
		const Buckets* parent;

        BucketLoc() : segm_id(-1), b(-1), parent(nullptr) {}
        BucketLoc(segm_t segm_id, rpos_t b, const Buckets* parent) : segm_id(segm_id), b(b), parent(parent) {}

        bool operator==(const BucketLoc &other) const {
            return segm_id == other.segm_id && b == other.b;
        }
        friend std::ostream& operator<<(std::ostream& os, const BucketLoc& b) {
            os << "(" << b.segm_id << "," << b.b << ")";
            return os;
        }

		rpos_t begin() const {
			return b * parent->len;
		}

		rpos_t end() const {
			return (b + 2) * parent->len;
		}
	};
	struct BucketHash {
		std::size_t operator()(const BucketLoc& b) const {
			return std::hash<segm_t>()(b.segm_id) ^ (std::hash<rpos_t>()(b.b) << 1);
		}
	};

	struct BucketContent {
		qpos_t matches;
		int codirection;
		rpos_t r_min, r_max;

		BucketContent() : matches(0), codirection(0), r_min(std::numeric_limits<rpos_t>::max()), r_max(-1) {}
		BucketContent(qpos_t matches, int codirection, rpos_t r_min, rpos_t r_max) : matches(matches), codirection(codirection), r_min(r_min), r_max(r_max) {}

		friend std::ostream& operator<<(std::ostream& os, const BucketContent& b) {
			os << "(" << b.matches << "," << b.codirection << "," << b.r_min << "," << b.r_max << ")";
			return os;
		}

		BucketContent operator+=(const BucketContent &other) {
			matches += other.matches;
			codirection += other.codirection;
			r_min = std::min(r_min, other.r_min);
			r_max = std::max(r_max, other.r_max);

			return *this;
		}
	};

	using bucket_map_t = ankerl::unordered_dense::map<BucketLoc, BucketContent, BucketHash, std::equal_to<BucketLoc> >;

	struct Pos {
		segm_t segm_id;
		qpos_t r;

		Pos() : segm_id(-1), r(-1) {}
		Pos(segm_t segm_id, qpos_t r)
			: segm_id(segm_id), r(r) {}
	};

	//Pos get_pos(segm_t segm_id, qpos_t r) const {
	//	return Pos(segm_id, r, this);
	//}

    void add_to_pos(const Pos& p, const BucketContent &content) {
        qpos_t b = p.r/len;
        _B[ BucketLoc(p.segm_id, b, this) ] += content;
        if (b>0) _B[ BucketLoc(p.segm_id, b-1, this) ] += content;
    }

    void assign_to_pos(const Pos& p, const BucketContent &content) {
        qpos_t b = p.r/len;
        _B[ BucketLoc(p.segm_id, b, this) ] = content;
        if (b>0) _B[ BucketLoc(p.segm_id, b-1, this) ] = content;
    }

    void add_to_bucket(BucketLoc b, const BucketContent &content) {
		b.parent = this;
        _B[b] += content;
    }

	bool delete_bucket(BucketLoc b) {
		return _B.erase(b);
	}

	qpos_t get_bucket_len() const {
		return len;
	}

	qpos_t get_matches(const BucketLoc& b) const {
		auto it = _B.find(b);
		if (it == _B.end()) return 0;
		return it->second.matches;
	}

	int size() const {
		return _B.size();
	}

	class UnorderedIterator {
	public:
		using unordered_iterator = bucket_map_t::const_iterator;
		using iterator_category = std::forward_iterator_tag;
		using value_type = Pos;
		using difference_type = std::ptrdiff_t;
		using pointer = Pos*;
		using reference = Pos&;

		UnorderedIterator(unordered_iterator it) : it_(it) {}

	    const bucket_map_t::value_type& operator*() const {
			return *it_;
		}

        const bucket_map_t::value_type* operator->() const {
            return &(*it_);
        }

		UnorderedIterator& operator++() {
			++it_;
			return *this;
		}

		bool operator!=(const UnorderedIterator& other) const {
			return it_ != other.it_;
		}

	private:
		unordered_iterator it_;
		const Buckets* parent_;
	};

	UnorderedIterator unordered_begin() const {
		return UnorderedIterator(_B.begin());
	}

	UnorderedIterator unordered_end() const {
		return UnorderedIterator(_B.end());
	}

	class OrderedIterator {
	public:
		using sorted_vector_t = std::vector<std::reference_wrapper<bucket_map_t::value_type>>;
		using sorted_iterator = typename sorted_vector_t::iterator;
		using iterator_category = std::forward_iterator_tag;
		using value_type = Pos;
		using difference_type = std::ptrdiff_t;
		using pointer = Pos*;
		using reference = Pos&;

		OrderedIterator(sorted_iterator it) : it_(it) {}

		bucket_map_t::value_type& operator*() {
			return it_->get();
		}

        bucket_map_t::value_type* operator->() {
            return &(it_->get());
        }

		OrderedIterator& operator++() {
			++it_;
			return *this;
		}

		bool operator!=(const OrderedIterator& other) const {
			return it_ != other.it_;
		}

	private:
		sorted_iterator it_;
	};

    // TODO: sort once, not every time; make sure nothing changes
	// Helper to create a sorted vector
	std::vector<std::reference_wrapper<bucket_map_t::value_type>> get_sorted_buckets() {
		std::vector<std::reference_wrapper<bucket_map_t::value_type>> sorted_buckets;
		sorted_buckets.reserve(_B.size());
		for (auto& bucket : _B) {
			sorted_buckets.push_back(std::ref(bucket));
		}
		std::sort(sorted_buckets.begin(), sorted_buckets.end(),
				  [](const auto& a, const auto& b) {
					  // Define sorting criteria, e.g., by segm_id and then r
					//  if (a.get().first.segm_id != b.get().first.segm_id)
					//	  return a.get().first.segm_id < b.get().first.segm_id;
                    // Sort by decreasing number of matches
					  return a.get().second.matches > b.get().second.matches; 
				  });
		return sorted_buckets;
	}

	OrderedIterator ordered_begin() {
		sorted_buckets_ = get_sorted_buckets();
		return OrderedIterator(sorted_buckets_.begin());
	}

	OrderedIterator ordered_end() {
		return OrderedIterator(sorted_buckets_.end());
	}

private:
	std::vector<std::reference_wrapper<bucket_map_t::value_type>> sorted_buckets_;

	qpos_t len;  // half-length; the bucket spans [r*len, (r+2)*len)
	bucket_map_t _B;
};

} // namespace sweepmap
