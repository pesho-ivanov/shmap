#pragma once

#include <ankerl/unordered_dense.h>
#include <functional>  // for std::hash

#include "utils.h"

namespace sweepmap {

class Buckets {
public:
	Buckets(qpos_t len) : len(len) {}

    struct Bucket {
        segm_t segm_id;		// segm_id refers to tidx.T[segm_id]
        rpos_t b;   		// b refers to interval [lmax*b, lmax*(b+2)) in tidx.T[segm_id].kmers
		const Buckets* parent;

        Bucket() : segm_id(-1), b(-1), parent(nullptr) {}
        Bucket(segm_t segm_id, rpos_t b, const Buckets* parent) : segm_id(segm_id), b(b), parent(parent) {}

        bool operator==(const Bucket &other) const {
            return segm_id == other.segm_id && b == other.b;
        }
        friend std::ostream& operator<<(std::ostream& os, const Bucket& b) {
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
		std::size_t operator()(const Bucket& b) const {
			return std::hash<segm_t>()(b.segm_id) ^ (std::hash<rpos_t>()(b.b) << 1);
		}
	};
	using bucket_map_t = ankerl::unordered_dense::map<Bucket, qpos_t, BucketHash, std::equal_to<Bucket> >;

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

    void add_to_pos(const Pos& p, qpos_t matches) {
        qpos_t b = p.r/len;
        _B[ Bucket(p.segm_id, b, this) ] += matches;
        if (b>0) _B[ Bucket(p.segm_id, b-1, this) ] += matches;
    }

    void assign_to_pos(const Pos& p, qpos_t matches) {
        qpos_t b = p.r/len;
        _B[ Bucket(p.segm_id, b, this) ] = matches;
        if (b>0) _B[ Bucket(p.segm_id, b-1, this) ] = matches;
    }

    void add_to_bucket(Bucket b, qpos_t matches) {
		b.parent = this;
        _B[b] += matches;
    }

	bool delete_bucket(Bucket b) {
		return _B.erase(b);
	}

	qpos_t get_bucket_len() const {
		return len;
	}

	qpos_t get_matches(const Bucket& b) const {
		auto it = _B.find(b);
		if (it == _B.end()) return 0;
		return it->second;
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
		using sorted_vector_t = std::vector<bucket_map_t::value_type>;
		using sorted_iterator = typename sorted_vector_t::const_iterator;
		using iterator_category = std::forward_iterator_tag;
		using value_type = Pos;
		using difference_type = std::ptrdiff_t;
		using pointer = Pos*;
		using reference = Pos&;

		OrderedIterator(sorted_iterator it) : it_(it) {}

		const bucket_map_t::value_type& operator*() const {
			return *it_;
		}

        const bucket_map_t::value_type* operator->() const {
            return &(*it_);
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
	std::vector<bucket_map_t::value_type> get_sorted_buckets() const {
		std::vector<bucket_map_t::value_type> sorted_buckets(_B.begin(), _B.end());
		std::sort(sorted_buckets.begin(), sorted_buckets.end(),
				  [](const bucket_map_t::value_type& a, const bucket_map_t::value_type& b) {
					  // Define sorting criteria, e.g., by segm_id and then r
					//  if (a.first.segm_id != b.first.segm_id)
					//	  return a.first.segm_id < b.first.segm_id;
                    // Sort by decreasing number of matches
					  return a.second > b.second; 
				  });
		return sorted_buckets;
	}

	OrderedIterator ordered_begin() const {
		sorted_buckets_ = get_sorted_buckets();
		return OrderedIterator(sorted_buckets_.cbegin());
	}

	OrderedIterator ordered_end() const {
		return OrderedIterator(sorted_buckets_.cend());
	}

private:
	mutable std::vector<bucket_map_t::value_type> sorted_buckets_;

	qpos_t len;  // half-length; the bucket spans [r*len, (r+2)*len)
	bucket_map_t _B;
};

} // namespace sweepmap
