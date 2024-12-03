#pragma once

#include <ankerl/unordered_dense.h>
#include "utils.h"
#include "index.h"
//#include "edlib.h"

namespace sweepmap {

using Hist = ankerl::unordered_dense::map<hash_t, qpos_t>;
//typedef pair<bucket_t, qpos_t> Bucket;
//using Buckets = std::unordered_map<bucket_t, qpos_t>;
//using Buckets = gtl::flat_hash_map<bucket_t, qpos_t>;
//using Buckets = ankerl::unordered_dense::map<bucket_t, qpos_t>;
//using Hist = unordered_map<hash_t, qpos_t>;
//using Hist = gtl::flat_hash_map<hash_t, qpos_t>;

enum class Metric {
	JACCARD,
	CONTAINMENT_INDEX
};

class Matcher {
	const SketchIndex &tidx;
	Hist diff_hist;

public:
	Matcher(const SketchIndex &tidx) : tidx(tidx) {}
	Matcher(const SketchIndex &tidx, const Hist &diff_hist) : tidx(tidx), diff_hist(diff_hist) {}

	void update(const Hist &diff_hist) {
		this->diff_hist = diff_hist;
	}

    bool do_overlap(const string& query_id, const Buckets::Bucket& b) {
        auto parsed = ParsedQueryId::parse(query_id);
        if (!parsed.valid) return false;
		//cerr << query_id << ", " << b.segm_id << ", " << b.begin() << ", " << b.end() << endl;
		//cerr << "name: " << tidx.T[b.first.segm_id].name<< ", parsed.segm_id: " << parsed.segm_id << endl;
		assert(b.segm_id >= 0 && b.segm_id < (rpos_t)tidx.T.size());
		const auto &segm = tidx.T[b.segm_id];
		if (segm.name == parsed.segm_id) {
			//rpos_t bucket_start_t = b.b * bucket_l;
			//rpos_t bucket_end_t = bucket_start_t + 2 * bucket_l;
			assert(b.begin() >= 0 && b.begin() < (rpos_t)segm.kmers.size());
			rpos_t bucket_start_T = segm.kmers[b.begin()].r;
			//assert(bucket_end_t >= 0 && bucket_end_t < (rpos_t)segm.kmers.size());
			rpos_t bucket_end_T = b.end() < (rpos_t)segm.kmers.size() ? segm.kmers[b.end()].r : segm.sz;  // b.end() may not be a valid index

//			cerr << "bucket_start_t: " << bucket_start_t << ", bucket_end_t: " << bucket_end_t << ", bucket_start_T: " << bucket_start_T << ", bucket_end_T: " << bucket_end_T << endl;
			if (bucket_start_T < parsed.end_pos && bucket_end_T >= parsed.start_pos) {  // if overlap
				return true;
			}
		}
		return false;
	}

    bool lost_correct_mapping(const string& query_id, const Buckets& B) {
		bool res = true;
		for (auto it = B.unordered_begin(); it != B.unordered_end(); ++it)
				if (do_overlap(query_id, it->first)) {
					res = false;
					break;
				}
        return res;
    }

	static double ContainmentIndex(qpos_t intersection, qpos_t m) {
		//assert(intersection <= s_sz);
		assert(intersection <= m);
		assert(1.0*intersection / m <= 1.0);
		return 1.0*intersection / m;
	}

	static double Jaccard(qpos_t intersection, qpos_t m, qpos_t s_sz) {
		assert(intersection <= s_sz);
		assert(intersection <= m);
		assert(1.0*intersection / (m + s_sz - intersection) <= 1.0);
		return 1.0*intersection / (m + s_sz - intersection);
	}

	Matches collect_matches(const Buckets::Bucket &b, const unordered_map<hash_t, Seed> &p_ht) {
		Matches M;
		//cerr << b << endl;
		assert(b.parent != nullptr);
		assert(b.segm_id >= 0 && b.segm_id < (rpos_t)tidx.T.size());
		for (rpos_t i = b.begin(); i < std::min(b.end(), (rpos_t)tidx.T[b.segm_id].kmers.size()); i++) {
			assert(i < (rpos_t)tidx.T[b.segm_id].kmers.size());
			//cerr << "parent" << b.parent << " " << b << " " << i << " " << tidx.T.size() << " " << tidx.T[b.segm_id].kmers.size() << endl;
			const auto kmer = tidx.T[b.segm_id].kmers[i];
			const auto seed_it = p_ht.find(kmer.h);
			if (seed_it != p_ht.end()) {
				M.push_back(Match(seed_it->second, Hit(kmer, i, b.segm_id)));
				assert(M[ M.size()-1 ].hit.tpos == M[ M.size()-1 ].hit.tpos);
			}
		}
		return M;
	}

	Mapping bestIncludedJaccard(const vector<Match> &M, const qpos_t P_sz, qpos_t lmin, qpos_t lmax, const qpos_t m) {
		qpos_t intersection = 0;
		qpos_t same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise
		Mapping best;
		
		for (int i=1; i<(int)M.size(); i++)
			assert(M[i-1].hit.tpos < M[i].hit.tpos);
		assert(M.size() == 0 || (int)M.size() <= M.back().hit.tpos - M.front().hit.tpos + 1);

		auto l = M.begin(), r = M.begin();
		for(; l != M.end(); ++l) {
			// move r to the left until l+lmin<r
			for(; l + 1 < r && l->hit.tpos + lmin < r->hit.tpos; ) {
				--r;
				if (++diff_hist[r->seed.el.h] >= 1) {
					--intersection;
					same_strand_seeds -= r->is_same_strand() ? +1 : -1;
				}
			}

			// move r to the right until l+lmax<=r
			for(;  r != M.end() && r->hit.tpos <= l->hit.tpos + lmax; ++r) {
				assert(diff_hist.contains(r->seed.el.h));
				if (--diff_hist[r->seed.el.h] >= 0) {
					++intersection;
					same_strand_seeds += r->is_same_strand() ? +1 : -1;
				}
				
				assert (l->hit.r <= r->hit.r);
				if (l < r) {
					int s_sz = prev(r)->hit.tpos - l->hit.tpos + 1;
					double J = Jaccard(intersection, m, s_sz);
					best.update(0, P_sz-1, l->hit.r, prev(r)->hit.r, tidx.T[l->hit.segm_id], intersection, J, same_strand_seeds, l, prev(r)); // TODO: should it be prev(r) instead?
					//best = Mapping(H->params.k, P_sz, m, l->hit.r, r->hit.r, l->hit.segm_id, intersection, J, same_strand_seeds, l, r);
				}
			}

			assert(diff_hist.contains(l->seed.el.h));
			if (++diff_hist[l->seed.el.h] >= 1) {
				--intersection;
				same_strand_seeds -= l->is_same_strand() ? +1 : -1;
			}

			assert(intersection >= 0);
		}
		assert(l == M.end());
		assert(r == M.end());
		assert(intersection == 0);
		
		return best;
	}

	Mapping bestFixedLength(const vector<Match> &M, const qpos_t P_sz, qpos_t lmax, const qpos_t m, Metric metric) {
		qpos_t intersection = 0;
		qpos_t same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise
		Mapping best;
		
		for (int i=1; i<(int)M.size(); i++)
			assert(M[i-1].hit.tpos < M[i].hit.tpos);
		assert(M.size() == 0 || (int)M.size() <= M.back().hit.tpos - M.front().hit.tpos + 1);

		auto l = M.begin(), r = M.begin();
		for(; l != M.end(); ++l) {
			for(;  r != M.end()  && r->hit.tpos <= l->hit.tpos + lmax; ++r) {
				//&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				//&& r->hit.tpos + H->params.k <= l->hit.tpos + P_sz
				//&& r->hit.r + H->params.k <= l->hit.r + P_sz
				assert(diff_hist.contains(r->seed.el.h));
				if (--diff_hist[r->seed.el.h] >= 0) {
					++intersection;
					same_strand_seeds += r->is_same_strand() ? +1 : -1;
				}
				
				assert (l->hit.r <= r->hit.r);
			}

			double J;
			if (metric == Metric::JACCARD) {
				int s_sz = prev(r)->hit.tpos - l->hit.tpos + 1;
				J = Jaccard(intersection, m, s_sz);
			} else if (metric == Metric::CONTAINMENT_INDEX) {
				J = ContainmentIndex(intersection, m);
			} else {
				J = -1.0;
				assert(false);
			}

			assert(J >= -0.0);
			assert(J <= 1.0);
			if (J > best.score())
				best.update(0, P_sz-1, l->hit.r, prev(r)->hit.r, tidx.T[l->hit.segm_id], intersection, J, same_strand_seeds, l, prev(r)); // TODO: should it be prev(r) instead?
				//best = Mapping(H->params.k, P_sz, m, l->hit.r, prev(r)->hit.r, l->hit.segm_id, intersection, J, same_strand_seeds, l, prev(r));

			assert(diff_hist.contains(l->seed.el.h));
			if (++diff_hist[l->seed.el.h] >= 1) {
				--intersection;
				same_strand_seeds -= l->is_same_strand() ? +1 : -1;
			}

			assert(intersection >= 0);
		}
		assert(l == M.end());
		assert(r == M.end());
		assert(intersection == 0);
		
		return best;
	}
};

}  // namespace sweepmap