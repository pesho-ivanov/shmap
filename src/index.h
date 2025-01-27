#pragma once

#include <string>
#include <vector>

#include <ankerl/unordered_dense.h>
#include "../ext/gtl/phmap.hpp"
#include "../ext/gtl/vector.hpp"

#include "sketch.h"
#include "utils.h"
#include "io.h"
#include "handler.h"
//#include "rmq.h"

namespace sweepmap {

class SketchIndex {

public:
	std::vector<RefSegment> T;
	//ankerl::unordered_dense::map<hash_t, Hit> h2single;               // all sketched kmers with =1 hit
	//ankerl::unordered_dense::map<hash_t, gtl::vector<Hit>> h2multi;   // all sketched kmers with >1 hits
	gtl::flat_hash_map<hash_t, Hit> h2single;               // all sketched kmers with =1 hit
	gtl::flat_hash_map<hash_t, gtl::vector<Hit>> h2multi;   // all sketched kmers with >1 hits

	Handler *H;

	const RefSegment& get_segment(const string &segm_id) const {
		for (const auto &segm: T)
			if (segm.name == segm_id)
				return segm;
		throw std::runtime_error("Segment not found: " + segm_id);
	}

	void get_kmer_stats() const {
		std::vector<rpos_t> hist(10, 0);
		rpos_t max_occ = 0;
        H->C.inc("indexed_hits", h2single.size());
        H->C.inc("indexed_kmers", h2single.size());
		for (const auto& h2p : h2multi) {
			rpos_t occ = h2p.second.size();
			H->C.inc("indexed_hits", occ);
			H->C.inc("indexed_kmers");
			if (occ >= rpos_t(hist.size())-1) {
				hist.back() += occ;
				if (occ > max_occ)
					max_occ = occ;
			} else 
				hist[occ] += occ;
		}
		H->C.inc("indexed_highest_freq_kmer", max_occ);
	}

//	rpos_t count(hash_t h) const {
//		if (h2single.contains(h)) return 1;
//		else if (h2multi.contains(h)) return h2multi.at(h).size();
//		else return 0;
//	}

	inline rpos_t count(hash_t h) const {
		if (auto it = h2single.find(h); it != h2single.end())
			return 1;
		if (auto it = h2multi.find(h); it != h2multi.end())
			return it->second.size();
		return 0;
	}


	bool intersect(const Hit &hit, const rpos_t from, const rpos_t to) const {
		// TODO: account for segments
		return from <= hit.r && hit.r < to;
	}

//	bool get_matches_in_t_interval(std::vector<Match> *matches, const Seed &s, const rpos_t from, const rpos_t to) const {
//		// assume matches of each seed are sorted by position in T
//		// TODO: account for segments
//		if (s.hits_in_T == 1) {
//			auto &hit = h2single.at(s.kmer.h);
//			if (from <= hit.tpos && hit.tpos < to) {
//				matches->push_back(Match(s, hit));
//				return true;
//			}
//		} else if (s.hits_in_T > 1) {
//			// TODO: account for segments
//			// TODO: careful with left and right ends
//			const vector<Hit> &hits = h2multi.at(s.kmer.h);
//				auto it = lower_bound(hits.begin(), hits.end(), from, [](const Hit &hit, rpos_t pos) { return hit.tpos < pos; });
//				//auto it_r = lower_bound(hits.rbegin(), hits.rend(), to, [](const Hit &hit, rpos_t pos) { return hit.r > pos; });
//
//				//matches->push_back(Match(s, *it));
//				//matches->push_back(Match(s, *it_r));
//
//				for (; it != hits.end() && it->tpos < to; ++it)
//					matches->push_back(Match(s, *it));
//
//				return true;
////			}
//		}	
//		return false;
//	}

	Buckets::BucketContent matches_in_bucket(const Seed &s, const Buckets::BucketLoc b) const {
		Buckets::BucketContent bucket;
		bucket.seeds = s.occs_in_p;
		if (s.hits_in_T == 0) {
			return bucket;
		} else if (s.hits_in_T == 1) {
			auto &hit = h2single.at(s.kmer.h);
#ifdef SHMAP_ABS_POS
			if (b.begin() <= hit.r && hit.r < b.end()) {
#else
			if (b.begin() <= hit.tpos && hit.tpos < b.end()) {
#endif
				bucket.matches = 1;
				bucket.codirection = hit.strand == s.kmer.strand ? 1 : -1;
				bucket.r_min = hit.r;
				bucket.r_max = hit.r;
			}
		} else if (s.hits_in_T > 1) {
			const auto &hits = h2multi.at(s.kmer.h);
			auto it = lower_bound(hits.begin(), hits.end(), b.begin(), [&b](const Hit &hit, rpos_t pos) {
				if (hit.segm_id != b.segm_id)
					return hit.segm_id < b.segm_id;
#ifdef SHMAP_ABS_POS
				return hit.r < pos;
#else
				return hit.tpos < pos;
#endif
			});
			rpos_t matches = 0;
#ifdef SHMAP_ABS_POS	
			for (; it != hits.end() && it->segm_id == b.segm_id && it->r < b.end(); ++it) {
#else
			for (; it != hits.end() && it->segm_id == b.segm_id && it->tpos < b.end(); ++it) {
#endif
				matches += 1;
				bucket.codirection += it->strand == s.kmer.strand ? 1 : -1;
				bucket.r_min = std::min(bucket.r_min, it->r);
				bucket.r_max = std::max(bucket.r_max, it->r);
			}
			bucket.matches = min(matches, s.occs_in_p);
			//auto it_prev = it; if (it_prev != hits.begin()) { --it_prev; assert(it_prev->r < from); }
			//if (it != hits.end()) assert(from <= it->r);
			//cout << (it != hits.end() && it->r < to) << ", [" << from << ", " << to << ")" << ", it: " << (it != hits.end() ? it->r : -1) << endl;
			
			//return it != hits.end() && it->tpos < b.end();
		}
		return bucket;
	}

	void get_matches(gtl::vector<Match> *matches, const Seed &s) const {
		matches->reserve(matches->size() + s.hits_in_T);
		if (s.hits_in_T == 1) {
			matches->push_back(Match(s, h2single.at(s.kmer.h)));
		} else if (s.hits_in_T > 1) {
			for (const auto &hit: h2multi.at(s.kmer.h))
				matches->push_back(Match(s, hit));
		}	
	}

	void erase_frequent_kmers() {
		assert(H->params.max_matches != -1);
		std::vector<hash_t> blacklisted_h;
		for (const auto &[h, hits]: h2multi)
			if (rpos_t(hits.size()) > H->params.max_matches) {
				blacklisted_h.push_back(h);
				H->C.inc("blacklisted_kmers");
				H->C.inc("blacklisted_hits", hits.size());
			}

		for (auto h: blacklisted_h)
			h2multi.erase(h);
	}

	void populate_h2pos(const sketch_t& sketch, segm_t segm_id) {
		// TODO: skip creating the sketch structure
		// TODO: prepare the hits and initialize hashmaps in batch
		for (size_t tpos = 0; tpos < sketch.size(); ++tpos) {
			const Kmer& kmer = sketch[tpos];
			const auto hit = Hit(kmer, tpos, segm_id);
			if (!h2single.contains(kmer.h))
				h2single[kmer.h] = hit; 
            else if (H->params.max_matches == -1 || (rpos_t)h2multi[kmer.h].size() < H->params.max_matches + 1)
                h2multi[kmer.h].push_back(hit);
		}
	}

	void add_segment(const string &segm_name, const string &segm_seq, const sketch_t& sketch) {
		T.push_back(RefSegment(sketch, segm_name, segm_seq, segm_seq.size(), T.size()));
		H->C.inc("segments");
		H->C.inc("total_nucls", segm_seq.size());
		populate_h2pos(sketch, T.size()-1);
	}

	SketchIndex(Handler *H) : H(H) {
		H->C.init("segments");
		H->C.init("total_nucls");
	}

	void build_index(const std::string &tFile) {
		H->T.start("indexing");
		cerr << "Indexing " << tFile << "..." << endl;
		H->T.start("index_reading");
		read_fasta_klib(tFile, [this](const string &segm_name, const string &T, float progress) {
			H->T.stop("index_reading");
			H->T.start("index_sketching");
			sketch_t t = H->sketcher.sketch(T);
			H->T.stop("index_sketching");

			H->T.start("index_initializing");
			add_segment(segm_name, T, t);
			H->T.stop("index_initializing");

			string msg = std::to_string(H->C.count("segments")) + " indexed segment" + (H->C.count("segments") == 1 ? "" : "s");
			printProgress(std::cerr, progress, msg);

			H->T.start("index_reading");
		});
		cerr << std::endl;

		// if a kmer is present in both single and multi, we move it out of single to multi
		for (auto &[h, hits] : h2multi) {
			if (h2single.contains(h)) {
				hits.push_back(h2single.at(h));
				h2single.erase(h);
			}
			// sorting each list of matches allows for binary search
			sort(hits.begin(), hits.end(), [](const Hit &a, const Hit &b) {
				if (a.segm_id != b.segm_id)
					return a.segm_id < b.segm_id;
				return a.r < b.r;  // equivalent to a.tpos < b.tpos
			});
		}
		H->T.stop("index_reading");
		H->T.stop("indexing");

		get_kmer_stats();
        H->C.inc("blacklisted_kmers", 0);
        H->C.inc("blacklisted_hits", 0);
		if (H->params.max_matches != -1)
			erase_frequent_kmers();
		print_stats();
	}

	void print_stats() const {
		cerr << std::fixed << std::setprecision(1);
        //printMemoryUsage();
		cerr << " | total nucleotides:     " << H->C.count("total_nucls") << endl;
		cerr << " | index segments:        " << H->C.count("segments") << " (~" << 1.0*H->C.count("total_nucls") / H->C.count("segments") << " nb per segment)" << endl;
		for (const auto &segm: T)
			cerr << " | | " << segm.name << " (" << segm.sz << " nb)" << endl;
		cerr << " | indexed kmers:         " << H->C.count("indexed_kmers") << endl;
		cerr << " | indexed hits:          " << H->C.count("indexed_hits") << " ("
												<< double(H->params.k)*H->C.perc("indexed_hits", "total_nucls") << "\% of the index, "
												<< "~" << H->C.frac("indexed_hits", "indexed_kmers") << " per kmer)" << endl;
		cerr << " | | most frequent kmer:      " << H->C.count("indexed_highest_freq_kmer") << " times." << endl;
		cerr << " | | blacklisted kmers:       " << H->C.count("blacklisted_kmers") << " (" << H->C.perc("blacklisted_kmers", "indexed_kmers") << "\%)" << endl;
		cerr << " | | blacklisted hits:        " << H->C.count("blacklisted_hits") << " (" << H->C.perc("blacklisted_hits", "indexed_hits") << "\%)" << endl;
		cerr << " | indexing time:        "    << H->T.secs("indexing") << "s" << endl;
		cerr << " | | reading time:          " << H->T.secs("index_reading") << "s" << endl;
		cerr << " | | sketching time:        " << H->T.secs("index_sketching") << "s" << endl;
		cerr << " | | initializing time:     " << H->T.secs("index_initializing") << "s" << endl;
	}
};

} // namespace sweepmap

	//bool is_kmer_in_interval(const Seed &s, int from, int to) const {
	//	assert(s.hits_in_T > 0);
	//	if (s.hits_in_T == 1) {
	//		auto &hit = h2single.at(s.kmer.h);
	//		return intersect(hit, from, to);
	//	} else if (s.hits_in_T > 1) {
	//		const vector<Hit> &hits = h2multi.at(s.kmer.h);
	//		//for (int i=1; i<(int)hits.size(); i++) assert(hits[i-1].r <= hits[i].r);
	//		auto it = lower_bound(hits.begin(), hits.end(), from, [](const Hit &hit, int pos) { return hit.r < pos; });
	//		//auto it_prev = it; if (it_prev != hits.begin()) { --it_prev; assert(it_prev->r < from); }
	//		//if (it != hits.end()) assert(from <= it->r);
	//		//cout << (it != hits.end() && it->r < to) << ", [" << from << ", " << to << ")" << ", it: " << (it != hits.end() ? it->r : -1) << endl;
	//		return it != hits.end() && it->r < to;
	//	}
	//	return false;
	//}

	// returns the number of hits in the interval from, to
//	void match_kmer_in_interval(SegmentTree *hist, const Seed &s, const int from, const int to, vector<Match> *matches_freq) const {
//		// assume matches of each seed are sorted by position in T
//		if (s.hits_in_T == 1) {
//			auto &hit = h2single.at(s.kmer.h);
//			int curr_l = hist->from(hit);
//			int curr_r = hist->to(hit);
//			if (from <= curr_r && curr_l <= to) {
//				matches_freq->push_back(Match(s, hit));
//				hist->incRange(max(curr_l, from), min(curr_r, to));
//			}
//		} else if (s.hits_in_T > 1) {
//			const vector<Hit> &hits = h2multi.at(s.kmer.h);
//			auto it = lower_bound(hits.begin(), hits.end(), from, [&hist](const Hit &hit, int pos) {
//				return hist->to(hit) < pos;
//			});
//			assert(it==hits.begin() || hist->to(*(it-1)) < from);
//			assert(it==hits.end() || from <= hist->to(*it));
//
//			int prev_l = -1;  //hist->from(*it);
//			int prev_r = -1;  // hist->to(*it);
//			for (; it != hits.end() && prev_l <= to; ++it) {
//				int curr_l = hist->from(*it);
//				int curr_r = hist->to(*it);
//				assert(curr_l <= curr_r);
//				if (from <= curr_r && curr_l <= to)  // intersect?
//					matches_freq->push_back(Match(s, *it));
//				if (prev_r < curr_l) {  // if there will be a gap, push the current range
//					if (prev_l != -1) {
//#ifdef DEBUG
//						cerr << "  match_kmer_in_interval(" << s << ", " << from << ", " << to << ")" << endl;
//#endif
//						hist->incRange(max(prev_l, from), min(prev_r, to));
//					}
//					prev_l = curr_l;
//					prev_r = curr_r;
//				} else {
//					assert(prev_l <= curr_l);
//					prev_r = curr_r;
//				}
//			}
//			if (prev_l != -1 && prev_l <= to) {
//				assert (prev_l <= prev_r);
//#ifdef DEBUG
//				cerr << "  match_kmer_in_interval(" << s << ", " << from << ", " << to << ")" << endl;
//#endif
//				hist->incRange(max(prev_l, from), min(prev_r, to));
//			}	
//		}
//	}

			//if (hits.size() < 10) {
			//	for (int i=0; i<(int)hits.size(); i++)
			//		if (hits[i].r >= from) {
			//			matches->push_back(Match(s, hits[i]));
			//			for (int j=(int)hits.size()-1; j>i; j--)
			//				if (hits[j].r < to) {
			//					matches->push_back(Match(s, hits[j]));
			//					break;
			//				}
			//			return true;
			//		}
			//} else {
//			if (hits.size() < 10) {
//				for (int i=0; i<(int)hits.size(); i++)
//					if (hits[i].tpos >= from) {
//						matches->push_back(Match(s, hits[i]));
//						for (int j=(int)hits.size()-1; j>i; j--)
//							if (hits[j].tpos < to) {
//								matches->push_back(Match(s, hits[j]));
//								break;
//							}
//						return true;
//					}
//			} else {

//	bool get_matches_in_interval(std::vector<Match> *matches, const Seed &s, const bucket_t b, qpos_t lmax) const {
//		// assume matches of each seed are sorted by position in T
//		// TODO: account for segments
//		rpos_t from = b.b * lmax;
//		rpos_t to = (b.b + 2) * lmax;
//		if (s.hits_in_T == 1) {
//			auto &hit = h2single.at(s.kmer.h);
//			if (intersect(hit, from, to)) {
//				matches->push_back(Match(s, hit));
//				return true;
//			}
//		} else if (s.hits_in_T > 1) {
//			// TODO: account for segments
//			// TODO: careful with left and right ends
//			const vector<Hit> &hits = h2multi.at(s.kmer.h);
//				auto it = lower_bound(hits.begin(), hits.end(), from, [](const Hit &hit, rpos_t pos) {
////					if (b.segm_id != hit.segm_id)
////						return b.segm_id < hit.segm_id;
//					return hit.r < pos;  // hit.tpos!!
//				});
//				//auto it_r = lower_bound(hits.rbegin(), hits.rend(), to, [](const Hit &hit, rpos_t pos) { return hit.r > pos; });
//
//				//matches->push_back(Match(s, *it));
//				//matches->push_back(Match(s, *it_r));
//
//				for (; it != hits.end() && it->r < to; ++it)
//					matches->push_back(Match(s, *it));
//
//				return true;
//			//}
//		}	
//		return false;
//	}