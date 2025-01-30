#pragma once

#include "index.h"
#include "io.h"
#include "sketch.h"
#include "utils.h"
#include "handler.h"
#include "mapper.h"
#include "refine.h"

#include "../ext/gtl/vector.hpp"

namespace sweepmap {

template<bool abs_pos>
class AnalyseSimulatedReads {
public:
	const string &query_id;
	const string &P;
	const int P_sz;
	const h2cnt &diff_hist;
	const int m;
	const h2seed_t &p_ht;
	const SketchIndex &tidx;
	const Buckets<abs_pos> &B;
	const double theta;

	Matcher<abs_pos> matcher;

	rpos_t gt_start_nucl, gt_end_nucl;
	// calculated
	int bucket_l;
	rpos_t start, end;
	int segm_id;
	string segm_name;
	vector<BucketLoc> J_buckets;
	vector<BucketLoc> C_buckets;
	BucketLoc gt_b_l;
	BucketLoc gt_b_r; 
	BucketLoc gt_b_next;
	Matches gt_M_l;
	Matches gt_M_r;
	Matches gt_M_next;

	Mapping gt_mapping;
	Mapping gt_J_l;
	Mapping gt_J_r;
	Mapping gt_J_next;
	Mapping gt_C_l;
	Mapping gt_C_r;
	Mapping gt_C_next;
	Mapping gt_C_l_lmax;
	Mapping gt_C_r_lmax;

	tuple<rpos_t, rpos_t, rpos_t, rpos_t, int, string> GT_start_end(const string& query_id) {
		auto parsed = ParsedQueryId::parse(query_id);
		assert(parsed.valid);
		
		// Find index i where tidx.T[i].seq matches parsed.segm_id
		int segm_id = -1;
		for (int i = 0; i < (int)tidx.T.size(); i++)
			if (tidx.T[i].name == parsed.segm_id) {
				segm_id = i;
				break;
			}
		assert(segm_id >= 0 && segm_id < (rpos_t)tidx.T.size());

		if (abs_pos)
			return {parsed.start_pos, parsed.end_pos, parsed.start_pos, parsed.end_pos, segm_id, parsed.segm_id};
		else {
			const auto &segm = tidx.T[segm_id];
			qpos_t start = lower_bound(segm.kmers.begin(), segm.kmers.end(), parsed.start_pos, [](const auto &kmer, const auto &pos) { return kmer.r < pos; }) - segm.kmers.begin();
			qpos_t end  = lower_bound(segm.kmers.begin(), segm.kmers.end(), parsed.end_pos, [](const auto &kmer, const auto &pos) { return kmer.r < pos; }) - segm.kmers.begin();
			return {parsed.start_pos, parsed.end_pos, start, end, segm_id, parsed.segm_id};
		}
	}

	string vec2str(vector<BucketLoc> buckets, int bucket_l) {
		stringstream res;
		res << std::fixed << std::setprecision(4);
		res << "{";
		double max_J = 0.0, max_C = 0.0;
		for (auto &b: buckets) {
			auto M = matcher.collect_matches(b, p_ht);
			auto best_J_mapping	= matcher.bestIncludedJaccard(M, P_sz, m-2, m+2, m);
			auto best_C_mapping = matcher.bestFixedLength(M, P_sz, 2*bucket_l, m, MappingMetric::CONTAINMENT_INDEX);
			res << "(" << tidx.T[b.segm_id].name << ", " << b.b << ", " << best_J_mapping.score() << ", " << best_C_mapping.score() << "),";

			max_J = max(max_J, best_J_mapping.score());
			max_C = max(max_C, best_C_mapping.score());	
		}
		res << "}";
		res << "\tmax_J: " << max_J << "\tmax_C: " << max_C;
		return res.str();
	}

public:
	AnalyseSimulatedReads(const string& query_id, const string &P, int P_sz, const h2cnt &diff_hist, int m, const h2seed_t &p_ht, const SketchIndex &tidx, Buckets<abs_pos> &B, const double theta)
	 : query_id(query_id), P(P), P_sz(P_sz), diff_hist(diff_hist), m(m), p_ht(p_ht), tidx(tidx), B(B), theta(theta), matcher(tidx, B, diff_hist) {
		bucket_l = B.get_bucket_halflen();

		// Ground-truth
		auto parsed_orig = ParsedQueryId::parse(query_id);
		const auto &segm = tidx.get_segment(parsed_orig.segm_id);
		gt_mapping.paf = MappingPAF(0, P_sz, parsed_orig.strand, segm.name.c_str(), segm.sz, segm.id, parsed_orig.start_pos, parsed_orig.end_pos);

		tie(gt_start_nucl, gt_end_nucl, start, end, segm_id, segm_name) = GT_start_end(query_id);
		//update(0, P_sz-1, start, end, tidx.T[segm_id], -1, -1, -1, nullptr, nullptr); // TODO: should it be prev(r) instead?
		gt_b_l 		= BucketLoc(segm_id, max(0, start/bucket_l-1));
		gt_b_r 		= BucketLoc(segm_id, start/bucket_l);
		gt_b_next 	= BucketLoc(segm_id, start/bucket_l+1);
		gt_M_l 		= matcher.collect_matches(gt_b_l, p_ht);
		gt_M_r 		= matcher.collect_matches(gt_b_r, p_ht);
		gt_M_next 	= matcher.collect_matches(gt_b_next, p_ht);
		gt_J_l 		= matcher.bestIncludedJaccard(gt_M_l, P_sz, end-start-2, end-start+2, m);
		gt_J_r 		= matcher.bestIncludedJaccard(gt_M_r, P_sz, end-start-2, end-start+2, m);
		gt_J_next 	= matcher.bestIncludedJaccard(gt_M_next, P_sz, end-start-2, end-start+2, m);
		gt_C_l 		= matcher.bestFixedLength(gt_M_l, P_sz, 2*bucket_l, m, MappingMetric::CONTAINMENT_INDEX);
		gt_C_r 		= matcher.bestFixedLength(gt_M_r, P_sz, 2*bucket_l, m, MappingMetric::CONTAINMENT_INDEX);
		gt_C_next 	= matcher.bestFixedLength(gt_M_next, P_sz, 2*bucket_l, m, MappingMetric::CONTAINMENT_INDEX);
		gt_C_l_lmax = matcher.bestFixedLength(gt_M_l, P_sz, bucket_l, m, MappingMetric::CONTAINMENT_INDEX);
		
		// Buckets with Jaccard and Containment index >= theta
		vector<BucketLoc> J_buckets;
		vector<BucketLoc> C_buckets;
		auto sorted_buckets = B.get_sorted_buckets();
		for (auto &[b, content] : sorted_buckets) {
			auto M = matcher.collect_matches(b, p_ht);
			auto mapping_J = matcher.bestIncludedJaccard(M, P_sz, end-start-2, end-start+2, m);
			auto mapping_C = matcher.bestFixedLength(M, P_sz, 2*bucket_l, m, MappingMetric::CONTAINMENT_INDEX);
			mapping_J.set_bucket(b);
			if (mapping_J.score() >= theta) J_buckets.push_back(b);
			if (mapping_C.score() >= theta) C_buckets.push_back(b);
		}

		//static int not_overlapping_with_gt = 0;
		//bool J_overlaps_with_gt = false;
		//for (const auto &b_J: J_buckets)
		//	if (do_overlap(query_id, b_J, tidx)) {
		//		J_overlaps_with_gt = true;
		//		break;
		//	}
		//if (!J_overlaps_with_gt)
		//	++not_overlapping_with_gt;	
		//cerr << "not_overlapping_with_gt: " << not_overlapping_with_gt << endl;

		//TODO: make sure these are the same
		//gt_C_l.bucket.b
		//gt_C_right.bucket.b

	}

 	void print_tsv(ostream &paulout) {
		// open a new file stream
		static bool is_first_row = true;
		if (is_first_row) {
			paulout << "query_id"		// read name including the ground-truth
					   "\tm"			// sketch size of the read
			           "\ttheta"		// theta threshold
					   "\thl"           // bucket's half-length
					   "\tsegm"         // GT chromosome name
					   "\tgt_l_bucket"  // GT left bucket
					   "\tgt_r_bucket"  // GT right bucket
					   "\tgt_next_bucket"  // GT right bucket
					   "\tgt_J_l"       // GT left bucket: max included Jaccard
					   "\tgt_J_r"       // GT right bucket: max included Jaccard
					   "\tgt_J_next"    // GT next bucket: max included Jaccard
					   "\tgt_C_l"       // GT left bucket: max containment index of mapping of bucket length (2*lmax)
					   "\tgt_C_r"       // GT right bucket: max containment index of mapping of bucket length (2*lmax)
					   "\tgt_C_next"    // GT next bucket: max containment index of mapping of bucket length (2*lmax)
					   "\tgt_C_l_lmax"  // GT left bucket: max containment index of mapping of length lmax
					   "\tgt_C_r_lmax"  // GT right bucket: max containment index of mapping of length lmax
					   "\t#J>theta"     // number of buckets that contain a mapping with Jaccard >= theta
					   "\t#C>theta"     // number of buckets with Containment index >= theta
					   "\tJ>theta"      // buckets that contain a mapping with Jaccard >= theta
					   "\tC>theta"      // buckets with Containment index >= theta                                     
					   "\tmaxJ"         // maximal Jaccard of a reported bucket
					   "\tmaxC"         // maximal Containment index of a reported bucket
					   "\tP"      		// the read sequence
					   << endl;
			is_first_row = false;
		}
		paulout 	<< query_id
			<< "\t" << m						
			<< "\t" << theta			
			<< "\t" << bucket_l       			
			<< "\t" << segm_name       			
			<< "\t" << gt_b_l.b					
			<< "\t" << gt_b_r.b					
			<< "\t" << gt_b_next.b					
			<< "\t" << gt_J_l.score()					
			<< "\t" << gt_J_r.score()					
			<< "\t" << gt_J_next.score()					
			<< "\t" << gt_C_l.score()					
			<< "\t" << gt_C_r.score()					
			<< "\t" << gt_C_next.score()					
			<< "\t" << gt_C_l_lmax.score()			
			<< "\t" << gt_C_r_lmax.score()			
			<< "\t" << J_buckets.size()			
			<< "\t" << C_buckets.size()         
			<< "\t" << vec2str(J_buckets, bucket_l)
			<< "\t" << vec2str(C_buckets, bucket_l)
			<< "\t" << P
			<< endl;	
	}

	void print_paf(ostream &out) {
		out 
			<< "\tgt_segm:s:" << segm_name
			<< "\tgt_start_nucl:i:" << gt_start_nucl
			<< "\tgt_end_nucl:i:" << gt_end_nucl
			<< "\tbucket_l:i:" << bucket_l
			<< "\tgt_b_l:s:" << gt_b_l
			<< "\tgt_b_r:s:" << gt_b_r
			<< "\tgt_b_next:s:" << gt_b_next
			<< "\tgt_M_l:i:" << gt_M_l.size()					
			<< "\tgt_M_r:i:" << gt_M_r.size()					
			<< "\tgt_M_next:i:" << gt_M_next.size()					
			<< "\tgt_J_l:f:" << gt_J_l.score()
			<< "\tgt_J_r:f:" << gt_J_r.score()
			<< "\tgt_J_next:f:" << gt_J_next.score()
			<< "\tgt_C_l:f:" << gt_C_l.score()
			<< "\tgt_C_r:f:" << gt_C_r.score()
			<< "\tgt_C_next:f:" << gt_C_next.score()
			<< "\tgt_C_l_lmax:f:" << gt_C_l_lmax.score()
			<< "\tgt_C_r_lmax:f:" << gt_C_r_lmax.score();
			//<< "\tgt_C_l_overlap:f:" << Mapping::overlap(m, gt_C_l)
			//<< "\tgt_C_r_overlap:f:" << Mapping::overlap(m, gt_C_r)
			//<< "\tgt_C_next_overlap:f:" << Mapping::overlap(m, gt_C_next)
			//<< endl;	
	}

};

} // namespace sweepmap