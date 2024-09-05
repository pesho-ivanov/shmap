#pragma once

#include <algorithm>
#include <cassert>
#include <deque>
#include <iomanip>
#include <set>

//#include "../ext/edlib.h"
#include "../ext/pdqsort.h"

#include "index.h"
#include "io.h"
#include "sketch.h"
#include "utils.h"
#include "handler.h"
#include "mapper.h"

namespace sweepmap {

using std::setw;
using std::right;
using std::vector;

class BucketMapper : public Mapper {
	const SketchIndex &tidx;
	Handler *H;

	using hist_t = vector<int>;
	using bucket_t = int;
	using Buckets = unordered_map<bucket_t, vector<Match>>;

	vector<Seed> select_seeds(sketch_t& p) {
		H->T.start("unique_kmers");
		pdqsort_branchless(p.begin(), p.end(), [](const Kmer &a, const Kmer &b) {
            return a.h < b.h;
		});
		p.erase(std::unique(p.begin(), p.end()), p.end());
		H->T.stop("unique_kmers");

		H->T.start("collect_seed_info");
		vector<Seed> seeds;
		seeds.reserve(p.size());

		// TODO: limit The number of kmers in the pattern p
		for (int ppos = 0; ppos < (int)p.size(); ++ppos) {
			const auto &kmer = p[ppos];
			const auto count = tidx.count(kmer.h);
			if (count > 0)
				seeds.push_back(Seed(kmer, p[ppos].r, p[ppos].r, count));
		}
		H->T.stop("collect_seed_info");
        H->C.inc("collected_seeds", seeds.size());

		H->T.start("sort_seeds");
		pdqsort_branchless(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b) {
            return a.hits_in_T < b.hits_in_T;
		});
		H->T.stop("sort_seeds");

		return seeds;
	}

	bucket_t pos2bucket(pos_t pos, pos_t half_bucket, bool which_bucket) {
		return pos/half_bucket + (int)which_bucket;
	}

	pair<pos_t, pos_t> bucket2interval(bucket_t bucket, pos_t half_bucket) {
		return {bucket*half_bucket, (bucket+2)*half_bucket};
	}

	// Initializes the histogram of the pattern and the list of matches
	Buckets match_seeds(pos_t P_sz, pos_t p_sz, const vector<Seed> &seeds, double t_frac) {
		pos_t h = P_sz;
		int t_abs = t_frac * p_sz;
		Buckets M(h);  // M[bucket] - number of unique kmers from p starting at T[h*bucket, h*(bucket+2))
		size_t seed;

		H->T.start("match_infrequent");
		// match infrequent seeds to all buckets
		for (seed=0; seed<seeds.size(); seed++) {
			if ((int)seeds.size() - (int)seed < t_abs)	
				break;
			vector<Match> seed_matches;
			tidx.get_matches(&seed_matches, seeds[seed], seed);
			for (const auto &match: seed_matches)
				for (int which_bucket=0; which_bucket<2; ++which_bucket) {
					auto bucket = pos2bucket(match.hit.r, h, bool(which_bucket));
					M[bucket].push_back(match);
				}
		}
		H->T.stop("match_infrequent");

		H->T.start("match_frequent");
		// match frequent seeds to all buckets with an infrequent seed
		for (auto &[bucket, matches]: M) {
			auto bucket_interval = bucket2interval(bucket, h);
			for (auto seed_it = seeds.begin()+seed; seed_it != seeds.end(); ++seed_it)
				for (const auto &match: tidx.get_matches_in_interval(*seed_it, bucket_interval.first, bucket_interval.second))
					matches.push_back(match);
		}
		H->T.stop("match_frequent");
		
		// filter the promising buckets
		H->T.start("filter_promising_buckets");
		Buckets promising_buckets;
		for (const auto &[bucket, matches]: M)
			if ((int)matches.size() >= t_abs)
				promising_buckets[bucket] = matches;
		H->T.stop("filter_promising_buckets");
		return promising_buckets;
	}

	vector<Mapping> sweep(Buckets &promising_buckets, pos_t P_len, const vector<Seed> &seeds) {
		vector<Mapping> mappings;	// List of tripples <i, j, score> of matches
		if (promising_buckets.empty())
			return mappings;

		Mapping best(H->params.k, P_len, 0, -1, -1, -1, -1, -1, 0, promising_buckets.begin()->second.begin(), promising_buckets.begin()->second.end());
		Mapping second = best;

		for (auto &[bucket, matches]: promising_buckets) {
			sort(matches.begin(), matches.end(), [](const Match &a, const Match &b) {
				return a.hit.r < b.hit.r;
			}); 
			int same_strand_seeds = 0;
			for(const auto &m: matches)
				same_strand_seeds += m.is_same_strand() ? +1 : -1;
			auto thin_seeds_cnt = seeds.size();
			auto l = matches.front().hit.r;
			auto r = matches.back().hit.r;
			auto segm_id = matches.front().hit.segm_id;
			auto intersection = matches.size();

			auto m = Mapping(H->params.k, P_len, thin_seeds_cnt, l, r, segm_id, matches.size(), intersection, same_strand_seeds, matches.begin(), matches.end());
			//Mapping(       int k,    pos_t P_sz, int seeds, pos_t T_l, pos_t T_r, segm_t segm_id, pos_t s_sz, int intersection, int same_strand_seeds, std::vector<Match>::const_iterator l, std::vector<Match>::const_iterator r)

			//cerr << "Bucket: " << bucket << " " << m << endl;

			if (H->params.onlybest) {
				if (m.intersection > best.intersection) {  // if (intersection > best.intersection)
					if (best.T_l < m.T_l - 0.9*P_len)
						second = best;
					best = m;
				} else if (m.intersection > second.intersection && m.T_l > best.T_l + 0.9*P_len) {
					second = m;
				}
			} else {
				if (m.J > H->params.tThres) {
					mappings.push_back(m);
				}
			}
		}

		if (H->params.onlybest && best.intersection != -1) { // && best.J > H->params.tThres)
			best.mapq = (best.intersection > 5 && best.J > 0.1 && best.J > 1.2*second.J) ? 60 : 0;
			best.J2 = second.J;
			mappings.push_back(best);
		}

		return mappings;
	}

  public:
	BucketMapper(const SketchIndex &tidx, Handler *H) : tidx(tidx), H(H) {
        H->C.inc("seeds_limit_reached", 0);
        H->C.inc("unmapped_reads", 0);
        if (H->params.tThres < 0.0 || H->params.tThres > 1.0) {
            cerr << "tThres = " << H->params.tThres << " outside of [0,1]." << endl;
            exit(1);
        }
    }

	void map(const string &pFile) {
		cerr << "Mapping reads using BucketMap " << pFile << "..." << endl;

		H->C.inc("spurious_matches", 0);
		H->C.inc("J", 0);
		H->C.inc("mappings", 0);
		H->C.inc("sketched_kmers", 0);
		H->C.inc("total_edit_distance", 0);

		H->T.start("mapping");
		H->T.start("query_reading");
		read_fasta_klib(pFile, [this](kseq_t *seq) {
			H->T.stop("query_reading");
			H->T.start("query_mapping");
			H->T.start("sketching");
			sketch_t p = H->sketcher.sketch(seq->seq.s);
			H->T.stop("sketching");

			string query_id = seq->name.s;
			pos_t P_sz = (pos_t)seq->seq.l;

			H->C.inc("read_len", P_sz);

			Timer read_mapping_time;
			read_mapping_time.start();
			H->T.start("seeding");
			vector<Seed> seeds = select_seeds(p);
			H->T.stop("seeding");

			H->T.start("matching");
			Buckets promising_buckets = match_seeds(P_sz, p.size(), seeds, H->params.tThres);
			H->T.stop("matching");

			H->T.start("sweep");
			vector<Mapping> mappings = sweep(promising_buckets, P_sz, seeds);
			H->T.stop("sweep");

			std::vector<Match> all_matches;
			for (const auto &[bucket, matches]: promising_buckets)
				all_matches.insert(all_matches.end(), matches.begin(), matches.end());

			for (auto &m: mappings) {
				m.map_time = read_mapping_time.secs() / (double)mappings.size();
				const auto &segm = tidx.T[m.segm_id];
//				if (H->params.sam) {
//					auto ed = m.print_sam(query_id, segm, (int)matches.size(), seq->seq.s, seq->seq.l);
//					H->C.inc("total_edit_distance", ed);
//				}
//				else

					m.print_paf(query_id, segm, all_matches);
				//  H->C.inc("spurious_matches", spurious_matches(m, matches));
				H->C.inc("J", int(10000.0*m.J));
				H->C.inc("mappings");
				H->C.inc("sketched_kmers", m.p_sz);
			}
//			H->C.inc("matches", matches.size());
			H->C.inc("reads");
			if (mappings.empty())
				H->C.inc("unmapped_reads");
			//H->T.stop("postproc");

			H->T.stop("query_mapping");
			H->T.start("query_reading");
		});
		H->T.stop("query_reading");
		H->T.stop("mapping");

		print_stats();
		print_time_stats();
	}
	
	void print_stats() {
		cerr << std::fixed << std::setprecision(1);
		cerr << "Mapping stats:" << endl;
		cerr << " | Total reads:           " << H->C.count("reads") << " (~" << 1.0*H->C.count("read_len") / H->C.count("reads") << " nb per read)" << endl;
		cerr << " |  | Unmapped reads:        " << H->C.count("unmapped_reads") << " (" << H->C.perc("unmapped_reads", "reads") << "%)" << endl;
		cerr << " | Sketched read kmers:   " << H->C.count("sketched_kmers") << " (" << H->C.frac("sketched_kmers", "reads") << " per read)" << endl;
		//cerr << " | Kmer matches:          " << H->C.count("matches") << " (" << H->C.frac("matches", "reads") << " per read)" << endl;
		//cerr << " | Seed limit reached:    " << H->C.count("seeds_limit_reached") << " (" << H->C.perc("seeds_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Matches limit reached: " << H->C.count("matches_limit_reached") << " (" << H->C.perc("matches_limit_reached", "reads") << "%)" << endl;
		//cerr << " | Spurious matches:      " << H->C.count("spurious_matches") << " (" << H->C.perc("spurious_matches", "matches") << "%)" << endl;
		//cerr << " | Discarded seeds:       " << H->C.count("discarded_seeds") << " (" << H->C.perc("discarded_seeds", "collected_seeds") << "%)" << endl;
		//cerr << " | Average Jaccard:       " << H->C.frac("J", "mappings") / 10000.0 << endl;
		//cerr << " | Average edit dist:     " << H->C.frac("total_edit_distance", "mappings") << endl;
		//print_time_stats();
        //printMemoryUsage();
	}

    void print_time_stats() {
        cerr << std::fixed << std::setprecision(1);
        cerr << " | Map:                   "         << setw(5) << right << H->T.secs("mapping")           << " (" << setw(5) << right << H->T.range_ratio("query_mapping") << "x)" << endl; //setw(4) << right << H->C.count("reads") / H->T.secs("total") << " reads per sec)" << endl;
        cerr << " |  | load queries:           "     << setw(5) << right << H->T.secs("query_reading")     << " (" << setw(4) << right << H->T.perc("query_reading", "mapping")       << "\%, " << setw(5) << right << H->T.range_ratio("query_reading") << "x)" << endl;
        cerr << " |  | sketch reads:           "     << setw(5) << right << H->T.secs("sketching")         << " (" << setw(4) << right << H->T.perc("sketching", "mapping")           << "\%, " << setw(5) << right << H->T.range_ratio("sketching") << "x)" << endl;
        cerr << " |  | seeding:                "     << setw(5) << right << H->T.secs("seeding")           << " (" << setw(4) << right << H->T.perc("seeding", "mapping")             << "\%, " << setw(5) << right << H->T.range_ratio("seeding") << "x)" << endl;
        cerr << " |  |  | collect seed info:       " << setw(5) << right << H->T.secs("collect_seed_info") << " (" << setw(4) << right << H->T.perc("collect_seed_info", "seeding")   << "\%, " << setw(5) << right << H->T.range_ratio("collect_seed_info") << "x)" << endl;
//        cerr << " |  |  | thin sketch:             " << setw(5) << right << H->T.secs("thin_sketch")       << " (" << setw(4) << right << H->T.perc("thin_sketch", "seeding")         << "\%, " << setw(5) << right << H->T.range_ratio("thin_sketch") << "x)" << endl;
//        cerr << " |  |  | unique seeds:            " << setw(5) << right << H->T.secs("unique_seeds")      << " (" << setw(4) << right << H->T.perc("unique_seeds", "seeding")        << "\%, " << setw(5) << right << H->T.range_ratio("unique_seeds") << "x)" << endl;
        cerr << " |  |  | sort seeds:              " << setw(5) << right << H->T.secs("sort_seeds")        << " (" << setw(4) << right << H->T.perc("sort_seeds", "seeding")          << "\%, " << setw(5) << right << H->T.range_ratio("sort_seeds") << "x)" << endl;
        cerr << " |  | matching seeds:         "     << setw(5) << right << H->T.secs("matching")          << " (" << setw(4) << right << H->T.perc("matching", "mapping")            << "\%, " << setw(5) << right << H->T.range_ratio("matching") << "x)" << endl;
        cerr << " |  |  | match infrequent:        " << setw(5) << right << H->T.secs("match_infrequent")  << " (" << setw(4) << right << H->T.perc("match_infrequent", "matching")   << "\%, " << setw(5) << right << H->T.range_ratio("match_infrequent") << "x)" << endl;
        cerr << " |  |  | match frequent:          " << setw(5) << right << H->T.secs("match_frequent")    << " (" << setw(4) << right << H->T.perc("match_frequent", "matching")     << "\%, " << setw(5) << right << H->T.range_ratio("match_frequent") << "x)" << endl;
        cerr << " |  |  | filter promising buckets:" << setw(5) << right << H->T.secs("filter_promising_buckets")    << " (" << setw(4) << right << H->T.perc("filter_promising_buckets", "matching")     << "\%, " << setw(5) << right << H->T.range_ratio("filter_promising_buckets") << "x)" << endl;
        cerr << " |  | sweep:                  "     << setw(5) << right << H->T.secs("sweep")             << " (" << setw(4) << right << H->T.perc("sweep", "mapping")               << "\%, " << setw(5) << right << H->T.range_ratio("sweep") << "x)" << endl;
//        cerr << " |  | post proc:              "     << setw(5) << right << H->T.secs("postproc")          << " (" << setw(4) << right << H->T.perc("postproc", "mapping")            << "\%, " << setw(5) << right << H->T.range_ratio("postproc") << "x)" << endl;
    }
};

}  // namespace sweepmap