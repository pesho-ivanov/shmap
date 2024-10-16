#pragma once

#include <set>
#include <map>

#include "index.h"
#include "io.h"
#include "sketch.h"
#include "utils.h"
#include "handler.h"
#include "mapper.h"

namespace sweepmap {
	
class JaccMapper : public Mapper {
	const SketchIndex &tidx;
	Handler *H;

    using Kmer = Seed;
	using hist_t = vector<int>;
	using Kmers = vector<Kmer>;
	using Matches = vector<Match>;
	using Intervals = vector<pair<int, int>>;

	Kmers select_seeds(sketch_t& p) {
		H->T.start("collect_kmer_info");
		Kmers kmers;
		kmers.reserve(p.size());

		for (int ppos = 0; ppos < (int)p.size(); ++ppos) {
			const auto &kmer = p[ppos];
			const auto count = tidx.count(kmer.h);
			if (count > 0)  // TODO: comment out!
				kmers.push_back(Kmer(kmer, p[ppos].r, p[ppos].r, count, kmers.size()));
		}
		H->T.stop("collect_kmer_info");
        H->C.inc("collected_kmers", kmers.size());

		H->T.start("sort_kmers");
		pdqsort_branchless(kmers.begin(), kmers.end(), [](const Kmer &a, const Kmer &b) {
            return a.hits_in_T < b.hits_in_T;
		});
		H->T.stop("sort_kmers");

		return kmers;
	}

	void sweep(vector<Match> &M, const pos_t P_sz, const Kmers &kmers, vector<Mapping> *mappings) {
		unordered_map<int, int> diff_hist;
        //unordered_multiset<int> diff_hist;

		for (auto &kmer: kmers)
			++diff_hist[kmer.kmer.h];

		int intersection = 0;
		int same_strand_seeds = 0;  // positive for more overlapping strands (fw/fw or bw/bw); negative otherwise

		H->T.start("sweep-sort");
		sort(M.begin(), M.end(), [](const Match &a, const Match &b) { return a.hit.r < b.hit.r; }); 
		H->T.stop("sweep-sort");

#ifdef DEBUG
		for (auto &kmer: kmers) {
			assert(diff_hist[kmer.kmer.h] > 0);
			cerr << "kmer: " << kmer.kmer.h << " " << kmer.kmer.r << " " << kmer.hits_in_T << endl;
		}

		for (int i=1; i<(int)M.size(); i++) {
			assert(M[i-1].hit.r <= M[i].hit.r);
			cerr << "match: " << M[i].hit.segm_id << " " << M[i].hit.r << " " << M[i].seed.kmer.h << " " << M[i].seed.kmer.r << endl;
		}
#endif

		H->T.start("sweep-main");
		// Increase the left point end of the window [l,r) one by one. O(matches)
		for(auto l = M.begin(), r = M.begin(); l != M.end(); ++l) {
			// Increase the right end of the window [l,r) until it gets out.
			//cerr << "l: " << l->hit.segm_id << " " << l->hit.r << " " << l->seed.kmer.h << " " << l->seed.kmer.r << endl;
			for(;  r != M.end()
				&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				&& r->hit.r + H->params.k <= l->hit.r + P_sz
				; ++r) {
				same_strand_seeds += r->is_same_strand() ? +1 : -1;  // change to r inside the loop
				if (--diff_hist[r->seed.kmer.h] >= 0)
					++intersection;
				assert (l->hit.r <= r->hit.r);
			}

			auto m = Mapping(H->params.k, P_sz, kmers.size(), l->hit.r, prev(r)->hit.r, l->hit.segm_id, pos_t(r-l), intersection, same_strand_seeds, l, r);
			if (m.J >= H->params.tThres)
				mappings->push_back(m);

			same_strand_seeds -= l->is_same_strand() ? +1 : -1;
			if (++diff_hist[l->seed.kmer.h] >= 1)
				--intersection;

			assert(intersection >= 0);
		}
		H->T.stop("sweep-main");
		assert(intersection == 0);
		assert(same_strand_seeds == 0);
	}

  public:
	JaccMapper(const SketchIndex &tidx, Handler *H) : tidx(tidx), H(H) {
        H->C.inc("seeds_limit_reached", 0);
        H->C.inc("mapped_reads", 0);
        if (H->params.tThres < 0.0 || H->params.tThres > 1.0) {
            cerr << "tThres = " << H->params.tThres << " outside of [0,1]." << endl;
            exit(1);
        }
    }

	void map(const string &pFile) {
		cerr << "Mapping reads using JaccMap " << "..." << endl;

		H->C.inc("kmers", 0);
		H->C.inc("seeds", 0);
		H->C.inc("matches", 0);
		H->C.inc("matches_infreq", 0);
		H->C.inc("matches_freq", 0);
		H->C.inc("spurious_matches", 0);
		H->C.inc("J", 0);
		H->C.inc("mappings", 0);
		H->C.inc("sketched_kmers", 0);
		H->C.inc("total_edit_distance", 0);
		H->C.inc("intersection_diff", 0);

		H->T.start("mapping");
		H->T.start("query_reading");

		read_fasta_klib(pFile, [this](kseq_t *seq) {
			H->C.inc("reads");
			H->T.stop("query_reading");
			H->T.start("query_mapping");
				H->T.start("sketching");
					sketch_t p = H->sketcher.sketch(seq->seq.s);
					H->T.stop("sketching");

				H->T.start("prepare");
					string query_id = seq->name.s;
					pos_t P_sz = (pos_t)seq->seq.l;

					H->C.inc("kmers", p.size());
					H->C.inc("read_len", P_sz);
					Timer read_mapping_time;  // TODO: change to H->T.start
					read_mapping_time.start();
				H->T.stop("prepare");

				H->T.start("seeding");
					Kmers kmers = select_seeds(p);
					//assert(kmers.size() == p.size());
					H->T.stop("seeding");
					H->C.inc("kmers", kmers.size());

				int lmax = int(kmers.size() / H->params.tThres);			// maximum length of a similar mapping
				int S = int((1.0 - H->params.tThres) * kmers.size()) + 1;	// any similar mapping includes at least 1 seed match
				std::unordered_map<int, int> M;  							// M[b] -- #matched kmers[0...i] in [bl, (b+2)l)

				H->T.start("match_infrequent");
				for (int i = 0; i < S; i++) {
					Kmer seed = kmers[i];
					if (seed.hits_in_T > 0) {
						Matches matches_infreq;
						std::unordered_set<int> matched_buckets;
						tidx.get_matches(&matches_infreq, seed);
						for (auto &m: matches_infreq) {
							matched_buckets.insert(int(m.hit.tpos / lmax));
							matched_buckets.insert(int(m.hit.tpos / lmax) + 1);
						}
						for (int b: matched_buckets)
							++M[b];
					}
				}
				H->T.stop("match_infrequent");

				//{
				//	int total = 0; for (auto &[b, cnt]: M) total += cnt;	
				//	cerr << "Before. buckets: " << M.size() << ", total: " << total << endl;
				//}

				H->T.start("match_frequent");
				for (int i = S; i < (int)kmers.size(); i++) {
					vector<int> to_erase;

					for (auto &[b, cnt]: M) {
						if (tidx.is_kmer_in_t_interval(kmers[i], b*lmax, (b+2)*lmax))
							++cnt;
						if (cnt <= i-S+1) {
							to_erase.push_back(b);
						}
					}
							
					for (auto b: to_erase)
						M.erase(b);
				}
				H->T.stop("match_frequent");

				//{
				//	int total = 0; for (auto &[b, cnt]: M) total += cnt;	
				//	cerr << "After. buckets: " << M.size() << ", total: " << total << endl;
				//}

				vector<Mapping> mappings;
				int total_matches = 0; //matches.size();
				for (const auto &[b, cnt]: M) {
					H->T.start("match_collect");
					Matches matches;
					for (const auto &kmer: kmers)
						tidx.get_matches_in_t_interval(&matches, kmer, b*lmax, (b+2)*lmax);
					H->T.stop("match_collect");
					total_matches += matches.size();

					H->T.start("sweep");
						sweep(matches, P_sz, kmers, &mappings);
					H->T.stop("sweep");
				}

				if (mappings.size() > 0) {
					H->C.inc("mapped_reads");
					//H->C.inc("intersection_diff", t_abs - mappings[0].intersection);
				}

				if (H->params.onlybest && mappings.size() > 0) {
					Mapping best = *mappings.begin(); 
					for (auto &m: mappings) {
						//cerr << m << endl;
						if (m.J > best.J)
							best = m;
					}
					
					mappings.clear();
					mappings.push_back(best);
				}

				H->T.start("output");
					read_mapping_time.stop();

					for (auto &m: mappings) {
						m.map_time = read_mapping_time.secs() / (double)mappings.size();
						const auto &segm = tidx.T[m.segm_id];
		//				if (H->params.sam) {
		//					auto ed = m.print_sam(query_id, segm, (int)matches.size(), seq->seq.s, seq->seq.l);
		//					H->C.inc("total_edit_distance", ed);
		//				}
		//				else
							m.print_paf(query_id, segm, total_matches);
						//  H->C.inc("spurious_matches", spurious_matches(m, matches));
						H->C.inc("J", int(10000.0*m.J));
						H->C.inc("mappings");
					}
		//			H->C.inc("matches", matches.size());
					//H->T.stop("postproc");

				H->T.stop("output");
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
		cerr << "Mapping:" << endl;
		cerr << " | Total reads:           " << H->C.count("reads") << " (~" << 1.0*H->C.count("read_len") / H->C.count("reads") << " nb per read)" << endl;
		cerr << " |  | mapped:               " << H->C.count("mapped_reads") << " (" << H->C.perc("mapped_reads", "reads") << "%)" << endl;
		cerr << " |  |  | intersect. diff:     " << H->C.frac("intersection_diff", "mapped_reads") << " per mapped read" << endl;
		cerr << " | Read kmers (total):    " << H->C.count("kmers") << " (" << H->C.frac("kmers", "reads") << " per read)" << endl;
		cerr << " |  | unique:                 " << H->C.count("seeds") << " (" << H->C.frac("seeds", "kmers") << ")" << endl;
		cerr << " | Matches:               " << H->C.count("matches") << " (" << H->C.frac("matches", "reads") << " per read)" << endl;
		cerr << " |  | infrequent:             " << H->C.count("matches_infreq") << " (" << H->C.perc("matches_infreq", "matches") << "%)" << endl;
		cerr << " |  | frequent:               " << H->C.count("matches_freq") << " (" << H->C.perc("matches_freq", "matches") << "%)" << endl;
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
        cerr << " | Runtime:                   "     << setw(5) << right << H->T.secs("mapping")       << " sec (" << setw(5) << right << H->T.range_ratio("query_mapping") << "x)" << endl; //setw(4) << right << H->C.count("reads") / H->T.secs("total") << " reads per sec)" << endl;
        cerr << " |  | load queries:           "     << setw(5) << right << H->T.secs("query_reading")     << " (" << setw(4) << right << H->T.perc("query_reading", "mapping")       << "\%, " << setw(5) << right << H->T.range_ratio("query_reading") << "x)" << endl;
        cerr << " |  | prepare:                "     << setw(5) << right << H->T.secs("prepare")           << " (" << setw(4) << right << H->T.perc("prepare", "mapping")             << "\%, " << setw(5) << right << H->T.range_ratio("prepare") << "x)" << endl;
        cerr << " |  | sketch reads:           "     << setw(5) << right << H->T.secs("sketching")         << " (" << setw(4) << right << H->T.perc("sketching", "mapping")           << "\%, " << setw(5) << right << H->T.range_ratio("sketching") << "x)" << endl;
        cerr << " |  | seeding:                "     << setw(5) << right << H->T.secs("seeding")           << " (" << setw(4) << right << H->T.perc("seeding", "mapping")             << "\%, " << setw(5) << right << H->T.range_ratio("seeding") << "x)" << endl;
        cerr << " |  |  | collect seed info:       " << setw(5) << right << H->T.secs("collect_kmer_info") << " (" << setw(4) << right << H->T.perc("collect_kmer_info", "seeding")   << "\%, " << setw(5) << right << H->T.range_ratio("collect_kmer_info") << "x)" << endl;
//        cerr << " |  |  | thin sketch:             " << setw(5) << right << H->T.secs("thin_sketch")       << " (" << setw(4) << right << H->T.perc("thin_sketch", "seeding")         << "\%, " << setw(5) << right << H->T.range_ratio("thin_sketch") << "x)" << endl;
//        cerr << " |  |  | unique seeds:            " << setw(5) << right << H->T.secs("unique_seeds")      << " (" << setw(4) << right << H->T.perc("unique_seeds", "seeding")        << "\%, " << setw(5) << right << H->T.range_ratio("unique_seeds") << "x)" << endl;
        cerr << " |  |  | sort seeds:              " << setw(5) << right << H->T.secs("sort_kmers")        << " (" << setw(4) << right << H->T.perc("sort_kmers", "seeding")          << "\%, " << setw(5) << right << H->T.range_ratio("sort_kmers") << "x)" << endl;
        cerr << " |  | match infrequent:        " << setw(5) << right << H->T.secs("match_infrequent")  << " (" << setw(4) << right << H->T.perc("match_infrequent", "mapping")   << "\%, " << setw(5) << right << H->T.range_ratio("match_infrequent") << "x)" << endl;
        cerr << " |  | match frequent:          " << setw(5) << right << H->T.secs("match_frequent")    << " (" << setw(4) << right << H->T.perc("match_frequent", "mapping")     << "\%, " << setw(5) << right << H->T.range_ratio("match_frequent") << "x)" << endl;
        cerr << " |  | matches collect:         " << setw(5) << right << H->T.secs("match_collect")     << " (" << setw(4) << right << H->T.perc("match_collect", "mapping")      << "\%, " << setw(5) << right << H->T.range_ratio("match_collect") << "x)" << endl;
//        cerr << " |  |  | get intervals:           " << setw(5) << right << H->T.secs("get_intervals")     << " (" << setw(4) << right << H->T.perc("get_intervals", "mapping")      << "\%, " << setw(5) << right << H->T.range_ratio("get_intervals") << "x)" << endl;
//        cerr << " |  |  | filter promising buckets:" << setw(5) << right << H->T.secs("filter_promising_buckets")    << " (" << setw(4) << right << H->T.perc("filter_promising_buckets", "mapping")     << "\%, " << setw(5) << right << H->T.range_ratio("filter_promising_buckets") << "x)" << endl;
        cerr << " |  | sweep:                  "     << setw(5) << right << H->T.secs("sweep")             << " (" << setw(4) << right << H->T.perc("sweep", "mapping")               << "\%, " << setw(5) << right << H->T.range_ratio("sweep") << "x)" << endl;
        cerr << " |  |  | sort matches:            " << setw(5) << right << H->T.secs("sweep-sort")        << " (" << setw(4) << right << H->T.perc("sweep-sort", "sweep")            << "\%, " << setw(5) << right << H->T.range_ratio("sweep-sort") << "x)" << endl;
        cerr << " |  |  | main:                    " << setw(5) << right << H->T.secs("sweep-main")        << " (" << setw(4) << right << H->T.perc("sweep-main", "sweep")            << "\%, " << setw(5) << right << H->T.range_ratio("sweep-main") << "x)" << endl;
        cerr << " |  | prepare:                "     << setw(5) << right << H->T.secs("output")            << " (" << setw(4) << right << H->T.perc("output", "mapping")              << "\%, " << setw(5) << right << H->T.range_ratio("output") << "x)" << endl;
//        cerr << " |  | post proc:              "     << setw(5) << right << H->T.secs("postproc")          << " (" << setw(4) << right << H->T.perc("postproc", "mapping")            << "\%, " << setw(5) << right << H->T.range_ratio("postproc") << "x)" << endl;
    }
};

}  // namespace sweepmap