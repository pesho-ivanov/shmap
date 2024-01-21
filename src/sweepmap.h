#ifndef SWEEP_HPP
#define SWEEP_HPP

#include <algorithm>
#include <cassert>
#include <deque>
#include <iomanip>
#include <set>

#include "io.h"
#include "sketch.h"

struct Seed {
	pos_t P_r;
	hash_t kmer;
	const vector<Hit> *hits_in_T;

	Seed() {}
	Seed(pos_t P_r, hash_t kmer, const vector<Hit> *hits_in_T) :
		P_r(P_r), kmer(kmer), hits_in_T(hits_in_T) {}
};

struct Match {
	hash_t kmer;
	pos_t P_r;
	Hit hit;
	//pos_t T_r;
	//pos_t t_pos;
	//bool strand;  // 0 - forward, 1 - reverse

	//char get_strand() const {
	//	return strand ? '-' : '+';
	//}

	//Match() {}
	Match(hash_t kmer, pos_t P_r, const Hit &hit) :
		kmer(kmer), P_r(P_r), hit(hit) {}

	//void print() const {
	//	cout << "Match: " << kmer << " " << P_r << " " << T_r << " " << t_pos << endl;
	//}
};

struct Mapping {
	int k; 	   // kmer size
	pos_t P_sz;     // pattern size |P| bp 
	pos_t seeds;     // number of seeds (subset of the sketch kmers)
	int matches;  // M.size() -- total number of matches in `t' 
	pos_t T_l;      // the position of the leftmost nucleotide of the mapping
	pos_t T_r;      // the position of the rightmost nucleotide of the mapping
	int segm_id;
	pos_t s_sz;      // the position of the rightmost nucleotide of the mapping
	int xmin;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	double sim;     // the number of kmers in the intersection between the pattern and its mapping in `t'
	double map_time;
	char strand;    // '+' or '-'

    Mapping() {}
	Mapping(int k, pos_t P_sz, int seeds, pos_t T_l, pos_t T_r, int segm_id, pos_t s_sz, int xmin)
		: k(k), P_sz(P_sz), seeds(seeds), T_l(T_l), T_r(T_r), segm_id(segm_id), s_sz(s_sz), xmin(xmin), sim(double(xmin) / std::max(seeds, s_sz)) {}

	// --- https://github.com/lh3/miniasm/blob/master/PAF.md ---
    void print_paf(const string &query_id, const RefSegment &segm, const int matches) const {
		cout << query_id  			// Query sequence name
			<< "\t" << P_sz     // query sequence length
			<< "\t" << 0   // query start (0-based; closed)
			<< "\t" << P_sz  // query end (0-based; open)
			<< "\t" << strand   // '+' or '-'
			<< "\t" << segm.name // reference segment name
			<< "\t" << segm.seq.size() // T_sz -- target sequence length
			<< "\t" << T_l  // target start on original strand (0-based)
			<< "\t" << T_r  // target start on original strand (0-based)
			<< "\t" << P_sz  // TODO: fix; Number of residue matches (number of nucleotide matches)
			<< "\t" << P_sz  // TODO: fix; Alignment block length: total number of sequence matches, mismatches, insertions and deletions in the alignment
			<< "\t" << 60  // Mapping quality (0-255; 255 for missing)
// ----- end of required PAF fields -----
			<< "\t" << "k:i:" << k
			<< "\t" << "p:i:" << seeds  // sketches
			<< "\t" << "M:i:" << matches // kmer matches in T
			<< "\t" << "s:i:" << s_sz
			<< "\t" << "I:i:" << xmin  // intersection of `p` and `s` [kmers]
			<< "\t" << "J:f:" << sim   // similarity [0; 1]
			<< "\t" << "t:f:" << map_time
			<< endl;
	}
};

class SweepMap {
	const SketchIndex &tidx;
	const params_t &params;
	Timers *T;
	Counters *C;

	using hist_t = ankerl::unordered_dense::map<hash_t, int>;

	vector<Seed> choose_seeds(const Sketch& p, int max_seeds, hist_t *hist) {
		vector<Seed> seeds;
		seeds.reserve(p.size());

		// TODO: limit the number of kmers in the pattern p
		// TODO: Add pattern kmers from the pattern to p_hist 
		// TODO: make p_hist smaller since we don't use all seeds
		T->start("collect_seed_info");
		for (const auto &curr: p) {
			auto hist_it = hist->find(curr.kmer);
			if (hist_it != hist->end()) {
				++(hist_it->second);
			} else {
				(*hist)[curr.kmer] = 1;
				auto t_it = tidx.h2pos.find(curr.kmer);
				if (t_it != tidx.h2pos.end())
					seeds.push_back(Seed(curr.r, curr.kmer, &t_it->second));
			}
		}
		T->stop("collect_seed_info");

		T->start("sort_seeds");
		sort(seeds.begin(), seeds.end(), [](const Seed &a, const Seed &b) {
			// Sort the matching kmers by number of hits in the reference
			return a.hits_in_T->size() < b.hits_in_T->size();
		});
		T->stop("sort_seeds");

		if ((int)seeds.size() > max_seeds)
			seeds.resize(max_seeds);

		return seeds;
	}

	// Initializes the histogram of the pattern and the list of matches
	vector<Match> match_seeds(pos_t p_sz, const vector<Seed> &seeds) {
		T->start("collect_matches");
		// Get MAX_SEEDS of kmers with the lowest number of hits.
		vector<Match> matches;
		matches.reserve(P_MULTIPLICITY * p_sz);

		int total_hits = 0;
		for (const auto &seed: seeds) {
			for (const auto &hit: *(seed.hits_in_T)) {
				assert(hit.segm_id >= 0 && (size_t)hit.segm_id < tidx.T.size());
				matches.push_back(Match(seed.kmer, seed.P_r, hit));
			}
			if ((total_hits += (int)seed.hits_in_T->size()) > (int)params.max_matches) {
				C->inc("matches_limit_reached");
				break;
			}
		}
		T->stop("collect_matches");

		T->start("sort_matches");
		sort(matches.begin(), matches.end(), [](const Match &a, const Match &b) {
			// Preparation for sweeping: sort M by ascending positions within one reference segment.
			if (a.hit.segm_id != b.hit.segm_id)
				return a.hit.segm_id < b.hit.segm_id;
			return a.hit.r < b.hit.r;
		});
		T->stop("sort_matches");

		return matches;
	}

	// vector<hash_t> diff_hist;  // diff_hist[kmer_hash] = #occurences in `p` - #occurences in `s`
	// vector<Match> M;   	   // for all kmers from P in T: <kmer_hash, last_kmer_pos_in_T> * |P| sorted by second
	const vector<Mapping> sweep(hist_t &diff_hist, const vector<Match> &M, const pos_t P_len, const int seeds) {
		vector<Mapping> mappings;	// List of tripples <i, j, score> of matches

		int xmin = 0;
		Mapping best(params.k, P_len, 0, -1, -1, -1, -1, -1);
		//best.k = params.k;
		//best.P_sz = P_len;
		//best.p_sz = p_len;
		//best.xmin = 0;
		//best.sim = 0.0;

		// Increase the left point end of the window [l,r) one by one. O(matches)
		int i = 0, j = 0;
		for(auto l = M.begin(), r = M.begin(); l != M.end(); ++l, ++i) {
			// Increase the right end of the window [l,r) until it gets out.
			for(;  r != M.end()
				&& l->hit.segm_id == r->hit.segm_id   // make sure they are in the same segment since we sweep over all matches
				&& r->hit.r + params.k <= l->hit.r + P_len
				; ++r, ++j) {
				// If taking this kmer from T increases the intersection with P.
				if (--diff_hist[r->kmer] >= 0)
					++xmin;
				assert (l->hit.r <= r->hit.r);
			}

			auto m = Mapping(params.k, P_len, seeds, l->hit.r, prev(r)->hit.r, l->hit.segm_id, pos_t(r-l), xmin);
			if (params.onlybest) {
				//if (xmin > best.xmin) {
				if (m.sim > best.sim) {
					best = m;
				}
			} else {
				if (m.sim > params.tThres) {
					mappings.push_back(m);
				}
			}

			// Prepare for the next step by moving `l` to the right.
			if (++diff_hist[l->kmer] > 0)
				--xmin;

			assert(xmin >= 0);
		}
		assert (xmin == 0);

		if (params.onlybest && best.xmin != -1) { // && best.J > params.tThres)
			mappings.push_back(best);
		}

		return mappings;
	}

	// Return only reasonable matches (i.e. those that are not J-dominated by
	// another overlapping match). Runs in O(|all|).
	vector<Mapping> filter_reasonable(const vector<Mapping> &all, const pos_t P_len) {
		vector<Mapping> reasonable;
		deque<Mapping> recent;

		// Minimal separation between mappings to be considered reasonable
		pos_t sep = pos_t((1.0 - params.tThres) * double(P_len));

		// The deque `recent' is sorted decreasingly by J
		//					  _________`recent'_________
		//                   /                          \
		// ---------------- | High J ... Mid J ... Low J | current J
		// already removed    deque.back ... deque.front    to add next
		for (const auto &next: all) {
			// 1. Prepare for adding `curr' by removing from the deque back all
			//    mappings that are too far to the left. This keeps the deque
			//    within |P| from back to front. A mapping can become reasonable
			//    only after getting removed.
			while(!recent.empty() && next.T_l - recent.back().T_l > sep) {
				// If the mapping is not marked as unreasonable (coverted by a preivous better mapping)
				if (recent.back().matches != -1) {
					// Take the leftmost mapping.
					reasonable.push_back(recent.back());
					// Mark the next closeby mappings as not reasonable
					for (auto it=recent.rbegin(); it!=recent.rend() && it->T_l - recent.back().T_l < sep; ++it)
						it->matches = -1;
				}
				// Remove the mapping that is already too much behind.
				recent.pop_back();
			}
			assert(recent.empty() || (recent.back().T_r <= recent.front().T_r && recent.front().T_r <= next.T_r));

			// Now all the mappings in `recent' are close to `curr' 
			// 2. Remove from the deque front all mappings that are strictly
			//    less similar than the current J. This keeps the deque sorted
			//    descending in J from left to right
			while(!recent.empty() && recent.front().xmin < next.xmin)
				recent.pop_front();
			assert(recent.empty() || (recent.back().xmin >= recent.front().xmin && recent.front().xmin >= next.xmin));

			// 3. Add the next mapping to the front
			recent.push_front(next);	 

			// 4. If there is another mapping in the deque, it is near and better.
			if (recent.size() > 1)
				recent.front().matches = -1;
		}

		// 5. Add the last mapping if it is reasonable
		if (!recent.empty() && recent.back().matches != -1)
			reasonable.push_back(recent.back());

		return reasonable;
	}

  public:
	SweepMap(const SketchIndex &tidx, const params_t &params, Timers *T, Counters *C)
		: tidx(tidx), params(params), T(T), C(C) {
			if (params.tThres < 0.0 || params.tThres > 1.0) {
				cerr << "tThres = " << params.tThres << " outside of [0,1]." << endl;
				exit(1);
			}
		}

	bool is_same_strand(const string &P, const string &S) {
		auto p_fw = buildFMHSketch_onlyfw(P, params.k, params.hFrac);
		auto s_fw = buildFMHSketch_onlyfw(S, params.k, params.hFrac);
		auto s_rc = buildFMHSketch_onlyfw(revCompl(S), params.k, params.hFrac);

		auto p_fw_ = set<kmer_with_pos_t>(p_fw.begin(), p_fw.end());
		auto s_fw_ = set<kmer_with_pos_t>(s_fw.begin(), s_fw.end());
		auto s_rc_ = set<kmer_with_pos_t>(s_rc.begin(), s_rc.end());

		std::set<kmer_with_pos_t> sect_fw, sect_rc;
		std::set_intersection(p_fw_.begin(), p_fw_.end(), s_fw_.begin(), s_fw_.end(), std::inserter(sect_fw, sect_fw.begin()));
		std::set_intersection(p_fw_.begin(), p_fw_.end(), s_rc_.begin(), s_rc_.end(), std::inserter(sect_rc, sect_rc.begin()));
		//cout << "intersection sizes: " << " " << s_fw.size() << " " << s_rc.size() << " " << sect_fw.size() << " " << sect_rc.size() << endl;
		return sect_fw.size() > sect_rc.size();
	}

	void map(const string &pFile, const blmers_t &bLstmers) {
		T->start("mapping");
		T->start("query_reading");
		read_fasta_klib(pFile, [this, &bLstmers](kseq_t *seq) {
			T->stop("query_reading");
			T->start("query_mapping");
			T->start("sketching");
			Sketch p = buildFMHSketch(seq->seq.s, params.k, params.hFrac, bLstmers);
			T->stop("sketching");

			string query_id = seq->name.s;
			pos_t P_sz = (pos_t)seq->seq.l;

			C->inc("read_len", P_sz);
			hist_t p_hist;

			T->start("seeding");
			vector<Seed> seeds = choose_seeds(p, params.max_seeds, &p_hist);
			T->stop("seeding");

			T->start("matching");
			vector<Match> matches = match_seeds(p.size(), seeds);
			T->stop("matching");

			//cout << "matches: " << matches.size() << endl;

//			for(auto m: matches) {
//				//m.print();
//				string P_kmer = string(seq->seq.s).substr(m.P_r-params.k, params.k);
//				string T_kmer = tidx.T.substr(m.T_r-params.k, params.k);
//				cout << P_kmer << " " << T_kmer << " " << revCompl(T_kmer) << endl;
//			}

			T->start("sweep");
			vector<Mapping> mappings = sweep(p_hist, matches, P_sz, seeds.size());
			T->stop("sweep");

			T->start("postproc");
			auto reasonable_mappings = params.overlaps ? mappings : filter_reasonable(mappings, P_sz);
			for (auto &m: reasonable_mappings) {
				const auto &segm = tidx.T[m.segm_id];
				//cout << m.T_l << " " << m.T_r << " " << m.segm_id << " " << segm.seq.size() << endl;
				m.strand = is_same_strand(seq->seq.s, segm.seq.substr(m.T_l, m.T_r-m.T_l)) ? '+' : '-';
				m.map_time = (T->secs("sketching") + T->secs("seeding") + T->secs("matching") + T->secs("sweep")) / (double)reasonable_mappings.size();
				m.print_paf(query_id, segm, (int)matches.size());
				C->inc("similarity", int(10000.0*m.sim));
				C->inc("mappings");
				C->inc("sketched_kmers", m.seeds);
				C->inc("matches", m.matches);
			}
			C->inc("reads");
			if (reasonable_mappings.empty())
				C->inc("unmapped_reads");

			T->stop("postproc");
			T->stop("query_mapping");
			T->start("query_reading");
		});
		T->stop("query_reading");
		T->stop("mapping");
	}

	static void print_params(const params_t &params) {
		cerr << "Params:" << endl;
		cerr << " | reference:             " << params.tFile << endl;
		cerr << " | queries:               " << params.pFile << endl;
		cerr << " | k:                     " << params.k << endl;
		//cerr << " | w:                     " << params.w << endl;
		//cerr << " | blacklist file:        " << (params.bLstFl.size() ? params.bLstFl : "-") << endl;
		cerr << " | hFrac:                 " << params.hFrac << endl;
		cerr << " | max_seeds (S):         " << params.max_seeds << endl;
		cerr << " | max_matches (M):       " << params.max_matches << endl;
		cerr << " | onlybest:              " << params.onlybest << endl;
		cerr << " | tThres:                " << params.tThres << endl;
	}

	void print_report(const Counters &C, const Timers &T) {
		cerr << fixed << setprecision(1);
		// TODO: report average Jaccard similarity
		cerr << "Stats:" << endl;
		cerr << " | Total reads:           " << C.count("reads") << " (~" << 1.0*C.count("read_len") / C.count("reads") << " nb per read)" << endl;
		cerr << " | Sketched read kmers:   " << C.count("sketched_kmers") << " (" << C.frac("sketched_kmers", "reads") << " per read)" << endl;
		cerr << " | Kmer matches:          " << C.count("matches") << " (" << C.frac("matches", "reads") << " per read)" << endl;
		//cerr << " | Seed limit reached:    " << C.count("seeds_limit_reached") << " (" << C.perc("seeds_limit_reached", "reads") << "%)" << endl;
		cerr << " | Matches limit reached: " << C.count("matches_limit_reached") << " (" << C.perc("matches_limit_reached", "reads") << "%)" << endl;
		cerr << " | Unmapped reads:        " << C.count("unmapped_reads") << " (" << C.perc("unmapped_reads", "reads") << "%)" << endl;
		cerr << " | Average similarity:    " << C.frac("similarity", "mappings") / 10000.0 << endl;
		cerr << "Total time [sec]:         "         << setw(5) << right << T.secs("total")             << " (" << setw(4) << right << C.count("reads") / T.secs("total")       << " reads per sec)" << endl;
		cerr << " | Index:                 "         << setw(5) << right << T.secs("indexing")          << " (" << setw(4) << right << T.perc("indexing", "total")              << "\%)" << endl;
		cerr << " |  | loading:                "     << setw(5) << right << T.secs("index_reading")     << " (" << setw(4) << right << T.perc("index_reading", "indexing")      << "\%)" << endl;
		cerr << " |  | sketch:                 "     << setw(5) << right << T.secs("index_sketching")   << " (" << setw(4) << right << T.perc("index_sketching", "indexing")    << "\%)" << endl;
		cerr << " |  | initialize:             "     << setw(5) << right << T.secs("index_initializing")<< " (" << setw(4) << right << T.perc("index_initializing", "indexing") << "\%)" << endl;
		cerr << " | Map:                   "         << setw(5) << right << T.secs("mapping")           << " (" << setw(4) << right << T.perc("mapping", "total")               << "\%, " << setw(5) << right << T.range_ratio("query_mapping") << "x)" << endl;
		cerr << " |  | load queries:           "     << setw(5) << right << T.secs("query_reading")     << " (" << setw(4) << right << T.perc("query_reading", "mapping")       << "\%, " << setw(5) << right << T.range_ratio("query_reading") << "x)" << endl;
		cerr << " |  | sketch reads:           "     << setw(5) << right << T.secs("sketching")         << " (" << setw(4) << right << T.perc("sketching", "mapping")           << "\%, " << setw(5) << right << T.range_ratio("sketching") << "x)" << endl;
		cerr << " |  | seeding:                "     << setw(5) << right << T.secs("seeding")           << " (" << setw(4) << right << T.perc("seeding", "mapping")             << "\%, " << setw(5) << right << T.range_ratio("seeding") << "x)" << endl;
		cerr << " |  |  | collect seed info:       " << setw(5) << right << T.secs("collect_seed_info") << " (" << setw(4) << right << T.perc("collect_seed_info", "seeding")   << "\%, " << setw(5) << right << T.range_ratio("collect_seed_info") << "x)" << endl;
		cerr << " |  |  | sort seeds by #matches:  " << setw(5) << right << T.secs("sort_seeds")        << " (" << setw(4) << right << T.perc("sort_seeds", "seeding")          << "\%, " << setw(5) << right << T.range_ratio("sort_seeds") << "x)" << endl;
		cerr << " |  | matching seeds:         "     << setw(5) << right << T.secs("matching")          << " (" << setw(4) << right << T.perc("matching", "mapping")            << "\%, " << setw(5) << right << T.range_ratio("matching") << "x)" << endl;
		cerr << " |  |  | collect matches:         " << setw(5) << right << T.secs("collect_matches")   << " (" << setw(4) << right << T.perc("collect_matches", "matching")    << "\%, " << setw(5) << right << T.range_ratio("collect_matches") << "x)" << endl;
		cerr << " |  |  | sort matches:            " << setw(5) << right << T.secs("sort_matches")      << " (" << setw(4) << right << T.perc("sort_matches", "matching")       << "\%, " << setw(5) << right << T.range_ratio("sort_matches") << "x)" << endl;
		cerr << " |  | sweep:                  "     << setw(5) << right << T.secs("sweep")             << " (" << setw(4) << right << T.perc("sweep", "mapping")               << "\%, " << setw(5) << right << T.range_ratio("sweep") << "x)" << endl;
		cerr << " |  | post proc:              "     << setw(5) << right << T.secs("postproc")          << " (" << setw(4) << right << T.perc("postproc", "mapping")            << "\%, " << setw(5) << right << T.range_ratio("postproc") << "x)" << endl;
	}
};

#endif