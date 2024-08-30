#pragma once

#include <string>
#include <vector>

#include <ankerl/unordered_dense.h>

#include "sketch.h"
#include "utils.h"
#include "io.h"
#include "handler.h"

namespace sweepmap {

class SketchIndex {

public:
	std::vector<RefSegment> T;
	ankerl::unordered_dense::map<hash_t, Hit> h2single;               // all sketched kmers with =1 hit
	ankerl::unordered_dense::map<hash_t, std::vector<Hit>> h2multi;   // all sketched kmers with >1 hits

	Handler *H;

	void get_kmer_stats() const {
		std::vector<int> hist(10, 0);
		int max_occ = 0;
        H->C.inc("indexed_hits", h2single.size());
        H->C.inc("indexed_kmers", h2single.size());
		for (const auto& h2p : h2multi) {
			int occ = h2p.second.size();
			H->C.inc("indexed_hits", occ);
			H->C.inc("indexed_kmers");
			if (occ >= (int)hist.size()-1) {
				hist.back() += occ;
				if (occ > max_occ)
					max_occ = occ;
			} else 
				hist[occ] += occ;
		}
		H->C.inc("indexed_highest_freq_kmer", max_occ);
	}

	int count(hash_t h) const {
		if (h2single.contains(h)) return 1;
		else if (h2multi.contains(h)) return h2multi.at(h).size();
		else return 0;
	}

	bool intersect(const Hit &hit, const int from, const int to) const {
		// TODO: account for segments
		return hit.r >= from && hit.r < to;
	}

	std::vector<Match> get_matches_in_interval(const Seed &s, const int from, const int to) const {
		// assume matches of each seed are sorted by position in T
		// TODO: account for segments
		std::vector<Match> matches;
		if (s.hits_in_T == 1) {
			auto &hit = h2single.at(s.kmer.h);
			if (intersect(hit, from, to))
				matches.push_back(Match(s, hit, -1));
		} else {
			const vector<Hit> &hits = h2multi.at(s.kmer.h);
			// TODO: account for segments
			// TODO: careful with left and right ends
			auto it = lower_bound(hits.begin(), hits.end(), from, [](const Hit &hit, int pos) { return hit.r < pos; });
			for (; it != hits.end() && it->r < to; ++it)
				matches.push_back(Match(s, *it, -1));
		}	
		return matches;
	}

	void get_matches(std::vector<Match> *matches, const Seed &s, int seed_num) const {
		matches->reserve(s.hits_in_T);
		if (s.hits_in_T == 1) {
			matches->push_back(Match(s, h2single.at(s.kmer.h), seed_num));
		} else {
			for (const auto &hit: h2multi.at(s.kmer.h))
				matches->push_back(Match(s, hit, seed_num));
		}	
	}

	void erase_frequent_kmers() {
		assert(H->params.max_matches != -1);
		std::vector<hash_t> blacklisted_h;
		for (const auto &[h, hits]: h2multi)
			if ((int)hits.size() > H->params.max_matches) {
				blacklisted_h.push_back(h);
				H->C.inc("blacklisted_kmers");
				H->C.inc("blacklisted_hits", hits.size());
			}

		for (auto h: blacklisted_h)
			h2multi.erase(h);
	}

	void populate_h2pos(const sketch_t& sketch, int segm_id) {
		// TODO: skip creating the sketch structure
		for (size_t tpos = 0; tpos < sketch.size(); ++tpos) {
			const Kmer& kmer = sketch[tpos];
			const auto hit = Hit(kmer, tpos, segm_id);
			if (!h2single.contains(kmer.h))
				h2single[kmer.h] = hit; 
            else if (H->params.max_matches == -1 || (int)h2multi[kmer.h].size() < H->params.max_matches + 1)
                h2multi[kmer.h].push_back(hit);
		}
	}

	void add_segment(const kseq_t *seq, const sketch_t& sketch) {
		string segm_name = seq->name.s;
		string segm_seq = seq->seq.s;
		int segm_size = seq->seq.l;
		T.push_back(RefSegment(sketch, segm_name, segm_seq, segm_size));
		H->C.inc("segments");
		H->C.inc("total_nucls", segm_size);
		populate_h2pos(sketch, T.size()-1);
	}

	SketchIndex(Handler *H) : H(H) {}

	void build_index(const std::string &tFile) {
		H->T.start("indexing");
		cerr << "Indexing " << H->params.tFile << "..." << endl;
		H->T.start("index_reading");
		read_fasta_klib(H->params.tFile, [this](kseq_t *seq) {
			H->T.stop("index_reading");
			H->T.start("index_sketching");
			sketch_t t = H->sketcher.sketch(seq->seq.s);
			H->T.stop("index_sketching");

			H->T.start("index_initializing");
			add_segment(seq, t);
			H->T.stop("index_initializing");

			H->T.start("index_reading");
		});

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
				return a.r < b.r;
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
		cerr << "Index stats:" << endl;
        //printMemoryUsage();
		cerr << " | total nucleotides:     " << H->C.count("total_nucls") << endl;
		cerr << " | index segments:        " << H->C.count("segments") << " (~" << 1.0*H->C.count("total_nucls") / H->C.count("segments") << " nb per segment)" << endl;
		cerr << " | indexed kmers:         " << H->C.count("indexed_kmers") << endl;
		cerr << " | indexed hits:          " << H->C.count("indexed_hits") << " ("
												<< double(H->params.k)*H->C.perc("indexed_hits", "total_nucls") << "\% of the index, "
												<< "~" << H->C.frac("indexed_hits", "indexed_kmers") << " per kmer)" << endl;
		cerr << " | | most frequent kmer:      " << H->C.count("indexed_highest_freq_kmer") << " times." << endl;
		cerr << " | | blacklisted kmers:       " << H->C.count("blacklisted_kmers") << " (" << H->C.perc("blacklisted_kmers", "indexed_kmers") << "\%)" << endl;
		cerr << " | | blacklisted hits:        " << H->C.count("blacklisted_hits") << " (" << H->C.perc("blacklisted_hits", "indexed_hits") << "\%)" << endl;
	}
};

} // namespace sweepmap