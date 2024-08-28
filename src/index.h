#pragma once

#include <string>
#include <vector>

#include <ankerl/unordered_dense.h>

#include "sketch.h"
#include "utils.h"
#include "io.h"

namespace sweepmap {

class SketchIndex {

public:
	std::vector<RefSegment> T;
	ankerl::unordered_dense::map<hash_t, Hit> h2single;               // all sketched kmers with =1 hit
	ankerl::unordered_dense::map<hash_t, std::vector<Hit>> h2multi;   // all sketched kmers with >1 hits

	const params_t &params;
	Counters *C;
	Timers *timer;
	const FracMinHash &sketcher;

	void get_kmer_stats() {
		std::vector<int> hist(10, 0);
		int max_occ = 0;
        C->inc("indexed_hits", h2single.size());
        C->inc("indexed_kmers", h2single.size());
		for (const auto& h2p : h2multi) {
			int occ = h2p.second.size();
			C->inc("indexed_hits", occ);
			C->inc("indexed_kmers");
			if (occ >= (int)hist.size()-1) {
				hist.back() += occ;
				if (occ > max_occ)
					max_occ = occ;
			} else 
				hist[occ] += occ;
		}
		C->inc("indexed_highest_freq_kmer", max_occ);
	}

	int count(hash_t h) const {
		if (h2single.contains(h)) return 1;
		else if (h2multi.contains(h)) return h2multi.at(h).size();
		else return 0;
	}

	void add_matches(std::vector<Match> *matches, const Seed &s, int seed_num) const {
		if (s.hits_in_T == 1) {
			assert(h2single.contains(s.kmer.h));
			matches->push_back(Match(s, h2single.at(s.kmer.h), seed_num));
				
		} else {
			assert(s.hits_in_T > 1);
			assert(h2multi.contains(s.kmer.h));
			for (const auto &hit: h2multi.at(s.kmer.h))
				matches->push_back(Match(s, hit, seed_num));
		}	
	}

	void erase_frequent_kmers() {
		std::vector<hash_t> blacklisted_h;
		for (const auto &[h, hits]: h2multi)
			if (hits.size() > (size_t)params.max_matches) {
				blacklisted_h.push_back(h);
				C->inc("blacklisted_kmers");
				C->inc("blacklisted_hits", hits.size());
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
            else if (h2multi[kmer.h].size() < (size_t)params.max_matches + 1)
                h2multi[kmer.h].push_back(hit);
		}
	}

	void add_segment(const kseq_t *seq, const sketch_t& sketch) {
		string segm_name = seq->name.s;
		string segm_seq = seq->seq.s;
		int segm_size = seq->seq.l;
		T.push_back(RefSegment(sketch, segm_name, segm_seq, segm_size));
		C->inc("segments");
		C->inc("total_nucls", segm_size);
		populate_h2pos(sketch, T.size()-1);
	}

	SketchIndex(const params_t &params, Counters *C, Timers *timer, FracMinHash &sketcher)
		: params(params), C(C), timer(timer), sketcher(sketcher) {}

	void build_index(const std::string &tFile) {
		timer->start("indexing");
		cerr << "Indexing " << params.tFile << "..." << endl;
		timer->start("index_reading");
		read_fasta_klib(params.tFile, [this](kseq_t *seq) {
			timer->stop("index_reading");
			timer->start("index_sketching");
			sketch_t t = sketcher.sketch(seq->seq.s);
			timer->stop("index_sketching");

			timer->start("index_initializing");
			add_segment(seq, t);
			timer->stop("index_initializing");

			timer->start("index_reading");
		});

		// if a kmer is present in both single and multi, we move it out of single to multi
		for (auto &[h, hits] : h2multi) {
			if (h2single.contains(h)) {
				hits.push_back(h2single.at(h));
				h2single.erase(h);
			}
		}
		timer->stop("index_reading");
		timer->stop("indexing");

		get_kmer_stats();
        C->inc("blacklisted_kmers", 0);
        C->inc("blacklisted_hits", 0);
		erase_frequent_kmers();
		print_stats();
	}

	void print_stats() {
		cerr << std::fixed << std::setprecision(1);
		cerr << "Index stats:" << endl;
        printMemoryUsage();
		cerr << " | total nucleotides:     " << C->count("total_nucls") << endl;
		cerr << " | index segments:        " << C->count("segments") << " (~" << 1.0*C->count("total_nucls") / C->count("segments") << " nb per segment)" << endl;
		cerr << " | indexed kmers:         " << C->count("indexed_kmers") << endl;
		cerr << " | indexed hits:          " << C->count("indexed_hits") << " ("
												<< double(params.k)*C->perc("indexed_hits", "total_nucls") << "\% of the index, "
												<< "~" << C->frac("indexed_hits", "indexed_kmers") << " per kmer)" << endl;
		cerr << " | | most frequent kmer:      " << C->count("indexed_highest_freq_kmer") << " times." << endl;
		cerr << " | | blacklisted kmers:       " << C->count("blacklisted_kmers") << " (" << C->perc("blacklisted_kmers", "indexed_kmers") << "\%)" << endl;
		cerr << " | | blacklisted hits:        " << C->count("blacklisted_hits") << " (" << C->perc("blacklisted_hits", "indexed_hits") << "\%)" << endl;
	}
};

} // namespace sweepmap