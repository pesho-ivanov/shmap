#ifndef SWEEPMAP_INDEX_HPP
#define SWEEPMAP_INDEX_HPP

#include <iostream>
#include <string>
#include <vector>

#include <ankerl/unordered_dense.h>

#include "sketch.h"
#include "utils.h"
#include "io.h"

namespace sweepmap {

struct RefSegment {
	std::string name;
	int sz;
	RefSegment(const std::string &name, const int sz) : name(name), sz(sz) {}
};

struct Hit {  // TODO: compress all in 32bit
	pos_t r;
	bool strand;
	segm_t segm_id;
	Hit(const Kmer &kmer, pos_t pos, segm_t segm_id)
		: r(kmer.r), strand(kmer.strand),
		 segm_id(segm_id) {}
};

class SketchIndex {

public:
	std::vector<RefSegment> T;
	const params_t &params;
	ankerl::unordered_dense::map<hash_t, std::vector<Hit>> h2pos;
	Timers *timer;
	Counters *C;

	void print_kmer_hist() {
		std::vector<int> hist(10, 0);
		int max_occ = 0;
		for (const auto& h2p : h2pos) {
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

	void apply_blacklist() {
		for (auto hits: h2pos)
			if (hits.second.size() > (size_t)params.max_matches) {
				C->inc("blacklisted_kmers");
				C->inc("blacklisted_hits", hits.second.size());
				h2pos.erase(hits.first);  // TODO: use the iterator instead
			}
	}

	void populate_h2pos(const Sketch& sketch, int segm_id) {
		h2pos.reserve(sketch.size());
		for (size_t tpos = 0; tpos < sketch.size(); ++tpos) {
			const Kmer& kmer = sketch[tpos];
			h2pos[kmer.h].push_back(Hit(kmer, tpos, segm_id));
		}
	}

	void add_segment(const Sketch& sketch, const std::string &name, int T_sz) {
		T.push_back(RefSegment(name, T_sz));
		C->inc("segments");
		C->inc("total_nucls", T_sz);
		populate_h2pos(sketch, T.size()-1);
	}

	void print_stats() {
		cerr << std::fixed << std::setprecision(1);
		cerr << "Index stats:" << endl;
		cerr << " | total nucleotides:     " << C->count("total_nucls") << endl;
		cerr << " | index segments:        " << C->count("segments") << " (~" << 1.0*C->count("total_nucls") / C->count("segments") << " nb per segment)" << endl;
		cerr << " | indexed kmers:         " << C->count("indexed_kmers") << endl;
		cerr << " | indexed hits:          " << C->count("indexed_hits") << " ("
												<< double(params.k)*C->perc("indexed_hits", "total_nucls") << "\% of the index, "
												<< "~" << C->frac("indexed_hits", "indexed_kmers") << " per kmer)" << endl;
		cerr << " | | most frequent kmer:      " << C->count("indexed_highest_freq_kmer") << " times." << endl;
		cerr << " | | blacklisted kmers:       " << C->count("blacklisted_kmers") << " (" << C->perc("blacklisted_kmers", "indexed_kmers") << "\%)" << endl;
		cerr << " | | blacklisted hits:        " << C->count("blacklisted_hits") << " (" << C->perc("blacklisted_hits", "indexed_hits") << "\%)" << endl;
//		cerr << std::fixed << std::setprecision(2);
//		cerr << "Histogram of " << C->count("indexed_total_kmers") << " kmers ("
//			<< 100.0*double(C->count("indexed_different_kmers"))/double(C->count("indexed_total_kmers")) << "\% different) covering "
//			<< 100.0*double(params.k)*double(C->count("indexed_total_kmers"))/double(C->count("total_nucls")) << "\% of the " << C->count("total_nucls") << "nb index" << endl;
//		for (size_t i=0; i<hist.size(); ++i)
//			if (hist[i] > 0)
//				cerr << std::setw(5) << std::right << i << (i<hist.size()-1?" ":"+") << "occ: " << std::setw(9) << std::right << hist[i] << " kmers (" << 100.0*double(hist[i])/double(kmers) << "\%)" << endl;
	}

	SketchIndex(const params_t &params, Timers *timer, Counters *C)
		: params(params), timer(timer), C(C) {}

	void index(const std::string &tFile) {
		timer->start("indexing");
		cerr << "Indexing " << params.tFile << "..." << endl;
		timer->start("index_reading");
		read_fasta_klib(params.tFile, [this](kseq_t *seq) {
			timer->stop("index_reading");
			timer->start("index_sketching");
			Sketch t = buildFMHSketch(seq->seq.s, params.k, params.hFrac);
			timer->stop("index_sketching");

			timer->start("index_initializing");
			add_segment(t, seq->name.s, seq->seq.l);
			timer->stop("index_initializing");

			timer->start("index_reading");
		});
		timer->stop("index_reading");
		timer->stop("indexing");

		print_kmer_hist();
		apply_blacklist();

		//C->inc("index_memory_MB", get_current_memory_MB());

		print_stats();
	}
};

} // namespace sweepmap

#endif