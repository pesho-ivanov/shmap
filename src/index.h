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

class SketchIndex {

public:
	std::vector<RefSegment> T;
	int total_size;  // total size of all segments // TODO: change type to size_t
	const params_t &params;
	ankerl::unordered_dense::map<hash_t, std::vector<Hit>> h2pos;
	Timers *timer;
	Counters *C;

	void print_kmer_hist() {
		std::vector<int> hist(10, 0);
		int kmers = 0, different_kmers = 0, max_occ = 0;
		for (const auto& h2p : h2pos) {
			int occ = h2p.second.size();
			kmers += occ;
			++different_kmers;
			if (occ >= (int)hist.size()-1) {
				hist.back() += occ;
				if (occ > max_occ)
					max_occ = occ;
			} else 
				hist[occ] += occ;
		}

		cerr << std::fixed << std::setprecision(2);
		cerr << "Histogram of " << kmers << " kmers ("
			<< 100.0*double(different_kmers)/double(kmers) << "\% different) covering "
			<< 100.0*double(params.k)*double(kmers)/double(total_size) << "\% of the " << total_size << "nb index" << endl;
		for (size_t i=0; i<hist.size(); ++i)
			if (hist[i] > 0)
				cerr << std::setw(5) << std::right << i << (i<hist.size()-1?" ":"+") << "occ: " << std::setw(9) << std::right << hist[i] << " kmers (" << 100.0*double(hist[i])/double(kmers) << "\%)" << endl;
		cerr << "The most frequent kmer occurs " << max_occ << " times." << endl;
	}

	void apply_blacklist() {
		for (auto hits: h2pos)
			if (hits.second.size() > (size_t)params.max_matches) {
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
		total_size += T_sz;
		populate_h2pos(sketch, T.size()-1);
	}

	SketchIndex(const params_t &params, Timers *timer, Counters *C)
		: total_size(0), params(params), timer(timer), C(C) {}

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

		apply_blacklist();
		print_kmer_hist();
		
		C->inc("T_sz", total_size);
		//C->inc("index_memory_MB", get_current_memory_MB());
	}
};

} // namespace sweepmap

#endif