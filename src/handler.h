#pragma once

#include <iostream>
using namespace std;

#include "io.h"
#include "utils.h"
#include "sketch.h"

namespace sweepmap {

class Handler {
public:
	Counters C;
	Timers T;
	params_t params;
    FracMinHash sketcher;

    Handler(const params_t &params)
            : params(params),
              sketcher(params.k, params.hFrac, &C, &T) {
        T.start("total");
        if (!params.paramsFile.empty()) {
            cerr << "Writing parameters to " << params.paramsFile << "..." << endl;
            auto fout = std::ofstream(params.paramsFile);
            params.print(fout, false);
        } else {
            params.print(cerr, true);
        }
        this->params.print_display(std::cerr);
    }

    sketch_t sketch(const std::string& s) {
        return sketcher.sketch(s);
    }

	void print_sketching_stats() {
        T.stop("total");

		cerr << std::fixed << std::setprecision(1);
		cerr << "Sketching:" << endl;
		cerr << " | Sketched sequences:    " << C.count("sketched_seqs") << " (" << C.count("sketched_len") << " nb)" << endl;
		cerr << " | Kmers:                 " << C.count("sketched_kmers") << endl;
		//cerr << " |  | original:               " << C.count("original_kmers") << " (" << C.perc("original_kmers", "sketched_kmers") << "%)" << endl;
	}

    void printMemoryUsage() {
        std::ifstream statusFile("/proc/self/status");
        std::string line;

        if (statusFile.is_open()) {
            while (getline(statusFile, line)) {
                // Check for memory usage lines
                if (line.find("VmSize:") != std::string::npos || line.find("VmRSS:") != std::string::npos) {
                    std::cerr << line << std::endl;
                }
            }
            statusFile.close();
        } else {
            std::cerr << "Unable to open /proc/self/status" << std::endl;
        }
    }
};

} // namespace sweepmap