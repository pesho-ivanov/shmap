#pragma once

#include <fstream>
#include <functional>
#include <iomanip>
#include <string>
#include <vector>
#include <iostream>
#include <sstream>
#include <cassert>

#include <zlib.h>
#include "klib/kseq.h"
#include "utils.h"
#include "cmd_line_parser/include/parser.hpp"  // use Giulio's custom command line parser
#include "../ext/tracy/public/tracy/Tracy.hpp"

namespace sweepmap {

KSEQ_INIT(gzFile, gzread)

using std::cerr;
using std::string;
using std::pair;
using std::ifstream;
using std::endl;

// The options supported by sweepmap are now described by params_t.
struct params_t {
    // Required file names.
    string pFile, tFile;

    // Options with an argument:
    int k;             // K-mer length
    double hFrac;      // FracMinHash ratio
    int max_seeds;     // Maximum seeds in a sketch
    int max_matches;   // Maximum seed matches in a sketch
    double theta;      // Homology threshold
    double min_diff;   // Minimum difference between the best and second best mapping
    double max_overlap;// Maximum overlap between the best and second best mapping
    string paramsFile; // File to output parameters (TSV)
    Metric metric;     // Optimization metric (see utils.h)
    int verbose;       // Verbosity level

    // Boolean options (default false)
    bool sam;              // Output in SAM format (default: PAF)
    bool normalize;        // Normalize scores by length
    bool no_bucket_pruning;// Disables bucket pruning
    bool one_sweep;        // Run one sweep on all matches
    bool abs_pos;          // Use absolute positions instead of kmer positions

    params_t() :
        k(15), hFrac(0.05), max_seeds(-1), max_matches(-1),
        theta(0.9), min_diff(0.02), max_overlap(0.5),
        metric(Metric::Containment), verbose(0),
        sam(false), normalize(false),
        no_bucket_pruning(false), one_sweep(false), abs_pos(false)
    {}

    // Prints the parameters in either human-readable or TSV format.
    void print(std::ostream& out, bool human) const {
        std::vector<pair<string, string>> m;
        m.push_back({"pFile", pFile});
        m.push_back({"tFile", tFile});
        m.push_back({"k", std::to_string(k)});
        m.push_back({"hFrac", std::to_string(hFrac)});
        m.push_back({"max_seeds", std::to_string(max_seeds)});
        m.push_back({"max_matches", std::to_string(max_matches)});
        m.push_back({"tThres", std::to_string(theta)});
        m.push_back({"min_diff", std::to_string(min_diff)});
        m.push_back({"max_overlap", std::to_string(max_overlap)});
        m.push_back({"paramsFile", paramsFile});
        m.push_back({"metric", mapping_metric_str(metric)});

        m.push_back({"sam", std::to_string(sam)});
        m.push_back({"normalize", std::to_string(normalize)});
        m.push_back({"verbose", std::to_string(verbose)});

        m.push_back({"no-bucket-pruning", std::to_string(no_bucket_pruning)});
        m.push_back({"one-sweep", std::to_string(one_sweep)});
        m.push_back({"abs-pos", std::to_string(abs_pos)});

        if (human) {
            out << "Parameters:" << endl;
            for (auto& p : m)
                out << std::setw(20) << std::right << p.first << ": " << p.second << endl;
        } else {
            for (auto& p : m)
                out << p.first << "\t";
            out << endl;
            for (auto& p : m)
                out << p.second << "\t";
        }
    }

    // A short display of key parameters.
    void print_display(std::ostream& out) {
        out << "Params:" << endl;
        out << " | reference:             " << tFile << endl;
        out << " | reads:                 " << pFile << endl;
        out << " | metric:                " << mapping_metric_str(metric) << endl;
        out << " | k:                     " << k << endl;
        out << " | hFrac:                 " << hFrac << endl;
        out << " | verbose:               " << verbose << endl;
        out << " | tThres:                " << theta << endl;
        out << " | min_diff:              " << min_diff << endl;
        out << " | max_overlap:           " << max_overlap << endl;
        out << " | no-bucket-pruning:     " << no_bucket_pruning << endl;
        out << " | abs-pos:               " << abs_pos << endl;
    }

    // This method replaces the getopt-based parsing with the parser.hpp based parser.
    bool prsArgs(int argc, char** argv) {
        using namespace cmd_line_parser;
        parser p(argc, argv);

        // Required options.
        // name, description, shorthand, required argument, boolean option (default is false)
        p.add("pattern",    "Pattern sequences file (FASTA format)",                            "-p", true, false);
        p.add("text",       "Text sequence file (FASTA format)",                                "-s", true, false);

        // Options with an argument.
        p.add("ksize",      "K-mer length to be used for sketches (positive integer)",          "-k", false, false);
        p.add("hashratio",  "FracMinHash ratio in (0,1]",                                       "-r", false, false);
        p.add("max_seeds",  "Maximum seeds in a sketch (positive integer)",                     "-S", false, false);
        p.add("max_matches","Maximum seed matches in a sketch (positive integer)",              "-M", false, false);
        p.add("threshold",  "Homology threshold in [0,1]",                                      "-t", false, false);
        p.add("min_diff",   "Minimum difference between best and second best mapping",          "-d", false, false);
        p.add("max_overlap","Maximum overlap between best and second best mapping (0,1]",       "-o", false, false);
        p.add("metric",     "Optimization metric: bucket_SH, bucket_LCS, Containment, Jaccard", "-m", false, false);
        p.add("params",     "Output file with parameters (tsv)",                                "-z", false, false);
        p.add("verbose",    "Verbosity level: 0 for none, 1 for some, 2 for all",               "-v", false, false);

        // Boolean (flag) options.
        p.add("sam",        "Output in SAM format (PAF by default)",                            "-a", false, true);
        p.add("normalize",  "Normalize scores by length",                                       "-n", false, true);
        p.add("no_bucket_pruning", "Disables bucket pruning",                                   "-P", false, true);
        p.add("one_sweep",  "Disregards seed heuristic and runs one sweep on all matches",      "-B", false, true);
        p.add("abs_pos",    "Use absolute positions instead of kmer positions",                 "-F", false, true);

        if (!p.parse())
            return false;

        try {
            // Retrieve required parameters.
            pFile = p.get<std::string>("pattern");
            tFile = p.get<std::string>("text");

            // Retrieve optional options if provided (otherwise the defaults remain).
            if (p.parsed("ksize")) {
                k = p.get<int>("ksize");
                if (k <= 0) {
                    std::cerr << "ERROR: K-mer length (-k) should be positive. You provided " << k << "." << std::endl;
                    return false;
                }
            }
            if (p.parsed("hashratio")) {
                hFrac = p.get<double>("hashratio");
                if (hFrac <= 0 || hFrac > 1.0) {
                    std::cerr << "ERROR: Given hash ratio (-r) should be between 0 and 1. You provided " << hFrac << "." << std::endl;
                    return false;
                }
            }
            if (p.parsed("max_seeds")) {
                max_seeds = p.get<int>("max_seeds");
                if (max_seeds <= 0) {
                    std::cerr << "ERROR: The number of seeds (-S) should be positive. You provided " << max_seeds << "." << std::endl;
                    return false;
                }
            }
            if (p.parsed("max_matches")) {
                max_matches = p.get<int>("max_matches");
                if (max_matches <= 0) {
                    std::cerr << "ERROR: The number of seed matches (-M) should be positive. You provided " << max_matches << "." << std::endl;
                    return false;
                }
            }
            if (p.parsed("threshold")) {
                theta = p.get<double>("threshold");
                if (theta < 0.0 || theta > 1.0) {
                    std::cerr << "ERROR: The threshold (-t) should be between 0 and 1. You provided " << theta << "." << std::endl;
                    return false;
                }
            }
            if (p.parsed("min_diff")) {
                min_diff = p.get<double>("min_diff");
                if (min_diff < 0.0) {
                    std::cerr << "ERROR: The minimum difference (-d) should be >= 0. You provided " << min_diff << "." << std::endl;
                    return false;
                }
            }
            if (p.parsed("max_overlap")) {
                max_overlap = p.get<double>("max_overlap");
                if (max_overlap < 0.0 || max_overlap > 1.0) {
                    std::cerr << "ERROR: The maximum overlap (-o) should be between 0 and 1. You provided " << max_overlap << "." << std::endl;
                    return false;
                }
            }
            if (p.parsed("verbose")) {
                verbose = p.get<int>("verbose");
                if (verbose < 0 || verbose > 2) {
                    std::cerr << "ERROR: --verbose (-v) should be 0, 1, or 2. You provided " << verbose << "." << std::endl;
                    return false;
                }
            }
            if (p.parsed("params")) {
                paramsFile = p.get<std::string>("params");
            }
            if (p.parsed("metric")) {
                std::string metric_str = p.get<std::string>("metric");
                metric = mapping_metric_from_str(metric_str);
            }

            // For boolean options, simply retrieve the flags.
            sam              = p.get<bool>("sam");
            normalize        = p.get<bool>("normalize");
            no_bucket_pruning= p.get<bool>("no_bucket_pruning");
            one_sweep        = p.get<bool>("one_sweep");
            abs_pos          = p.get<bool>("abs_pos");
        } catch (std::runtime_error const& e) {
            std::cerr << "ERROR: " << e.what() << std::endl;
            return false;
        }
        return true;
    }
};

// This function remains unchanged: it reads a FASTA file using kseq.
inline void read_fasta_klib(const std::string& filename,
    std::function<void(const std::string&, const std::string&, float)> callback)
{
    ZoneScoped;
    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int l;

    std::ifstream in_for_size(filename, std::ifstream::ate | std::ifstream::binary);
    auto total_bytes = in_for_size.tellg(); 

    while ((l = kseq_read(seq)) >= 0) {
        string query_id = seq->name.s;
        string P = seq->seq.s;
        float percentage = 1.0 * gztell(fp) / total_bytes;
        callback(query_id, P, percentage);
    }
    if (l == -2)
        cerr << "ERROR: truncated quality string" << endl;
    else if (l == -3)
        cerr << "ERROR: error reading stream" << endl;

    kseq_destroy(seq);
    gzclose(fp);
}

} // namespace sweepmap