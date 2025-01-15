#include "io.h"

namespace sweepmap {


void read_fasta_klib(const std::string& filename, std::function<void(const std::string&, const std::string&, float)> callback) {
    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int l;

    std::ifstream in_for_size(filename, std::ifstream::ate | std::ifstream::binary);
    auto total_bytes = in_for_size.tellg(); 

    while ((l = kseq_read(seq)) >= 0) {
        string query_id = seq->name.s;
        string P = seq->seq.s;
        //print_progress_bar(cerr, 1.0 * gztell(fp) / total_bytes);
        //std::cerr << "Mapping progress:" << 100.0 * gztell(fp) / total_bytes << "%\r";
        float percentage = 1.0 * gztell(fp) / total_bytes;
        callback(query_id, P, percentage);
    }

    //if (l != -1)  // -1  end of file; do nothing
    if (l == -2) cerr << "ERROR: truncated quality string" << endl;  // -2   truncated quality string
    else if (l == -3) cerr << "ERROR: error reading stream" << endl;  // -3   error reading stream

    kseq_destroy(seq);
    gzclose(fp);
}

}