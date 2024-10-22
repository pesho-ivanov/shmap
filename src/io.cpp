#include "io.h"

namespace sweepmap {

void read_fasta_klib(const std::string& filename, std::function<void(kseq_t*)> callback) {
    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int l;

    while ((l = kseq_read(seq)) >= 0) {
        callback(seq);
    }

    //if (l != -1)  // -1  end of file; do nothing
    if (l == -2) cerr << "ERROR: truncated quality string" << endl;  // -2   truncated quality string
    else if (l == -3) cerr << "ERROR: error reading stream" << endl;  // -3   error reading stream

    kseq_destroy(seq);
    gzclose(fp);
}

}