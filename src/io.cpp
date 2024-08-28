#include "io.h"

namespace sweepmap {

void read_fasta_klib(const std::string& filename, std::function<void(kseq_t*)> callback) {
    gzFile fp = gzopen(filename.c_str(), "r");
    kseq_t *seq = kseq_init(fp);
    int l;

    while ((l = kseq_read(seq)) >= 0) {
        callback(seq);
    }

    kseq_destroy(seq);
    gzclose(fp);
}

}