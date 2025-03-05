#include "klib/kseq.h"
#include <zlib.h>
#include <string>
#include <iostream>

KSEQ_INIT(gzFile, gzread)

int main(int argc, char *argv[]) {
    if (argc != 2) {
        std::cerr << "Usage: " << argv[0] << " <input.fa>" << std::endl;
        return 1;
    }

    // Open input file
    gzFile fp = gzopen(argv[1], "r");
    if (!fp) {
        std::cerr << "Failed to open " << argv[1] << std::endl;
        return 1;
    }

    // Initialize parser
    kseq_t *seq = kseq_init(fp);

    // Print TSV header
    std::cout << "query_id\tsequence" << std::endl;

    // Read sequences
    while (kseq_read(seq) >= 0) {
        // Print sequence name (without '>')
        std::cout << seq->name.s << "\t";
        
        // Print sequence (already comes without newlines from kseq)
        std::cout << seq->seq.s << std::endl;
    }

    // Clean up
    kseq_destroy(seq);
    gzclose(fp);

    return 0;
}
