CC = g++
CXX_STANDARD = -std=c++2a
CFLAGS = -g -O3 -march=native -lm -lpthread -Wno-unused-result
LIBS = minimap2/libminimap2.a -lz
DEPFLAGS = -MMD -MP

SRCS = src/sweep.cpp src/IO.h src/Index.h src/Measures.h src/Sketch.h
SWEEP_BIN = ./sweep

CHR_Y = simulations/genomes/t2thumanChrY.fasta
READS_ESKEMAP = simulations/reads/t2thumanChrY_sr0.0001090909090909091_dr0.0009818181818181818_i0.0009090909090909091_sd7361077429744071834_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20.fasta

all: sweep

$(SWEEP_BIN): $(SRCS)
	$(CC) $(CXX_STANDARD) $(CFLAGS) $< -o $@ $(LIBS)

run: sweep
	mkdir -p ../out
#	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 >out/sweep.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -z out/sweep-b.params >out/sweep-b.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z out/sweep-b-a.params >out/sweep-b-a.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -z out/sweep-b-e.params >out/sweep-b-e.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -a extend_equally -z out/sweep-b-e-a.params >out/sweep-b-e-a.out

clean:
	rm -f sweep