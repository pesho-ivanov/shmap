CC = g++
CXX_STANDARD = -std=c++2a
CFLAGS = -g -O3 -march=native -lm -lpthread -Wall -Wextra -Wno-unused-parameter  -Wno-unused-result 
LIBS = minimap2/libminimap2.a -lz
DEPFLAGS = -MMD -MP

SRCS = src/sweep.cpp src/IO.h src/Sketch.h
SWEEP_BIN = ./sweep
MINIMAP_BIN = minimap2
MAPQUIK_BIN = mapquik
ESKEMAP_BIN = eskemap
EDLIB_BIN = ~/libs/edlib/build/bin/edlib-aligner

CHR_Y = simulations/genomes/t2thumanChrY.fasta
READS_ESKEMAP = simulations/reads/t2thumanChrY_sr0.0001090909090909091_dr0.0009818181818181818_i0.0009090909090909091_sd7361077429744071834_lmn100_lmx1000000_lavg9000_ls7000_dp10_rm20.fasta

all: sweep

$(SWEEP_BIN): $(SRCS)
	$(CC) $(CXX_STANDARD) $(CFLAGS) $< -o $@ $(LIBS)

eval_abe: sweep
	mkdir -p ../out
#	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 >out/sweep.out

	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -z out/sweep-b.params >out/sweep-b.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z out/sweep-b-a.params >out/sweep-b-a.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -z out/sweep-b-e.params >out/sweep-b-e.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -a extend_equally -z out/sweep-b-e-a.params >out/sweep-b-e-a.out

	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -z out/sweep.params >out/sweep.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -a extend_equally -z out/sweep-a.params >out/sweep-a.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -e consecutive -z out/sweep-e.params >out/sweep-e.out
	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -e consecutive -a extend_equally -z out/sweep-e-a.params >out/sweep-e-a.out

eval_multiple_alignments: sweep
#	time $(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -a fine 		     -z out/sweep.params >out/sweep.out
	time $(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -a fine -x		     -z out/sweep-x.params >out/sweep-x.out
#	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -x -z out/sweep-b-a.params >out/sweep-b-a-x.out
#	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z out/sweep-b-a.params >out/sweep-b-a.out 2>out/sweep-b-a.cerr
#	$(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -o -z out/sweep-b-a.params >out/sweep-b-a-o.out 2>out/sweep-b-a-o.cerr

debug: sweep
	$(SWEEP_BIN) -s $(CHR_Y) -p simulations/reads/1.fasta $(READS_ESKEMAP) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z out/sweep-1.params >out/sweep-1.out 2>out/sweep-1.cerr

run_sweep: sweep
	time $(SWEEP_BIN) -s $(CHR_Y) -p $(READS_ESKEMAP) -k 15 -a fine -x -z out/sweep-Y-x.params >out/sweep-Y-x.out

run_minimap:
	$(MINIMAP_BIN) -x map-hifi -t 1 $(CHR_Y) $(READS_ESKEMAP) >out/minimap-Y.paf

run_mapquik:
	$(MAPQUIK_BIN) $(READS_ESKEMAP) --reference $(CHR_Y) --threads 1 -p out/mapquik-Y

run_eskemap:
	time $(ESKEMAP_BIN) -p $(READS_ESKEMAP) -s $(CHR_Y) -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -N >out/eskemap-Y.out

#run_edlib:
#	time $(EDLIB_BIN) $(READS_ESKEMAP) $(CHR_Y) >out/edlib-Y.out

run: run_sweep run_minimap run_mapquik

clean:
	rm -f sweep