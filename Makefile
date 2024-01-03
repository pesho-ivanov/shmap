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

DIR = newevals
REFNAME = chm13
ACCURACY = 0.99
DEPTH = 0.01
MEANLEN = 10000  # hifi
REF = $(DIR)/$(REFNAME).fa
READS_PREFIX = reads-$(REFNAME)
READS = $(DIR)/$(READS_PREFIX).fa
OUTDIR = $(DIR)/$(REFNAME)-a$(ACCURACY)-d$(DEPTH)

MAX_SEEDS = 10000
MAX_MATCHES = 100 300 1000 3000 10000 30000 100000
#10000000

all: sweep

$(SWEEP_BIN): $(SRCS)
	$(CC) $(CXX_STANDARD) $(CFLAGS) $< -o $@ $(LIBS)

gen_reads:
	pbsim \
		   $(REF) \
		   --model_qc $(DIR)/model_qc_clr \
		   --accuracy-mean $(ACCURACY)\
		   --accuracy-sd 0\
		   --depth $(DEPTH)\
		   --prefix $(READS_PREFIX)\
		   --length-mean $(MEANLEN)

	samtools faidx $(REF)
	paftools.js pbsim2fq $(REF).fai "$(READS_PREFIX)"_*.maf >$(READS)_
	rm -f "$(READS_PREFIX)"_*.maf "$(READS_PREFIX)"_*.ref "$(READS_PREFIX)"_*.fastq
	awk '/^>/ {printit = /+$$/} printit' $(READS)_ >$(READS)

eval_abe: sweep
	mkdir -p ../out
#	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 >out/sweep.out

	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -z $(DIR)/out/sweep-b.params >out/sweep-b.out
	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a.out
	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -z $(DIR)/out/sweep-b-e.params >out/sweep-b-e.out
	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -a extend_equally -z $(DIR)/out/sweep-b-e-a.params >out/sweep-b-e-a.out

	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -z $(DIR)/out/sweep.params >out/sweep.out
	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -a extend_equally -z $(DIR)/out/sweep-a.params >out/sweep-a.out
	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -e consecutive -z $(DIR)/out/sweep-e.params >out/sweep-e.out
	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -e consecutive -a extend_equally -z $(DIR)/out/sweep-e-a.params >out/sweep-e-a.out

eval_multiple_alignments: sweep
#	time $(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -a fine 		     -z $(DIR)/out/sweep.params >out/sweep.out
	time $(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -a fine -x		     -z $(DIR)/out/sweep-x.params >out/sweep-x.out
#	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -x -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a-x.out
#	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a.out 2>out/sweep-b-a.cerr
#	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -o -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a-o.out 2>out/sweep-b-a-o.cerr

debug: sweep
	$(SWEEP_BIN) -s $(REF) -p simulations/reads/1.fasta $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-1.params >out/sweep-1.out 2>out/sweep-1.cerr

eval_sweep_max_seeds: sweep
	@DIR=newevals/eval_sweep-K22; \
	mkdir -p $${DIR}; \
	for maxseeds in $(MAX_SEEDS); do \
		for maxmatches in $(MAX_MATCHES); do \
			f=$${DIR}/"sweep-S$${maxseeds}-M$${maxmatches}"; \
			echo "Processing $${f}"; \
			time $(SWEEP_BIN) -s $(REF) -p $(READS) -t 0.0 -x -k 22 -z sweep.params -S $${maxseeds} -M $${maxmatches} $${f}.params >$${f}.paf 2>$${f}.log; \
			paftools.js mapeval $${f}.paf | tee $${f}.eval; \
		done \
    done

run_sweep: sweep
	mkdir -p $(OUTDIR)
	time -o $(OUTDIR)/sweep.time $(SWEEP_BIN) -s $(REF) -p $(READS) -t 0.0 -x -k 20 -z $(OUTDIR)/sweep.params -S 300 -M 10000 $(OUTDIR)/sweep.params >$(OUTDIR)/sweep.paf 2>$(OUTDIR)/sweep.log
	paftools.js mapeval $(OUTDIR)/sweep.paf | tee $(OUTDIR)/sweep.eval

run_minimap:
	mkdir -p $(OUTDIR)
	time -o $(OUTDIR)/minimap.time $(MINIMAP_BIN) -x map-hifi -t 1 --secondary=no $(REF) $(READS) >$(OUTDIR)/minimap.paf 2>$(OUTDIR)/minimap.log
	paftools.js mapeval $(OUTDIR)/minimap.paf | tee $(OUTDIR)/minimap.eval

run_mapquik:
	mkdir -p $(OUTDIR)
	time -o $(OUTDIR)/mapquik.time $(MAPQUIK_BIN) $(READS) --reference $(REF) --threads 1 -p $(OUTDIR)/mapquik 2>$(OUTDIR)/mapquik.log
	paftools.js mapeval $(OUTDIR)/mapquik.paf | tee $(OUTDIR)/mapquik.eval

run_eskemap:
	mkdir -p $(OUTDIR)
	time -o $(OUTDIR)/eskemap.time $(ESKEMAP_BIN) -p $(READS) -s $(REF) -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -N >$(OUTDIR)/eskemap.out 2>$(OUTDIR)/eskemap.log

run: run_sweep run_minimap run_mapquik

clean:
	rm -f sweep
