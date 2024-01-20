SHELL := /bin/bash
CC = g++
CXX_STANDARD = -std=c++20
CFLAGS = -O2 -g -march=native -lm -lpthread -Igtl/ -Wall -Wextra -Wno-unused-parameter -Wno-unused-result -Wno-comment -fpermissive #-Wconversion 
LIBS = minimap2/libminimap2.a -lz
DEPFLAGS = -MMD -MP

TIME_CMD = /usr/bin/time -f "%U\t%M"

SRCS = src/sweepmap.cpp src/sweepmap.h src/io.h src/sketch.h src/utils.h src/kseq.h
SWEEP_BIN = sweepmap
MINIMAP_BIN = minimap2
BLEND_BIN = ~/libs/blend/bin/blend
VERITYMAP_BIN = ~/libs/VerityMap/veritymap/main.py
MAPQUIK_BIN = mapquik
ESKEMAP_BIN = eskemap
MERYL_BIN = ~/libs/Winnowmap/bin/meryl
WINNOWMAP_BIN = ~/libs/Winnowmap/bin/winnowmap
EDLIB_BIN = ~/libs/edlib/build/bin/edlib-aligner

REFNAME ?= chm13-1B
ACCURACY ?= 0.99
DEPTH ?= 1
MEANLEN ?= 10000

K ?= 22
R ?= 0.1 
#W ?= 20
S ?= 200
M ?= 2000
T ?= 0.0

K_SLOW ?= 22
R_SLOW ?= 0.1
S_SLOW ?= 1000
M_SLOW ?= 10000
T_SLOW ?= 0.0

DIR = evals
REF = $(DIR)/refs/$(REFNAME).fa
READS_PARAMS = $(REFNAME)-a$(ACCURACY)-d$(DEPTH)-l$(MEANLEN)
READS_PREFIX = reads-$(READS_PARAMS)
ALLREADS_DIR = $(DIR)/reads
READS = $(ALLREADS_DIR)/$(READS_PREFIX).fa
ALLOUT_DIR = $(DIR)/out
OUTDIR = $(ALLOUT_DIR)/$(READS_PARAMS)

MAX_SEEDS = 10000
MAX_MATCHES = 100 300 1000 3000 10000 30000 100000

all: sweepmap

$(SWEEP_BIN): $(SRCS)
	$(CC) $(CXX_STANDARD) $(CFLAGS) $< -o $@ $(LIBS)

gen_reads:
ifeq ($(wildcard $(READS)),)
	mkdir -p $(ALLREADS_DIR)
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
#	mv $(READS)_ $(READS)
	awk '/^>/ {printit = /+$$/} printit' $(READS)_ >$(READS)
endif

eval_abe: sweep gen_reads
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

eval_multiple_alignments: sweep gen_reads
#	time $(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -a fine 		     -z $(DIR)/out/sweep.params >out/sweep.out
	time $(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -a fine -x		     -z $(DIR)/out/sweep-x.params >out/sweep-x.out
#	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -x -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a-x.out
#	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a.out 2>out/sweep-b-a.cerr
#	$(SWEEP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -o -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a-o.out 2>out/sweep-b-a-o.cerr

#debug: sweep gen_reads
#	$(SWEEP_BIN) -s $(REF) -p simulations/reads/1.fasta $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-1.params >out/sweep-1.out 2>out/sweep-1.cerr

eval_sweep_max_seeds: sweep gen_reads
	@DIR=evals/eval_sweep-K22; \
	mkdir -p $${DIR}; \
	for maxseeds in $(MAX_SEEDS); do \
		for maxmatches in $(MAX_MATCHES); do \
			f=$${DIR}/"sweep-S$${maxseeds}-M$${maxmatches}"; \
			echo "Processing $${f}"; \
			time $(SWEEP_BIN) -s $(REF) -p $(READS) -t 0.0 -x -k 22 -z sweep.params -S $${maxseeds} -M $${maxmatches} $${f}.params >$${f}.paf 2>$${f}.log; \
			paftools.js mapeval $${f}.paf | tee $${f}.eval; \
		done \
    done

eval_sweep: sweepmap gen_reads
	mkdir -p $(OUTDIR)
	$(TIME_CMD) -o $(OUTDIR)/sweep.time ./$(SWEEP_BIN) -s $(REF) -p $(READS)_ -z $(OUTDIR)/sweep.params -x -t $(T) -k $(K) -r $(R) -S $(S) -M $(M) 2> >(tee $(OUTDIR)/sweep.log) >$(OUTDIR)/sweep.paf 
	paftools.js mapeval $(OUTDIR)/sweep.paf | tee $(OUTDIR)/sweep.eval
	paftools.js mapeval -Q 0 $(OUTDIR)/sweep.paf >$(OUTDIR)/sweep.wrong

eval_sweep_slow: sweepmap gen_reads
	mkdir -p $(OUTDIR)
	$(TIME_CMD) -o $(OUTDIR)/sweep-slow.time ./$(SWEEP_BIN) -s $(REF) -p $(READS) -z $(OUTDIR)/sweep-slow.params -x -t $(T_SLOW) -k $(K_SLOW) -r $(R_SLOW) -S $(S_SLOW) -M $(M_SLOW) 2> >(tee $(OUTDIR)/sweep-slow.log) >$(OUTDIR)/sweep-slow.paf 
	paftools.js mapeval $(OUTDIR)/sweep-slow.paf | tee $(OUTDIR)/sweep-slow.eval
	paftools.js mapeval -Q 0 $(OUTDIR)/sweep-slow.paf >$(OUTDIR)/sweep-slow.wrong

eval_winnowmap: gen_reads
	mkdir -p $(OUTDIR)
	if [ ! -f winnowmap_$(REF)_repetitive_k15.txt ]; then \
		$(MERYL_BIN) count k=15 output merylDB $(REF); \
		$(MERYL_BIN) print greater-than distinct=0.9998 merylDB > winnowmap_$(REF)_repetitive_k15.txt; \
	fi
	$(TIME_CMD) -o $(OUTDIR)/winnowmap.time $(WINNOWMAP_BIN) -W winnowmap_$(REF)_repetitive_k15.txt -x map-pb -t 1 $(REF) $(READS) 2> >(tee $(OUTDIR)/winnowmap.log) >$(OUTDIR)/winnowmap.paf 
	paftools.js mapeval $(OUTDIR)/winnowmap.paf | tee $(OUTDIR)/winnowmap.eval

eval_minimap: gen_reads
	mkdir -p $(OUTDIR)
	$(TIME_CMD) -o $(OUTDIR)/minimap.time $(MINIMAP_BIN) -x map-hifi -t 1 --secondary=no $(REF) $(READS) 2> >(tee $(OUTDIR)/minimap.log) >$(OUTDIR)/minimap.paf 
	paftools.js mapeval $(OUTDIR)/minimap.paf | tee $(OUTDIR)/minimap.eval

eval_blend: gen_reads
	mkdir -p $(OUTDIR)
	$(TIME_CMD) -o $(OUTDIR)/blend.time $(BLEND_BIN) -x map-hifi -t 1 -N 0 $(REF) $(READS) 2> >(tee $(OUTDIR)/blend.log) >$(OUTDIR)/blend.paf 
	paftools.js mapeval $(OUTDIR)/blend.paf | tee $(OUTDIR)/blend.eval

eval_veritymap:
	mkdir -p $(OUTDIR)
	$(TIME_CMD) -o $(OUTDIR)/veritymap.time $(VERITYMAP_BIN) -d hifi -t 1 --reads $(READS) $(REF) -o $(OUTDIR)/veritymap 2> >(tee $(OUTDIR)/veritymap.log) >$(OUTDIR)/veritymap.paf 
	paftools.js mapeval $(OUTDIR)/veritymap.paf | tee $(OUTDIR)/veritymap.eval


eval_mapquik: gen_reads
	mkdir -p $(OUTDIR)
	$(TIME_CMD) -o $(OUTDIR)/mapquik.time $(MAPQUIK_BIN) $(READS) --reference $(REF) --threads 1 -p $(OUTDIR)/mapquik | tee $(OUTDIR)/mapquik.log
	paftools.js mapeval $(OUTDIR)/mapquik.paf | tee $(OUTDIR)/mapquik.eval

eval_eskemap: gen_reads
	mkdir -p $(OUTDIR)
	$(TIME_CMD) -o $(OUTDIR)/eskemap.time $(ESKEMAP_BIN) -p $(READS) -s $(REF) -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -N >$(OUTDIR)/eskemap.out 2>$(OUTDIR)/eskemap.log

eval: eval_sweep eval_minimap eval_mapquik eval_winnowmap

clean:
	rm -r $(SWEEP_BIN)

clean_evals:
	rm -r $(ALLREADS_DIR)
	rm -r $(ALLOUT_DIR)
