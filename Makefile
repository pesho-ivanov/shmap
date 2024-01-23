SHELL := /bin/bash
CC = g++
CXX_STANDARD = -std=c++20
CFLAGS = -O2 -g -march=native -lm -lpthread -Igtl/ -Wall -Wextra -Wno-unused-parameter -Wno-unused-result -Wno-comment -fpermissive #-Wconversion 
LIBS = minimap2/libminimap2.a -lz
DEPFLAGS = -MMD -MP

TIME_CMD = /usr/bin/time -f "%U\t%M"

SRCS = src/sweepmap.cpp src/sweepmap.h src/io.h src/sketch.h src/utils.h src/kseq.h
SWEEPMAP_BIN = ./sweepmap
MINIMAP_BIN = minimap2
BLEND_BIN = ~/libs/blend/bin/blend
MAPQUIK_BIN = mapquik
MERYL_BIN = ~/libs/Winnowmap/bin/meryl
WINNOWMAP_BIN = ~/libs/Winnowmap/bin/winnowmap

REFNAME ?= chm13
ACCURACY ?= 0.99
DEPTH ?= 1
MEANLEN ?= 10000

K ?= 22
R ?= 0.1 
S ?= 300
M ?= 3000
T ?= 0.0

K_SLOW ?= 22
R_SLOW ?= 0.1
S_SLOW ?= 1000
M_SLOW ?= 30000
T_SLOW ?= 0.0

DIR = evals
REF_DIR = $(DIR)/refs
REF = $(REF_DIR)/$(REFNAME).fa
READS_PARAMS = $(REFNAME)-a$(ACCURACY)-d$(DEPTH)-l$(MEANLEN)
READS_PREFIX = reads-$(READS_PARAMS)
ALLREADS_DIR = $(DIR)/reads
READS = $(ALLREADS_DIR)/$(READS_PREFIX).fa
ONE_READ = $(ALLREADS_DIR)/$(READS_PREFIX).oneread.fa
ALLOUT_DIR = $(DIR)/out
OUTDIR = $(ALLOUT_DIR)/$(READS_PARAMS)

SWEEPMAP_PREF = $(OUTDIR)/sweepmap/sweepmap
SWEEPMAP_SLOW_PREF = $(OUTDIR)/sweepmap-slow/sweepmap-slow
MINIMAP_PREF = $(OUTDIR)/minimap/minimap
BLEND_PREF = $(OUTDIR)/blend/blend
MAPQUIK_PREF = $(OUTDIR)/mapquik/mapquik
WINNOWMAP_PREF = $(OUTDIR)/winnowmap/winnowmap

MAX_SEEDS = 10 30 100 300 1000 3000 10000
MAX_MATCHES = 100 300 1000 3000 10000 30000 100000 300000

Ks = 14 16 18 20 22 24 26
Rs = 0.01 0.05 0.1 0.15 0.2

all: sweepmap

$(SWEEPMAP_BIN): $(SRCS)
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
	paftools.js pbsim2fq $(REF).fai "$(READS_PREFIX)"_*.maf >$(READS)
	rm -f "$(READS_PREFIX)"_*.maf "$(READS_PREFIX)"_*.ref "$(READS_PREFIX)"_*.fastq
#	mv $(READS)_ $(READS)
	head -n 2 $(READS) >$(ONE_READ)
#	awk '/^>/ {printit = /+$$/} printit' $(READS)_ >$(READS)
endif

eval_sketching: sweepmap gen_reads
	@DIR=$(OUTDIR)/sketching; \
	mkdir -p $${DIR}; \
	for k in $(Ks); do \
		for r in $(Rs); do \
			f=$${DIR}/"sweepmap-K$${k}-R$${r}"; \
			echo "Processing $${f}"; \
			$(TIME_CMD) -o $${f}.index.time $(SWEEPMAP_BIN) -s $(REF) -p $(ONE_READ) -x -t $(T) -k $${k} -r $${r} -S $(S) -M $(M) 2>&1 >/dev/null; \
			$(TIME_CMD) -o $${f}.time $(SWEEPMAP_BIN) -s $(REF) -p $(READS) -z $${f}.params -x -t $(T) -k $${k} -r $${r} -S $(S) -M $(M) 2> >(tee $${f}.log) >$${f}.paf; \
			paftools.js mapeval $${f}.paf | tee $${f}.eval; \
		done \
    done

eval_thinning: sweepmap gen_reads
	@DIR=$(OUTDIR)/thinning; \
	mkdir -p $${DIR}; \
	for s in $(MAX_SEEDS); do \
		for m in $(MAX_MATCHES); do \
			f=$${DIR}/"sweepmap-S$${s}-M$${m}"; \
			echo "Processing $${f}"; \
			$(TIME_CMD) -o $${f}.index.time $(SWEEPMAP_BIN) -s $(REF) -p $(ONE_READ) -x -t $(T) -k $(K) -r $(R) -S $${s} -M $${m} 2>&1 >/dev/null; \
			$(TIME_CMD) -o $${f}.time $(SWEEPMAP_BIN) -s $(REF) -p $(READS) -z $${f}.params -x -t $(T) -k $(K) -r $(R) -S $${s} -M $${m} 2> >(tee $${f}.log) >$${f}.paf; \
			paftools.js mapeval $${f}.paf | tee $${f}.eval; \
		done \
    done

eval_sweepmap: sweepmap gen_reads
	@mkdir -p $(shell dirname $(SWEEPMAP_PREF))
	$(TIME_CMD) -o $(SWEEPMAP_PREF).index.time $(SWEEPMAP_BIN) -s $(REF) -p $(ONE_READ) -x -t $(T) -k $(K) -r $(R) -S $(S) -M $(M) 2>&1 >/dev/null
	$(TIME_CMD) -o $(SWEEPMAP_PREF).time $(SWEEPMAP_BIN) -s $(REF) -p $(READS) -z $(SWEEPMAP_PREF).params -x -t $(T) -k $(K) -r $(R) -S $(S) -M $(M) 2> >(tee $(SWEEPMAP_PREF).log) >$(SWEEPMAP_PREF).paf
	paftools.js mapeval $(SWEEPMAP_PREF).paf | tee $(SWEEPMAP_PREF).eval
	paftools.js mapeval -Q 60 $(SWEEPMAP_PREF).paf >$(SWEEPMAP_PREF).wrong

eval_sweepmap_slow: sweepmap gen_reads
	@mkdir -p $(shell dirname $(SWEEPMAP_SLOW_PREF))
	$(TIME_CMD) -o $(SWEEPMAP_SLOW_PREF).index.time $(SWEEPMAP_BIN) -s $(REF) -p $(ONE_READ) -z $(SWEEPMAP_SLOW_PREF).params -x -t $(T_SLOW) -k $(K_SLOW) -r $(R_SLOW) -S $(S_SLOW) -M $(M_SLOW) >/dev/null 2>&1
	$(TIME_CMD) -o $(SWEEPMAP_SLOW_PREF).time $(SWEEPMAP_BIN) -s $(REF) -p $(READS) -z $(SWEEPMAP_SLOW_PREF).params -x -t $(T_SLOW) -k $(K_SLOW) -r $(R_SLOW) -S $(S_SLOW) -M $(M_SLOW) 2> >(tee $(SWEEPMAP_SLOW_PREF).log) >$(SWEEPMAP_SLOW_PREF).paf 
	paftools.js mapeval $(SWEEPMAP_SLOW_PREF).paf | tee $(SWEEPMAP_SLOW_PREF).eval
	paftools.js mapeval -Q 0 $(SWEEPMAP_SLOW_PREF).paf >$(SWEEPMAP_SLOW_PREF).wrong

eval_winnowmap: gen_reads
	@mkdir -p $(shell dirname $(WINNOWMAP_PREF))
	if [ ! -f $(REF_DIR)/winnowmap_$(REF)_repetitive_k15.txt ]; then \
		$(MERYL_BIN) count k=15 output merylDB $(REF); \
		$(MERYL_BIN) print greater-than distinct=0.9998 merylDB > $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt; \
	fi
	$(TIME_CMD) -o $(WINNOWMAP_PREF).index.time $(WINNOWMAP_BIN) -W $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt -x map-pb -t 1 $(REF) $(ONE_READ) >/dev/null 2>&1 
	$(TIME_CMD) -o $(WINNOWMAP_PREF).time $(WINNOWMAP_BIN) -W $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt -x map-pb -t 1 $(REF) $(READS) 2> >(tee $(WINNOWMAP_PREF).log) >$(WINNOWMAP_PREF).paf 
	paftools.js mapeval $(WINNOWMAP_PREF).paf | tee $(WINNOWMAP_PREF).eval

eval_minimap: gen_reads
	@mkdir -p $(shell dirname $(MINIMAP_PREF))
	$(TIME_CMD) -o $(MINIMAP_PREF).index.time $(MINIMAP_BIN) -x map-hifi -t 1 --secondary=no $(REF) $(ONE_READ) >/dev/null 2>&1
	$(TIME_CMD) -o $(MINIMAP_PREF).time $(MINIMAP_BIN) -x map-hifi -t 1 --secondary=no $(REF) $(READS) 2> >(tee $(MINIMAP_PREF).log) >$(MINIMAP_PREF).paf 
	paftools.js mapeval $(MINIMAP_PREF).paf | tee $(MINIMAP_PREF).eval

eval_blend: gen_reads
	@mkdir -p $(shell dirname $(BLEND_PREF))
	$(TIME_CMD) -o $(BLEND_PREF).index.time $(BLEND_BIN) -x map-hifi -t 1 -N 0 $(REF) $(ONE_READ) >/dev/null 2>&1
	$(TIME_CMD) -o $(BLEND_PREF).time $(BLEND_BIN) -x map-hifi -t 1 -N 0 $(REF) $(READS) 2> >(tee $(BLEND_PREF).log) >$(BLEND_PREF).paf 
	paftools.js mapeval $(BLEND_PREF).paf | tee $(BLEND_PREF).eval

eval_mapquik: gen_reads
	@mkdir -p $(shell dirname $(MAPQUIK_PREF))
	$(TIME_CMD) -o $(MAPQUIK_PREF).index.time $(MAPQUIK_BIN) $(ONE_READ) --reference $(REF) --threads 1  >/dev/null
	$(TIME_CMD) -o $(MAPQUIK_PREF).time $(MAPQUIK_BIN) $(READS) --reference $(REF) --threads 1 -p $(MAPQUIK_PREF) | tee $(MAPQUIK_PREF).log
	paftools.js mapeval $(MAPQUIK_PREF).paf | tee $(MAPQUIK_PREF).eval

eval_all: eval_sweepmap eval_sweepmap_slow eval_mapquik eval_blend eval_minimap #eval_winnowmap

clean:
	rm -r $(SWEEPMAP_BIN)

clean_evals:
	rm -r $(ALLREADS_DIR)
	rm -r $(ALLOUT_DIR)

#EDLIB_BIN = ~/libs/edlib/build/bin/edlib-aligner
#VERITYMAP_BIN = ~/libs/VerityMap/veritymap/main.py
#ESKEMAP_BIN = eskemap

#eval_veritymap:
#	mkdir -p $(OUTDIR)
#	$(TIME_CMD) -o $(OUTDIR)/veritymap.time $(VERITYMAP_BIN) -d hifi -t 1 --reads $(READS) $(REF) -o $(OUTDIR)/veritymap 2> >(tee $(OUTDIR)/veritymap.log) >$(OUTDIR)/veritymap.paf 
#	paftools.js mapeval $(OUTDIR)/veritymap.paf | tee $(OUTDIR)/veritymap.eval

#eval_eskemap: gen_reads
#	mkdir -p $(OUTDIR)
#	$(TIME_CMD) -o $(OUTDIR)/eskemap.time $(ESKEMAP_BIN) -p $(READS) -s $(REF) -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -N >$(OUTDIR)/eskemap.out 2>$(OUTDIR)/eskemap.log

#eval_abe: sweep gen_reads
#	mkdir -p ../out
##	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 >out/sweep.out
#
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -z $(DIR)/out/sweep-b.params >out/sweep-b.out
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a.out
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -z $(DIR)/out/sweep-b-e.params >out/sweep-b-e.out
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -a extend_equally -z $(DIR)/out/sweep-b-e-a.params >out/sweep-b-e-a.out
#
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -z $(DIR)/out/sweep.params >out/sweep.out
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -a extend_equally -z $(DIR)/out/sweep-a.params >out/sweep-a.out
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -e consecutive -z $(DIR)/out/sweep-e.params >out/sweep-e.out
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -e consecutive -a extend_equally -z $(DIR)/out/sweep-e-a.params >out/sweep-e-a.out

#eval_multiple_alignments: sweep gen_reads
#	time $(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -a fine 		     -z $(DIR)/out/sweep.params >out/sweep.out
#	time $(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -a fine -x		     -z $(DIR)/out/sweep-x.params >out/sweep-x.out
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -x -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a-x.out
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a.out 2>out/sweep-b-a.cerr
#	$(SWEEPMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -o -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a-o.out 2>out/sweep-b-a-o.cerr

#debug: sweep gen_reads
#	$(SWEEPMAP_BIN) -s $(REF) -p simulations/reads/1.fasta $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-1.params >out/sweep-1.out 2>out/sweep-1.cerr