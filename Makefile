SHELL := /bin/bash
CC = g++
CXX_STANDARD = -std=c++20
DEBUG_FLAGS = -g -DDEBUG
RELEASE_FLAGS = -O3 -DNDEBUG -flto
CFLAGS = -march=native -lm -lpthread -Igtl/ -isystem ext/ -Wall -Wextra -Wno-unused-parameter -Wno-unused-result -Wno-comment -fpermissive -flto -fopenmp #-Wconversion 
ifeq ($(DEBUG), 1)
    CFLAGS += $(DEBUG_FLAGS)
	SHMAP_BIN = ./debug/shmap
	FLAG_FILE = debug_flag
else
    CFLAGS += $(RELEASE_FLAGS)
	SHMAP_BIN = ./release/shmap
	FLAG_FILE = release_flag
endif
LIBS = -lz
DEPFLAGS = -MMD -MP

TIME_CMD = /usr/bin/time -f "%U\t%M"

SRCS = src/map.cpp src/mapper.cpp src/io.cpp ext/edlib.cpp
OBJS = $(SRCS:.cpp=.o)
DEPS = $(OBJS:.o=.d)

MINIMAP_BIN = minimap2
MM2_BIN = ~/libs/mm2-fast/minimap2
BLEND_BIN = ~/libs/blend/bin/blend
MAPQUIK_BIN = mapquik
MERYL_BIN = ~/libs/Winnowmap/bin/meryl
WINNOWMAP_BIN = ~/libs/Winnowmap/bin/winnowmap
SURVIVOR_BIN = ~/libs/SURVIVOR/Debug/SURVIVOR

INF = 999999

ARGS_SHMAP ?= 

REFNAME ?= t2tChrY
ACCURACY ?= 0.99
DEPTH ?= 1
MEANLEN ?= 10000
READSIM_REFNAME ?= $(REFNAME)

K ?= 25
R ?= 0.05
T ?= 0.75
MIN_DIFF ?= 0.02
MAX_OVERLAP ?= 0.5

THETAS = 0.95 0.9 0.85 0.8 0.75 0.7 0.65 0.6 0.55 0.5  
PAUL_THETAS = 1.0 0.9 0.8 0.7 0.6 0.5 0.4 0.3

Ks = 14 16 18 20 22 24 26
Rs = 0.01 0.05 0.1 0.15 0.2

M ?= 100000
S ?= 30000

MAX_SEEDS = 10 30 100 300 1000 3000 10000
MAX_MATCHES = 100 300 1000 3000 10000 30000 100000 300000

K_SLOW ?= $(K)
R_SLOW ?= $(R)
T_SLOW ?= $(T)

DIR = evals
REF_DIR = $(DIR)/refs
REF = $(REF_DIR)/$(REFNAME).fa
READSIM_REF = $(REF_DIR)/$(READSIM_REFNAME).fa
READS_PREFIX ?= $(REFNAME)-reads$(READSIM_REFNAME)-a$(ACCURACY)-d$(DEPTH)-l$(MEANLEN)
ALLREADS_DIR = $(DIR)/reads
READS = $(ALLREADS_DIR)/$(READS_PREFIX).fa
#ifeq ($(wildcard $(READS)),)
#	READS = $(foreach ext,fasta fq fastq,$(if $(wildcard $(ALLREADS_DIR)/$(READS_PREFIX).$(ext)),$(ALLREADS_DIR)/$(READS_PREFIX).$(ext),))
#endif
#READS = $(ALLREADS_DIR)/$(READS_PREFIX).fa
#ifeq ($(wildcard $(READS)),)
#    READS = $(ALLREADS_DIR)/$(READS_PREFIX).fq
#endif
ONE_READ = $(ALLREADS_DIR)/$(READS_PREFIX).oneread.fa
ALLOUT_DIR = $(DIR)/out
#OUTDIR = $(ALLOUT_DIR)/$(READS_PREFIX)/k$(K)-r$(R)-s$(S)-m$(M)-t$(T)
OUTDIR = $(ALLOUT_DIR)/$(READS_PREFIX)

PAFTOOLS = ./ext/paftools.js
LIFT_FASTA = $(DIR)/convert_fasta_with_args.py
CHAIN_FILE = $(DIR)/refs/hg002v1.1_to_CHM13v2.0.chain

SHMAP_PREF         = $(ALLOUT_DIR)/shmap/$(READS_PREFIX)/shmap
SHMAP_NOPRUNE_PREF = $(ALLOUT_DIR)/shmap-noprune/$(READS_PREFIX)/shmap-noprune
SHMAP_ONESWEEP_PREF= $(ALLOUT_DIR)/shmap-onesweep/$(READS_PREFIX)/shmap-onesweep # NOT USED
MINIMAP_PREF       = $(ALLOUT_DIR)/minimap/$(READS_PREFIX)/minimap
MM2_PREF           = $(ALLOUT_DIR)/mm2-fast/$(READS_PREFIX)/mm2-fast   # NOT USED
BLEND_PREF         = $(ALLOUT_DIR)/blend/$(READS_PREFIX)/blend
MAPQUIK_PREF       = $(ALLOUT_DIR)/mapquik/$(READS_PREFIX)/mapquik
WINNOWMAP_PREF     = $(ALLOUT_DIR)/winnowmap/$(READS_PREFIX)/winnowmap

all: $(SHMAP_BIN)

# Add the flag file as a prerequisite and create it if needed
$(SHMAP_BIN): $(OBJS) | $(FLAG_FILE)
	mkdir -p $(shell dirname $@)
	$(CXX) $(CXX_STANDARD) $(CFLAGS) -o $@ $(OBJS) $(LIBS)

# Create flag files if they don't exist
$(FLAG_FILE):
	@rm -f debug_flag release_flag
	@touch $@

%.o: %.cpp
	$(CXX) $(CXX_STANDARD) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

-include $(DEPS)

clean:
	rm -f $(OBJS) debug_flag release_flag
	rm -rf debug release
	$(MAKE) clean_test

simulate_SVs:
	cd $(REF_DIR);\
	$(SURVIVOR_BIN) simSV $(REFNAME).fa SURVIVOR.params 0 0 $(REFNAME)-SVs;\
	mv $(REFNAME)-SVs.fasta $(REFNAME)-SVs.fa

gen_reads:
ifeq ($(wildcard $(READS)),)
	echo $(READS)
	mkdir -p $(ALLREADS_DIR)
	pbsim \
		   $(READSIM_REF) \
		   --model_qc $(DIR)/model_qc_clr \
		   --accuracy-mean $(ACCURACY)\
		   --accuracy-sd 0\
		   --depth $(DEPTH)\
		   --prefix $(READS_PREFIX)\
		   --length-mean $(MEANLEN)

	samtools faidx $(READSIM_REF)
	-paftools.js pbsim2fq $(READSIM_REF).fai "$(READS_PREFIX)"_*.maf >$(READS).unshuf
	~/miniconda3/bin/seqkit shuffle $(READS).unshuf -o $(READS)
	rm -f "$(READS_PREFIX)"_*.maf "$(READS_PREFIX)"_*.ref "$(READS_PREFIX)"_*.fastq

	if [ "$(READSIM_REFNAME)" != "$(REFNAME)" ]; then \
		echo "Lifting over reads from $(READSIM_REFNAME) to $(REFNAME)"; \
		python $(LIFT_FASTA) -f $(READS) -c $(CHAIN_FILE) -o $(READS).lifted; \
		mv $(READS).lifted $(READS); \
	else \
		echo "READSIM_REFNAME and REFNAME are the same. No liftover required."; \
	fi

# take only positive strand reads
#	mv $(READS)_ $(READS)
#	awk '/^>/ {printit = /+$$/} printit' $(READS)_ >$(READS)
endif
ifeq ($(wildcard $(ONE_READ)),)
	head -n 2 $(READS) >$(ONE_READ)
endif

eval_sketching: sweepmap gen_reads
	@DIR=$(OUTDIR)/sketching; \
	mkdir -p $${DIR}; \
	for k in $(Ks); do \
		for r in $(Rs); do \
			f=$${DIR}/"sweepmap-K$${k}-R$${r}"; \
			echo "Processing $${f}"; \
			$(TIME_CMD) -o $${f}.index.time $(SHMAP_BIN) -s $(REF) -p $(ONE_READ) -x -t $(T) -k $${k} -r $${r} -S $(S) -M $(M) 2>&1 >/dev/null; \
			$(TIME_CMD) -o $${f}.time $(SHMAP_BIN) -s $(REF) -p $(READS) -z $${f}.params -x -t $(T) -k $${k} -r $${r} -S $(S) -M $(M) 2> >(tee $${f}.log) >$${f}.paf; \
			-paftools.js mapeval $${f}.paf | tee $${f}.eval; \
		done \
    done

eval_thinning: sweepmap gen_reads
	@DIR=$(OUTDIR)/thinning; \
	mkdir -p $${DIR}; \
	for s in $(MAX_SEEDS); do \
		for m in $(MAX_MATCHES); do \
			f=$${DIR}/"sweepmap-S$${s}-M$${m}"; \
			echo "Processing $${f}"; \
			$(TIME_CMD) -o $${f}.index.time $(SHMAP_BIN) -s $(REF) -p $(ONE_READ) -x -t $(T) -k $(K) -r $(R) -S $${s} -M $${m} 2>&1 >/dev/null; \
			$(TIME_CMD) -o $${f}.time $(SHMAP_BIN) -s $(REF) -p $(READS) -z $${f}.params -x -t $(T) -k $(K) -r $(R) -S $${s} -M $${m} 2> >(tee $${f}.log) >$${f}.paf; \
			-paftools.js mapeval $${f}.paf | tee $${f}.eval; \
		done \
    done

fdr_per_theta: $(SHMAP_BIN) gen_reads
	@DIR=$(OUTDIR)/fdr_per_theta; \
	mkdir -p $${DIR}; \
	for t in $(THETAS); do \
		f=$${DIR}/"shmap-T$${t}"; \
		echo "Processing $${f}"; \
		$(TIME_CMD) -o $${f}.time $(SHMAP_BIN) -s $(REF) -p $(READS) -z $${f}.params -x -t $${t} -k $(K) -r $(R) 2> >(tee $${f}.log) >$${f}.paf; \
		-paftools.js mapeval $${f}.paf | tee $${f}.eval; \
	done

eval_shmap: $(SHMAP_BIN) gen_reads
	@mkdir -p $(shell dirname $(SHMAP_PREF))
#	$(TIME_CMD) -o $(SHMAP_PREF).index.time $(SHMAP_BIN) -s $(REF) -p $(ONE_READ) -k $(K) -r $(R) -t $(T) $(ARGS_SHMAP) 2>/dev/null >/dev/null
	$(TIME_CMD) -o $(SHMAP_PREF).time $(SHMAP_BIN) -s $(REF) -p $(READS) -z $(SHMAP_PREF).params -k $(K) -r $(R) -t $(T) -d $(MIN_DIFF) -o $(MAX_OVERLAP) $(ARGS_SHMAP)    2> >(tee $(SHMAP_PREF).log) > $(SHMAP_PREF).paf
	$(PAFTOOLS) mapeval -r 0.1 $(SHMAP_PREF).paf 2>/dev/null | tee $(SHMAP_PREF).eval || true
	$(PAFTOOLS) mapeval -r 0.1 -Q 0 $(SHMAP_PREF).paf >$(SHMAP_PREF).wrong 2>/dev/null || true

eval_shmap_noprune: $(SHMAP_BIN) gen_reads
	@mkdir -p $(shell dirname $(SHMAP_NOPRUNE_PREF))
	$(TIME_CMD) -o $(SHMAP_NOPRUNE_PREF).index.time $(SHMAP_BIN) -s $(REF) -p $(ONE_READ) -k $(K) -r $(R) -t $(T) -M $(M) -x -b 2>/dev/null >/dev/null
	$(TIME_CMD) -o $(SHMAP_NOPRUNE_PREF).time $(SHMAP_BIN) -s $(REF) -p $(READS) -z $(SHMAP_NOPRUNE_PREF).params -k $(K) -r $(R) -t $(T) -M $(M) -x -b    2> >(tee $(SHMAP_NOPRUNE_PREF).log) > $(SHMAP_NOPRUNE_PREF).paf
	-paftools.js mapeval -r 0.1 $(SHMAP_NOPRUNE_PREF).paf | tee $(SHMAP_NOPRUNE_PREF).eval
	-$(PAFTOOLS) mapeval -r 0.1 -Q 0 $(SHMAP_NOPRUNE_PREF).paf >$(SHMAP_NOPRUNE_PREF).wrong

eval_shmap_onesweep: $(SHMAP_BIN) gen_reads
	@mkdir -p $(shell dirname $(SHMAP_ONESWEEP_PREF))
	$(TIME_CMD) -o $(SHMAP_ONESWEEP_PREF).index.time $(SHMAP_BIN) -s $(REF) -p $(ONE_READ) -k $(K) -r $(R) -t $(T) -M $(M) -x -B 2>/dev/null >/dev/null
	$(TIME_CMD) -o $(SHMAP_ONESWEEP_PREF).time $(SHMAP_BIN) -s $(REF) -p $(READS) -z $(SHMAP_ONESWEEP_PREF).params -k $(K) -r $(R) -t $(T) -M $(M) -x -B    2> >(tee $(SHMAP_ONESWEEP_PREF).log) > $(SHMAP_ONESWEEP_PREF).paf
	-paftools.js mapeval -r 0.1 $(SHMAP_ONESWEEP_PREF).paf | tee $(SHMAP_ONESWEEP_PREF).eval
	-$(PAFTOOLS) mapeval -r 0.1 -Q 0 $(SHMAP_ONESWEEP_PREF).paf >$(SHMAP_ONESWEEP_PREF).wrong

eval_winnowmap: gen_reads
	@mkdir -p $(shell dirname $(WINNOWMAP_PREF))
	if [ ! -f $(REF_DIR)/winnowmap_$(REF)_repetitive_k15.txt ]; then \
		$(MERYL_BIN) count k=15 output merylDB $(REF); \
		$(MERYL_BIN) print greater-than distinct=0.9998 merylDB > $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt; \
	fi
	$(TIME_CMD) -o $(WINNOWMAP_PREF).index.time $(WINNOWMAP_BIN) -W $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt -x map-pb -t 1 --secondary=no --sv-off $(REF) $(ONE_READ) >/dev/null 2>/dev/null
	$(TIME_CMD) -o $(WINNOWMAP_PREF).time $(WINNOWMAP_BIN) -W $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt -x map-pb -t 1 --secondary=no --sv-off -M 0 --hard-mask-level $(REF) $(READS) 2> >(tee $(WINNOWMAP_PREF).log) >$(WINNOWMAP_PREF).paf 
	-paftools.js mapeval $(WINNOWMAP_PREF).paf | tee $(WINNOWMAP_PREF).eval

eval_minimap: gen_reads
	@mkdir -p $(shell dirname $(MINIMAP_PREF))
	$(TIME_CMD) -o $(MINIMAP_PREF).index.time $(MINIMAP_BIN) -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level -a $(REF) $(ONE_READ) >/dev/null 2>/dev/null
#	$(TIME_CMD) -o $(MINIMAP_PREF).time       $(MINIMAP_BIN) -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level -a $(REF) $(READS) 2> >(tee $(MINIMAP_PREF).log) >$(MINIMAP_PREF).bam
	$(TIME_CMD) -o $(MINIMAP_PREF).time       $(MINIMAP_BIN) -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level $(REF) $(READS) 2> >(tee $(MINIMAP_PREF).log) >$(MINIMAP_PREF).paf
	-paftools.js mapeval $(MINIMAP_PREF).paf | tee $(MINIMAP_PREF).eval

eval_mm2: gen_reads
	@mkdir -p $(shell dirname $(MM2_PREF))
	$(TIME_CMD) -o $(MM2_PREF).index.time $(MM2_BIN) -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level -a $(REF) $(ONE_READ) >/dev/null 2>/dev/null
	$(TIME_CMD) -o $(MM2_PREF).time       $(MM2_BIN) -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level $(REF) $(READS) 2> >(tee $(MM2_PREF).log) >$(MM2_PREF).paf
	-paftools.js mapeval $(MM2_PREF).paf | tee $(MM2_PREF).eval

eval_blend: gen_reads
	@mkdir -p $(shell dirname $(BLEND_PREF))
	$(TIME_CMD) -o $(BLEND_PREF).index.time $(BLEND_BIN) -x map-hifi -t 1 -N 0 $(REF) $(ONE_READ) >/dev/null 2>/dev/null
	$(TIME_CMD) -o $(BLEND_PREF).time $(BLEND_BIN) -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level $(REF) $(READS) 2> >(tee $(BLEND_PREF).log) >$(BLEND_PREF).paf 
	-paftools.js mapeval $(BLEND_PREF).paf | tee $(BLEND_PREF).eval

eval_mapquik: gen_reads
	@mkdir -p $(shell dirname $(MAPQUIK_PREF))
	$(TIME_CMD) -o $(MAPQUIK_PREF).index.time $(MAPQUIK_BIN) $(ONE_READ) --reference $(REF) --threads 1  >/dev/null
	$(TIME_CMD) -o $(MAPQUIK_PREF).time $(MAPQUIK_BIN) $(READS) --reference $(REF) --threads 1 -p $(MAPQUIK_PREF) | tee $(MAPQUIK_PREF).log
	-paftools.js mapeval $(MAPQUIK_PREF).paf | tee $(MAPQUIK_PREF).eval

#eval_tools: eval_shmap eval_shmap_noprune eval_minimap eval_mapquik eval_blend #eval_winnowmap #eval_shmap_onesweep #eval_mm2 
eval_tools: eval_minimap eval_mapquik eval_blend eval_winnowmap # eval_shmap

eval_tools_on_datasets:
#	make eval_tools REFNAME=t2tChrY DEPTH=10
#	make eval_tools REFNAME=chm13   DEPTH=1
#	make eval_tools REFNAME=t2tChrY DEPTH=10  MEANLEN=24000
#	make eval_tools REFNAME=chm13   READS_PREFIX=HG002_24kb_10G
	make eval_tools REFNAME=chm13   READS_PREFIX=hg002_hifi_10000

eval_shmap_on_datasets:
	make eval_shmap REFNAME=t2tChrY DEPTH=10
	make eval_shmap REFNAME=chm13   DEPTH=1
	make eval_shmap REFNAME=t2tChrY DEPTH=10 	MEANLEN=24000
	make eval_shmap REFNAME=chm13   READS_PREFIX=HG002_24kb_10G

eval_winnowmap_on_datasets:
	make eval_winnowmap REFNAME=t2tChrY DEPTH=10
	make eval_winnowmap REFNAME=chm13   DEPTH=1
	make eval_winnowmap REFNAME=t2tChrY DEPTH=10  MEANLEN=24000
	make eval_winnowmap REFNAME=chm13   READS_PREFIX=HG002_24kb_10G

#SHMAP_PREF         = $(ALLOUT_DIR)/pauls_experiment/$(READS_PREFIX)/
pauls_experiment: $(SHMAP_BIN) gen_reads
	@DIR=$(ALLOUT_DIR)/pauls_experiment; \
	for t in $(PAUL_THETAS); do \
		pref=$${DIR}/"T$${t}"; \
		make eval_shmap ALLOUT_DIR=$${pref} T=$${t} REFNAME=t2tChrY DEPTH=0.1; \
		make eval_shmap ALLOUT_DIR=$${pref} T=$${t} REFNAME=chm13   DEPTH=0.01; \
		make eval_shmap ALLOUT_DIR=$${pref} T=$${t} REFNAME=t2tChrY DEPTH=0.1  MEANLEN=24000; \
	done

#eval_tools_on_SV:
#	make eval_tools REFNAME=t2tChrY READSIM_REFNAME=t2tChrY-SVs DEPTH=0.1

clean_evals:
	rm -r $(ALLREADS_DIR)
	rm -r $(ALLOUT_DIR)

# 10GB reads: (10737418240 bytes):
# curl -r 0-10737418239  https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/deepconsensus_predictions/hg002_24kb/two_smrt_cells/HG002_24kb_132517_061936_2fl_DC_hifi_reads.fastq >HG002_24kb_132517_061936_2fl_DC_hifi_reads.10GB.fastq

# TODO: add curl -r 0-1073741823 https://storage.googleapis.com/brain-genomics-public/research/deepconsensus/publication/deepconsensus_predictions/hg002_24kb/two_smrt_cells/HG002_24kb_132517_061936_2fl_DC_hifi_reads.fastq >HG002_24kb_132517_061936_2fl_DC_hifi_reads.1GB.fastq
# TODO: `seqtk seq -AU` for mapquik
# TODO: `seqtk seq -A HG002_24kb.fastq >HG002_24kb.fa`

# new HG002 dataset
# https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifirevio/m84005_220827_014912_s1.hifi_reads.fastq.gz

#EDLIB_BIN = ~/libs/edlib/build/bin/edlib-aligner
#VERITYMAP_BIN = ~/libs/VerityMap/veritymap/main.py
#ESKEMAP_BIN = eskemap

#eval_veritymap:
#	mkdir -p $(OUTDIR)
#	$(TIME_CMD) -o $(OUTDIR)/veritymap.time $(VERITYMAP_BIN) -d hifi -t 1 --reads $(READS) $(REF) -o $(OUTDIR)/veritymap 2> >(tee $(OUTDIR)/veritymap.log) >$(OUTDIR)/veritymap.paf 
#	-paftools.js mapeval $(OUTDIR)/veritymap.paf | tee $(OUTDIR)/veritymap.eval

#eval_eskemap: gen_reads
#	mkdir -p $(OUTDIR)
#	$(TIME_CMD) -o $(OUTDIR)/eskemap.time $(ESKEMAP_BIN) -p $(READS) -s $(REF) -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -N >$(OUTDIR)/eskemap.out 2>$(OUTDIR)/eskemap.log

#eval_abe: sweep gen_reads
#	mkdir -p ../out
##	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 >out/sweep.out
#
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -z $(DIR)/out/sweep-b.params >out/sweep-b.out
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a.out
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -z $(DIR)/out/sweep-b-e.params >out/sweep-b-e.out
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -e consecutive -a extend_equally -z $(DIR)/out/sweep-b-e-a.params >out/sweep-b-e-a.out

#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -z $(DIR)/out/sweep.params >out/sweep.out
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -a extend_equally -z $(DIR)/out/sweep-a.params >out/sweep-a.out
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -e consecutive -z $(DIR)/out/sweep-e.params >out/sweep-e.out
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -e consecutive -a extend_equally -z $(DIR)/out/sweep-e-a.params >out/sweep-e-a.out

#eval_multiple_alignments: sweep gen_reads
#	time $(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -a fine 		     -z $(DIR)/out/sweep.params >out/sweep.out
#	time $(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -a fine -x		     -z $(DIR)/out/sweep-x.params >out/sweep-x.out
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -x -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a-x.out
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a.out 2>out/sweep-b-a.cerr
#	$(SHMAP_BIN) -s $(REF) -p $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -o -z $(DIR)/out/sweep-b-a.params >out/sweep-b-a-o.out 2>out/sweep-b-a-o.cerr

#debug: sweep gen_reads
#	$(SHMAP_BIN) -s $(REF) -p simulations/reads/1.fasta $(READS) -k 15 -b highAbundKmersMiniK15w10Lrgr100BtStrnds.txt -a extend_equally -z $(DIR)/out/sweep-1.params >out/sweep-1.out 2>out/sweep-1.cerr


#eval_sweepmap_sam: $(SHMAP_BIN) gen_reads
#	@mkdir -p $(shell dirname $(SWEEPMAP_PREF))
#	$(TIME_CMD) -o $(SWEEPMAP_PREF).index.time $(SHMAP_BIN) -s $(REF) -p $(ONE_READ) -x -t $(T) -k $(K) -r $(R) -S $(S) -M $(M) 2>/dev/null >/dev/null
#	$(TIME_CMD) -o $(SWEEPMAP_PREF).time $(SHMAP_BIN) -s $(REF) -p $(READS) -z $(SWEEPMAP_PREF).params -x -t $(T) -k $(K) -r $(R) -S $(S) -M $(M) -a 2> >(tee $(SWEEPMAP_PREF).log) >$(SWEEPMAP_PREF).sam
#	-paftools.js mapeval $(SWEEPMAP_PREF).sam | tee $(SWEEPMAP_PREF).eval
#	@-paftools.js mapeval -Q 60 $(SWEEPMAP_PREF).sam >$(SWEEPMAP_PREF).wrong
#

# Test-specific variables
TEST_SRCS = src/test_shmap.cpp
TEST_OBJS = src/mapper.o src/io.o ext/edlib.o src/test_shmap.o  # Remove map.o as it contains main()
TEST_EXEC = ./test_shmap

# Test target
test: $(TEST_OBJS)
	$(CXX) $(CXX_STANDARD) $(CFLAGS) -o $(TEST_EXEC) $(TEST_OBJS) $(LIBS)
	./$(TEST_EXEC) --no-intro --no-version 2>/dev/null

# Compile test objects with TESTING defined
src/test_shmap.o: src/test_shmap.cpp
	$(CXX) $(CXX_STANDARD) $(CFLAGS) $(DEPFLAGS) -c $< -o $@

clean_test:
	rm -f $(TEST_OBJS) $(TEST_EXEC)