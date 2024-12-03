# Tool commands
SHMAP_CMD := SHMAP_BIN + " -s " + REF + " -p " + READS + " -k " + K + " -r " + R + " -t " + T + " -x"
MINIMAP_INDEX_CMD := "minimap2 -x map-pb -d " + MINIMAP_PREF + ".mmi " + REF
MINIMAP_MAP_CMD := "minimap2 -x map-pb " + MINIMAP_PREF + ".mmi " + READS

# Variables
SHMAP_BIN := "./release/shmap"
PAFTOOLS := "./ext/paftools.js"
TIME_CMD := "/usr/bin/time -f \"%U\t%M\""
DIR := "evals"
REF_DIR := DIR + "/refs"
ALLREADS_DIR := DIR + "/reads"
ALLOUT_DIR := DIR + "/out"

# Default values for parameters
K := "25"
R := "0.05"
T := "0.7"
REFNAME := "t2tChrY"
ACCURACY := "0.99"
DEPTH := "1"
MEANLEN := "10000"
READSIM_REFNAME := REFNAME

# Computed paths
REF := REF_DIR + "/" + REFNAME + ".fa"
READSIM_REF := REF_DIR + "/" + READSIM_REFNAME + ".fa"
READS_PREFIX := REFNAME + "-reads" + READSIM_REFNAME + "-a" + ACCURACY + "-d" + DEPTH + "-l" + MEANLEN
READS := ALLREADS_DIR + "/" + READS_PREFIX + ".fa"
ONE_READ := ALLREADS_DIR + "/" + READS_PREFIX + ".oneread.fa"
SHMAP_PREF := ALLOUT_DIR + "/shmap/" + READS_PREFIX + "/shmap"
MINIMAP_PREF := ALLOUT_DIR + "/minimap/" + READS_PREFIX + "/minimap"

# Tool-specific configurations
shmap:
    just eval "shmap" {{SHMAP_PREF}} "{{SHMAP_CMD}}" "{{SHMAP_CMD}} -z {{SHMAP_PREF}}.params"

minimap:
    just eval "minimap2" {{MINIMAP_PREF}} "{{MINIMAP_INDEX_CMD}}" "{{MINIMAP_MAP_CMD}}"

# Generate simulated reads if they don't exist
gen_reads:
    #!/usr/bin/env sh
    if [ ! -f {{READS}} ]; then
        echo {{READS}}
        mkdir -p {{ALLREADS_DIR}}
        pbsim \
            {{READSIM_REF}} \
            --model_qc {{DIR}}/model_qc_clr \
            --accuracy-mean {{ACCURACY}} \
            --accuracy-sd 0 \
            --depth {{DEPTH}} \
            --prefix {{READS_PREFIX}} \
            --length-mean {{MEANLEN}}

        samtools faidx {{READSIM_REF}}
        paftools.js pbsim2fq {{READSIM_REF}}.fai "{{READS_PREFIX}}"_*.maf >{{READS}}.unshuf
        ~/miniconda3/bin/seqkit shuffle {{READS}}.unshuf -o {{READS}}
        rm -f "{{READS_PREFIX}}"_*.maf "{{READS_PREFIX}}"_*.ref "{{READS_PREFIX}}"_*.fastq
    fi
    if [ ! -f {{ONE_READ}} ]; then
        head -n 2 {{READS}} >{{ONE_READ}}
    fi

# Generic evaluation logic
eval name prefix index_cmd map_cmd: gen_reads
    #!/usr/bin/env bash
    echo "Evaluating {{name}}..."
    mkdir -p $(dirname {{prefix}})
    {{TIME_CMD}} -o {{prefix}}.index.time {{index_cmd}} 2>/dev/null
    {{TIME_CMD}} -o {{prefix}}.time {{map_cmd}} 2> >(tee {{prefix}}.log) > {{prefix}}.paf
    {{PAFTOOLS}} mapeval -r 0.1 {{prefix}}.paf 2>/dev/null | tee {{prefix}}.eval || true
    {{PAFTOOLS}} mapeval -r 0.1 -Q 0 {{prefix}}.paf >{{prefix}}.wrong 2>/dev/null || true
