# Snakefile
import os
from collections import OrderedDict

# Directories and Files
OUTDIR    = "evals/snake_out"
REF_DIR   = "evals/refs"
READS_DIR = "evals/reads"
REFNAME   = "chm13-chr1"
REF       = "{REF_DIR}/{REFNAME}.fa"

#READS_PREFIX ?= $(REFNAME)-reads$(READSIM_REFNAME)-a$(ACCURACY)-d$(DEPTH)-l$(MEANLEN)
ACCURACY  = '?'
DEPTH     = 1
MEANLEN   = '?'
READSIM_REFNAME = {REFNAME}
ONE_READ  = "{READS_DIR}/{REFNAME}-reads{READSIM_REFNAME}-a{ACCURACY}-d{DEPTH}-l{MEANLEN}.oneread.fa"
READS     = "{READS_DIR}/{REFNAME}-reads{READSIM_REFNAME}-a{ACCURACY}-d{DEPTH}-l{MEANLEN}.fa"

# Timing and evaluation tool
TIME_CMD  = "/usr/bin/time -f '%U\t%M'"
PAFTTOOLS = "evals/ext/paftools.js"

# Mapping parameters
K  = 25
R  = 0.01
T  = 0.4

# Executable paths (adjust as needed)
SHMAP     = "./release/shmap"
MINIMAP   = "~/libs/minimap2/minimap2"
MM2       = "~/libs/mm2-fast/minimap2"
BLEND     = "~/libs/blend/bin/blend"
MAPQUIK   = "~/libs/mapquik/target/release/mapquik"
MERYL     = "~/libs/Winnowmap/bin/meryl"
WINNOWMAP = "~/libs/Winnowmap/bin/winnowmap"
MASHMAP1  = "~/libs/mashmap"
MASHMAP3  = "~/libs/MashMap/build/bin/mashmap"
ASTARIX   = "~/libs/astarix/release/astarix"

tools = {
    "shmap":    f"{SHMAP} -s {{REF}} -p {{INPUT}} -k {{K}} -r {{R}} -t {{T}} -z {{outdir}}/shmap.params",
    "minimap":  f"{MINIMAP} -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level {{REF}} {{INPUT}}",
    "mm2":      f"{MM2} -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level {{REF}} {{INPUT}}",
    "blend":    f"{BLEND} -x map-hifi -t 1 -N 0 {{REF}} {{INPUT}}",
    "mapquik":  f"{MAPQUIK} {{INPUT}} --reference {{REF}} --threads 1",
    "winnowmap":f"{MERYL} count k=15 output merylDB {{REF}}; {MERYL} print greater-than distinct=0.9998 merylDB > {{REF_DIR}}/winnowmap_{{REFNAME}}_repetitive_k15.txt; {WINNOWMAP} -W {{REF_DIR}}/winnowmap_{{REFNAME}}_repetitive_k15.txt -x map-pb -t 1 --secondary=no --sv-off {{REF}} {{INPUT}}",
    "mashmap1": f"{MASHMAP1} -s {{REF}} -q {{INPUT}} -o {{outdir}}/{{tool}}.paf",
    "mashmap3": f"{MASHMAP3} -t 1 --noSplit --pi 90 -r {{REF}} -q {{INPUT}} -o {{outdir}}/{{tool}}.paf",
    "astarix":  f"{ASTARIX} align-optimal -g {{REF}} -q {{INPUT}} -o {{outdir}}/{{tool}}.paf",
}

tool_list = list(tools.keys())

all_eval_files = expand("{outdir}/{tool}/{tool}.eval", outdir=OUTDIR, tool=tool_list)

rule all:
    input:
        all_eval_files

rule eval_tool:
    input:
    wildcard_constraints:
        tool = "|".join(tool_list)
    output:
        eval = "{outdir}/{tool}/{tool}.eval"
    run:
        outdir = os.path.join(OUTDIR, wildcards.tool)
        shell(f"mkdir -p {outdir}")
        context_index = {
            "REF": REF, "INPUT": ONE_READ, "K": K, "R": R, "T": T,
            "outdir": outdir, "tool": wildcards.tool, "REF_DIR": REF_DIR, "REFNAME": REFNAME
        }
        context_map = context_index.copy()
        context_map["INPUT"] = READS
        tool_cmd = tools[wildcards.tool].format(**context_map)
        shell(f"{TIME_CMD} -o {outdir}/{wildcards.tool}.index.time " + tool_cmd + " > /dev/null 2>&1")
        shell(f"{TIME_CMD} -o {outdir}/{wildcards.tool}.time " + tool_cmd)
        shell(f"{PAFTTOOLS} mapeval -r 0.1 {outdir}/{wildcards.tool}.paf | tee {outdir}/{wildcards.tool}.eval")
        shell(f"{PAFTTOOLS} mapeval -r 0.1 -Q 0 {outdir}/{wildcards.tool}.paf > {outdir}/{wildcards.tool}.wrong")