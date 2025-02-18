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
ALL_READS = "{READS_DIR}/{REFNAME}-reads{READSIM_REFNAME}-a{ACCURACY}-d{DEPTH}-l{MEANLEN}.fa"

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
    "shmap":    f"{SHMAP} -s {{ref}} -p {{reads}} -k {{k}} -r {{r}} -t {{t}} -z {{out_pref}}.params > {{out_pref}}.paf",
    "mm2":      f"{MM2} -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level {{ref}} {{reads}} > {{out_pref}}.paf",
    "blend":    f"{BLEND} -x map-hifi -t 1 -N 0 {{ref}} {{reads}} > {{out_pref}}.paf",
    "mashmap3": f"{MASHMAP3} -t 1 --noSplit --pi 90 -r {{ref}} -q {{reads}} -o {{out_pref}}.paf",
    "mapquik":  f"{MAPQUIK} {{reads}} --reference {{ref}} --threads 1 -p {{out_pref}}",
    "minimap":  f"{MINIMAP} -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level {{ref}} {{reads}} > {{out_pref}}.paf",
#    "winnowmap":f"{MERYL} count k=15 output merylDB {{ref}}; {MERYL} print greater-than distinct=0.9998 merylDB > {{ref_dir}}/winnowmap_{{REFNAME}}_repetitive_k15.txt; {WINNOWMAP} -W {{ref_dir}}/winnowmap_{{REFNAME}}_repetitive_k15.txt -x map-pb -t 1 --secondary=no --sv-off {{ref}} {{reads}} > {{out_pref}}.paf",
#    "mashmap1": f"{MASHMAP1} -s {{ref}} -q {{reads}} -o {{out_pref}}.paf",
#    "astarix":  f"{ASTARIX} align-optimal -g {{ref}} -q {{reads}} -o {{out_pref}}.paf",
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
        out_pref = os.path.join(outdir, wildcards.tool)
        shell(f"mkdir -p {outdir}")
        context_index = {
            "ref": REF, "reads": ONE_READ, "out_pref": out_pref,
            "outdir": outdir, "tool": wildcards.tool, "ref_dir": REF_DIR, "REFNAME": REFNAME,
            "k": K, "r": R, "t": T,
        }
        tool_index_cmd = tools[wildcards.tool].format(**context_index)
        context_map = context_index.copy()
        context_map["reads"] = ALL_READS
        tool_map_cmd = tools[wildcards.tool].format(**context_map)
        shell(f"{TIME_CMD} -o {outdir}/{wildcards.tool}.index.time " + tool_index_cmd + " > /dev/null 2>&1")
        shell(f"{TIME_CMD} -o {outdir}/{wildcards.tool}.time " + tool_map_cmd + " 2> {out_pref}.log")
        shell(f"{PAFTTOOLS} mapeval -r 0.1 {out_pref}.paf | tee {out_pref}.eval")
        shell(f"{PAFTTOOLS} mapeval -r 0.1 -Q 0 {out_pref}.paf > {out_pref}.wrong")