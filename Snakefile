# Snakefile
import os
from collections import OrderedDict

# Configuration and constants
config = {
    # Directories
    "outdir": "evals/snake_out",
    "ref_dir": "evals/refs",
    "reads_dir": "evals/reads",
    
    # Reference settings
    "refname": "chm13-chr1",

    # Read simulation parameters
    "readsim_refname": "chm13-chr1",
    "accuracy": "?",
    "depth": 1,
    "meanlen": "?",
    
    # Mapping parameters
    "k": 25,
    "r": 0.01,
    "t": 0.4,
    
    # Tools and paths
    "time_cmd": "/usr/bin/time -f '%U\t%M'",
    "paftools": "evals/ext/paftools.js",
    
    # Tool executables
    "tools": {
        "shmap": "./release/shmap",
        "minimap": "~/libs/minimap2/minimap2",
        "mm2": "~/libs/mm2-fast/minimap2",
        "blend": "~/libs/blend/bin/blend",
        "mapquik": "~/libs/mapquik/target/release/mapquik",
        "meryl": "~/libs/Winnowmap/bin/meryl",
        "winnowmap": "~/libs/Winnowmap/bin/winnowmap", 
        "mashmap3": "~/libs/MashMap/build/bin/mashmap",
    }
}

# Tool command templates
tool_commands = {
    "shmap": "{tool} -s {ref} -p {reads} -k {k} -r {r} -t {t} -z {out_pref}.params > {out_pref}.paf",
    "mm2": "{tool} -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level {ref} {reads} > {out_pref}.paf",
    "blend": "{tool} -x map-hifi -t 1 -N 0 {ref} {reads} > {out_pref}.paf",
    "mashmap3": "{tool} -t 1 --noSplit --pi 90 -r {ref} -q {reads} -o {out_pref}.paf",
    "mapquik": "{tool} {reads} --reference {ref} --threads 1 -p {out_pref}",
    "minimap": "{tool} -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level {ref} {reads} > {out_pref}.paf",
}

run_tools = ["shmap"] #, "mm2", "blend", "mashmap3", "mapquik", "minimap"]

# Rules
rule all:
    input:
        expand("{outdir}/{tool}/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}/{tool}.eval", 
               outdir=config["outdir"], 
#               reads_dir=config["reads_dir"],
               refname=config["refname"],
               readsim_refname=config["readsim_refname"],
               accuracy=config["accuracy"],
               depth=config["depth"],
               meanlen=config["meanlen"],
               tool=run_tools)

rule eval_tool:
    input:
        reads = multiext(config["reads_dir"] + "/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}", ".oneread.fa", ".fa")
    output:
        eval = "{outdir}/{tool}/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}/{tool}.eval"
    wildcard_constraints:
        tool = '|'.join(tool_commands.keys())
    run:
        outdir = os.path.dirname(output.eval)
        out_pref = os.path.splitext(output.eval)[0]
        shell(f"mkdir -p {outdir}")

        # Prepare context for command formatting
        context = {
            "tool": config["tools"][wildcards.tool],
            "ref": f"{config['ref_dir']}/{config['refname']}.fa",
            "out_pref": out_pref,
            "k": config["k"],
            "r": config["r"],
            "t": config["t"],
        }

        # Index phase
        context["reads"] = input.reads[0]
        index_cmd = tool_commands[wildcards.tool].format(**context)
        shell(f"{config['time_cmd']} -o {outdir}/{wildcards.tool}.index.time {index_cmd} > /dev/null 2>&1")

        # Mapping phase
        context["reads"] = input.reads[1]
        map_cmd = tool_commands[wildcards.tool].format(**context)
        shell(f"{config['time_cmd']} -o {outdir}/{wildcards.tool}.time {map_cmd} 2> {out_pref}.log")

        # Evaluation
        shell(f"{config['paftools']} mapeval -r 0.1 {out_pref}.paf | tee {out_pref}.eval")
        shell(f"{config['paftools']} mapeval -r 0.1 -Q 0 {out_pref}.paf > {out_pref}.wrong")