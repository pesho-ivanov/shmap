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
    "depth": 0.2,
    "meanlen": "?",
    
    # Mapping parameters
    "k": 25,
    "r": 0.01,
    "t": 0.4,
    
    # Tools and paths
    "time_cmd": "/usr/bin/time -f '%U\t%M'",
    "paftools": "evals/ext/paftools.js",
    
    # Tool executables
    "tool_bin": {
        "shmap": "./release/shmap",
        "minimap": "~/libs/minimap2/minimap2",
        "mm2": "~/libs/mm2-fast/minimap2",
        "blend": "~/libs/blend/bin/blend",
        "mapquik": "~/libs/mapquik/target/release/mapquik",
        "meryl": "~/libs/Winnowmap/bin/meryl",
        "winnowmap": "~/libs/Winnowmap/bin/winnowmap", 
        "mashmap3": "~/libs/MashMap/build/bin/mashmap",
    },

    # Add these to the config section
    "pbsim3": "~/libs/pbsim3/src/pbsim",
#    "seqkit": "~/miniconda3/bin/seqkit",
    "stream_lift_fasta": "evals/stream_fasta.py",
    "chain_file": "evals/refs/hg002v1.1_to_CHM13v2.0.chain",
#    "seqtk": "~/miniconda3/bin/seqtk",
    "hg002": {
        "fq": "evals/reads/HG002.fq",
        "fa": "evals/reads/HG002.fa",
        "url": "https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifirevio/m84005_220827_014912_s1.hifi_reads.fastq.gz"
    }
}

# Tool command templates
tool_pattern = {
    "shmap": "{tool_bin} -s {ref} -p {reads} -k {k} -r {r} -t {t} -z {out_pref}.params > {out_pref}.paf",
    "mm2": "{tool_bin} -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level {ref} {reads} > {out_pref}.paf",
    "blend": "{tool_bin} -x map-hifi -t 1 -N 0 {ref} {reads} > {out_pref}.paf",
    "mashmap3": "{tool_bin} -t 1 --noSplit --pi 90 -r {ref} -q {reads} -o {out_pref}.paf",
    "mapquik": "{tool_bin} {reads} --reference {ref} --threads 1 -p {out_pref}",
    "minimap": "{tool_bin} -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level {ref} {reads} > {out_pref}.paf",
}

run_tools = [
    ("shmap", "k25-r0.01-t0.4-d0.075-o0.3-mJaccard"),
#    ("mm2", "default"),
    ("blend", "default"),
#    ("mashmap3", "default"),
#    ("mapquik", "default"),
#    ("minimap", "default"),
]

# Rules
rule all:
    input:
        expand("{outdir}/{tool_params[0]}/{tool_params[1]}/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}/{tool_params[0]}.eval", 
               outdir=config["outdir"], 
#               reads_dir=config["reads_dir"],
               refname=config["refname"],
               readsim_refname=config["readsim_refname"],
               accuracy=config["accuracy"],
               depth=config["depth"],
               meanlen=config["meanlen"],
               tool_params=run_tools),

rule eval_tool:
    input:
        reads = multiext(config["reads_dir"] + "/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}", ".oneread.fa", ".fa")
    output:
        eval = "{outdir}/{tool}/{params}/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}/{tool}.eval"
    #wildcard_constraints:
    #    tool = '|'.join([k.split("-")[0] for k in tool_pattern.keys()])
    #benchmark:
    #    "benchmarks/{sample}.bwa.benchmark.txt"
    run:
        outdir = os.path.dirname(output.eval)
        out_pref = os.path.splitext(output.eval)[0]
        shell(f"mkdir -p {outdir}")

        # Prepare context for command formatting
        context = {
            "tool": wildcards.tool,
            "tool_bin": config["tool_bin"][wildcards.tool],
            "ref": f"{config['ref_dir']}/{config['refname']}.fa",
            "out_pref": out_pref,
        }
        
        # Tool parameters, e.g. "k25-r0.01-t0.4-d0.075-o0.3-mJaccard"
        for part in wildcards.params.split("-"):
            param = part[0]
            value = part[1:]
            context[param] = value

        # Index phase
        context["reads"] = input.reads[0]
        print('context: ', context)
        index_cmd = tool_pattern[context['tool']].format(**context)
        shell(f"{config['time_cmd']} -o {outdir}/{context['tool']}.index.time {index_cmd} > /dev/null 2>&1")

        # Mapping phase
        context["reads"] = input.reads[1]
        map_cmd = tool_pattern[context['tool']].format(**context)
        shell(f"{config['time_cmd']} -o {outdir}/{context['tool']}.time {map_cmd} 2> {out_pref}.log")

        # Evaluation
        shell(f"{config['paftools']} mapeval -r 0.1 {out_pref}.paf | tee {out_pref}.eval")
        shell(f"{config['paftools']} mapeval -r 0.1 -Q 0 {out_pref}.paf > {out_pref}.wrong")

# Add this rule after the existing rules
rule gen_reads:
    input:
        ref = config["ref_dir"] + "/{readsim_refname}.fa",
        real_reads = config["hg002"]["fq"]
    output:
        reads = "evals/reads/{readsim_refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}.fa",
        one_read = "evals/reads/{readsim_refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}.oneread.fa"
    params:
        reads_prefix = lambda w: f"{w.readsim_refname}-reads{w.readsim_refname}-a{w.accuracy}-d{w.depth}-l{w.meanlen}",
        real_reads_tmp = "evals/reads/HG002.fq.tmp"
    conda:
        "envs/simulate.yaml"
    shell:
        """
        # Clean up reference headers
        sed -i 's/^\\(>[^[:space:]]*\\).*/\\1/' {input.ref}
        
        # Prepare subset of real reads
        head -n 400000 {input.real_reads} > {params.real_reads_tmp}
        
        # Run PBSIM3
        {config[pbsim3]} --strategy wgs \
            --method sample \
            --sample {params.real_reads_tmp} \
            --genome {input.ref} \
            --depth {wildcards.depth} \
            --prefix {params.reads_prefix} \
            --no-fastq 1

        # Index reference
        seqkit faidx {input.ref}
        
        # Convert PBSIM output to FASTQ
        {config[paftools]} pbsim2fq {input.ref}.fai "{params.reads_prefix}"_*.maf > {output.reads}.unshuf
        
        # Shuffle reads
        seqkit shuffle -2 {output.reads}.unshuf -o {output.reads}
        
        # Clean up intermediate files
        rm -f "{params.reads_prefix}"_*.maf "{params.reads_prefix}"_*.ref "{params.reads_prefix}"_*.fastq {params.real_reads_tmp}
        
        # Perform liftover if needed
        if [ "{wildcards.readsim_refname}" != "{config[refname]}" ]; then
            echo "Lifting over reads from {wildcards.readsim_refname} to {config[refname]}"
            python {config[stream_lift_fasta]} {output.reads} {config[chain_file]} > {output.reads}.lifted
            mv {output.reads}.lifted {output.reads}
        fi
        
        # Create one read file
        head -n 2 {output.reads} > {output.one_read}
        """

# Add this rule after existing rules
rule get_hg002_reads:
    output:
        fq = config["hg002"]["fq"],
        fa = config["hg002"]["fa"]
    conda:
        "envs/download.yaml"
    message:
        "Downloading and preparing HG002 reads (~84GB). This may take a while..."
    shell:
        """
        # Create output directory if it doesn't exist
        mkdir -p $(dirname {output.fq})
        
        # Download and decompress with proper pipe handling
        set -o pipefail  # Make pipe failures propagate
        curl -f -s {config[hg002][url]} | gunzip | head -n 40000 > {output.fq} || true
        
        # Convert fastq to fasta with all capital letters (some tools require all caps)
        seqtk seq -AU {output.fq} > {output.fa}
        """