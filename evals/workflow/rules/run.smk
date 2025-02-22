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

best_shmap_params = "k25-r0.01-t0.4-d0.075-o0.3-mJaccard"
robustness_config = {
    "variations": { 
        "k": [15, 17, 19, 21, 23, 25, 27, 29, 31, 33, 35],
#        "r": [0.005, 0.01, 0.015, 0.02],
#        "t": [0.2, 0.3, 0.4, 0.5],
#        "d": [0.05, 0.075, 0.1, 0.15],
#        "o": [0.0, 0.1, 0.2, 0.3, 0.4, 0.5, 0.6, 0.7, 0.8, 0.9, 1.0],
#        "m": ["Jaccard", "Hamming", "Levenshtein"],
    },
}

def substitute(params, var):
    pattern = var + r"([^-]+)"
    return [re.sub(pattern, var + str(val), params) for val in robustness_config["variations"][var]]

# For a given {var} (e.g. k), run shmap for each value (e.g. 15, 17, etc.), and aggregate the robustness.
rule robustness_eval:
    input:
        lambda wildcards: expand("{outdir}-robustness-{var}/{tool}/{dataset}/{params}/{tool}.eval",
           outdir=config["outdir"],
           refname=config["refname"],
           readsim_refname=config["readsim_refname"],
           accuracy=config["accuracy"],
           depth=config["depth"],
           meanlen=config["meanlen"],
           tool=wildcards.tool,
           dataset=wildcards.dataset,
           params=substitute(wildcards.params, var=wildcards.var),
           var=wildcards.var)
    output:
        "{outdir}-robustness-{var}/{tool}/{dataset}/{params}/summary.eval"
    shell:
        "cat {input} > {output}"

# Execute a tool for a given set of parameters.
rule run_tool:
    input:
        reads = multiext(config["reads_dir"] + "/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}", ".oneread.fa", ".fa")
    output:
        paf = "{outdir}/{tool}/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}/{params}/{tool}.paf"
    #wildcard_constraints:
    #    tool = '|'.join([k.split("-")[0] for k in tool_pattern.keys()])
    run:
        outdir = os.path.dirname(output.paf)
        out_pref = os.path.splitext(output.paf)[0]
        shell(f"mkdir -p {outdir}")

        # Prepare context for command formatting
        context = {
            "tool": wildcards.tool,
            "tool_bin": config["tool_bin"][wildcards.tool],
            "ref": f"{config['ref_dir']}/{config['refname']}.fa",
            "out_pref": out_pref,
        }
        
        # Tool parameters, e.g. "k25-r0.01-t0.4-d0.075-o0.3-mJaccard"; works only for single-letter parameters
        [context.update({part[0]: part[1:]}) for part in wildcards.params.split("-")]

        # Index phase
        context["reads"] = input.reads[0]
        #print('context: ', context)
        index_cmd = tool_pattern[context['tool']].format(**context)
        shell(f"{config['time_cmd']} -o {outdir}/{context['tool']}.index.time {index_cmd} > /dev/null 2>&1")

        # Mapping phase
        context["reads"] = input.reads[1]
        map_cmd = tool_pattern[context['tool']].format(**context)
        shell(f"{config['time_cmd']} -o {outdir}/{context['tool']}.time {map_cmd} 2> {out_pref}.log")

# Evaluate the mapping results of a tool execution.
rule mapeval:
    input:
        paf = "{prefix}.paf"
    output:
        eval = "{prefix}.eval",
        wrong = "{prefix}.wrong"
    shell:
        """
        {config[paftools]} mapeval -r 0.1 {input.paf} | tee {output.eval}
        {config[paftools]} mapeval -r 0.1 -Q 0 {input.paf} > {output.wrong}
        """