# Generate reads for a given reference.
rule gen_reads:
    input:
        ref = config["ref_dir"] + "/{readsim_refname}.fa",
        real_reads = config["hg002"]["fq"]
    output:
        reads = "reads/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}.fa",
        one_read = "reads/{refname}-reads{readsim_refname}-a{accuracy}-d{depth}-l{meanlen}.oneread.fa"
    params:
        reads_prefix = lambda w: f"{w.readsim_refname}-reads{w.readsim_refname}-a{w.accuracy}-d{w.depth}-l{w.meanlen}",
        real_reads_tmp = "reads/HG002.fq.tmp"
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
        if [ "{wildcards.readsim_refname}" != "{wildcards.refname}" ]; then
            echo "Lifting over reads from {wildcards.readsim_refname} to {wildcards.refname}"
            python {config[stream_lift_fasta]} {output.reads} {config[chain_file]} > {output.reads}.lifted
            mv {output.reads}.lifted {output.reads}
        fi
        
        # Create one read file
        head -n 2 {output.reads} > {output.one_read}
        """

# Get HG002 reads.
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