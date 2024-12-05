This is a description of the evaluation of Map-SHmap, needed to reproduce the results in our paper.

### Generate simulated reads using PBSIM3
For our simulated experiments, we use PBSIM3 to generated reads from the CHM13 reference: the full reference for S2, and only from the Y chromosome for S3.
```bash
pbsim --strategy wgs \
	--method sample \
	--sample hg002-hifi.fq \
	--depth 30 \
	--genome chm13.fa  # or chm13_Y.fa
```
### Datasets
The R1 dataset was downloaded from
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/scratch/HG002/sequencing/hifirevio/m84005_220827_014912_s1.hifi_reads.fastq.gz
and more information about this read set can be found in 
https://github.com/marbl/HG002/blob/main/Sequencing_data.md

### Reference
Details about the CHM13 reference can be found at 
https://github.com/marbl/CHM13
and it can be downloaded using the link
https://s3-us-west-2.amazonaws.com/human-pangenomics/T2T/CHM13/assemblies/analysis_set/chm13v2.0.fa.gz.

### Running tools
The tools produce `.paf` files with mappings results. We measure their running time and memory using the `time` command. In order to separate the indexing time from the mapping time, we run the tools twice for each dataset: once to align only a single read (effectively measuring the time and memory of the indexing), and once on all reads (subtracting the indexing time gives us the mapping time).

We run minimap2 with HiFi configuration and no secondary chains:

Winnwomap requires the frequent kmers to be found before running the mapping. We did not find a way to run `Winnowmap` on a single thread without changing its source code. Therefore, we accepted that it uses several threds:
```bash
$(MERYL_BIN) count k=15 output merylDB $(REF)
$(MERYL_BIN) print greater-than distinct=0.9998 merylDB > $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt
```

# Tool Versions and Commands

This document provides the commands used to run the mappers and simulate reads.

## Table of Tools and Parameters

| **Tool**       | **Version**     | **Arguments**                                                                                       | **Sketch Parameters**       |
|-----------------|-----------------|-----------------------------------------------------------------------------------------------------|-----------------------------|
| `minimap2`     | 2.26-r1175      | `-x map-hifi -t 1 --secondary=no --hard-mask-level [ref.fa] [reads.fa]`                                              | `k=19, skip=19`            |
| `winnowmap`    | 2.03            | `-W rep_[ref]_k15.txt -x map-pb -t 1 --secondary=no --sv-off -M 0 --hard-mask-level [ref.fa] [reads.fa]`                                 | `k=15, w=50`               |
| `blend`        | 1.0             | `-x map-hifi -t 1 -N 0 [ref.fa] [reads.fa]`                                                       | `k=19, skip=50`            |
| `mapquik`      | 0.1.0           | `[reads.fa] --reference [ref.fa] --threads 1`                                                     | `k=5, l=31`                |
| `shmap`        | 0.1             | `-s [ref.fa] -p [reads.fa] -k 25 -r 0.05 -t 0.9`                                                  | `k=25, r=0.05`             |
| `pbsim`        |                 | `--strategy wgs --method sample --sample hg002-hifi.fq --depth 30`                                |                             |
|                 |                 | `--genome chm13.fa [or chm13_Y.fa]`                                                              |                             |

### Additional Notes
- `rep_[ref]_k15.txt` is a list of frequent kmers generated using the following commands:
  ```bash
  meryl count k=15 output merylDB [ref.fa]
  meryl print greater-than distinct=0.9998 merylDB
