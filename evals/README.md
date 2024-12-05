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

### Running tools
The tools produce `.paf` files with mappings results. We measure their running time and memory using the `time` command. In order to separate the indexing time from the mapping time, we run the tools twice for each dataset: once to align only a single read (effectively measuring the time and memory of the indexing), and once on all reads (subtracting the indexing time gives us the mapping time).

We run minimap2 with HiFi configuration and no secondary chains:
```bash
minimap2 -x map-hifi -t 1 --secondary=no -M 0 --hard-mask-level $(REF) $(READS) 2> >(tee $(MINIMAP_PREF).log) >$(MINIMAP_PREF).paf
```

We run mapquick with default parameters on a single thread:
```bash
mapquik $(READS) --reference $(REF) --threads 1 -p $(MAPQUIK_PREF)
```

Winnwomap requires the frequent kmers to be found before running the mapping. We did not find a way to run `Winnowmap` on a single thread without changing its source code. Therefore, we accepted that it uses several threds:
```bash
$(MERYL_BIN) count k=15 output merylDB $(REF)
$(MERYL_BIN) print greater-than distinct=0.9998 merylDB > $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt

$(WINNOWMAP_BIN) -W $(REF_DIR)/winnowmap_$(REFNAME)_repetitive_k15.txt -x map-pb -t 1 --secondary=no --sv-off -M 0 --hard-mask-level $(REF) $(READS) >$(WINNOWMAP_PREF).paf 
```