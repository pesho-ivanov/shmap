# pbsim (version 1) always adds some errors, so it can't simulate exactly perfect reads, but it can simulate with variable error rates
# note: pbsim (version 1) only has a 1% precision for error rate (can't simulate 0.5% error rate for example)

#ref=chm13.genome.fa
#reads=nearperfect-chm13.10X.24kb
ref=t2tChrY.fa
reads=out/reads-ChrY

pbsim \
       $ref \
       --model_qc  model_qc_clr \
       --accuracy-mean 0.85\
       --accuracy-sd 0\
       --depth 1\
       --prefix $reads\
       --length-mean 10000 #hifi

#mason_simulator -ir $ref -n 100 -o reads_tmp.fa -oa mason_ali.sam

samtools faidx $ref
#paftools.js mason2fq mason_ali.sam > reads_paf.fa
paftools.js pbsim2fq $ref.fai "$reads"_*.maf > $reads.fa
rm -f "$reads"_*.maf "$reads"_*.ref "$reads"_*.fastq
