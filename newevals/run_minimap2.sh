echo "minimap2 -------------"

ref=t2tChrY.fa
reads=out/reads-ChrY.fa

# minimap2
time minimap2 -x map-hifi -t 1 $ref $reads > out/minimap2.paf
paftools.js mapeval minimap2.paf >minimap2.evals

