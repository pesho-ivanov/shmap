echo "------------- mapquik -------------"

ref=t2tChrY.fa
reads=out/reads-ChrY.fa

time mapquik $reads --reference $ref -p mapquik --threads 1 > out/mapquik.out
paftools.js mapeval mapquik.paf
#tail -4 mapquik.out
