SIM_DATA=../../ngsTools/ngsSim/examples
ANGSD=../../ngsTools/angsd


##### Clean up
rm -f testA_*


N_IND=24
N_SITES=10000


##### TRUE genotypes
cat $SIM_DATA/testA.geno | perl -s -p -e 's/0 0/0/g; s/(\w) \1/2/g; s/\w \w/1/g; $n=s/2/2/g; tr/02/20/ if($n>$n_ind/2)' -- -n_ind=$N_IND | awk '{print "chrSIM\t"NR"\t"$0}' | gzip -cfn --best > testA_T.geno.gz
../ngsDist -n_threads 10 -seed 12345 -verbose 0 -geno testA_T.geno.gz -n_ind $N_IND -n_sites $N_SITES -labels testA.labels -n_boot_rep 5                     -out_prefix testA_T
../ngsDist -n_threads 10 -seed 12345 -verbose 0 -geno testA_T.geno.gz -n_ind $N_IND -n_sites $N_SITES -labels testA.labels -n_boot_rep 5 -boot_block_size 10 -out_prefix testA_T-10


##### Called genotypes
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf -1 -doGeno 2 -out testA_2
../ngsDist -n_threads 10 -seed 12345 -verbose 0 -geno testA_2.geno.gz -n_ind $N_IND -n_sites $N_SITES -labels testA.labels -n_boot_rep 5                     -out_prefix testA_2
../ngsDist -n_threads 10 -seed 12345 -verbose 0 -geno testA_2.geno.gz -n_ind $N_IND -n_sites $N_SITES -labels testA.labels -n_boot_rep 5 -boot_block_size 10 -out_prefix testA_2-10

##### Genotypes' likelihood and posterior probabilities
$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf -1 -doGeno 32 -out testA_32
gunzip testA_32.geno.gz
../ngsDist -n_threads 10 -seed 12345 -verbose 0 -geno testA_32.geno -probs -log_scale -n_ind $N_IND -n_sites $N_SITES -labels testA.labels -n_boot_rep 5                     -out_prefix testA_32
../ngsDist -n_threads 10 -seed 12345 -verbose 0 -geno testA_32.geno -probs -log_scale -n_ind $N_IND -n_sites $N_SITES -labels testA.labels -n_boot_rep 5 -boot_block_size 10 -out_prefix testA_32-10

$ANGSD/angsd -glf $SIM_DATA/testA.glf.gz -fai $SIM_DATA/testAF.ANC.fai -nInd $N_IND -doMajorMinor 1 -doPost 1 -doMaf -1 -doGeno 8 -out testA_8
../ngsDist -n_threads 10 -seed 12345 -verbose 0 -geno testA_8.geno.gz -probs -n_ind $N_IND -n_sites $N_SITES -labels testA.labels -n_boot_rep 5                     -out_prefix testA_8
../ngsDist -n_threads 10 -seed 12345 -verbose 0 -geno testA_8.geno.gz -probs -n_ind $N_IND -n_sites $N_SITES -labels testA.labels -n_boot_rep 5 -boot_block_size 10 -out_prefix testA_8-10


##### Check MD5
rm -f *.arg
md5sum testA_* | sort -k 2,2 > /tmp/test.md5
if diff /tmp/test.md5 test.md5 > /dev/null
then
    echo "ngsDist: All tests OK!"
else
    echo "ngsDist: test(s) failed!"
fi
