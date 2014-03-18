SIM_DATA=../../ngsTools/ngsSim/examples
ANGSD=../../ngsTools/angsd


##### Clean up
rm testA_*


##### Genotypes' likelihood and posterior probabilities
$ANGSD/angsd -sim1 $SIM_DATA/testA.glf.gz -nInd 24 -doMajorMinor 1 -doPost 1 -doMaf 2 -doGeno 32 -out testA_32
gunzip testA_32.geno.gz
../ngsDist -g testA_32.geno -n 24 -s 10000 -o testA_32 -lkl -n_threads 10 -verbose 0 -labels testA.labels

$ANGSD/angsd -sim1 $SIM_DATA/testA.glf.gz -nInd 24 -doMajorMinor 1 -doPost 1 -doMaf 2 -doGeno 8 -out testA_8
../ngsDist -g testA_8.geno.gz -n 24 -s 10000 -o testA_8 -lkl -n_threads 10 -verbose 0 -labels testA.labels



##### Check MD5
rm -f *.mafs.gz *.arg
md5sum testA_* | sort -k 2,2 > /tmp/test.md5
if diff /tmp/test.md5 test.md5 > /dev/null
then
    echo "ngsDist: All tests OK!"
else
    echo "ngsDist: test(s) failed!"
fi
