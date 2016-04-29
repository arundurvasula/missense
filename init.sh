set -e

#1. make dirs
mkdir -p data scripts results software
mkdir -p scripts/data-munge

#2. download some conversion scripts
cd scripts/data-munge

if [ ! -f eigenstrat2vcf.py ]; then
    wget https://raw.githubusercontent.com/mathii/gdc/master/eigenstrat2vcf.py
fi
if [ ! -f pyEigenstrat.py ]; then
    wget https://raw.githubusercontent.com/mathii/pyEigenstrat/master/pyEigenstrat.py
fi 
if [ ! -f gdc.py ]; then
    wget https://raw.githubusercontent.com/mathii/gdc/master/gdc.py
fi

#3. get some data
cd ../../data/
if [ ! -f MathiesonEtAl_genotypes_April2016.tar.gz ]; then
    wget http://genetics.med.harvard.edu/reich/Reich_Lab/Datasets_files/MathiesonEtAl_genotypes_April2016.tar.gz
    tar -xvzf MathiesonEtAl_genotypes_April2016.tar.gz
fi

if [ ! -f ExAC.r0.3.1.sites.vep.vcf.gz ]; then
    wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz
    wget ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz.tbi
fi

#4. convert data
cd MathiesonEtAl_genotypes/
if [ ! -f full230.csv ]; then
    wget https://gist.githubusercontent.com/arundurvasula/62a1cd1884d0f3c8931e60a81546939c/raw/e4372e99741dab8f82e281273e74cadd77cb9b41/full230.csv
fi
if [ ! -f full230.vcf.gz ]; then
    /usr/bin/python ../../scripts/data-munge/eigenstrat2vcf.py -r ./full230 | bgzip -c > ./full230.vcf.gz
fi

# 4a. get some software
cd ../../software
git clone https://github.com/Schraiber/selection.git
cd selection
g++ -O3 *.cpp -lgsl -lgslcblas -lm -o sr


#5. done

echo "vcf for ancient genomes created in data/MathiesonEtAl_genotypes/full230.vcf.gz"
echo "ExAC data downloaded into data/ExAC.r0.3.1.sites.vep.vcf.gz"

#6. prepare for some analyses

echo "----"
echo "Subsetting ExAC data for only SNPs in full 230 that are also 1) missense mutations, 2) biallelic, 3) MAF > 0.1"

cd ../../
bcftools query -f '%CHROM\t%POS\n' data/MathiesonEtAl_genotypes/full230.vcf.gz > data/MathiesonEtAl_genotypes/full230.sites.txt
bcftools view data/ExAC.r0.3.1.sites.vep.vcf.gz -R data/MathiesonEtAl_genotypes/full230.sites.txt -q 0.1 -i INFO/CSQ "~" "*missense_variant*" -m2 -M2 -v snps| bgzip -c > data/ExAC.biallelic.missense.maf-0.1.vcf.gz
