set -e

#1. make dirs
mkdir -p data scripts results software
mkdir -p scripts/data-munge

#2. download some conversion scripts
cd scripts/data-munge
if [ ! -f eigenstrat2vcf.py ]; then
    wget --no-check-certificate https://raw.githubusercontent.com/mathii/gdc/master/eigenstrat2vcf.py
fi
if [ ! -f pyEigenstrat.py ]; then
    wget --no-check-certificate https://raw.githubusercontent.com/mathii/pyEigenstrat/master/pyEigenstrat.py
fi 
if [ ! -f gdc.py ]; then
    wget --no-check-certificate https://raw.githubusercontent.com/mathii/gdc/master/gdc.py
fi

#3. get some data
cd ../../data/
if [ ! -f MathiesonEtAl_genotypes_April2016.tar.gz ]; then
    wget --no-check-certificate http://genetics.med.harvard.edu/reich/Reich_Lab/Datasets_files/MathiesonEtAl_genotypes_April2016.tar.gz
    tar -xvzf MathiesonEtAl_genotypes_April2016.tar.gz
fi

if [ ! -f ExAC.r0.3.1.sites.vep.vcf.gz ]; then
    wget --no-check-certificate ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz
    wget --no-check-certificate ftp://ftp.broadinstitute.org/pub/ExAC_release/release0.3.1/ExAC.r0.3.1.sites.vep.vcf.gz.tbi
fi

mkdir -p pop_ids
if [ ! -f pop_ids/AEN ]; then
    cd pop_ids
    wget --no-check-certificate https://gist.github.com/arundurvasula/62a1cd1884d0f3c8931e60a81546939c/raw/bb8c2be022ccd00e2b984bdcd9e0a96730c27cd9/AEN
    wget --no-check-certificate https://gist.github.com/arundurvasula/62a1cd1884d0f3c8931e60a81546939c/raw/bb8c2be022ccd00e2b984bdcd9e0a96730c27cd9/CEM
    wget --no-check-certificate https://gist.github.com/arundurvasula/62a1cd1884d0f3c8931e60a81546939c/raw/bb8c2be022ccd00e2b984bdcd9e0a96730c27cd9/CLB
    wget --no-check-certificate https://gist.github.com/arundurvasula/62a1cd1884d0f3c8931e60a81546939c/raw/bb8c2be022ccd00e2b984bdcd9e0a96730c27cd9/HG
    wget --no-check-certificate https://gist.github.com/arundurvasula/62a1cd1884d0f3c8931e60a81546939c/raw/bb8c2be022ccd00e2b984bdcd9e0a96730c27cd9/INC
    wget --no-check-certificate https://gist.github.com/arundurvasula/62a1cd1884d0f3c8931e60a81546939c/raw/bb8c2be022ccd00e2b984bdcd9e0a96730c27cd9/STP
    cd ..
fi

#4. convert data
cd MathiesonEtAl_genotypes/
if [ ! -f full230.csv ]; then
    wget --no-check-certificate https://gist.github.com/arundurvasula/62a1cd1884d0f3c8931e60a81546939c/raw/bb8c2be022ccd00e2b984bdcd9e0a96730c27cd9/full230.csv
fi

if [ ! -f full230.vcf.gz ]; then
    /usr/bin/python ../../scripts/data-munge/eigenstrat2vcf.py -r ./full230 | bgzip -c > ./full230.vcf.gz
    bcftools index ./full230.vcf.gz
fi

# 4a. get some software
cd ../../software
if [ ! -d selection ]; then
    git clone https://github.com/Schraiber/selection.git
fi

cd selection
if [ ! -f sr ]; then 
    g++ -O3 *.cpp -lgsl -lgslcblas -lm -o sr
fi

#5. done with prep

echo "----"
echo "--> VCF for ancient genomes created in data/MathiesonEtAl_genotypes/full230.vcf.gz"
echo "--> ExAC data downloaded into data/ExAC.r0.3.1.sites.vep.vcf.gz"

#6. prepare for some analyses

echo "----"
echo "--> Subsetting ExAC data for only SNPs in full 230 that are also 1) missense mutations, 2) biallelic, 3) MAF > 0.1"

cd ../../
if [ ! -f data/ExAC.biallelic.missense.maf-0.1.vcf.gz ]; then
    bcftools query -f '%CHROM\t%POS\n' data/MathiesonEtAl_genotypes/full230.vcf.gz > data/MathiesonEtAl_genotypes/full230.sites.txt
    bcftools view data/ExAC.r0.3.1.sites.vep.vcf.gz -R data/MathiesonEtAl_genotypes/full230.sites.txt -q 0.1 -i INFO/CSQ "~" "*missense_variant*" -m2 -M2 -v snps| bgzip -c > data/ExAC.biallelic.missense.maf-0.1.vcf.gz
fi
echo "--> Subsetting ExAC data for only SNPs in full 230 that are also 1) missense mutations, 2) biallelic, 3) MAF < 0.1"

if [ ! -f data/ExAC.biallelic.missense.maf-lt-0.1.vcf.gz ]; then
    bcftools view data/ExAC.r0.3.1.sites.vep.vcf.gz -R data/MathiesonEtAl_genotypes/full230.sites.txt -Q 0.1 -i INFO/CSQ "~" "*missense_variant*" -m2 -M2 -v snps| bgzip -c > data/ExAC.biallelic.missense.maf-lt-0.1.vcf.gz
fi

echo "----"
echo "--> Subsetting Mathieson genotypes for only missense mutations and MAF > 0.1"
if [ ! -f data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-0.1.vcf.gz ]; then
    bcftools query -f '%CHROM\t%POS\n' data/ExAC.biallelic.missense.maf-0.1.vcf.gz > data/ExAC.biallelic.missense.maf-0.1.sites.txt
    bcftools view data/MathiesonEtAl_genotypes/full230.vcf.gz -R data/ExAC.biallelic.missense.maf-0.1.sites.txt | bgzip -c > data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-0.1.vcf.gz
fi

echo "--> Subsetting Mathieson genotypes for only missense mutations and MAF < 0.1"
if [ ! -f data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-lt-0.1.vcf.gz ]; then
    bcftools query -f '%CHROM\t%POS\n' data/ExAC.biallelic.missense.maf-lt-0.1.vcf.gz > data/ExAC.biallelic.missense.maf-lt-0.1.sites.txt
    bcftools view data/MathiesonEtAl_genotypes/full230.vcf.gz -R data/ExAC.biallelic.missense.maf-lt-0.1.sites.txt | bgzip -c > data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-lt-0.1.vcf.gz
fi

echo "----"
echo "--> Calculating allele frequencies for 6 ancient populations (MAF > 0.1)"
if [ ! -f results/AEN-maf-0.1.frq ]; then
    for pop in `ls data/pop_ids/*`;
    do
        p=`basename ${pop}`
        vcftools --gzvcf data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-0.1.vcf.gz --freq2 --out results/${p}-maf-0.1 --keep ${pop}
    done 
fi

echo "--> Calculating allele frequencies for 6 ancient populations (MAF < 0.1)"
if [ ! -f results/AEN-maf-lt-0.1.frq ]; then
    for pop in `ls data/pop_ids/*`;
    do
        p=`basename ${pop}`
        vcftools --gzvcf data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-lt-0.1.vcf.gz --freq2 --out results/${p}-maf-lt-0.1 --keep ${pop}
    done
fi

echo "--> Calculating Fst between ancient populationss."
if [ ! -f results/AEN.CLB.weir.fst ]; then
    for pop1 in `ls data/pop_ids/*`;
    do
	for pop2 in `ls data/pop_ids/*`;
	do
            p1=`basename ${pop1}`
            p2=`basename ${pop2}`
	    vcftools --gzvcf data/MathiesonEtAl_genotypes/full230.vcf.gz --weir-fst-pop ${pop1} --weir-fst-pop ${pop2} --out results/${p1}.${p2}
	done
    done
fi

echo "----"
echo "--> Done. Check out scripts/missense.Rmd to continue analysis"
