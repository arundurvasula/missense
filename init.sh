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
echo "--> Getting missense sites in the ExAC data"
cd ../../
if [ ! -f data/ExAC.missense.txt ]; then
    awk '/missense_variant/ {print $1"\t"$2}' <(gunzip -c data/ExAC.r0.3.1.sites.vep.vcf.gz) > data/ExAC.missense.txt
fi
if [ ! -f data/ExAC.missense.vcf.gz ]; then
    bcftools view data/ExAC.r0.3.1.sites.vep.vcf.gz -R data/ExAC.missense.txt -Oz > data/ExAC.missense.vcf.gz
    bcftools index data/ExAC.missense.vcf.gz
fi

echo "--> Getting synonymous sites in the ExAC data"
if [ ! -f data/ExAC.syn.txt ]; then
    awk '/synonymous_variant/ {print $1"\t"$2}' <(gunzip -c data/ExAC.r0.3.1.sites.vep.vcf.gz) > data/ExAC.syn.txt
fi
if [ ! -f data/ExAC.synonymous.vcf.gz ]; then
    bcftools view data/ExAC.r0.3.1.sites.vep.vcf.gz -R data/ExAC.syn.txt -Oz > data/ExAC.synonymous.vcf.gz
    bcftools index data/ExAC.synonymous.vcf.gz
fi

echo "--> Subsetting ExAC data for only SNPs in full 230 that are also 1) missense mutations, 2) biallelic, 3) MAF > 0.1"
if [ ! -f data/ExAC.biallelic.missense.maf-0.1.vcf.gz ]; then
    bcftools query -f '%CHROM\t%POS\n' data/MathiesonEtAl_genotypes/full230.vcf.gz > data/MathiesonEtAl_genotypes/full230.sites.txt
    bcftools view data/ExAC.missense.vcf.gz -R data/MathiesonEtAl_genotypes/full230.sites.txt -q 0.1 -m2 -M2 -v snps| bgzip -c > data/ExAC.biallelic.missense.maf-0.1.vcf.gz
fi
echo "--> Subsetting ExAC data for only SNPs in full 230 that are also 1) missense mutations, 2) biallelic, 3) MAF < 0.0001"
if [ ! -f data/ExAC.biallelic.missense.maf-lt-0.1.vcf.gz ]; then
    bcftools view data/ExAC.missense.vcf.gz -R data/MathiesonEtAl_genotypes/full230.sites.txt -Q 0.0001 -m2 -M2 -v snps| bgzip -c > data/ExAC.biallelic.missense.maf-lt-0.1.vcf.gz
fi
echo "--> Subsetting ExAC data for only SNPs in full 230 that are also 1) synonymous mutations, 2) biallelic (don't care about frequency)"
if [ ! -f data/ExAC.biallelic.synonymous.vcf.gz ]; then
    bcftools view data/ExAC.synonymous.vcf.gz -R data/MathiesonEtAl_genotypes/full230.sites.txt -m2 -M2 -v snps| bgzip -c > data/ExAC.biallelic.synonymous.vcf.gz
fi

echo "----"
echo "--> Subsetting Mathieson genotypes for only missense mutations and MAF > 0.1"
if [ ! -f data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-0.1.vcf.gz ]; then
    bcftools query -f '%CHROM\t%POS\n' data/ExAC.biallelic.missense.maf-0.1.vcf.gz > data/ExAC.biallelic.missense.maf-0.1.sites.txt
    bcftools view data/MathiesonEtAl_genotypes/full230.vcf.gz -R data/ExAC.biallelic.missense.maf-0.1.sites.txt | bgzip -c > data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-0.1.vcf.gz
fi
echo "--> Subsetting Mathieson genotypes for only missense mutations and MAF < 0.0001"
if [ ! -f data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-lt-0.1.vcf.gz ]; then
    bcftools query -f '%CHROM\t%POS\n' data/ExAC.biallelic.missense.maf-lt-0.1.vcf.gz > data/ExAC.biallelic.missense.maf-lt-0.1.sites.txt
    bcftools view data/MathiesonEtAl_genotypes/full230.vcf.gz -R data/ExAC.biallelic.missense.maf-lt-0.1.sites.txt | bgzip -c > data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-lt-0.1.vcf.gz
fi
echo "--> Subsetting Mathieson genotypes for only synonymous mutations"
if [ ! -f data/MathiesonEtAl_genotypes/full230.biallelic.synonymous.vcf.gz ]; then
    bcftools query -f '%CHROM\t%POS\n' data/ExAC.biallelic.synonymous.vcf.gz > data/ExAC.biallelic.synonymous.sites.txt
    bcftools view data/MathiesonEtAl_genotypes/full230.vcf.gz -R data/ExAC.biallelic.synonymous.sites.txt | bgzip -c > data/MathiesonEtAl_genotypes/full230.biallelic.synonymous.vcf.gz
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
echo "--> Calculating allele frequencies for 6 ancient populations (MAF < 0.0001)"
if [ ! -f results/AEN-maf-lt-0.1.frq ]; then
    for pop in `ls data/pop_ids/*`;
    do
        p=`basename ${pop}`
        vcftools --gzvcf data/MathiesonEtAl_genotypes/full230.biallelic.missense.maf-lt-0.1.vcf.gz --freq2 --out results/${p}-maf-lt-0.1 --keep ${pop}
    done
fi
echo "--> Calculating allele frequencies for 6 ancient populations (synonymous)"
if [ ! -f results/AEN-syn.frq ]; then
    for pop in `ls data/pop_ids/*`;
    do
        p=`basename ${pop}`
        vcftools --gzvcf data/MathiesonEtAl_genotypes/full230.biallelic.synonymous.vcf.gz --freq2 --out results/${p}-syn --keep ${pop}
    done
fi
echo "--> Calculating Fst between ancient populations."
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


# 7. write ancient allele frequencies in order (AEN, CEM, STP, CLB)
echo "----"
echo "--> Gathering ancient allele frequencies"
if [ ! -f results/ancient.derived-maf-0.1.csv ]; then
    Rscript scripts/freq-table.R
fi

echo "--> Appending present day allele counts"
if [ ! -f results/nchr-maf-0.1.csv ]; then
    paste -d "," results/ancient.nchr-maf-0.1.csv <(bcftools query -f "%INFO/AN_NFE\n" data/ExAC.biallelic.missense.maf-0.1.vcf.gz) > results/nchr-maf-0.1.csv
fi
if [ ! -f results/nchr-maf-lt-0.1.csv ]; then
    paste -d "," results/ancient.nchr-maf-lt-0.1.csv <(bcftools query -f "%INFO/AN_NFE\n" data/ExAC.biallelic.missense.maf-lt-0.1.vcf.gz) > results/nchr-maf-lt-0.1.csv
fi
if [ ! -f results/nchr-syn.csv ]; then
    paste -d "," results/ancient.nchr-syn.csv <(bcftools query -f "%INFO/AN_NFE\n" data/ExAC.biallelic.synonymous.vcf.gz) > results/nchr-syn.csv
fi

echo "--> Appending present day sample size"
## USE AC_EUR and AN_EUR for nchr and derived. Maybe also remove frequencies if I don't really need them (the stuff above)
if [ ! -f results/derived-maf-0.1.csv ]; then
    paste -d "," results/ancient.derived-maf-0.1.csv <(bcftools query -f "%INFO/AC_NFE\n" data/ExAC.biallelic.missense.maf-0.1.vcf.gz) > results/derived-maf-0.1.csv
fi
if [ ! -f results/derived-maf-lt-0.1.csv ]; then
    paste -d "," results/ancient.derived-maf-lt-0.1.csv <(bcftools query -f "%INFO/AC_NFE\n" data/ExAC.biallelic.missense.maf-lt-0.1.vcf.gz) > results/derived-maf-lt-0.1.csv
fi
if [ ! -f results/derived-syn.csv ]; then
    paste -d "," results/ancient.derived-syn.csv <(bcftools query -f "%INFO/AC_NFE\n" data/ExAC.biallelic.synonymous.vcf.gz) > results/derived-syn.csv
fi


# 8. use the allele frequencies to figure out selection coefficients and allele ages
mkdir -p results/selection
if [ ! -f data/maf-0.1.sites.txt ]; then
    tail -n +2 results/AEN-maf-0.1.frq | awk '{print $1"-"$2}' > data/maf-0.1.sites.txt
fi
if [ ! -f data/maf-lt-0.1.sites.txt ]; then
    tail -n +2 results/AEN-maf-lt-0.1.frq | awk '{print $1"-"$2}' > data/maf-lt-0.1.sites.txt
fi
if [ ! -f data/syn.sites.txt ]; then
    tail -n +2 results/AEN-syn.frq | awk '{print $1"-"$2}' > data/syn.sites.txt
fi
echo "--> Estimating allele frequency trajectories for MAF > 0.1"
# time is based on Ne = 3500, y=8400,7700,5400,4900,0. g=y/2N
if [ ! -f results/selection/maf-0.1-1-914852.param ]; then
    while read der <&3 && read nchr <&4 && read sites <&5; do
	if [[ ${der} == "0,0,0,0,0" ]]; then
            continue
        fi
       	software/selection/sr -X ${der} -N ${nchr} -T -1.2,-1.1,-0.77,-0.7,0 -n 100000 -f 2000 -s 100 -P software/selection/constant.pop -a -o results/selection/maf-0.1-${sites}
	rm results/selection/maf-0.1-${sites}.time
	rm results/selection/maf-0.1-${sites}.traj
    done 3< results/derived-maf-0.1.csv 4< results/nchr-maf-0.1.csv 5< data/maf-0.1.sites.txt
fi
echo "--> Estimating allele frequency trajectories for MAF < 0.0001"
if [ ! -f results/selection/maf-lt-0.1-1-2489248.param ]; then
    while read der <&3 && read nchr <&4 && read sites <&5; do
        if [[ ${der} == "0,0,0,0,0" ]]; then
	    continue
	fi
	software/selection/sr -X ${der} -N ${nchr} -T -1.2,-1.1,-0.77,-0.7,0 -n 100000 -f 2000 -s 100 -P software/selection/constant.pop -a -o results/selection/maf-lt-0.1-${sites}
	rm results/selection/maf-lt-0.1-${sites}.time
	rm results/selection/maf-lt-0.1-${sites}.traj
    done 3< results/derived-maf-lt-0.1.csv 4< results/nchr-maf-lt-0.1.csv 5< data/maf-lt-0.1.sites.txt
fi
echo "--> Estimating allele frequency trajectories for synonymous mutations"
if [ ! -f results/selection/syn-1-914852.param ]; then
    while read der <&3 && read nchr <&4 && read sites <&5; do
        if [[ ${der} == "0,0,0,0,0" ]]; then
            continue
        fi
	software/selection/sr -X ${der} -N ${nchr} -T -1.2,-1.1,-0.77,-0.7,0 -n 100000 -f 2000 -s 100 -P software/selection/constant.pop -a -o results/selection/syn-${sites}
	rm results/selection/syn-${sites}.time
	rm results/selection/syn-${sites}.traj
    done 3< results/derived-syn.csv 4< results/nchr-syn.csv 5< data/syn.sites.txt
fi

echo "--> Plotting allele ages and selection coefficients."
Rscript scripts/plot-s.R