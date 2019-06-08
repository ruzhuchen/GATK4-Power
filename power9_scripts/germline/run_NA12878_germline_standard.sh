#!/bin/bash

################################
# Power9 tests - WGS pipeline for Power9 system
################################

#Change this path for your test
export GATK_HOME=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG
export BMK_HOME=
export REF_HOME=

export PATH=$GATK_HOME/bin:$PATH
export GATK_LOCAL_JAR=$GATK_HOME/gatk-4.1.0.0/libs/gatk.jar
export GATK_SPARK_JAR=$GATK_HOME/gatk-4.1.0.0/libs/gatk-spark.jar
export LD_LIBRARY_PATH=$GATK_HOME/gatk-4.1.0.0/libs:$LD_LIBRARY_PATH

workPath=$BMK_HOME/work_dir
ref=$REF_HOME/Homo_sapiens_assembly38.fasta
ref_dir=$REF_HOME

if [ -d "$workPath" ]; then
   echo "workPath already created"
else
   mkdir -p $workPath
fi
vcfFile=${workPath}/NA12878_merged.vcf

knownSites=($ref_dir/Homo_sapiens_assembly38.dbsnp138.sort.vcf
           $ref_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
           $ref_dir/Homo_sapiens_assembly38.known_indels.vcf.gz)

vcfHapmap=$ref_dir/hapmap_3.3.hg38.vcf.gz
vcfOmni=$ref_dir/1000G_omni2.5.hg38.vcf.gz
vcfGlk=$ref_dir/1000G_phase1.snps.high_confidence.hg38.vcf
vcfMills=$ref_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
vcfDbsnp=$ref_dir/dbsnp_138.hg38.vcf.gz

cd $workPath
# BWA MEM MAPPING and SAMTOOLS SORTING ############

input1=./input/NA12878_1.fastq.gz
input2=./input/NA12878_2.fastq.gz
output=$workPath/NA12878_hg38.bwa.bam

/usr/bin/time -v -o time_bwa.log taskset -c 0-159 bwa mem -t 160 -Ma \
     -R '@RG\tID:sample_lane\tSM:sample\tPL:illumina\tLB:sample\tPU:lane' \
     $ref $input1 $input2 | samtools view -bS - -@ 20 | samtools sort - -@ 40 -m 8G -T $input1 -o $output

# Markduplicates ############

input=$workPath/NA12878_hg38.bwa.bam
output=$workPath/NA12878_hg38.md.bam
/usr/bin/time -v -o time_markDup.log gatk --java-options "-Xmx30G -XX:+UseParallelGC -XX:ParallelGCThreads=4" MarkDuplicates \
      -I $input -O $output \
      -M $workPath/NA12878_bwa_metrics.txt --MAX_RECORDS_IN_RAM 5000000 -MAX_SEQS 5000000 \
      --OPTICAL_DUPLICATE_PIXEL_DISTANCE 2500 -MAX_FILE_HANDLES 1000

#create index ## Maybe no need to index
samtools index -@ 40 $output

# BASE QUALITY SCORE RECALIBRATION ##########
#Setup knownSites
for i in ${!knownSites[*]}
  do
    if [ $i == 0 ]
    then
      knownSiteArg="--known-sites ${knownSites[i]}"
    else
      knownSiteArg="${knownSiteArg} --known-sites ${knownSites[i]}"
    fi
done

input=$workPath/NA12878_hg38.md.bam
chr=0
chr2=3
for i in `seq -f '%04g' 0 39`
do
outfile=$workPath/NA12878_hg38.md.bam.br_$i.table
/usr/bin/time -v -o time_gatkBaseRecalibrator.log taskset -c $chr-$chr2 gatk \
         --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=20" BaseRecalibrator \
         -L $GATK_HOME/benchmarks/intervals/40c/$i-scattered.interval_list \
         -R $ref -I $input $knownSiteArg -O $outfile &
    chr=$(($chr+4))
    chr2=$(($chr2+4))
done
wait

## ApplyBQSR (note replacing PrintReads of GATK3)
chr=0
chr2=3
for i in `seq -f '%04g' 0 39`
do
bqfile=$workPath/NA12878_hg38.md.bam.br_$i.table
output=$workPath/NA12878_hg38.br.recal_$i.bam
/usr/bin/time -v -o time_gatkApplyBQSR.log taskset -c $chr-$chr2 gatk \
       --java-options "-Xmx4G -Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSR -R $ref -I $input \
       -L $GATK_HOME/benchmarks/intervals/40c/$i-scattered.interval_list  -bqsr $bqfile \
       --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O $output &
    chr=$(($chr+4))
    chr2=$(($chr2+4))
done
wait
	
# VARIANT CALLING optimized for 36 cores Power9 ##########
# Used Onle native pairhmm  # Run Haplotypecaller with split bam file
chr=0
chr2=3
for i in `seq -f '%04g' 0 39`
do
infile=$workPath/NA12878_hg38.br.recal_$i.bam
outfile=$workPath/NA12878_hg38.br.recal_$i.g.vcf
/usr/bin/time -v -o time_gatkHaplotypeCaller.log taskset -c $chr-$chr2 gatk \
      --java-options "-Xmx4G -Djava.library.path=$GATK_HOME/gatk-4.1.0.0/libs" HaplotypeCaller \
      -R ${ref} -I $infile \
      -L $GATK_HOME/benchmarks/intervals/40c/$i-scattered.interval_list \
      --native-pair-hmm-threads 4 --smith-waterman FASTEST_AVAILABLE \
      -O $outfile -ERC GVCF -stand-call-conf 10 &
   chr=$(($chr+4))
   chr2=$(($chr2+4))
done
wait

# combinge gvcf
gvcfFiles=()
for i in NA12878_hg38.br.recal_*.g.vcf
do
gvcfFiles[${#gvcfFiles[*]}]=$i
done

for i in ${!gvcfFiles[*]}
do
  if [ $i == 0 ]
  then
    gvcfFilesArg="--variant ${gvcfFiles[i]}"
  else
    gvcfFilesArg="${gvcfFilesArg} --variant ${gvcfFiles[i]}"
  fi
done

/usr/bin/time -v -o time_gatkCombineGVCFs.log gatk --java-options "-Xmx4G" CombineGVCFs -R $ref $gvcfFilesArg -O ${vcfFile%vcf}g.vcf
  
# genotype gvcf files
input=${vcfFile%vcf}g.vcf
/usr/bin/time -v -o time_gatkGenotypeGVCFs.log gatk --java-options "-Xmx4G" GenotypeGVCFs -R $ref -V $input -O $vcfFile
 
# VARIANT QUALITY SCORE RECALIBRATION (VQSR) ########
/usr/bin/time -v -o time_gatkVariantRecalibratorSNP.log gatk --java-options "-Xmx4G" VariantRecalibrator \
     -V $vcfFile -O ${workPath}/NA12878_recalibrate_SNP.recal \
     -mode SNP --tranches-file ${workPath}/NA12878_recalibrate_SNP.tranches \
     -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR -an MQ \
     --max-gaussians 6 \
     -resource:hapmap,known=false,training=true,truth=true,prior=15.0 $vcfHapmap \
     -resource:omni,known=false,training=true,truth=true,prior=12.0 $vcfOmni \
     -resource:1000G,known=false,training=true,truth=false,prior=10.0 $vcfGlk \
     -resource:dbsnp,known=true,training=false,truth=false,prior=7.0 $vcfDbsnp
  
  # Apply recalibration to SNPs
/usr/bin/time -v -o time_gatkApplyVQSRSNP.log gatk --java-options "-Xmx4G" ApplyVQSR \
     -V $vcfFile -O ${workPath}/NA12878_recalibrated_snps_raw_indels.vcf \
     --recal-file ${workPath}/NA12878_recalibrate_SNP.recal \
     --tranches-file ${workPath}/NA12878_recalibrate_SNP.tranches \
     -truth-sensitivity-filter-level 99.5 --create-output-variant-index true -mode SNP
  
  # Run Variant Recalibrator on Indels
/usr/bin/time -v -o time_gatkVariantRecalibratorIndel.log gatk --java-options "-Xmx4G" VariantRecalibrator \
     -V ${workPath}/NA12878_recalibrated_snps_raw_indels.vcf -O ${workPath}/NA12878_recalibrate_INDEL.recal \
     -mode INDEL --tranches-file ${workPath}/NA12878_recalibrate_INDEL.tranches \
     -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR \
     --max-gaussians 4 \
     -resource:mills,known=false,training=true,truth=true,prior=12.0 $vcfMills \
     -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 $vcfDbsnp \
       
   # Apply recalibration to INDELs  
/usr/bin/time -v -o time_gatkApplyVQSRIndel.log gatk --java-options "-Xmx4G" ApplyVQSR \
    -V ${workPath}/NA12878_recalibrated_snps_raw_indels.vcf -O ${vcfFile%.vcf}.recal_snps_indels.vcf \
    --recal-file ${workPath}/NA12878_recalibrate_INDEL.recal --tranches-file ${workPath}/NA12878_recalibrate_INDEL.tranches  \
    -truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode INDEL
