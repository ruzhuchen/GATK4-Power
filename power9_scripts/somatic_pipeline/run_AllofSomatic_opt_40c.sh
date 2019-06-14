#!/bin/bash

################################
# Power9 tests - WGS pipeline for Power9 system
################################

#Change this path for your test
export GATK_HOME=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG
export BMK_HOME=
export REF_HOME=

export PATH=$GATK_HOME/bin:$PATH
export GATK_LOCAL_JAR=$GATK_HOME/gatk-4.1.2.0/libs/gatk.jar
export GATK_SPARK_JAR=$GATK_HOME/gatk-4.1.2.0/libs/gatk-spark.jar
export LD_LIBRARY_PATH=$GATK_HOME/gatk-4.1.2.0/libs:$LD_LIBRARY_PATH

workPath=$BMK_HOME/work_dir
ref=$REF_HOME/Homo_sapiens_assembly38.fasta
ref_dir=$REF_HOME

if [ -d "$workPath" ]; then
   echo "workPath already created"
else
   mkdir -p $workPath
fi
knownSites=($ref_dir/Homo_sapiens_assembly38.dbsnp138.sort.vcf
           $ref_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
           $ref_dir/Homo_sapiens_assembly38.known_indels.vcf.gz)

cd $workPath
echo "Start pipeline at $(date)"

# BWA MEM MAPPING and SAMTOOLS SORTING ############
# The two inputs YX160000128N/T were combined 
input1=$workPath/../input/YX160000128T_1.fq.gz
input2=$workPath/../input/YX160000128T_2.fq.gz
input3=$workPath/../input/YX160000128N_1.fq.gz
input4=$workPath/../input/YX160000128N_2.fq.gz
output1=$workPath/YX160000128T.bwa.bam
output2=$workPath/YX160000128N.bwa.bam

/usr/bin/time -v -o time_bwa_T.log taskset -c 0-159:4 bwa mem -t 40 -Ma \
     -R '@RG\tID:YX160000128T_lane\tSM:YX160000128T\tPL:illumina\tLB:YX160000128T\tPU:lane' \
     $ref $input1 $input2 \
     | samtools view -bS - -@ 20 | samtools sort - -@ 20 -n -m 16G -T YX160000128T -o $output1

/usr/bin/time -v -o time_bwa_N.log taskset -c 0-159:4 bwa mem -t 40 -Ma \
     -R '@RG\tID:YX160000128N_lane\tSM:YX160000128N\tPL:illumina\tLB:YX160000128N\tPU:lane' \
     $ref $input3 $input4 \
     | samtools view -bS - -@ 20 | samtools sort - -@ 20 -n -m 16G -T YX160000128N -o $output2
# Markduplicates ############

input1=$workPath/YX160000128T.bwa.bam
input2=$workPath/YX160000128N.bwa.bam
output1=$workPath/YX160000128T.md.bam
output2=$workPath/YX160000128N.md.bam

/usr/bin/time -v -o time_Fixmate_mk1.log taskset -c 0-159:4 samtools fixmate -m -@ 40 $input1 fixmate1.bam
/usr/bin/time -v -o time_Fixmate_mk2.log taskset -c 0-159:4 samtools fixmate -m -@ 40 $input2 fixmate2.bam
/usr/bin/time -v -o time_Sort_mk1.log samtools sort -@ 40 -o sorted1.bam fixmate1.bam
/usr/bin/time -v -o time_Sort_mk2.log samtools sort -@ 40 -o sorted2.bam fixmate2.bam
/usr/bin/time -v -o time_Markduplicates1.log taskset -c 0-159:4 samtools markdup -s -@ 40 sorted1.bam $output1
/usr/bin/time -v -o time_Markduplicates2.log taskset -c 0-159:4 samtools markdup -s -@ 40 sorted2.bam $output2

rm -f fixmate*.bam sorted*.bam

#create index
samtools index -@ 40 $output1
samtools index -@ 40 $output2

# BASE QUALITY SCORE RECALIBRATION ##########
	
input1=$workPath/YX160000128T.md.bam
input2=$workPath/YX160000128N.md.bam

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

# Tumor
outfile1=$workPath/YX160000128T_hg38.md.bam.br.table
/usr/bin/time -v -o time_gatkBaseRecalibratorT.log taskset -c 0-159:4 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
         BaseRecalibratorSpark \
         -R $ref -I $input1 $knownSiteArg -O $outfile1 \
         -- --spark-runner LOCAL --spark-master local[40] --conf spark.local.dir=$workPath

# Normal
outfile2=$workPath/YX160000128N_hg38.md.bam.br.table
/usr/bin/time -v -o time_gatkBaseRecalibratorN.log taskset -c 0-159:4 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
         BaseRecalibratorSpark \
         -R $ref -I $input2 $knownSiteArg -O $outfile2 \
         -- --spark-runner LOCAL --spark-master local[40] --conf spark.local.dir=$workPath

# ApplyBQSR Tumor
bqfile1=$workPath/YX160000128T_hg38.md.bam.br.table
output1=$workPath/YX160000128T_hg38.br.recal.bam
/usr/bin/time -v -o time_gatkApplyBQSR_T.log taskset -c 0-159:4 gatk \
       --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSRSpark -R $ref -I $input1 \
       -bqsr $bqfile1 \
       --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O $output1 \
       -- --spark-runner LOCAL --spark-master local[40] --conf spark.local.dir=$workPath
samtools index -@ 40 $output1

# ApplyBQSR for Normal
bqfile2=$workPath/YX160000128N_hg38.md.bam.br.table
output2=$workPath/YX160000128N_hg38.br.recal.bam
/usr/bin/time -v -o time_gatkApplyBQSR_N.log taskset -c 0-159:4 gatk \
       --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSRSpark -R $ref -I $input2 \
       -bqsr $bqfile2 \
       --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O $output2 \
       -- --spark-runner LOCAL --spark-master local[40] --conf spark.local.dir=$workPath
samtools index -@ 40 $output2

# SOMATIC PIPELINE: MUTECT2 ##########
# Based on blog https://software.broadinstitute.org/gatk/documentation/article?id=11136
# 1. Call somatic short variants and generate a bamout
Tumor_input=$workPath/YX160000128T_hg38.br.recal.bam
Normal_input=$workPath/YX160000128N_hg38.br.recal.bam
chr=0
chr2=3
for i in `seq -f '%04g' 0 39`
do
M2_output=$workPath/YX160000128_m2_$i.vcf.gz
/usr/bin/time -v -o time_gatkMutect2_m2.log taskset -c $chr-$chr2 gatk \
       --java-options "-Xmx4g -Djava.library.path=$GATK_HOME/gatk-4.1.2.0/libs" Mutect2 \
       -R $ref -I $Tumor_input -tumor YX160000128T \
       -I $Normal_input -normal YX160000128N \
       --native-pair-hmm-threads 4 \
       -pon $ref_dir/somatic/M2PoN_4.0_WGS_for_public.vcf \
       --germline-resource $ref_dir/somatic/af-only-gnomad.hg38.vcf.gz \
       -L $GATK_HOME/benchmarks/intervals/40c/$i-scattered.interval_list \
       -O $M2_output &
   chr=$(($chr+4))
   chr2=$(($chr2+4))
done
wait

# Merge VCF files
vFiles=()
for i in YX160000128_somatic_m2_*.vcf.gz
do
   vFiles[${#vFiles[*]}]=$i
done
for i in ${!vFiles[*]}
do
   if [ $i == 0 ]
   then
      vFilesArg="-I ${vFiles[i]}"
   else
      vFilesArg="${vFilesArg} -I ${vFiles[i]}"
   fi
done

/usr/bin/time -v -o time_gatkMergeVCF-m2.log gatk --java-options "-Xmx4g" MergeVcfs \
      $vFilesArg -R $ref \
      -O YX160000128_somatic_unfiltered.vcf.gz

rm -f YX160000128_somatic_m2_*.vcf.gz*

# 2. Estimate cross-sample contamination using GetPileupSummaries and CalculateContamination
/usr/bin/time -v -o time_gatkGetPileupSummaries.log gatk --java-options "-Xmx4g" GetPileupSummaries \
      -I $workPath/YX160000128T_hg38.br.recal.bam \
      -V $ref_dir/somatic/small_exac_common_3.hg38.vcf.gz \
      -L $ref_dir/wgs_calling_regions.hg38.interval_list \
      -O YX160000128T_getpileupsummaries.table

/usr/bin/time -v -o time_gatkCalculateContamination.log gatk --java-options "-Xmx4g" CalculateContamination \
      -I YX160000128T_getpileupsummaries.table \
      -O YX160000128T_calculatecontamination.table

# 3. Filter for confident somatic calls using FilterMutectCalls
/usr/bin/time -v -o time_gatkFilterMutectCalls.log gatk --java-options "-Xmx4g" FilterMutectCalls \
      -V YX160000128_somatic_unfiltered.vcf.gz \
      --contamination-table YX160000128T_calculatecontamination.table \
      -O YX160000128_somatic_filtered.vcf.gz

##### OPTINAL STEPS #####
# 4. Estimate artifacts with CollectSequencingArtifactMetrics and filter them with FilterByOrientationBias
/usr/bin/time -v -o time_gatkCollectSequencingArtifactMetrics.log gatk --java-options "-Xmx4g" CollectSequencingArtifactMetrics \
      -R $ref -I $workPath/YX160000128T_hg38.br.recal.bam \
      -O YX160000128_tumor_artifact

/usr/bin/time -v -o time_gatkFilterByOrientationBias.log gatk --java-options "-Xmx4g" FilterByOrientationBias \
      -AM 'G/T' -AM 'C/T' \
      -V YX160000128_somatic_filtered.vcf.gz \
      -P YX160000128_tumor_artifact.pre_adapter_detail_metrics \
      -O YX160000128_somatic_twicefiltered.vcf.gz

echo "Finish pipeline at $(date)"