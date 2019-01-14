#!/bin/bash

################################
# Power9 tests - WGS Somatic pipeline for 40 cores Power9 system
################################
export PATH=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/bin:$PATH
export GATK_LOCAL_JAR=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs/gatk.jar
export GATK_SPARK_JAR=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs/gatk-spark.jar
export LD_LIBRARY_PATH=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs:$LD_LIBRARY_PATH

workPath=/gpfs/gpfs_4mb/rchen/Power9/GATK4/benchmarks/YX160000128/work_opt
ref=/gpfs/gpfs_4mb/rchen/Power9/GATK4/benchmarks/Ref/Homo_sapiens_assembly38.fasta
ref_dir=/gpfs/gpfs_4mb/rchen/Power9/GATK4/benchmarks/Ref

knownSites=(/gpfs/gpfs_4mb/rchen/Power9/GATK4/benchmarks/Ref/dbsnp_146.hg38.vcf.gz
           /gpfs/gpfs_4mb/rchen/Power9/GATK4/benchmarks/Ref/Homo_sapiens_assembly38.known_indels.vcf.gz
           /gpfs/gpfs_4mb/rchen/Power9/GATK4/benchmarks/Ref/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz)

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

/usr/bin/time -v -o time_bwa_T.log taskset -c 0-143:4 bwa mem -t 36 -Ma \
     -R '@RG\tID:YX160000128T_lane\tSM:YX160000128T\tPL:illumina\tLB:YX160000128T\tPU:lane' \
     $ref $input1 $input2 \
     | samtools view -bS - -@ 18 | samtools sort - -@ 18 -n -m 16G -T YX160000128T -o $output1

/usr/bin/time -v -o time_bwa_N.log taskset -c 0-143:4 bwa mem -t 36 -Ma \
     -R '@RG\tID:YX160000128N_lane\tSM:YX160000128N\tPL:illumina\tLB:YX160000128N\tPU:lane' \
     $ref $input3 $input4 \
     | samtools view -bS - -@ 18 | samtools sort - -@ 18 -n -m 16G -T YX160000128N -o $output2
# Markduplicates ############

input1=$workPath/YX160000128T.bwa.bam
input2=$workPath/YX160000128N.bwa.bam
output1=$workPath/YX160000128T.md.bam
output2=$workPath/YX160000128N.md.bam

/usr/bin/time -v -o time_Fixmate_mk1.log taskset -c 0-143:4 samtools fixmate -m -@ 36 $input1 fixmate1.bam
/usr/bin/time -v -o time_Fixmate_mk2.log taskset -c 0-143:4 samtools fixmate -m -@ 36 $input2 fixmate2.bam
/usr/bin/time -v -o time_Sort_mk1.log samtools sort -@ 36 -o sorted1.bam fixmate1.bam
/usr/bin/time -v -o time_Sort_mk2.log samtools sort -@ 36 -o sorted2.bam fixmate2.bam
/usr/bin/time -v -o time_Markduplicates1.log taskset -c 0-143:4 samtools markdup -s -@ 36 sorted1.bam $output1
/usr/bin/time -v -o time_Markduplicates2.log taskset -c 0-159:4 samtools markdup -s -@ 36 sorted2.bam $output2

rm -f fixmate*.bam sorted*.bam

#create index
samtools index -@ 36 $output1
samtools index -@ 36 $output2

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
/usr/bin/time -v -o time_gatkBaseRecalibratorT.log taskset -c 0-143:4 gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
         BaseRecalibratorSpark \
         -R $ref -I $input1 $knownSiteArg -O $outfile1 \
         -- --spark-runner LOCAL --spark-master local[36] --conf spark.local.dir=$workPath

# Normal
outfile2=$workPath/YX160000128N_hg38.md.bam.br.table
/usr/bin/time -v -o time_gatkBaseRecalibratorN.log taskset -c 0-143:4 gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
         BaseRecalibratorSpark \
         -R $ref -I $input2 $knownSiteArg -O $outfile2 \
         -- --spark-runner LOCAL --spark-master local[36] --conf spark.local.dir=$workPath

# ApplyBQSR Tumor
bqfile1=$workPath/YX160000128T_hg38.md.bam.br.table
output1=$workPath/YX160000128T_hg38.br.recal.bam
/usr/bin/time -v -o time_gatkApplyBQSR_T.log taskset -c 0-143:4 gatk \
       --java-options "-Xmx4G -Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSRSpark -R $ref -I $input1 \
       -bqsr $bqfile1 \
       --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O $output1 \
       -- --spark-runner LOCAL --spark-master local[36] --conf spark.local.dir=$workPath
samtools index -@ 36 $output1

# ApplyBQSR for Normal
bqfile2=$workPath/YX160000128N_hg38.md.bam.br.table
output2=$workPath/YX160000128N_hg38.br.recal.bam
/usr/bin/time -v -o time_gatkApplyBQSR_N.log taskset -c 0-143:4 gatk \
       --java-options "-Xmx4G -Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSRSpark -R $ref -I $input2 \
       -bqsr $bqfile2 \
       --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O $output2 \
       -- --spark-runner LOCAL --spark-master local[36] --conf spark.local.dir=$workPath
samtools index -@ 36 $output2

# SOMATIC PIPELINE: MUTECT2 ##########
# Based on blog https://software.broadinstitute.org/gatk/documentation/article?id=11136
# 1. Call somatic short variants and generate a bamout
Tumor_input=$workPath/YX160000128T_hg38.br.recal.bam
Normal_input=$workPath/YX160000128N_hg38.br.recal.bam

chr=0
chr2=15
for i in `seq 1 22`
do
M2_output=$workPath/YX160000128_somatic_m2_$i.vcf.gz
if [ $i -lt 5 ]
then
/usr/bin/time -v -o time_gatkMutect2_step1.log taskset -c $chr-$chr2 gatk \
       --java-options "-Xmx4g -Djava.library.path=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs" Mutect2 \
       -R $ref -I $Tumor_input -I $Normal_input \
       -tumor YX160000128T \
       -normal YX160000128N \
       --native-pair-hmm-threads 16 \
       -pon $ref_dir/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz \
       --germline-resource $ref_dir/somatic/af-only-gnomad.hg38.vcf.gz \
       --af-of-alleles-not-in-resource 0.0000025 \
       --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
       -L chr$i \
       -O $M2_output &
   if [ $i -lt 4 ]
   then
      chr=$(($chr+16))
      chr2=$(($chr2+16))
   else
      chr=$(($chr+16))
      chr2=$(($chr2+4))
   fi
else
/usr/bin/time -v -o time_gatkMutect2_step1.log taskset -c $chr-$chr2 gatk \
       --java-options "-Xmx4g -Djava.library.path=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs" Mutect2 \
       -R $ref -I $Tumor_input -I $Normal_input \
       -tumor YX160000128T \
       -normal YX160000128N \
       --native-pair-hmm-threads 4 \
       -pon $ref_dir/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz \
       --germline-resource $ref_dir/somatic/af-only-gnomad.hg38.vcf.gz \
       --af-of-alleles-not-in-resource 0.0000025 \
       --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
       -L chr$i \
       -O $M2_output &
   chr=$(($chr+4))
   chr2=$(($chr2+4))
fi
done
for i in X Y M
do
M2_output=$workPath/YX160000128_somatic_m2_$i.vcf.gz
if [ $i == X ]
then
thd=4
else
thd=2
fi
/usr/bin/time -v -o time_gatkMutect2_step1.log taskset -c $chr-$chr2 gatk \
       --java-options "-Xmx4g -Djava.library.path=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs" Mutect2 \
       -R $ref -I $Tumor_input -I $Normal_input \
       -tumor YX160000128T \
       -normal YX160000128N \
       --native-pair-hmm-threads $thd \
       -pon $ref_dir/somatic/somatic-hg38_1000g_pon.hg38.vcf.gz \
       --germline-resource $ref_dir/somatic/af-only-gnomad.hg38.vcf.gz \
       --af-of-alleles-not-in-resource 0.0000025 \
       --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
       -L chr$i \
       -O $M2_output &
if [ $i == X ]
then
   chr=$(($chr+4))
   chr2=$(($chr2+2))
else
   chr=$(($chr+2))
   chr2=$(($chr2+2))
fi
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
      -O YX160000128_somatic_m2.vcf.gz

rm -f YX160000128_somatic_m2_*.vcf.gz*

# 2.1 Create a sites-only PoN in "tumor-only mode" on normal sample

chr=0
chr2=15
for i in `seq 1 22`
do
N_Output=YX160000128_PoN_$i.vcf.gz
if [ $i -lt 5 ]
then
/usr/bin/time -v -o time_gatkMutect2_step2-1.log taskset -c $chr-$chr2 gatk \
      --java-options "-Xmx4g -Djava.library.path=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs" Mutect2 \
      -R $ref -I $workPath/YX160000128N_hg38.br.recal.bam \
      -tumor YX160000128N \
       --native-pair-hmm-threads 16 \
      --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
      -L chr$i \
      -O $N_Output &
   if [ $i -lt 4 ]
   then
      chr=$(($chr+16))
      chr2=$(($chr2+16))
   else
      chr=$(($chr+16))
      chr2=$(($chr2+4))
   fi
else
/usr/bin/time -v -o time_gatkMutect2_step2-1.log taskset -c $chr-$chr2 gatk \
      --java-options "-Xmx4g -Djava.library.path=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs" Mutect2 \
      -R $ref -I $workPath/YX160000128N_hg38.br.recal.bam \
      -tumor YX160000128N \
       --native-pair-hmm-threads 4 \
      --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
      -L chr$i \
      -O $N_Output &
   chr=$(($chr+4))
   chr2=$(($chr2+4))
fi
done
for i in X Y M
do
N_Output=YX160000128_PoN_$i.vcf.gz
if [ $i == X ]
then
thd=4
else
thd=2
fi
/usr/bin/time -v -o time_gatkMutect2_step2-xym.log taskset -c $chr-$chr2 gatk \
      --java-options "-Xmx4g -Djava.library.path=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG/gatk-4.0.11.0/libs" Mutect2 \
      -R $ref -I $workPath/YX160000128N_hg38.br.recal.bam \
      -tumor YX160000128N \
       --native-pair-hmm-threads $thd \
      --disable-read-filter MateOnSameContigOrNoMappedMateReadFilter \
      -L chr$i \
      -O $N_Output &
if [ $i == X ]
then
   chr=$(($chr+4))
   chr2=$(($chr2+2))
else
   chr=$(($chr+2))
   chr2=$(($chr2+2))
fi
done
wait

vvFiles=()
for i in YX160000128_PoN_*.vcf.gz
do
   vvFiles[${#vvFiles[*]}]=$i
done
for i in ${!vvFiles[*]}
do
   if [ $i == 0 ]
   then
      vvFilesArg="-I ${vvFiles[i]}"
   else
      vvFilesArg="${vvFilesArg} -I ${vvFiles[i]}"
   fi
done
/usr/bin/time -v -o time_gatkMergeVCF-PonVCF.log gatk --java-options "-Xmx4g" MergeVcfs \
      $vvFilesArg -R $ref \
      -O YX160000128_PoN.vcf.gz
rm -f YX160000128_PoN_*.vcf.gz*

# Ruzhu: no need to run this step, since only one VCF file
# 2.2 Create a sites-only PoN with CreateSomaticPanelOfNormals
#/usr/bin/time -v -o time_gatkCreateSomaticPanelOfNormals.log gatk --java-options "-Xmx4g" CreateSomaticPanelOfNormals \
#      -vcfs YX160000128_PoN.vcf.gz \
#      -O YX160000128_hg38.pon.vcf.gz

# 3. Estimate cross-sample contamination using GetPileupSummaries and CalculateContamination
/usr/bin/time -v -o time_gatkGetPileupSummaries.log gatk --java-options "-Xmx4g" GetPileupSummaries \
      -I $workPath/YX160000128T_hg38.br.recal.bam \
      -V $ref_dir/somatic/small_exac_common_3.hg38.vcf.gz \
      -L $ref_dir/wgs_calling_regions.hg38.interval_list \
      -O YX160000128T_getpileupsummaries.table

/usr/bin/time -v -o time_gatkCalculateContamination.log gatk --java-options "-Xmx4g" CalculateContamination \
      -I YX160000128T_getpileupsummaries.table \
      -O YX160000128T_calculatecontamination.table

# 4. Filter for confident somatic calls using FilterMutectCalls
/usr/bin/time -v -o time_gatkFilterMutectCalls.log gatk --java-options "-Xmx4g" FilterMutectCalls \
      -V YX160000128_somatic_m2.vcf.gz \
      --contamination-table YX160000128T_calculatecontamination.table \
      -O YX160000128_somatic_oncefiltered.vcf.gz

##### OPTINAL STEPS #####
# 5. Estimate artifacts with CollectSequencingArtifactMetrics and filter them with FilterByOrientationBias
/usr/bin/time -v -o time_gatkCollectSequencingArtifactMetrics.log gatk --java-options "-Xmx4g" CollectSequencingArtifactMetrics \
      -R $ref -I $workPath/YX160000128T_hg38.br.recal.bam \
      -O YX160000128_tumor_artifact

/usr/bin/time -v -o time_gatkFilterByOrientationBias.log gatk --java-options "-Xmx4g" FilterByOrientationBias \
      -AM 'G/T' -AM 'C/T' \
      -V YX160000128_somatic_oncefiltered.vcf.gz \
      -P YX160000128_tumor_artifact.pre_adapter_detail_metrics \

echo "Finish pipeline at $(date)"
      -O YX160000128_somatic_twicefiltered.vcf.gz
