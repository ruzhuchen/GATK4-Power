#!/bin/bash

################################
# Power9 tests - WGS pipeline for Power9 system
################################

#Change this path for your test
export GATK_HOME=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG
export BMK_HOME=`pwd`
export RF_HOME=`pwd`/Ref

export PATH=$GATK_HOME/bin:$PATH
export GATK_LOCAL_JAR=$GATK_HOME/gatk-4.1.0.0/libs/gatk.jar
export GATK_SPARK_JAR=$GATK_HOME/gatk-4.1.0.0/libs/gatk-spark.jar
export LD_LIBRARY_PATH=$GATK_HOME/gatk-4.1.0.0/libs:$LD_LIBRARY_PATH

workPath=$BMK_HOME/work_dir
ref=$BMK_HOME/Ref/Homo_sapiens_assembly38.fasta
ref_dir=$BMK_HOME/Ref

if [ -d "$workPath" ]; then
   echo "workPath already created"
else
   mkdir -p $workPath
fi
knownSites=($ref_dir/Homo_sapiens_assembly38.dbsnp138.sort.vcf
           $ref_dir/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
           $ref_dir/Homo_sapiens_assembly38.known_indels.vcf.gz)

cd $workPath

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
/usr/bin/time -v -o time_gatkBaseRecalibratorT.log taskset -c 0-143:4 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
         BaseRecalibratorSpark \
         -R $ref -I $input1 $knownSiteArg -O $outfile1 \
         -- --spark-runner LOCAL --spark-master local[36] --conf spark.local.dir=$workPath

# Normal
outfile2=$workPath/YX160000128N_hg38.md.bam.br.table
/usr/bin/time -v -o time_gatkBaseRecalibratorN.log taskset -c 0-143:4 gatk --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" \
         BaseRecalibratorSpark \
         -R $ref -I $input2 $knownSiteArg -O $outfile2 \
         -- --spark-runner LOCAL --spark-master local[36] --conf spark.local.dir=$workPath

# ApplyBQSR Tumor
bqfile1=$workPath/YX160000128T_hg38.md.bam.br.table
output1=$workPath/YX160000128T_hg38.br.recal.bam
/usr/bin/time -v -o time_gatkApplyBQSR_T.log taskset -c 0-143:4 gatk \
       --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSRSpark -R $ref -I $input1 \
       -bqsr $bqfile1 \
       --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O $output1 \
       -- --spark-runner LOCAL --spark-master local[36] --conf spark.local.dir=$workPath
samtools index -@ 36 $output1

# ApplyBQSR for Normal
bqfile2=$workPath/YX160000128N_hg38.md.bam.br.table
output2=$workPath/YX160000128N_hg38.br.recal.bam
/usr/bin/time -v -o time_gatkApplyBQSR_N.log taskset -c 0-143:4 gatk \
       --java-options "-Xmx8G -XX:+UseParallelGC -XX:ParallelGCThreads=4" ApplyBQSRSpark -R $ref -I $input2 \
       -bqsr $bqfile2 \
       --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O $output2 \
       -- --spark-runner LOCAL --spark-master local[36] --conf spark.local.dir=$workPath
samtools index -@ 36 $output2
