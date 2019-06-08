#!/bin/bash

################################
# Power9 tests - WGS for 40 cores Power9 system
################################

#Change this path for your test
export GATK_HOME=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG
export BMK_HOME=`pwd`
export RF_HOME=`pwd`/Ref

export PATH=$GATK_HOME/bin:$PATH
export GATK_LOCAL_JAR=$GATK_HOME/gatk-4.1.0.0/libs/gatk.jar
export GATK_SPARK_JAR=$GATK_HOME/gatk-4.1.0.0/libs/gatk-spark.jar
export LD_LIBRARY_PATH=$GATK_HOME/gatk-4.1.0.0/libs:$LD_LIBRARY_PATH

workPath=$BMK_HOME/intervals
ref=$BMK_HOME/Ref/Homo_sapiens_assembly38.fasta
ref_dir=$BMK_HOME/Ref

mkdir -p $workPath
cd $workPath

gatk --java-options "-Xmx4G" SplitIntervals \
    -R $ref \
    -L $ref_dir/wgs_calling_regions.hg38.interval_list \
    -scatter 40 \
    -O $workPath/40c
