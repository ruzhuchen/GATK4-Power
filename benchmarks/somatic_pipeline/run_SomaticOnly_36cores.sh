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
       --java-options "-Xmx4g -Djava.library.path=$GATK_HOME/gatk-4.1.0.0/libs" Mutect2 \
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
       --java-options "-Xmx4g -Djava.library.path=$GATK_HOME/gatk-4.1.0.0/libs" Mutect2 \
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
       --java-options "-Xmx4g -Djava.library.path=$GATK_HOME/gatk-4.1.0.0/libs" Mutect2 \
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
      --java-options "-Xmx4g -Djava.library.path=$GATK_HOME/gatk-4.1.0.0/libs" Mutect2 \
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
      --java-options "-Xmx4g -Djava.library.path=$GATK_HOME/gatk-4.1.0.0/libs" Mutect2 \
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
      --java-options "-Xmx4g -Djava.library.path=$GATK_HOME/gatk-4.1.0.0/libs" Mutect2 \
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
# 5. (Optional) Estimate artifacts with CollectSequencingArtifactMetrics and filter them with FilterByOrientationBias
/usr/bin/time -v -o time_gatkCollectSequencingArtifactMetrics.log gatk --java-options "-Xmx4g" CollectSequencingArtifactMetrics \
      -R $ref -I $workPath/YX160000128T_hg38.br.recal.bam \
      -O YX160000128_tumor_artifact

/usr/bin/time -v -o time_gatkFilterByOrientationBias.log gatk --java-options "-Xmx4g" FilterByOrientationBias \
      -AM 'G/T' -AM 'C/T' \
      -V YX160000128_somatic_oncefiltered.vcf.gz \
      -P YX160000128_tumor_artifact.pre_adapter_detail_metrics \
      -O YX160000128_somatic_twicefiltered.vcf.gz
