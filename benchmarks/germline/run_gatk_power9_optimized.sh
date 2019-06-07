#!/bin/bash
#
# Power9 Germline Scripts (all systems) 
############################################
#Change this path for your test
GATK_HOME=/gpfs/gpfs_4mb/rchen/Power9/GATK4/P9_PKG
REF_HOME=/gpfs/gpfs_4mb/rchen/Power9/GATK4/benchmarks/Ref
BMK_HOME=/gpfs/gpfs_4mb/rchen/Power9/GATK4/benchmarks/NA12878
INPUT_FOLDER=$BMK_HOME/input

nprs=$(awk '/processor/' /proc/cpuinfo | wc -l)

usage()
{
   echo "Usage: -g gatk_home -r reference_dir -w benchmark_dir -i input_dir"
}
opt_error()
{ 
   echo "The $1 path does not exist "
   echo "Please check $1 folder path!"
   exit 1
}

if [ "$#" -eq 0 ]; then
   echo "The default paths will be used. If you like to set your own paths, see usage: "
   echo "GATK_HOME=$GATK_HOME"
   echo "REF_HOME=$REF_HOME"
   echo "BMK_HOME=$BMK_HOME"
   echo "INPUT_FOLDER=$INPUT_FOLDER"
fi
read -p "Are all path settings correct? [Y or N] " -n 1 -r
echo " "
if [[ ! $REPLY =~ ^[Yy]$ ]]; then 
   usage
   exit 1
fi
read -p "Using Power9 sam2bam? [Y | N]" -n 1 -r
echo " " 
if [[ $REPLY =~ ^[Yy]$ ]];then
   opt=yes
   filename=run_gatk_sam2bam-$((nprs/4))c.sh
else
   opt=no
   filename=run_gatk-$((nprs/4))c.sh
fi
echo "the choice is $opt and the script is $filename"

while [ "$1" != "" ]; do
   case $1 in
      -g | --gatk_home ) shift
                   GATK_HOME=$1
                   ;;
      -r | --REF_HOME )   shift
                   REF_HOME=$1
                   ;;
      -w | --work_dir ) shift
                   BMK_HOME=$1
                   ;;
      -i | --input_dir ) shift
                   INPUT_FOLDER=$1
                   ;;
      -h | --help)  usage
                   exit
                   ;;
      * ) usage
          exit 1
   esac
   shift
done
if [ ! -d "$GATK_HOME" ]; then
   opt_error $GATK_HOME 
fi
if [ ! -d "$REF_HOME" ]; then
   opt_error $REF_HOME
fi
# Below is the pipeline script
sam2bam()
{
if [ "$opt" = "yes" ];then
cat <<EOF
# MAPPING ############
# bwa:
  mapped=()
  for i in \${fastqFolder}/*_1.fastq.gz
  do
    filename=\$( basename \$i _1.fastq.gz )
    mapFile=\${filename}_bwa.bam
    vcfFile=\${filename}_merged.vcf
    mapped[\${#mapped[*]}]=\$mapFile
    /usr/bin/time -v -o timeBWA.log taskset -c 0-$((nprs-1)) bwa mem -t $nprs -Ma \\
        -R @RG\\tID:\${filename}\\tSM:\${filename}\\tPL:ILLUMINA\\tLB:\${filename} -Y \$ref \\
        \$i \${i%1.fastq.gz}2.fastq.gz >\${mapFile%bam}sam
 done

step1=\$(date -u +%s)
bwa_time=\$((\$step1 - \$start_time))
Elapsed[\${#Elapsed[*]}]="BWA+Samtool is \$bwa_time seconds"
# MARK DUPLICATES #######
  for mapFile in \${mapped[*]}
  do
    # mark duplicates
   export use_storage_mode=yes
   export BAM_PAGEFILE=\$workPath/pf
   /usr/bin/time -v -o time_Markduplicates.log sam2bam sam2bam -d -p -Fibm_markdup:r -o\${mapFile%.bam}_dedup.bam \${mapFile%bam}sam
  done
EOF
else
cat <<EOF
# MAPPING ############
# bwa:
  mapped=()
  for i in \${fastqFolder}/*_1.fastq.gz
  do
    filename=\$( basename \$i _1.fastq.gz )
    mapFile=\${filename}_bwa.bam
    vcfFile=\${filename}_merged.vcf
    mapped[\${#mapped[*]}]=\$mapFile
    /usr/bin/time -v -o timeBWA.log taskset -c 0-$((nprs-1)) bwa mem -t $nprs -Ma \\
        -R @RG\\tID:\${filename}\\tSM:\${filename}\\tPL:ILLUMINA\\tLB:\${filename} -Y \$ref \\
        \$i \${i%1.fastq.gz}2.fastq.gz | samtools sort - -@ $((nprs/4)) -n -m 6G -T \${i%R1.fastq.gz} -o \$mapFile
 done

step1=\$(date -u +%s)
bwa_time=\$((\$step1 - \$start_time))
Elapsed[\${#Elapsed[*]}]="BWA+Samtool is \$bwa_time seconds"

# MARK DUPLICATES #######
  for mapFile in \${mapped[*]}
  do
    # mark duplicates
   /usr/bin/time -v -o time_Fixmate_mk.log taskset -c 0-$((nprs-1)):4 samtools fixmate -m -@ $((nprs/4)) \${mapFile} fixmate.bam
   /usr/bin/time -v -o time_Sort_mk.log taskset -c 0-$((nprs-1)):4 samtools sort -@ $((nprs/4)) -m 8G -o sorted.bam fixmate.bam
   /usr/bin/time -v -o time_Markduplicates.log taskset -c 0-$((nprs-1)):4 samtools markdup -s -@ $((nprs/4)) sorted.bam \${mapFile%.bam}_dedup.bam
  done
  rm -f fixmate.bam sorted.bam
EOF
fi
}

cat >$filename<<EOF
#!/bin/bash

# GATK4 Script for running on $((nprs/4)) cores $(hostname)

GATK_HOME=$GATK_HOME
REF_HOME=$REF_HOME
BMK_HOME=$BMK_HOME

fastqFolder=$INPUT_FOLDER
ref=\$REF_HOME/Homo_sapiens_assembly38.fasta

export PATH=\$GATK_HOME/bin:\$PATH
export GATK_LOCAL_JAR=\$GATK_HOME/gatk-4.1.0.0/libs/gatk.jar
export GATK_SPARK_JAR=\$GATK_HOME/gatk-4.1.0.0/libs/gatk-spark.jar
export LD_LIBRARY_PATH=\$GATK_HOME/gatk-4.1.0.0/libs:\$LD_LIBRARY_PATH

knownSites=(\$REF_HOME/Homo_sapiens_assembly38.dbsnp138.sort.vcf
           \$REF_HOME/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
           \$REF_HOME/Homo_sapiens_assembly38.known_indels.vcf.gz)

vcfHapmap=\$REF_HOME/hapmap_3.3.hg38.vcf.gz
vcfOmni=\$REF_HOME/1000G_omni2.5.hg38.vcf.gz
vcfGlk=\$REF_HOME/1000G_phase1.snps.high_confidence.hg38.vcf
vcfMills=\$REF_HOME/Mills_and_1000G_gold_standard.indels.hg38.vcf.gz
vcfDbsnp=\$REF_HOME/dbsnp_138.hg38.vcf.gz

# Preparing Interval lists before starting pileline
if [ ! -d "\$REF_HOME/intervals/$((nprs/4))c" ]; then
   if [ ! "\$REF_HOME/intervals" ]; then
      mkdir -p \$REF_HOME/intervals
   fi
   echo "Interval_list file not existsi, creating one"
   gatk --java-options "-Xmx4G" ScatterIntervalsByNs \\
       -R \$ref -O \$REF_HOME/hg38_scatter.interval_list
   gatk --java-options "-Xmx4G" SplitIntervals \\
      -R \$ref -L \$REF_HOME/hg38_scatter.interval_list \\
      -scatter $((nprs/4)) \\
      -O \$REF_HOME/intervals/$((nprs/4))c
fi

workPath=\$BMK_HOME/work\$$
if [ -d "\$workPath" ]; then
   echo "workPath already created"
else
   mkdir -p \$workPath
fi
cd \$workPath
echo "Starting WES pipeline at \$(date)"
Elapsed=()
start_time=\$(date -u +%s)
$(sam2bam)

step2=\$(date -u +%s)
mkd_time=\$((\$step2 - \$step1))
Elapsed[\${#Elapsed[*]}]="MarkDuplicates is \$mkd_time seconds"

# delete SAM file
rm -f \${mapFile%bam}sam

for mapFile in \${mapped[*]}
do 
    # 1) Analyse covariation patterns in the dataset
    for i in \${!knownSites[*]}
    do
      if [ \$i == 0 ]
      then
        knownSiteArg="--known-sites \${knownSites[i]}"
      else
        knownSiteArg="\${knownSiteArg} --known-sites \${knownSites[i]}"
      fi      
    done

   for i in \$(seq -f '%04g' 0 $(((nprs/4)-1)))
   do
      outfile=\${mapFile%.bam}_dedup_recal_data_\$i.table
      /usr/bin/time -v -o time_gatkBaseRecalibrator.log gatk \\
         --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" BaseRecalibrator \\
         -L \$REF_HOME/intervals/$((nprs/4))c/\$i-scattered.interval_list \\
         -R \$ref -I \${mapFile%.bam}_dedup.bam \$knownSiteArg -O \$outfile &
   done
   wait

step3=\$(date -u +%s)
elp_time=\$((\$step3 - \$step2))
Elapsed[\${#Elapsed[*]}]="BaseRecalibration is \$elp_time seconds"

## ApplyBQSR (note replacing PrintReads of GATK3)
   for i in \$(seq -f '%04g' 0 $(((nprs/4)-1)))
   do
      bqfile=\${mapFile%.bam}_dedup_recal_data_\$i.table
      output=\${mapFile%.bam}_dedup_recal_\$i.bam
      /usr/bin/time -v -o time_gatkApplyBQSR.log gatk \\
         --java-options "-Xmx4G " ApplyBQSR -R \$ref -I \${mapFile%.bam}_dedup.bam \\
         -L \$REF_HOME/intervals/$((nprs/4))c/\$i-scattered.interval_list  -bqsr \$bqfile \\
         --static-quantized-quals 10 --static-quantized-quals 20 --static-quantized-quals 30 -O \$output &
   done
   wait
done
  
step4=\$(date -u +%s)
elp_time=\$((\$step4 - \$step3))
Elapsed[\${#Elapsed[*]}]="ApplyBQSR is \$elp_time seconds"

# Gather Bam files from scattered ApplyBQSR
   bamFiles=()
   for i in \${mapFile%.bam}_dedup_recal_*.bam
   do
      bamFiles[\${#bamFiles[*]}]=\$i
   done
   for i in \${!bamFiles[*]}
   do
     if [ \$i == 0 ]
     then
       bamFilesArg="-I \${bamFiles[i]}"
     else
       bamFilesArg="\${bamFilesArg} -I \${bamFiles[i]}"
     fi
   done
   /usr/bin/time -v -o time_gatkGatherBamFiles.log gatk --java-options "-Xmx4G" GatherBamFiles \\
      \$bamFilesArg --CREATE_INDEX true \\
      -O \${mapFile%.bam}_dedup_recal.bam
step44=\$(date -u +%s)
elp_time=\$((\$step44 - \$step4))
Elapsed[\${#Elapsed[*]}]="gatherBAM files is \$elp_time seconds"

# VARIANT CALLING optimized for Power9 ##########
# Used Onle native pairhmm  # Run Haplotypecaller with split bam file
   for i in \$(seq -f '%04g' 0 $(((nprs/4)-1)))
   do
     infile=\${mapFile%.bam}_dedup_recal_\$i.bam
     outfile=\${mapFile%.bam}_dedup_recal_\$i.g.vcf
     /usr/bin/time -v -o time_gatkHaplotypeCaller.log gatk \\
        --java-options "-Xmx4G -Djava.library.path=\$GATK_HOME/gatk-4.1.0.0/libs -XX:+UseParallelGC -XX:ParallelGCThreads=1" HaplotypeCaller \\
        -R \${ref} -I \$infile -pairHMM VSX_LOGLESS_CACHING \\
        -L \$REF_HOME/intervals/$((nprs/4))c/\$i-scattered.interval_list \\
        --native-pair-hmm-threads 8 --smith-waterman FASTEST_AVAILABLE \\
        -O \$outfile -ERC GVCF -stand-call-conf 10 &
   done
   wait

step5=\$(date -u +%s)
elp_time=\$((\$step5 - \$step44))
Elapsed[\${#Elapsed[*]}]="Haplotypecaller is \$elp_time seconds"

# combinge gvcf
   gvcfFiles=()
   for i in \${mapFile%.bam}_dedup_recal_*.g.vcf
   do
      gvcfFiles[\${#gvcfFiles[*]}]=\$i
   done

   for i in \${!gvcfFiles[*]}
   do
     if [ \$i == 0 ]
     then
        gvcfFilesArg="-I \${gvcfFiles[i]}"
     else
        gvcfFilesArg="\${gvcfFilesArg} -I \${gvcfFiles[i]}"
      fi
   done
   /usr/bin/time -v -o time_gatkGatherVCFs.log gatk --java-options "-Xmx4G" GatherVcfs -R \$ref \$gvcfFilesArg -O \${vcfFile%vcf}g.vcf
   
step6=\$(date -u +%s)
elp_time=\$((\$step6 - \$step5))
Elapsed[\${#Elapsed[*]}]="ComnbineGVCFs is \$elp_time seconds"

  # genotype gvcf files
   for i in \$(seq -f '%04g' 0 $(((nprs/4)-1)))
   do
     /usr/bin/time -v -o time_gatkGenotypeGVCFs.log gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=4" GenotypeGVCFs \\
        -R \$ref -V \${mapFile%.bam}_dedup_recal_\$i.g.vcf \\
        -L \$REF_HOME/intervals/$((nprs/4))c/\$i-scattered.interval_list \\
        -O \${vcfFile%.vcf}_\$i.vcf &
   done
   wait
   vcfFiles=()
   for i in \${vcfFile%.vcf}_*.vcf
   do
      vcfFiles[\${#vcfFiles[*]}]=\$i
   done

   for i in \${!vcfFiles[*]}
   do
     if [ \$i == 0 ]
     then
        vcfFilesArg="-I \${vcfFiles[i]}"
     else
        vcfFilesArg="\${vcfFilesArg} -I \${vcfFiles[i]}"
      fi
   done
   /usr/bin/time -v -o time_gatkGatherMergedVCFs.log gatk --java-options "-Xmx4G" GatherVcfs -R \$ref \$vcfFilesArg -O \$vcfFile

step7=\$(date -u +%s)
elp_time=\$((\$step7 - \$step6))
Elapsed[\${#Elapsed[*]}]="GenotypeGVCFs is \$elp_time seconds"

# Delete splited VCF files
   rm -f \${mapFile%.bam}_dedup_recal_*.g.vcf* \${mapFile%.bam}_dedup_recal_*.ba*
   rm -f \${mapFile%.bam}_dedup_recal_data_*.table
   rm -f \${vcfFile%.vcf}_*.vcf

# VARIANT QUALITY SCORE RECALIBRATION ########
  /usr/bin/time -v -o time_gatkVariantRecalibratorSNP.log taskset -c 0-$((nprs-1)):4 gatk \\
       --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$((nprs/4))" \\
        VariantRecalibrator -V \$vcfFile -O \${workPath}/recalibrate_SNP.recal \\
       -mode SNP --tranches-file \${workPath}/recalibrate_SNP.tranches \\
       -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR -an MQ \\
       --max-gaussians 6 \\
       -resource:hapmap,known=false,training=true,truth=true,prior=15.0 \$vcfHapmap \\
       -resource:omni,known=false,training=true,truth=true,prior=12.0 \$vcfOmni \\
       -resource:1000G,known=false,training=true,truth=false,prior=10.0 \$vcfGlk \\
       -resource:dbsnp,known=true,training=false,truth=false,prior=7.0 \$vcfDbsnp
  
step8=\$(date -u +%s)
elp_time=\$((\$step8 - \$step7))
Elapsed[\${#Elapsed[*]}]="VariantRecalibratorSNP is \$elp_time seconds"

  # 2. Apply recalibration to SNPs
  /usr/bin/time -v -o time_gatkApplyVQSRSNP.log taskset -c 0-$((nprs-1)):8 gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$((nprs/8))" \\
       ApplyVQSR -V \$vcfFile -O \${workPath}/recalibrated_snps_raw_indels.vcf \\
       --recal-file \${workPath}/recalibrate_SNP.recal --tranches-file \${workPath}/recalibrate_SNP.tranches \\
       -truth-sensitivity-filter-level 99.5 --create-output-variant-index true -mode SNP
  
step9=\$(date -u +%s)
elp_time=\$((\$step9 - \$step8))
Elapsed[\${#Elapsed[*]}]="ApplyVQSR-SNP is \$elp_time seconds"

  # 3. Run Variant Recalibrator on Indels
  /usr/bin/time -v -o time_gatkVariantRecalibratorIndel.log taskset -c 0-$((nprs-1)):4 gatk \\
       --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$((nprs/4))" \\
       VariantRecalibrator -V \${workPath}/recalibrated_snps_raw_indels.vcf -O \${workPath}/recalibrate_INDEL.recal \\
       -mode INDEL --tranches-file \${workPath}/recalibrate_INDEL.tranches \\
       -tranche 100.0 -tranche 99.9 -tranche 99.0 -tranche 90.0 -an QD -an FS -an MQRankSum -an ReadPosRankSum -an SOR \\
       --max-gaussians 4 \\
       -resource:mills,known=false,training=true,truth=true,prior=12.0 \$vcfMills \\
       -resource:dbsnp,known=true,training=false,truth=false,prior=2.0 \$vcfDbsnp 
       
step10=\$(date -u +%s)
elp_time=\$((\$step10 - \$step9))
Elapsed[\${#Elapsed[*]}]="VariantRecalibratorIndel is \$elp_time seconds"

   # 4. Apply recalibration to INDELs  
   /usr/bin/time -v -o time_gatkApplyVQSRIndel.log taskset -c 0-$((nprs-1)):8 gatk --java-options "-Xmx4G -XX:+UseParallelGC -XX:ParallelGCThreads=$((nprs/8))" \\
       ApplyVQSR -V \${workPath}/recalibrated_snps_raw_indels.vcf -O \${vcfFile%.vcf}.recal.vcf \\
       --recal-file \${workPath}/recalibrate_INDEL.recal --tranches-file \${workPath}/recalibrate_INDEL.tranches  \\
       -truth-sensitivity-filter-level 99.0 --create-output-variant-index true -mode INDEL

step11=\$(date -u +%s)
elp_time=\$((\$step11 - \$step10))
Elapsed[\${#Elapsed[*]}]="ApplyVQSR-Indels is \$elp_time seconds"
elp_time=\$((\$step11 - \$start_time))
Elapsed[\${#Elapsed[*]}]="Total time for this pipeline is \$elp_time seconds"
echo "Finished pipeline at \$(date)"

for i in \${!Elapsed[*]}
do
echo \${Elapsed[\$i]}
done

EOF

chmod +x $filename

