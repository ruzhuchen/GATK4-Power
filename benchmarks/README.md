# Use This Scripts to Run GATK4 on Power9
Note: you may need to have advance toolchains version 10.0 installed on your system (e.g., /opt/at10.0)

This script is optimized for Power9 only.
* Please download reference files directly from Broad Institute's ftp
* To use the script, please modify the path accordingly
   Please set your system parameters if running on Spark mode:
   ```
   ulimit -n 20480
   ulimit ps unlimited
   ```
   
## Setting PATH to `$GATK_HOME/bin`
   ```
   export GATK_HOME=your-path-to-this-directory
   export PATH=$GATK_HOME/bin:$PATH;
   export GATK_LOCAL_JAR=$GATK_HOME/gatk-4.0.11.0/libs/gatk.jar
   export GATK_SPARK_JAR=$GATK_HOME/gatk-4.0.11.0/libs/gatk-spark.jar
   ```
## JAVA Options
   `--java-options "-Xmx4G"` #increase as need! if running in Spark mode, use `"-Xmx8G"` and up
## Using IBM native pairhmm library
   `--java-options "-Xmx4G -Djava.library.path=$GATK_HOME/libs"`
   * For HaplotypeCaller, here is the example script:
   ```
   /usr/bin/time -v -o time_gatkHaplotypeCaller.txt gatk --java-options "-Xmx4G -Djava.library.path=$GATK_HOME/libs" \
         HaplotypeCaller -R ${ref} -I your_dedup_recal.bam \
         -O your_dedup_recal.g.vcf -ERC GVCF -stand-call-conf 10 --native-pair-hmm-threads 20
   ```
## Use `./gatk --list` for available functions

contact: ruzhuchen@us.ibm.com

