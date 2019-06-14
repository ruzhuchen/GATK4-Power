# GATK4 Version 4.1.2.0 on Power9 with Native PairHMM
The pipeline executes from BWA MEM to ApplyVQSR (SNPs and Indels). 
* BWA
* Samtools
* sam2bam (optional)
The above tools available at `$GATK_HOME/bin` dir

## To use native PairHMM:

add ```--native-pair-hmm-threads``` to HaplotypeCaller and Mutect2

use ```--java-options "-Djava.library.path=$GATK_HOME/gatk-4.1.2.0/libs" ``` to run HaplotypeCaller and Mutect2

## This script is optimized for Power9 only. 
* Please download reference files directly from Broad Institute's ftp
* To use the script, please modify the path accordingly

   Please set your system parameters if running on Spark mode:
   ```ulimit -n 20480```
   ```ulimit ps unlimited```
* Contact ruzhuchen@us.ibm.com if you have any question,
