[![GATK v4.1.4.1](https://img.shields.io/badge/gatk%20source-4.1.4.1-green.svg)](https://github.com/broadinstitute/gatk/archive/4.1.4.1.tar.gz)[![License (3-Clause BSD)](https://img.shields.io/badge/license-BSD%203--Clause-blue.svg)](https://opensource.org/licenses/BSD-3-Clause)[![IBM Power9 Download](https://img.shields.io/badge/power9-download-blue.svg)](https://ibm.box.com/v/gatk-power4141)
# GATK4-Power 
Supporting GATK4 on IBM Power architecture, including enablement, performance tuning and optimization.
* Latest version 
  v4.1.4.1 built on 06-15-2019
* Older version
  See other branches
## This includes:
* GATK4 package built on IBM Pwer9 
  GATK4 on Power9 full package (downloadable) is available upon request.
* Benchmarking scripts to run GATK4 pipelines on IBM Power system
## Install GATK4 on Power9
* Clone GATK4-Power
 ```git clone ```
* Run installation script to complete the installation
 ```./install_gatk4-power9.x```
* Generate runscript (Germline)
 ```./run_gatk-power.x```
## Create docker image
* Run docker script and upload to your system
 ``` docker run -t gatk:4.1.4.1 ``` 
## Performance
<img src="https://github.com/ruzhuchen/NGS/blob/master/images/p9_performance.png" alt="GATK4 performance on Power9" title="GATK4 performance on Power9" width="949px">

Contact: Ruzhu Chen (ruzhuchen@us.ibm.com) for support

* Note to user: The Power9 download is password protected. Please contact me for more info.
