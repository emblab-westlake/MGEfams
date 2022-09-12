

<!--
 * @Author: zhangguoqing and Feng Ju (supervisor) 
 * @Correspondence email: jufeng@westlake.edu.cn
 * @Affiliation: EMBLab, Westlake University
 * @Date: 2022-04-15 07:00:00 (by FJ)
 * @LastEditTime: 2021-04-22 20:50:49
-->


# MGEfams
MGEfams: Fast and robust identification of Mobile Genetic Element (MGE) genes from genomic and metagenomic assemblies of high-throghput DNA sequencing. The MGEfams is a high-quality, manually cruated and structured subdatabase of profile Hidden Markov Models (HMM) for MGE gene annotation


This subdatabase consists of HMM models for MGE genes extracted from Pfam and TIGRFAMs databases, based on string match and expert manual cruation of their functional annotations to one of the following keywords: transposase, transposon, conjugative, integrase, integron, recombinase, resolvase, conjugal, mobilization, recombination, and plasmid, as recommended previously (Science.2012;337(6098):1107-1111; Nature.2014; 509(7502): 612–616)

---
### Development Record

MGEfams_v0.1, v0.2, v0.3 and v0.4 - created by Feng Ju in Oct 2018 

MGEfams_v0.3 and v0.4 - created by Feng Ju before Oct 2019 

MGEfams_v0.5 - created by Guoqing Zhang in April 2021 based on Pfam (v34.0) and TIGRFAMs (v15.0)

MGEfams_v0.6 -  In development  

---

### Update Notes
MGEfams_v5.0 was constructed based on the lastest version of Pfam (v34.0) and TIGRFAMs (v15.0), and will be regularly updated.



### Dependence
python: >=3.6  
HMMER3: >=3.3.2



### Usage

DESCRIPTION
MGEgenetyper version: 0.5.0
Detailed introducion

**optional arguments:**  

  **-i INPUT_FILE, --input INPUT_FILE**  
    the input file ORFs.faa  

  **-o [OUTPUT_FILE_NAME], --output [OUTPUT_FILE_NAME]**  
    the outputfile prefix name: eg. PRIFEX.tlout  


Required arguments:  
  **-db Sub_Database**  
  Sub database; Default Antibiotic Resistance Genes  

  **-DB synthesis_database**  
    synthesis database, Default Antibiotic Resistance Genes database  

  **{--cut_ga,--cut_nc,--cut_tc}**  
    hmm type; chose from {--cut_ga, --cut_nc, --cut_tc default:cut_ga  

  **-n N, --nproc N**  
  The number of CPUs to use for parallelizing the mapping [default 1]  


**Other arguments:**   

  **--check**  
   Only checks if the Default ARG DB is installed and installs it if not.

  **-v, --version**  
    Prints the current MetaPhlAn version and exit

  **-h, --help**  
    show this help message and exit

**Example:**


```
MGEfams-v0.5.py -i <protein_c_ORFs>.fasta -o <OUTPUT_NAME> -db MGEs.v0.5.hmm -DB Pfam-TIGRfams.hmm -n 2
```


Citation: He L, Huang X, Zhang G, Yuan L, Shen E, Zhang L, et al. Distinctive signatures of pathogenic and antibiotic resistant potentials in the hadal microbiome. Environmental Microbiome. 2022;17(1):19.  https://environmentalmicrobiome.biomedcentral.com/articles/10.1186/s40793-022-00413-5
