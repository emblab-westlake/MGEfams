

<!--
 * @Author: zhangguoqing
 * @Date: 2021-04-21 14:37:48
 * @LastEditTime: 2021-04-22 19:23:16
-->

# MGEfams
Fast and robust identificaiton of mobile genetic elements (MGE) genes from genomes and metagenome assemblies using MGEfams, a high-quality and manually cruated subdatabase of profile hidden Markov models for MGE genes


This subdatabase of Hidden Markov Models (HMMs) consists of HMM models for MGE genes extracted from Pfam (v 30.0) 27 and TIGRFAMs databases, based on string match in their functional annotations to one of the following keywords: transposase, transposon, conjugative, integrase, integron, recombinase, resolvase, conjugal, mobilization, recombination, and plasmid, as recommended previously (Forsberg et al 2012, Forsberg et al 2014)

---
### Development Record  

MGEfams_v0.1, v0.2, v0.3 and v0.4 - created in Oct 2018 

MGEfams_v0.3 and v0.4 - created before Oct 2019 

MGEfams_v0.5 - created in April 2021 based on Pfam (v34.0) and TIGRFAMs (v15.0)

MGEfams_v0.6 in development  

---

### Update Notes
MGEfams_v5.0 18 modules were updated by lastest Pfam and TIGRfams



### Dependence
python: >=3.6  
HMMER3: >=3.3.2



### Usage
`
-i INPUT_FILE [-o [OUTPUT_FILE_NAME]] [-db Sub_Database]
                  [-DB synthesis_database] [-n N] [--check] [-v] [-h]
                  [{--cut_ga,--cut_nc,--cut_tc}]


DESCRIPTION
FunGeneTyper version: 1.0.0
Detailed introducion
Usage


optional arguments:
  -i INPUT_FILE, --input INPUT_FILE
                        the input file ORFs.faa
  -o [OUTPUT_FILE_NAME], --output [OUTPUT_FILE_NAME]
                        the outputfile prefix name: eg. PRIFEX.tlout


Required arguments:
  -db Sub_Database      Sub database; Default Antibiotic Resistance Genes
  -DB synthesis_database
                        synthesis database, Default Antibiotic Resistance Genes database
  {--cut_ga,--cut_nc,--cut_tc}
                        hmm type; chose from {--cut_ga, --cut_nc, --cut_tc}
                         default:cut_ga
  -n N, --nproc N       The number of CPUs to use for parallelizing the mapping [default 1]


Other arguments:
  --check               Only checks if the Default ARG DB is installed and installs it if not.
  -v, --version         Prints the current MetaPhlAn version and exit
  -h, --help            show this help message and exit
`
```


MGEsgenetyper-v0.5.py -i <protein_c_ORFs>.fasta -o <OUTPUT_NAME> -db MGEs.v0.5.hmm -DB Pfam-TIGRfams.hmm -n 2
```
