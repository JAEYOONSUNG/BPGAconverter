# Welcome to the BPGAconverter, a combined scripts-based pipeline for comparative genomics analysis.
Downstream processing of BPGA output files to generate user-friendly XLSX files and publication-ready figures.

<img width="1356" alt="table" src="https://github.com/user-attachments/assets/c69f0733-dcbf-46be-867c-6fd7502ef9c7">
<img width="1155" alt="plot" src="https://github.com/user-attachments/assets/a7e04b45-60c3-4d8b-9baa-0fd795580e54">


## Project Introduction

Welcome to our GitHub repository, where we're excited to share a series of workflows designed to streamline processes in systems biology. This repository is composed of various scripts, each tailored to specific tasks within our broader research framework. Additionally, we're providing access to a curated database to enhance your research capabilities.
- Project lead: Jae-Yoon Sung
- Maintainers: Jae-Yoon Sung
- Contributors: Jae-Yoon Sung

## Key Features
- The BPGA pipeline encounters challenges when retrieving coding sequences (CDS) using protein_id, making it difficult to map back to the original sequence.
- To address this, all CDS entries are extracted and matched using strain-specific locus_tag. Subsequently, the number of strains analyzed and their CDS are reassigned.
- This approach of pooling and redistributing allows for restructuring the results obtained through BPGA's usearch based on desired columns of comprehensive datasets. 
- By default, the table is organized with a combination of locus_tag, product (annotation), and sequence identity (%) columns.
- Additionally, BPGAconverter enhances analytical convenience by allowing the matching of various types of information beyond just the product of the original CDS.

Curated Database: 
User-Friendly Documentation: Detailed documentation is available to guide you through the installation, setup, and utilization of both the scripts and the database.


## Getting Started:
To begin using our resources, please follow the steps outlined in our documentation. 
Whether you're looking to integrate our scripts into your existing projects or explore our database for new insights, we've provided all the necessary instructions to get you started.

### Requirements

The BPGAconverter is supported for macOS, Linux and Windows machines, which can provide an environment for using R.
It requires R version >=4.2.1 for release, and R version >=4.3 for devel.

To download and install a Bacterial Pan-Genome Analysis pipeline (BPGA), see the [BPGA website]([https://iicb.res.in/bpga/]).
- note: BPGA is supported for only Windows machine 

To download and install R, see the [R-project website](https://www.r-project.org/).

### Installation
```r
setwd([GenBank directory])
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("JAEYOONSUNG/BPGAconverter")
```



## Anaylsis flow

#### Warning
- The basic file for genomic analysis, known as a GenBank file, requires both sequence and annotation in full-format files such as gbff, gb, or gbk. Additionally, GenBank prefers a format based on the GeneMarkS2+ pipeline, and using a different annotation pipeline to obtain GenBank files may lead to errors.
- Processing more than 50 strains may cause errors depending on your computer’s processing power and R’s data capacity.

# Quick start
## Run BPGA downstream analysis
```r
setwd([Analysis directory]) # Set the working directory to the directory containing folders named GenBank, BPGA Supporting_files, and Results. The GenBank folder must include files with extensions such as .gbff, .gb, or .gbk.
library(BPGAconverter)
run_BPGAconvert() # for generating reconstructed table
run_BPGAplot() # for generating figure
```

### BPGA
**Caution:**
- This package requires GenBank files annotated using **<mark>GeneMarkS-2+**</mark>. GenBank files from other annotation pipelines may cause unexpected errors due to differences in file structure.
- The BPGA result folder must be located in the same working directory as the GenBank files.
- The BPGA program extracts genome names from the ORGANISM line, located directly below the SOURCE line in the GenBank file (when viewed in text format). To avoid confusion when matching columns with the same name due to missing strain names, it is recommended to add strain names to the ORGANISM line before running BPGA.
- BPGA supports multiple aligners, including Usearch, CD-HIT, and OrthoMCL. However, **<mark>BPGAconverter is currently designed to process results obtained specifically through Usearch analysis</mark>**. Support for CD-HIT and OrthoMCL will be added in future updates.


### EggNOG [Optional]
```python
emapper.py --cpu 20 --mp_start_method forkserver --data_dir [eggnog_data directory] -o out --output_dir [eggnog_output] --temp_dir [eggnog_output] --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg -i [fasta] --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel

```
- **Note:** [Strain of interest].emapper.annotations.xlsx or [Strain of interest]emapper.annotations.csv output are used for merging data.
  
## Contributing
We welcome contributions from the community! If you have suggestions for improvements, additional scripts, or updates to the database, please see our contributing guidelines for more information on how to get involved.


## License
This project is released under MIT licence, which allows for both personal and commercial use, modification, and distribution of our work, provided that proper credit is given.

We hope our resources will prove invaluable to your research in systems biology. For any questions or feedback, please don't hesitate to reach out through our GitHub issues or contact section.

## Citation
If you use this piepline, please cite:
```
[DNMB] DNMB: Accelerating the Domestication of Non-model Thermophilic Microorganisms Geobacillus stearothermohpilus as a Thermophilic Platform Cell.
             Jae-Yoon Sung, Hyungbin Kim, Seong Do Kim, Sang Jae Lee, Seong Bo Kim, and Dong-Woo Lee. 2024.
             XXX, XXX, https://doi.org/XXX
```
Please, cite also the underlying algorithm if it was used for the search step of DNMB:
```
[BGPA] : BPGA- an ultra-fast pan-genome analysis pipeline
         Narendrakumar M. Chaudhari, Vinod Kumar Gupta & Chitra Dutta.
         Scientific Reports 6, 24373 (2016). https://doi.org/10.1038/srep24373

[eggNOG-mapper v2] eggNOG-mapper v2: functional annotation, orthology assignments, and domain 
                   prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
                   Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
                   Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293


```
