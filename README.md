# BPGA_downstream
Downstream process for BPGA output file. To make user-friendly xlsx file

# Welcome to the BPGA downstream, a combined scripts-based pipeline for comparative genomics analysis.


## Project Introduction

Welcome to our GitHub repository, where we're excited to share a series of workflows designed to streamline processes in systems biology. This repository is composed of various scripts, each tailored to specific tasks within our broader research framework. Additionally, we're providing access to a curated database to enhance your research capabilities.
- Project lead: Jae-Yoon Sung
- Maintainers: Jae-Yoon Sung
- Contributors: Jae-Yoon Sung

## Key Features
Diverse Scripts: Our collection includes a range of scripts, each developed to address unique challenges in systems biology research.

Curated Database: 
User-Friendly Documentation: Detailed documentation is available to guide you through the installation, setup, and utilization of both the scripts and the database.

### Algorithms for analysis


### Getting Started:
To begin using our resources, please follow the steps outlined in our documentation. 
Whether you're looking to integrate our scripts into your existing projects or explore our database for new insights, we've provided all the necessary instructions to get you started.

### Requirements

The BPGA downstream is supported for macOS, Linux and Windows machines, which can provide an environment for using R.
It requires R version >=4.2.1 for release, and R version >=4.3 for devel.

To download and install BPGA, see the [BPGA website]([https://iicb.res.in/bpga/]).

To download and install R, see the [R-project website](https://www.r-project.org/).


#### Warning
The basic file for genomic analysis, known as a GenBank file, requires both sequence and annotation in full-format files such as gbff, gb, or gbk. Additionally, GenBank prefers a format based on the GeneMarkS2+ pipeline, and using a different annotation pipeline to obtain GenBank files may lead to errors.



## Anaylsis flow


```r
setwd([GenBank directory])
if (!requireNamespace("devtools", quietly = TRUE))
    install.packages("devtools")
    
devtools::install_github("JAEYOONSUNG/BPGA_downstream")
```
   - **Note:** ............

**BPGA**

```python
emapper.py --cpu 20 --mp_start_method forkserver --data_dir [eggnog_data directory] -o out --output_dir [eggnog_output] --temp_dir [eggnog_output] --override -m diamond --dmnd_ignore_warnings --dmnd_algo ctg -i [fasta] --evalue 0.001 --score 60 --pident 40 --query_cover 20 --subject_cover 20 --itype proteins --tax_scope auto --target_orthologs all --go_evidence non-electronic --pfam_realign none --report_orthologs --decorate_gff yes --excel

```

- **Note:** BPGA result


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
[BGPA]  v2: functional annotation, orthology assignments, and domain 
                   prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
                   Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
                   Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293

[eggNOG-mapper v2] eggNOG-mapper v2: functional annotation, orthology assignments, and domain 
                   prediction at the metagenomic scale. Carlos P. Cantalapiedra, 
                   Ana Hernandez-Plaza, Ivica Letunic, Peer Bork, Jaime Huerta-Cepas. 2021.
                   Molecular Biology and Evolution, msab293, https://doi.org/10.1093/molbev/msab293


```
