##  snakemake pipeline for ScRNA-seq virus integration detection

#### This pipeline is modified from nextflow pipeline [virus integration detection](https://nf-co.re/viralintegration/0.1.1/)

### requirements
- conda
- snakemake
- cellranger

### installation
```bash
conda env create -f ScRNA_HPV_integration_detection.yml
conda activate ScRNA_HPV_integration_detection
```
### usage 
#### fill all missing parameters in the Snakefile and run it.

### citation
If you use this pipeline, please cite the following paper:
```bibtex
@article {Li2025.01.17.633627,
	author = {Li, Shiting and Xia, Shaomiao and Lawas, Maria and Kulshreshtha, Aishani and Garb, Bailey F. and Perera, AA Chamila and Li, Chen and Qin, Tingting and Welch, Joshua D. and D{\textquoteright}Silva, Nisha J. and Rozek, Laura S. and Sartor, Maureen A.},
	title = {HPV integration in head and neck cancer: downstream splicing events and expression ratios linked with poor outcomes},
	elocation-id = {2025.01.17.633627},
	year = {2025},
	doi = {10.1101/2025.01.17.633627},
	URL = {https://www.biorxiv.org/content/early/2025/01/20/2025.01.17.633627},
	eprint = {https://www.biorxiv.org/content/early/2025/01/20/2025.01.17.633627.full.pdf},
	journal = {bioRxiv}
}
```
