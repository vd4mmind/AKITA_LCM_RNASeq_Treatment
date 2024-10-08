# Manuscript Title : Transcriptomics of SGLT2-positive early proximal tubule segments in mice: response to type 1 diabetes, SGLT1/2 inhibition or GLP1 receptor agonism

# Authors
Young Chul Kim<sup>1</sup> ,<sup>2</sup> *, Vivek Das<sup>3</sup>, Sadhana Kanoo<sup>1</sup> ,2</sup>, Huazhen Yao<sup>4</sup>, Stephanie M. Stanford<sup>5</sup>, Nunzio Bottini<sup>5</sup>, Anil Karihaloo<sup>7</sup>, Volker Vallon<sup>1</sup>,<sup>2</sup>,<sup>6</sup>  *

<sup>*</sup> Contributed equally

<sup>1</sup> Division of Nephrology & Hypertension, Department of Medicine, University of California San Diego, La Jolla, CA, USA

<sup>2</sup>  VA San Diego Healthcare System, San Diego, CA, USA

<sup>3</sup>  Novo Nordisk A/S, Søborg, Denmark

<sup>4</sup>  Institute for Genomic Medicine, University of California San Diego, La Jolla, CA, USA 

<sup>5</sup>  Division of Rheumatology, Allergy & Immunology, Department of Medicine, University of California San Diego, La Jolla, CA, USA

<sup>6</sup>  Department of Pharmacology, University of California San Diego, La Jolla, CA, USA

<sup>7</sup>  Novo Nordisk 33 Hayden Ave, Lexington, MA 02421 USAUSA


# Data and code information

## AKITA_LCM_RNASeq_Treatment
This study includes 54 mouse kidney samples stratified into 9 groups based on disease, treatment, and genotype explanatory variables. The table below summarizes the experimental design.
&nbsp;

| **Disease** | **Treatment** | **Genotype** | **Samples** |
|---|---|---|---|
| Non-Diabetic | Vehicle | WT | 6 |
| Non-Diabetic | Dapagliflozin | WT | 6 |
| Non-Diabetic | Vehicle | KO | 6 |
| Non-Diabetic | Dapagliflozin | KO | 6 |
| Diabetic | Vehicle | WT | 6 |
| Diabetic | Dapagliflozin | WT | 6 |
| Diabetic | Semaglutide | WT | 6 |
| Diabetic | Vehicle | KO | 6 |
| Diabetic | Dapagliflozin | KO | 6 |

The full metadata is available in the [manifest](https://eu.sbgenomics.com/u/novo-nordisk/akita-ckd-rnaseq-with-treatment/files/633eeb4fb17faa45cbbfafcf/) file.

## Data

The raw sequencing data is available from the volume directory: ``Bioinformatics_dk/RNAseq/VVDA_mouse_CKD_AKITA_treatment``

## Results

The table below lists all of the results files

&nbsp;

| File | Description |
|---|---|
| Data/ | Input FASTQ files |
| Genome/ | Reference genome files |
| workspace_EDA/ | Exploratory data analysis and results |
| workspace_RING/ | RING Rdata object files for differential analysis using DESeq2 |
| Manifest.csv | Manifest file |

*Note: The differential expression analysis was first performed with all samples included in the DESeq model. After inspecting the exploratory data analysis plots, it was clear some groups had much higher within-group variation than others. This would affect the performance of the DESeq model, so subsequent comparisons were made using only the samples within the respective groups which were being compared. See the [FAQ section](https://bioconductor.org/packages/release/bioc/vignettes/DESeq2/inst/doc/DESeq2.html#if-i-have-multiple-groups-should-i-run-all-together-or-split-into-pairs-of-groups) of the DESeq2 manual for more information.*

## Support and Acknowledgement

Vivek Das from Novo Nordisk A/S executed the preprocessing steps with the consultancy Zifo RnD Solutions (Riya Saju, James Ashmore) using in-house Seven Bridges pipeline.

For queries Vivek Das: vvda@novonordisk.com
