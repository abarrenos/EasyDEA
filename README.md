# EasyDEA: An Automatized Tool for Differential Expression Analysis

### *Adrian Barreno Sánchez - MSc Computationa Biology (UPM)*

#### 1) What is _EasyDEA_?
_EasyDEA_ is a command-line tool that performs automated Differential Expression Analysis (DEA) from RNA sequencing data. 
_EasyDEA_ emerged as a result of regular colaboration with wet-lab scientists, who recognize the challenges faced by non-expert users in conducting accurate transcriptional analysis. By addressing these limitations, _EasyDEA_ provides an accessible solution for researchers without prior programming experience and statistical knowledge to perform DEA effectively.

#### 2) _EasyDEA_ workflow
_EasyDEA_ combines different statistical approaches and functionalitites from some of the most widespread tools for DEA, including _edgeR_, _limma_ and _DESeq2_. The user can therefore decide whether to follow either _edgeR-limma_ or _DESeq2_ pipelines, or both, to efficiently identify significant differentially expressed genes. To achieve this, _EasyDEA_ performs all possible pairwise comparisons of the samples based on the design variable defined by the user.

Besides this, _EasyDEA_ performs commonly used downstream analyses of the differentially expressed genes to obtain a broad understanding of the biological and molecular processes regulated under varying experimental conditions. With this aim, _EasyDEA_ incorporates several [Bioconductor](https://bioconductor.org/) packages to perform Gene Annotation and Gene Set Enrichment Analysis (GSEA), fetching up-to-date information from online databases. Alternatively, the user can provide GTF/GFF files of the target organism to generate in-home offline annotation.

*EasyDEA* generates multiple output files throughout the analysis, including raw and annotated lists of differentially expressed genes for each comparison, filtered by Log 2-fold change and significance thresholds. Moreover, several diagnostic plots are produced that allows the detection of potential issues in the data and the analysis workflow. _EasyDEA_ also generates an HTML report to facilitate the exploration and interpretability of the GSEA results.

<br/>
<img src="https://github.com/abarrenos/EasyDEA/assets/113832131/a90d21b3-ea40-48bd-89a6-0ebfbe510ab0" alt="EasyDEA Workflow" width="500"/>
<br/>

_Image generated with Biorender.com_

#### 3) How to run _EasyDEA_?

