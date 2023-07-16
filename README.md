# EasyDEA: An Automated Tool for Differential Expression Analysis

### *Adrian Barreno Sánchez - MSc Computationa Biology (UPM)*

#### 1) What is _EasyDEA_?
_EasyDEA_ is a command-line tool that performs automated Differential Expression Analysis (DEA) from RNA sequencing data. 
_EasyDEA_ emerged as a result of regular colaboration with wet-lab scientists, who recognize the challenges faced by non-expert users in conducting accurate transcriptional analysis. By addressing these limitations, _EasyDEA_ provides an accessible solution for researchers without prior programming experience and statistical knowledge to perform DEA effectively.

#### 2) _EasyDEA_ workflow
_EasyDEA_ combines different statistical approaches and functionalitites from some of the most widespread tools for DEA, including _edgeR_, _limma_ and _DESeq2_. The user can therefore decide whether to follow either _edgeR-limma_ or _DESeq2_ pipelines, or both, to efficiently identify significant differentially expressed genes. To achieve this, _EasyDEA_ performs all possible pairwise comparisons of the samples based on the design variable defined by the user.

Besides this, _EasyDEA_ performs commonly used downstream analyses of the differentially expressed genes to obtain a broad understanding of the biological and molecular processes regulated under varying experimental conditions. With this aim, _EasyDEA_ incorporates several [Bioconductor](https://bioconductor.org/) packages to perform Gene Annotation and Gene Set Enrichment Analysis (GSEA), fetching up-to-date information from online databases. Alternatively, the user can provide GTF/GFF files of the target organism to generate in-home offline annotation.

_EasyDEA_ generates multiple output files throughout the analysis, including raw and annotated lists of differentially expressed genes for each comparison, filtered by Log 2-fold change and significance thresholds. Moreover, several diagnostic plots are produced that allows the detection of potential issues in the data and the analysis workflow. _EasyDEA_ also generates an HTML report to facilitate the exploration and interpretability of the GSEA results.

<br/>
<img src="https://github.com/abarrenos/EasyDEA/assets/113832131/a90d21b3-ea40-48bd-89a6-0ebfbe510ab0" alt="EasyDEA Workflow" width="500"/>
<br/>

_Image generated with Biorender.com_

#### 3) How to run _EasyDEA_?
To run the program, clone or download the repository onto your local device, within the desired installation path. To initiate the analysis, execute the main script ```easyDEA.R``` followed by custom command-line options and parameter values. Alternatively, users can provide as input a configuration text file containing a list of options and parameter values for the analysis. For more information about the available program options and default parameters, run the program followed by the option ```--help```.

```
./EasyDEA.R [--OPTIONS] [VALUES] [-h/--help] --config-file [filepath]
```

Even though _EasyDEA_ contains a set of default parameters for the analysis pipeline, users can define custom default values in a text file and save it with the name ```.EasyDEA.rc```. This configuration file can be located either at the program installation folder or within the user's HOME directory ```~/```, allowing the program to automatically recognize it. Importantly, command-line options will take precedence over file-defined options.

An example configuration file containing all available command line options and parameter values can be found [here](data/config_file.txt).

<br/><hr/>

### Contact details
For more information regarding _EasyDEA_ functioning and  problem reporting, please do not hesitate to contact me.
- [adrian.barreno@alumnos.upm.es](mailto:adrian.barreno@alumnos.upm.es)
- [linkedIn](https://www.linkedin.com/in/adri%C3%A1n-barreno-s%C3%A1nchez-890a45206/)


<img src="https://github.com/abarrenos/EasyDEA/assets/113832131/de17de66-9cb2-418e-9b3a-8ccb34d78616" alt="Logo CNB" height="90"/>
<img src="https://github.com/abarrenos/EasyDEA/assets/113832131/add67521-e92e-4f57-90b5-06c6e6b6a29c" alt="Logo CNB" height="90"/>
