<p align="right">
  <img src="https://robert-koch-institut.github.io/SARS-CoV-2-Infektionen_in_Deutschland/assets/RKI_Logo.png" style="width: 204px; height: 60px;">
</p>
<!-- HEADER_START: {"lang": "en"} -->


Documentation  

# Dynamic Transmission Model of Meningococcal Adolescent Booster Vaccination

<br> 
<br> 
<br> 

[**Felix Günther**](https://orcid.org/0000-0001-6582-1174)&sup1;

<br> 



&emsp;&emsp;&sup1; [Robert Koch Institute](https://www.rki.de/en) | [Unit 33](https://www.rki.de/fg33-en)

<br> 

**Cite**  
Günther, F. (2026). Dynamic Transmission Model of Meningococcal Adolescent Booster Vaccination. Zenodo. [https://doi.org/10.5281/zenodo.18161773](https://doi.org/10.5281/zenodo.18161773)


<br>

**Abstract**    
Implementation of a dynamic transmission model for modelling meningococcal transmission in Germany at the carrier level, as well as estimating the effects of different vaccination strategies. The software makes it possible to reproduce the results of the manuscript "The impact of introducing meningococcal C/ACWY booster vaccination among adolescents in Germany: a dynamic transmission modelling study".

<br>

**Table of Content**  

<!-- TOC_START: {"heading_depth": 2} -->
  - [Project information](#project-information)  
  - [Installation](#installation)  
  - [Code and structure of repository](#code-and-structure-of-repository)  
  - [Collaborate](#collaborate)
<!-- TOC_END -->

<br>

<!-- HEADER_END -->
---  

This repository contains code for the manuscript Günther, Felix, et al. *"The impact of introducing meningococcal C/ACWY booster vaccination among adolescents in Germany: a dynamic transmission modelling study"* with which you can reproduce the manuscript's results and figures. A preprint of the manuscript is available on [medRxiv](https://www.medrxiv.org/content/10.1101/2024.12.19.24319393.abstract).

## Project information

This code was developed at the Robert Koch Institute as part of the project *AMSEC* funded by the G-BA. More information on the project (in German) can be found on the website of the [G-BA Innovationsfonds](https://innovationsfonds.g-ba.de/projekte/versorgungsforschung/amsec.214).

### Administrative and organizational information

This work was conducted by staff from [Unit 33 | Immunization/STIKO](https://www.rki.de/fg33-en). The publication of the code as well as the quality management of the metadata is done by department [MF 4 | Domain Specific Data and Research Data Management](https://www.rki.de/mf4-en). Questions regarding the publication infrastructure can be directed to the Open Data Team of the Department MF4 at [OpenData@rki.de](mailto:OpenData@rki.de).

### Motivation

In Germany, routine vaccination against invasive meningococcal disease (IMD) serogroup C and B is recommended for infants. Due to a second peak of incidence among adolescents, we developed a dynamic transmission model (DTM) for estimating the effectiveness and efficiency of introducing a meningococcal C or ACWY booster vaccination among adolescents.

## Installation

Our scripts are written in R. Please make sure you have the R programming language installed. The analysis is implemented based on a [{targets}](https://books.ropensci.org/targets/) pipeline and utilizes R-packages that are defined in the _targets.R file.

## Code and structure of repository

This repository contains an implementation of the developed dynamic transmission model (DTM) as well as code for the simulation of vaccination strategies conditional on estimated parameters from the model calibration. Simulation and estimation and summary of results are implemented based on a [{targets}](https://books.ropensci.org/targets/) pipeline. The repository is structured as follows: 

- *./_targets.R* is the target script file that configures and defines the analysis pipeline.
- Computation of the pipeline can be sourced by running the *./run_pipeline.R* script on a local computer or in a high-performance computing environment. Parallelization is implemented via base R functionality (within single targets, up to 10 cores per target) and the *crew* package (for parallel computation of multiple targets). The current specification expects availability of >30 cores for parallel computation, but this number of cores may be adjusted in the *_targets.R* file.
- The folder *./R/* contains multiple files of R-Code consisting of custom functions used for defining the single targets, e.g. functions for defining the DTM, data processing, simulation tasks, and plotting:
  - *ode_funs.R* contains code defining the differential equation-based DTM
  - *functions.R* contains custom functions for data processing or the definition of single analysis steps (targets) of the pipeline, like performing simulations or plotting
  - *ve_estimation.R* contains code for estimation of the VE parameters from external data
- The folder *./data/* contains input data for the simulations, mostly results of the model calibration extracted from a larger pipeline.
- The folder *./results/* contains output results of the simulations from running the pipeline, including Fig. 3 and 4 of the manuscript (results of the main simulations) and supplementary figures and tables related to the main analysis.
- The folder *./_targets/* is a storage for intermediate results from evaluating the analysis pipeline (following conventions of the targets package). These intermediate results are stored on the machine used to run the pipeline and are not committed to the github repository due to their file size (see the *.gitignore* file).


### Data

The code in this repository does not rely on any external data. All required data and input objects are contained in the ./data/ subfolder.

<!-- FOOTER_START: {"lang": "en"} -->
## Collaborate

If you want to participate in our project, feel free to fork this repo and send us pull requests. To make sure everything is working please use pre-commit. It will run a few tests and lints before a commit can be made. To install pre-commit, run

``pre-commit install``

### Metadata

To increase findability, the provided code is described with metadata. Versioning and DOI assignment are performed via [Zenodo.org](https://zenodo.org). The metadata prepared for import into Zenodo are stored in the [zenodo.json](https://github.com/robert-koch-institut/Meningococcal_BoosterVacc_DTM/blob/main/Metadata/zenodo.json). Documentation of the individual metadata variables can be found at https://developers.zenodo.org/representation.

> [Metadata/zenodo.json](https://github.com/robert-koch-institut/Meningococcal_BoosterVacc_DTM/blob/main/Metadata/zenodo.json)

### Publication platforms

This software publication is available on [Zenodo.org](http://Zenodo.org/), [GitHub.com](http://GitHub.com/) and [OpenCoDE](https://gitlab.opencode.de):  

- https://zenodo.org/communities/robertkochinstitut  
- https://github.com/robert-koch-institut  
- https://gitlab.opencode.de/robert-koch-institut


### License

The "Dynamic Transmission Model of Meningococcal Adolescent Booster Vaccination" code is licensed under the [MIT License](https://mit-license.org/). This means, that the code provided in the repository is freely available, with the condition of attributing the Robert Koch Institute as the source, for anyone to process and modify, create derivatives of the code and use it for commercial and non-commercial purposes. Further information about the license can be found in the [LICENSE](https://github.com/robert-koch-institut/Meningococcal_BoosterVacc_DTM/blob/main/LICENSE) file of the repository.  
<!-- FOOTER_END -->
