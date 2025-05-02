[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7827225.svg)](https://doi.org/10.5281/zenodo.7827225)

# flywalk
Code for a manuscipt of the flywalk.<br>
Takuto Okuno, Alexander Woodward, Hideyuki Okano, Junichi Hata (in submisstion)
["Functional connectivity, structural connectivity, and inter-individual variability in Drosophila melanogaster"](https://www.yahoo.com/)

## Requirements: Software
* MATLAB R2021b or later
* Parallel Computing Toolbox ver7.1 or later


## Installation
1. Download this code zip file.
2. Extract zip file under your working directory <work_path>.
3. Run the MATLAB software, and "Add Path" extracted directories (i.e. <work_path>/flywalk-main).
4. Move to <work_path>/flywalk-main directory and run the following demos.

## Demo Codes
<b>Demo</b><br>
The first demo shows the structural connectivity of hemibrain 52, CmKm and DistKm ROIs in Drosophila melanogaster (Fig.1a).<br>
Pre-processed connect list files should be downloaded from [zenodo](https://doi.org/10.5281/zenodo.7827225) and extracted under 'data' directory before running this code.
~~~
>> marmoAudGLMindividual
loading : data/s34waM3_1.nii.gz
apply mask atlas...
apply highpass filter (0.0078125 Hz) : tfMRI and design matrix...
process GLM with Tukey-Taper(8) estimation ...
done t=5.3043sec
P-value=0.05, T-value=1.6525
Tmax of GLM6 marmoAuCube1s34waM3_1CTukey8 : audio tmax=9.2829, tcnt=8111, mrv=1.6545
...
~~~
<div align="center">
<img src="data/demo1.jpg">
</div>

After calculation of GLM analysis for individual sessions, mixed-effects model (2nd analysis) could be applied (Fig.2c).


