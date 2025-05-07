[![License: MIT](https://img.shields.io/badge/License-MIT-success.svg)](https://opensource.org/licenses/MIT)
[![DOI](https://zenodo.org/badge/DOI/10.5281/zenodo.7827225.svg)](https://doi.org/10.5281/zenodo.7827225)

# flywalk
Code for a manuscipt of the flywalk.<br>
Takuto Okuno, Alexander Woodward, Hideyuki Okano, Junichi Hata (in submisstion)
["Functional connectivity, structural connectivity, and inter-individual variability in Drosophila melanogaster"](https://www.yahoo.com/)

<div align="center">
<img src="data/img/fig1.jpg" width="70%">
</div>

## Requirements: Software
* MATLAB R2021b or later
* Parallel Computing Toolbox ver7.1 or later


## Installation
1. Download this code zip file.
2. Extract zip file under your working directory <work_path>.
3. Run the MATLAB software, and "Add Path" extracted directories (i.e. <work_path>/flywalk-main).
4. Move to <work_path>/flywalk-main directory and run the following demos.

## Demo Codes
<b>Demo 1</b><br>
The first demo shows the structural connectivity (SC) of hemibrain primary 63, 52 ROIs, CmKm and DistKm ROIs in Drosophila melanogaster (Fig.1a, Fig.2a).<br>
Pre-processed connect list files should be downloaded from [zenodo](https://doi.org/10.5281/zenodo.7827225) and extracted under <work_path>/flywalk-main directory before running this code.<br>
(Please confirm to update data and results directory.)

~~~
>> makeStructConnectivity
...
~~~

<div align="center">
<img src="data/img/demo1.jpg" width="70%">
</div>
This demo shows several SC matrices of hemibrain, CmKm and DistKm ROIs.<br>
<br>
<b>Demo 2</b><br>
Second demo shows relation between spatial smoothing and FC-SC detection & correlation in DistKm and CmKm ROIs (Fig.2b).

~~~
>> plotFuncConnectivity
...
~~~

<div align="center">
<img src="data/img/demo2.jpg" width="100%">
</div>
From left to right, FC-SC correlation, FC-SC detection, and FC-SC detection & correlation, respectively. Smoothing size 0 to 300 with Polynomial & tCompCor methods. <br>

