# Muscle ATP Production Simulations (Complex I is bypassed during high intensity exercise)

- Models and scripts for simulations of ATP production in muscle, the repo contains: A) A small scale enzyme constrained stoichiometric model of metabolism of the single muscle fiber B) A whole body metabolic model incorporating complex I bypass.

- Abstract:
The exercising human must synthesize ATP at high rates. Humans produce ATP from glucose using either
an oxidative pathway or a fermentative pathway with lactate as byproduct. The oxidative pathway produces 
more ATP per glucose, but less ATP per gram protein, and there is hence a trade-off between endurance 
and power. To study how this affects single muscle fibers we integrated proteomics data with an enzyme 
constrained stoichiometric model of metabolism. We found that ATP synthesis is constrained by the 
specific activity of the metabolic enzymes. Additionally the amount of ATP that can be produced per gram 
of protein may be increased by bypassing complex I of the electron transport chain, suggesting a third
metabolic mode, distinct from the fully oxidative pathway. To test if the mode is of practical relevance 
we conducted a high resolved incremental exercise test. The gas exchange of the exercising subject was 
faithfully reproduced by a whole body metabolic model incorporating complex I bypass. We discuss how the
methods and results may be useful in the study of human endurance performance.




- KeyWords: Protein allocation, FBA, UCP3, vO2 max, lactate threshold, running economy, slow component

**Utilisation:** maximising ATP, predictive simulation, experimental data reconstruction; **Model Source:** HMR2.00; **Taxonomy:** Homo sapiens  **Condition:** Exercise 

- Reference: Submitted

- Pubmed ID: NaN

- Last update: 2018-09-05




This repository is administered by @avlant, Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology


## Installation
The files are scripts that run under Matlab, no installation is required. Open files in matlab and press the green play button, more details [here](https://se.mathworks.com/help/matlab/matlab_prog/create-scripts.html)

### Required Software:
Code has been tested under windows 7 using matlab and RAVEN.

* *_PROGRAMMING LANGUAGE/Version_*  (e.g.):
  * You need a functional Matlab installation of **Matlab_R_2018_b**
  * (optional, required for some functionality) The [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox for MATLAB. Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)
  
### Run
First run runMeFirst.m to set up the path to the embeded RAVEN functions  
To reproduce the figures of the paper run corresponding files (e.g. figure1b_XYZ.m for figure 1B)

* figure1b_paretoplot.m (runtime < 1 min) a pareto plot of protein and substrate efficency.
* figure1c_proteomics.m (runtime < 1 min) piechart of proteomics data and amount of protein per enzyme in the model.
* figure1d_fiberSimulation.m (runtime < 1 min) simulation results for muscle fiber constrained by proteomics data.
* figure1e_fiberVO2.m (runtime < 1 min) VO2 data of muscle fiber from literature.
* figure1f_bikeVO2.m (runtime < 1 min) comparison of whole body VO2 and Watt data from literature.
* figure2_connectedmodel.m (runtime < 5 min) multi tissue simulation of a subject compared with experimental measurments.
* figure3a_wattMaxSubstrate.m (runtime < 1 min) multi tissue simulation of maximum performance of a subject on different carbon sources.
* figure3b_wattMaxSubstrateDepletion (runtime < 10 min) multi tissue simulation of maximum performance over time with a finite substrate pool.
* figure4_slowComponent.m (runtime < 1 min) ODE simulation of lactate and oxygen dynamics.
* figureS1_enzymeUsage.m (runtime < 1 min) compare model flux with vmax
* figureS2_randomSample.m (runtime < 20 min) perturb parameters in figure 2.

### Addaption of the software to your data
The code is not intended to be used as a tool, but it is written in modular form, so minimal effort should in principle be required to replace parameter values, experimental data files, and stoichiometric models with custom versions. Benevolent behavior of the program cannot be guaranteed.

## FAQ
* I get the error message that "FaceAlpha" is not a property. Solution: update to the latest matlab version or  comment out the line or remove the "FaceAlpha", value pair.
