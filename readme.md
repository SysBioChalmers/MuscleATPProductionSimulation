# Muscle ATP Production Simulations (Complex I is bypassed during high intensity exercise)

- Models and scripts for simulations of ATP production in muscle, the repo contains: A) A small scale enzyme constrained stoichiometric model of metabolism of the single muscle fiber B) A whole body metabolic model incorporating complex I bypass.

- Abstract:
Human muscles are tailored towards ATP synthesis. The muscles convert glucose to lactate when exercising at high work rates. This is less nutrient efficient than respiration and there is hence a trade-off between endurance and power. To study this phenomenon, we integrated proteomics data from muscle fibers with an enzyme-constrained metabolic model. We found that ATP synthesis in muscle is constrained by the catalytic capacity of several enzymes. Additionally, bypassing complex I of the electron transport chain was found to increase the amount of ATP that can be synthesized per gram of protein, which suggested a third metabolic mode, distinct from the fully oxidative pathway. To test if this mode may occur in vivo, we conducted high resolved incremental exercise tests for five subjects. Their gas exchange at different work rates was accurately reproduced by a whole-body metabolic model incorporating complex I bypass. The study therefore provides new insight into metabolism during high intensity exercise.

- KeyWords: Protein allocation, FBA, UCP3, vO2max, lactate threshold, running economy, slow component

**Utilisation:** maximising ATP, predictive simulation, experimental data reconstruction; **Model Source:** HMR2.00; **Taxonomy:** Homo sapiens  **Condition:** Exercise 

- Reference: Submitted

- Pubmed ID: NaN

- Last update: 2019-08-29



This repository is administered by @avlant, Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology


## Installation
The repo contains scripts that run under Matlab, no installation is required. Open files in matlab and press the green play button, more details [here](https://se.mathworks.com/help/matlab/matlab_prog/create-scripts.html).

### Required Software:
Code has been tested under windows 10 using matlab and RAVEN.

* *_PROGRAMMING LANGUAGE/Version_*  (e.g.):
  * You need a functional Matlab installation of **Matlab_R_2018_b**
  * (optional, required for some functionality) The [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox for MATLAB. 
  
### Run
First run runMeFirst.m to set up the path to the embeded RAVEN functions  
To reproduce the figures of the paper run corresponding files (e.g. figure1b_XYZ.m for figure 1B)

* figure1b_paretoplot.m (runtime < 1 min) a pareto plot of protein and substrate efficency.
* figure2a_proteomics.m (runtime < 1 min) piechart of proteomics data and amount of protein per enzyme in the model.
* figure2b_fiberSimulation.m (runtime < 1 min) simulation results for muscle fiber constrained by proteomics data.
* figure2c_enzymeUsage.m (runtime < 1 min) compare model flux with vmax.
* figure2d_fiberVO2.m (runtime < 1 min) VO2 data of muscle fiber from literature compared with simulations.
* figure3a_bikeVO2literature.m (runtime < 1 min) comparison of whole body VO2 and Watt data from literature.
* figure3b_bikeVO2.m (runtime < 1 min) plot statistics
* figure4_connectedmodel.m (runtime < 5 min) multi tissue simulation of a subject compared with experimental measurments.
* figure5a_slowComponentData.m (runtime < 1 min) ODE simulation of lactate and oxygen dynamics compared with literature data.
* figure6a_wattMaxSubstrate.m (runtime < 1 min) multi tissue simulation of maximum performance of a subject on different carbon sources.
* figure6b_wattMaxSubstrateDepletion.m (runtime < 10 min) multi tissue simulation of maximum performance over time with a finite substrate pool.

* figureS1_enzymeUsage.m (runtime < 1 min) compare model flux with vmax
* figureS3_slope2y.m (runtime < 5 min) find best fit slopes for each subject.
* figureS4_parameteropimization.m (runtime < 20 min) find the best fit parameters for a subject.
* figureS5_randomSample.m (runtime < 20 min) perturb parameters in figure 2.
* figureS7_uncoupling.m (runtime < 5 min) simulate effects of lactate uptake in peripheral tissue.

### Addaption of the software to your data
The code is not intended to be used as a tool, but it is written in modular form, so minimal effort should in principle be required to replace parameter values, experimental data files, and stoichiometric models with custom versions. Benevolent behavior of the program cannot be guaranteed.

## FAQ
* I get the error message that "FaceAlpha" is not a property. Solution: update to the latest matlab version or comment out the line or remove the "FaceAlpha", value pair.