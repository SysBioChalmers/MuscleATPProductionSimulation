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

- KeyWords:

**Utilisation:** maximising ATP, predictive simulation, experimental data reconstruction; **Model Source:** HMR2.00; **Taxonomy:** Homo sapiens  **Condition:** Exercise 

- Reference: Submitted

- Pubmed ID: NaN

- Last update: 2018-09-05




This repository is administered by name @avlant, Division of Systems and Synthetic Biology, Department of Biology and Biological Engineering, Chalmers University of Technology


## Installation

### Required Software:

* *_PROGRAMMING LANGUAGE/Version_*  (e.g.):
  * You need a functional Matlab installation of **Matlab_R_2015_b** (MATLAB 7.3 and higher)
  * The [RAVEN](https://github.com/SysBioChalmers/RAVEN) toolbox for MATLAB. Add the directory to your Matlab path, instructions [here](https://se.mathworks.com/help/matlab/ref/addpath.html?requestedDomain=www.mathworks.com)

### Run
To reproduce the figures of the paper run corresponding files (e.g. figure1a.m for figure 1A)