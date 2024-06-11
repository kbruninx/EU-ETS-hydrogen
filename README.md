This code was developed to study the interaction between the power sector, the energy intensive industry and the hydrogen sector on the auctions of the European Emission Trading System, the energy-only electricity market, the hydrogen market. It calculates an equilibrium between a set of representative price-taking agents on these markets and allows enforcing renewable targets and different defnitions of the "additionality" principle (putting requirements on the electricity used in the electrolysis proces). The market stability reserve in EU ETS is fully considered, endogenously adapting the supply of emission allowances.  It employs an iterative price-search algorithm based on ADMM to calculate this equilibrium iteratively. The emissions of the energy-intensive industry are represented via marginal abatement cost curves, which can be automatically calibrated to reproduce allowance prices post MSR reform (2019). 

(Variations on) this research code documented below have been used in the following papers:

[1] K. Bruninx & Marten Ovaere, "COVID-19, Green Deal & the recovery plan permanently change emissions and prices in EU ETS Phase IV", Nature Communications, Volume 13, art. no. 1165, 2022.

[2] Kenneth Bruninx, Marten Ovaere, Erik Delarue, "The long-term impact of the market stability reserve on the EU emission trading system," Energy Economics, Volume 89, art. no. 104746, 2020.

[3] K. Bruninx & Marten Ovaere, “Estimating the impact of COVID-19 on emissions and emission allowance prices under EU ETS”, IAEE Energy Forum / Covid-19 Issue 2020, 2020. 


## Installation, hardware & software requirements
### Installation
After downloading the repository, no specific installation is required. The user should, however, have the following software installed to be able to execute the program:
- Julia (https://julialang.org/downloads/)
- Gurobi (https://www.gurobi.com/) and have a license for this solver. If the user doesn't have access to Gurobi, any alternative (open source) solver capable of solving quadratic programming problems can be used. 

The program can be executed in the Julia terminal. Any text editor can be used to modify the code. However, it may be more convenient to work with a tool such as Visual Studio Code (https://code.visualstudio.com/).

If the user does not have any of these programs installed, the installation should take approximately one hour.

### Software requirements: (to be updated) 
This code has been developed using Julia v1.5. The solver used is Gurobi v.9.0.

The following Julia packages are required:
- JuMP v.0.21.5
- Gurobi v.0.8.1
- DataFrames v.0.21.7
- CSV v.0.7.7
- YAML v0.4.2
- ProgressBars v0.7.1
- Printf

If the user does not have any of these programs installed, the installation should take less than one hour.

### Hardware requirements 
No specific hardware is required. Depending on the configuration (number of agents and markets considered), computational effort may significantly increase.

## License
The software is made available under the MIT license (https://opensource.org/licenses/MIT).
 
## Contributors
K. Bruninx (k.bruninx@tudelft.nl)
J. A. Moncada (jorgeandres.moncadaescudero@kuleuven.be)
A. Hoogsteyn (alexander.hoogsteyn@kuleuven.be)
