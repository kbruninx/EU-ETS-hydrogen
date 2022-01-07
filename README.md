THIS README IS OUTDATED! 

This code was developed to study the impact of the market stability reserve on the European Emission Trading System. It calculates an equilibrium between a representative price-taking agent on the ETS allowance auctions. The agent makes a trade-off between abatement and buying emission allowances, based on marginal abatement cost curves. It employs an iterative price-search algorithm based on ADMM to calculate this equilibrium iteratively. This allows considering a wide range of marginal abatement cost curves, which can be automatically calibrated to reproduce allowance prices post MSR reform (2019). The model allows studying the impact of emission allowance demand changes due to shocks, overlapping climate policies or EU ETS design changes as well as the impact of the convexity of the abatement cost curve.

The research code documented below was used in the following paper:

[1] K. Bruninx & Marten Ovaere, "COVID-19, Green Deal & the recovery plan permanently change emissions and prices in EU ETS Phase IV", Under review with Nature Communications, 2021.
Available online:

The numerical solution procedure based on ADMM has also been employed in the following paper, in which more details on its convergence can be found. Note, however, that the mathematical model in this paper is fundamentally different:

[2] Kenneth Bruninx, Marten Ovaere, Erik Delarue, "The long-term impact of the market stability reserve on the EU emission trading system," Energy Economics, Volume 89, 2020, art. no. 104746.

## Installation, hardware & software requirements
### Installation
After downloading the repository, no specific installation is required. The user should, however, have the following software installed to be able to execute the program:
- Julia (https://julialang.org/downloads/)
- Gurobi (https://www.gurobi.com/) and have a license for this solver. If the user doesn't have access to Gurobi, any alternative (open source) solver capable of solving quadratic programming problems can be used. This requires importing the relevant package (line 13), creating an environment (line 18) and creating an optimization model (line 53) in the "Main.jl" file.

The program can be executed in the Julia terminal. Any text editor can be used to modify the code. However, it may be more convenient to work with a tool such as Visual Studio Code (https://code.visualstudio.com/).

If the user does not have any of these programs installed, the installation should take approximately one hour.

### Software requirements:
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
No specific hardware is required. All cases in [1] required less than 10 minutes to be solved on a with a i7-6600U processor clocking at 2.6 Ghz and 8GB of RAM.

## Running the program
### Input
There are three places where the user can change the input.

1. On line 9 of the "Main.jl"-file, the user can specify the set of scenarios one would like to study, ranging from "start_scen" to "stop_scen". This determines which scenarios in the "overview_scenarios.csv"-file will be executed.
2. The specifications of each scenario can be found in the file "overview_scenarios.csv"
    - scen_number (Integer): the number of the scenario
    - LRF_2021 (Float): the Linear Reduction Factor applied in the period 2021-2030 (default: 0.022)
    - LRF_2031 (Float): the Linear Reduction Factor applied as of 2021 (default: 0.022)
    - COVID (0/1): consider the negative demand shock induced by COVID-19 yes (1) or no (0)
    - MSR (2018/2021): the MSR design considered (2018: current design, based on 2018 directive or 2021: as proposed in the Fit for 55 Package)
    - gamma (Float): exponent in the functional form of the abatement cost curve
    - op_dem (Float): magnitude of the negative demand shock (MtCO2/year)
    - op_supply (Float): magnitude of the negative supply shock (MtCO2/year)
    - start_op (Integer): starting year of overlapping policy
    - stop_op (Integer): end year of overlapping policy
    - ref_scen (Integer): reference scenario. If this is the same as the scenario number (scen_number), the beta-parameter in the marginal abatement cost curve will be calibrated to reproduce emission allowance prices in 2019. If it differs from the scenario number in the  same row, the beta-value from the reference scenario will be retrieved and used in the marginal abatement cost curve. Note that the user must ensure that the results of the reference scenario are available.
3. The file "overview_data.yaml" contains a number of input parameters that are common to all scenarios. Examples include the number of years the analysis considers, the historical emissions and prices that will be used to calibrate the marginal abatement cost curve and the parameters of the market stability reserve for each of the two designs. Note that the impact of COVID (240 MtCO2, linearly decreasing between 2020-2025) can only be switched on/off (see above).

### Running the code
The user should execute the "MAIN.jl"-file, after having specified the scenario one would like to study. The "MAIN.jl"-file is properly documented to maximize transparency.

If the scenario at hand is a scenario in which the marginal abatement cost needs to be calibrated, the ADMM procedure will be executed a number of times until the emission allowance price in 2019 is replicated with a user-specified tolerance (default: 0.1 €/MWh, see "overview_data.yaml"). The code will report (in the terminal) the progress of the ADMM procedure and the difference between the computed emission allowance price in 2019 and the historical value. Only the final result (i.e., based on the calibrated MACC) will be stored.

If the scenario at hand can use a calibrated marginal abatement cost curve from the specified reference scenario (see above), the ADMM procedure is executed once.

### Output & Postprocessing
Running the code will generate two output files:

1. "overview_results.csv", in which aggregate metrics of all simulations will be listed:
    -   n_iter: number of iterations required to compute equilibrium
    -	walltime: number seconds required to compute equilibrium
    -   PrimalResidual: Primal residual on solution (metric for quality of solution, lower is better, see [2])
    -   DualResidual: Dual residual on solution (metric for quality of solution, lower is better, see [2])		
    -   Beta: Employed beta-value in this scenario
    -	EUA_2021: Emission allowance price in 2021 (€/tCO2)
    -   CumulativeEmissions: Cumulative emissions 2020-end horizon (MtCO2)
    -   Cancellation: Cancellation volume over lifetime ETS (MtCO2)
    -   WBseal: Year in which waterbed is sealed (-)
    -   WBL: waterbed leakage. Non-zero if an overlapping policy was specified, NaN otherwise (-)
    -   dirWBL: waterbed leakage due to the direct effect. Non-zero if an overlapping policy was specified, NaN otherwise (-)
    -   indirWBL: waterbed leakage due to the indirect or price effect. Non-zero if an overlapping policy was specified, NaN otherwise (-)
2. "Scenario_X.csv" (X = scenario number, in the folder "Results"), a csv-file per simulation with more detailed results, per year (first column)
    -   CAP: annual cap on emissions (MtCO2)
    -   Supply: annual supply of emission allowances after corrections MSR (MtCO2)
    -   Cancellation: annual volume of invalidated or cancelled allowances (MtCO2)
    -   MSR: holdings of the MSR at the end of each year (MtCO2)
    -   TNAC: total number of allowances in circulation at the end of each year (MtCO2)
    -   Emissions: annual emissions (MtCO2)
    -   EUAprice: price of emission allowances (€/tCO2)
    -   EUAs: annually procured emission allowances (MtCO2)

### Demos & reproducing the results in [1]
The file "overview_scenarios.csv" contains seven example scenarios to illustrate the functionality of the code. The output associated with these example scenarios can be found in the relevant files and directories (see above). To recreate these results, the user can simply run the "Main.jl" file (indicating the relevant scenario numbers in line 9-10).

The procedure to reproduce the numerical results in [1] are documented in the Section "Methods" of that paper.

## Additional material

### Figures
The folder "Figures" contains the data in all figures in [1].

### Tables
The folder "Tables" contains the data in the Supplementary Tables in [1].

## License
The software is made available under the MIT license (https://opensource.org/licenses/MIT).

## DOI
10.5281/zenodo.4923069

## Contributors
K. Bruninx (kenneth.bruninx@kuleuven.be)
