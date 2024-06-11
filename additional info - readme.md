## Running the program
### Input (to be updated)
There are three places where the user can change the input.

1. On line 9 of the "Main.jl"-file, the user can specify the set of scenarios one would like to study, ranging from "start_scen" to "stop_scen". This determines which scenarios in the "overview_scenarios.csv"-file will be executed.
2. The specifications of each scenario can be found in the file "overview_scenarios.csv"
    - scen_number (Integer): the number of the scenario
    - LRF_2021 (Float): the Linear Reduction Factor applied in the period 2021-2030 (default: 0.022)
    - LRF_2031 (Float): the Linear Reduction Factor applied as of 2021 (default: 0.022)
    - COVID (0/1): consider the negative demand shock induced by COVID-19 yes (1) or no (0)
    - MSR (2018/2021): the MSR design considered (2018: current design, based on 2018 directive or 2021: as proposed in the Fit for 55 Package)
    - gamma (Float): exponent in the functional form of the abatement cost curve
    - ref_scen (Integer): reference scenario. If this is the same as the scenario number (scen_number), the beta-parameter in the marginal abatement cost curve will be calibrated to reproduce emission allowance prices in 2019. If it differs from the scenario number in the  same row, the beta-value from the reference scenario will be retrieved and used in the marginal abatement cost curve. Note that the user must ensure that the results of the reference scenario are available.
3. The file "overview_data.yaml" contains a number of input parameters that are common to all scenarios. Examples include the number of years the analysis considers, the historical emissions and prices that will be used to calibrate the marginal abatement cost curve and the parameters of the market stability reserve for each of the two designs. Note that the impact of COVID (240 MtCO2, linearly decreasing between 2020-2025) can only be switched on/off (see above).

### Running the code (to be updated)
The user should execute the "MAIN.jl"-file, after having specified the scenario one would like to study. The "MAIN.jl"-file is properly documented to maximize transparency.

If the scenario at hand is a scenario in which the marginal abatement cost needs to be calibrated, the ADMM procedure will be executed a number of times until the emission allowance price in 2019 is replicated with a user-specified tolerance (default: 0.1 €/MWh, see "overview_data.yaml"). The code will report (in the terminal) the progress of the ADMM procedure and the difference between the computed emission allowance price in 2019 and the historical value. Only the final result (i.e., based on the calibrated MACC) will be stored.

If the scenario at hand can use a calibrated marginal abatement cost curve from the specified reference scenario (see above), the ADMM procedure is executed once.

## Running the code on DelftBlue (TU Delft)

# Interactive jobs (to be updated for fair share accounting)
1. srun --job-name="your_job" --account=research-tpm-ess --partition=compute --time=00:30:00 --ntasks=1 --cpus-per-task=17 --mem-per-cpu=4GB --pty bash
2. cd /home/kbruninx/EU-ETS
3. module load 2022r2
4. module load julia 
5. julia --threads=16 MAIN.jl --start_scen 1 --stop_scen 15 > run.log

Tips: 
- The decision problems of the agents are solved in parallel in each iteration. Gurobi uses 4 threads to solve each problem.
- The number of cpus should be equal to the number of agents + 1 (1 master and 1 per agent)
- Set the number of cpu's equal to the number of threads used for julia x 4 (e.g., 17 x 4 = 68), as by default 4 threads are used by Gurobi
- Check resource use via seff [job_id]

Resources: 
- For basic info on how to use DelftBlue, see https://doc.dhpc.tudelft.nl/delftblue/crash-course/. 
- For interactive jobs, see https://doc.dhpc.tudelft.nl/delftblue/Slurm-scheduler/#interactive-use 


## Running the code on ThinKing (VSC) 

### Output & Postprocessing (to be updated)
Running the code will generate X output files:

1. "overview_results.csv", in which aggregate metrics of all simulations will be listed:
    -   n_iter: number of iterations required to compute equilibrium
    -	walltime: number seconds required to compute equilibrium
    -   PrimalResidual: Primal residual on solution, one for each considered market (metric for quality of solution, lower is better, see [2])
    -   DualResidual: Dual residual on solution, one for each considered market (metric for quality of solution, lower is better, see [2])		
    -   Beta: Employed beta-value in this scenario
    -	EUA_2021: Emission allowance price in 2021 (€/tCO2)
    -   CumulativeEmissions: Cumulative emissions 2020-end horizon (MtCO2)
2. "Scenario_X.csv" (X = scenario number, in the folder "Results"), a csv-file per simulation with more detailed results, per year (first column)
    -   CAP: annual cap on emissions (MtCO2)
    -   Supply: annual supply of emission allowances after corrections MSR (MtCO2)
    -   Cancellation: annual volume of invalidated or cancelled allowances (MtCO2)
    -   MSR: holdings of the MSR at the end of each year (MtCO2)
    -   TNAC: total number of allowances in circulation at the end of each year (MtCO2)
    -   Emissions: annual emissions (MtCO2)
    -   EUAprice: price of emission allowances (€/tCO2)
    -   EUAs: annually procured emission allowances (MtCO2)

### Demos & reproducing the results (to be updated)
## Some simulation
## Studying overlapping policies

### Extensions/future developments
- Allow user to study subset of markets, defined based on participating agents in yaml file (e.g. if "empty" power sector?)
- Additionality on monthly basis 
- Integrated selection of representative days?
