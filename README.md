# PEDso - Positive Energy District System Optimiser

Open-Source model for the optimisation of the (multi-) energy system of Positive Energy Districts (PEDs) in particular.

PEDso is created for an easy usage and without advanced preliminary knowledge in mathematical programming. Therefore, six spreadsheets are used to create and change the PED scenarios.

If the model is used for for academic work, please cite it accordingly:

https://doi.org/10.1016/j.energy.2022.124152, https://doi.org/10.3390/en14164864 and https://doi.org/10.5281/zenodo.5749268


## Usage of PEDso

1. Preparation of input files: Fill out the Excel input files according to your case study and if you want to run the annualised cost (AC) or net present value (NPV) optimisation.
2. Open Python Environment (PEDso is written using Spyder), open the *PEDso_AC* or *PEDso_NPV* python file and run it.
3. Save the directory link for the six spreadsheats in respective variables
4. Run: `m = build_model(in_tech, in_loc,in_timeseries, in_emisFules, in_tariffs, in_other, projectlife = 20, dt= 60, yeardivision = None, consHeat = "Yes", consCool = "Yes")"` and change the input according to your wishes. If the AC model is chosen, electric vehicles can be included and the file *Input_year_EV* and the selection *EV = Yes/No* needs to be given  **Attention:** Make sure that your function input is in line with you input in the spreadsheats!! --> Wait until the model is built with all its constraints. This might take a while ...
5. Solve the model: `SolverFactory('gurobi').solve(m,tee=True, keepfiles=True)` Instead of "gurobi" other, non-commercial solvers can be used here. Any solver used needs to be installed before
6. The model is solved. To make the most important results more readable and easier to plot, the helper file includes a function that creates a solution dictionary: run the *helpers* phython file.
7. Save the results in a variable of your choice --> Run `results = PED(m)` --> Now the results of the solved pyomo model m are stored in the dictionary results. 
8. To read the results from the dictionary use the the respective key `results["key"]`. The keys are the following: dict_keys(['obj', 'cost_I', 'cost_OMfix', 'cost_OMvar', 'cost_fuel', 'cost_CO2', 'cost_Ext', 'rev', 'capa_c', 'power', 'heat', 'cool', 'elEner', 'hEner', 'gridIn', 'gridInPE', 'gridEx', 'gridExPE', 'selfCons', 'selfCons_annual', 'gridInMax', 'emiss_CO2', 'emiss_CO2_annual', 'cons_fuel_vol', 'cons_fuel_e', 'emiss_PM', 'pv_groundPanelArea', 'batt_soc'])



This research has received funding from the European Union's Horizon 2020 research and innovation programme under the Marie Skłodowska-Curie Actions, Innovative Training Networks, Grant Agreement No. 812730.
