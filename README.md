# PEDso - Positive Energy District System Optimiser

Open-Source model for the optimisation of the (multi-) energy system of Positive Energy Districts (PEDs) in particular.

PEDso is created for an easy usage and without advanced preliminary knowledge in mathematical programming. Therefore, six spreadsheets are used to create and change the PED scenarios.

If the model is used for for academic work, please cite it accordingly:

https://doi.org/10.3390/en14164864 and https://doi.org/10.5281/zenodo.5749268


## Usage of PEDso

1. Preparation of input files: Fill out the six Excel input files according to your case study.
2. Open Python Environment (PEDso is written using Spyder), open the *main* python file and run it.
3. Save the directory link for the six spreadsheats in respective variables
4. Run: `m = build_model(in_tech, in_loc,in_timeseries, in_emisFules, in_tariffs, in_other, projectlife = 20, dt= 60, yeardivision = None, consHeat = "Yes", consCool = "Yes")"`
5. asd
