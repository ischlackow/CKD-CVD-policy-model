# Long-term health outcomes of people with reduced kidney function in the UK: a modelling study using population health data 

This release supplements the manuscript "Long-term health outcomes of people with reduced kidney function in the UK: a modelling study using population health data". It contains the following data:

* __Read codes__ folder: the list of Read and entity codes used for defining model variables;

* __R code__ folder: R scripts used to extract data from CPRD as well as generate Tables, Figures and other data for the manuscript  

Please note that the patient-level CPRD data used for analyses cannot be released and as such, the code is not executable and perhaps not complete. Rather, it is meant to be used for exploratory purposes.

As a general rules, R scripts with the _functions_ suffix contains functions that are called in _wrapper_ files. For example, the file `wrapper_intval.R` performs internal validation by calling functions from `extrapolation_functions_intval.R`

## Code for CPRD data extraction

This code is located in the __R code/CPRD data extraction__ folder. The following main rules of thumb, based on literature, were used in defining variables. Please see individual files for further detail and full information on how all variables were derived.

* Diabetes classification: it was assumed that diabetes of Type I if patient was <35 years old at the time it was diagnosed, and on insulins at cohort entry date

* Albuminuria classification: this was based on the [relationship between albuminuria and proteinuria](https://www.ncbi.nlm.nih.gov/pmc/articles/PMC4089693/table/tbl7/)

* Medication uptake: it was assumed that a patient was taking a medication when at least two prescriptions wree recorded, with the latest at most 28 days apart from the point of interest (eg cohort entry)

The data were derived in the following way

* Raw data on individual variables were extracted in the `extract_vars_all.R` file

* Creatinine tests were converted into eGFR values and CKD classification in the `eGFR.R` file. Albuminuria classification was performed in the `albuminuria.R` files, using functions from `albuminuria_functions.R`

* Flat file with all baseline characteristics was derived in `df_baseline_cohort.R` file.

* The file `functions.R` file contains generic functions used by other scripts in this folder

## Code for manuscript tables and figures

This code is located in the __R code/manuscript tables and figures__ folder. Specific outputs wree generated as follows:

* Demographics; adverse events; missing data (S1 Text: Table C & Table D; S1 Text: Table E; S1 Text: Table A.1): `wrapper_numbers.R`, `numbers_functions.R`

* Internal validation (Fig 2, S1 Text: Figure B & Figure C): `wrapper_intval.R`, `extrapolation_functions_intval.R`, `intval_results.R`

* Internal validation: c-statistic (S1 Text: Table H): `c-index.R`

* Application 1: Predicting remaining life expectancy of people with reduced kidney function (Fig 3; S1 Text: Table I): `wrapper_app_no_treat.R`, `extrapolation_functions_app_no_treat.R`, `applications_wrapper.R`

* Application 2: Quantifying the impact of partial and optimal use of cardiovascular prevention medications in people with reduced kidney function (Fig 4; S1 Text: Table J; S1 Text: Table A.2): `undertreatment_gap.R`, `wrapper_app_no_treat.R`, `wrapper_app_full_treat.R`, `extrapolation_functions_app_no_treat.R`, `extrapolation_functions_app_full_treat.R`, `applications_wrapper.R`

* Various auxiliary functions: `model_paper_functions.R`

## Contributors

* Iryna Schlackow, Health Economics Research Centre, University of Oxford, [ORCID iD 0000-0002-4154-1431](https://orcid.org/0000-0002-4154-1431)

* Claire Simons, Health Economics Research Centre, University of Oxford, [ORCID iD 0000-0003-0822-7261](https://orcid.org/0000-0003-0822-7261)

* Borislava Mihaylova, Health Economics Research Centre, University of Oxford,[ORCID iD 0000-0002-0951-1304](https://orcid.org/0000-0002-0951-1304)
