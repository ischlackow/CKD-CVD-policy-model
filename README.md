# Long-term health outcomes of people with reduced kidney function in the UK: a modelling study using population health data 

This release supplements the manuscript "Long-term health outcomes of people with reduced kidney function in the UK: a modelling study using population health data". It contains the list of Read codes used for defining model variables (`Read codes` folder) and R scripts used to produce Tables, Figures and other data for the manuscript (`R code` folder). 

Please note that the patient-level CPRD data used for analyses cannot be released and as such, the code is not executable. Rather, it is meant to be used for exploratory purposes.

As a general rules, R scripts with the _functions_ suffix contains functions that are called in _wrapper_ files. For example, the file `wrapper_intval.R` performs internal validation by calling functions from `extrapolation_functions_intval.R`

The following scripts were used for producing Tables, Figure and other data for the manuscript:

* Demographics; adverse events; missing data (S1 Text: Table C & Table D; S1 Text: Table E; S1 Text: Table A.1): `wrapper_numbers.R`, `numbers_functions.R`

* Internal validation (Fig 2, S1 Text: Figure B & Figure C): `wrapper_intval.R`, `extrapolation_functions_intval.R`, `intval_results.R`

* Internal validation: c-statistic (S1 Text: Table H): `c-index.R`

* Application 1: Predicting remaining life expectancy of people with reduced kidney function (Fig 3; S1 Text: Table I): `wrapper_app_no_treat.R`, `extrapolation_functions_app_no_treat.R`, `applications_wrapper.R`

* Application 2: Quantifying the impact of partial and optimal use of cardiovascular prevention medications in people with reduced kidney function (Fig 4; S1 Text: Table J; S1 Text: Table A.2): `undertreatment_gap.R`, `wrapper_app_no_treat.R`, `wrapper_app_full_treat.R`, `extrapolation_functions_app_no_treat.R`, `extrapolation_functions_app_full_treat.R`, `applications_wrapper.R`

* Various auxiliary functions: `model_paper_functions.R`

# Contributors

* Iryna Schlackow, Health Economics Research Centre, University of Oxford, [ORCID iD 0000-0002-4154-1431](https://orcid.org/0000-0002-4154-1431)

* Claire Simons, Health Economics Research Centre, University of Oxford, [ORCID iD 0000-0003-0822-7261](https://orcid.org/0000-0003-0822-7261)

* Borislava Mihaylova, Health Economics Research Centre, University of Oxford,[ORCID iD 0000-0002-0951-1304](https://orcid.org/0000-0002-0951-1304)
