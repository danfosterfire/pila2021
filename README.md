# pila2021

Code for data collection, processing, and analysis for the paper [Threats to the persistence of sugar pine (Pinus lambertiana) in the western USA](https://doi.org/10.1016/j.foreco.2023.121659).

## Abstract

This study assesses how numerous stressors shape the vital rates (survival, growth, and fecundity) of sugar pine across the vast majority of its range. Sugar pine (Pinus lambertiana Dougl.) is the largest Pinus species, an important timber species, and a component of several dry conifer forest types of western North America, in particular the extensive Sierra Nevada mixed conifer forest. The species faces several challenges in the Anthropocene, including a disrupted fire regime, an invasive pathogen, forest structure changes, and drought with ensuing bark beetle epidemics. Managers are concerned about the conservation outlook for sugar pine, but it is unclear where and how to best invest conservation resources. Using data from the US Forest Service’s Forest Inventory and Analysis program, we estimate the parameters for the vital rate functions and use these to construct an integral projection model which predicts the effects of various stressors on the asymptotic population growth rate. The asymptotic population growth rate is near or slightly below one even under reference conditions, and the actual abundance (in terms of both stem density and basal area) slightly declined over the duration of the study (2001–2019). The analysis reveals that wildfire, white pine blister rust, and forest density are key drivers of the demographic rates of sugar pine across its range. Drought and site dryness had lesser, but still meaningful, effects. Fire has strong negative effects on survival, resulting in a strongly negative population trajectory on burned sites. Conversely, lower than average forest density (plot-level basal area) results in a positive population growth rate via beneficial effects on individual growth. These results highlight the value of fire hazard mitigation, particularly where it also reduces forest density, in the conservation of this important species.

## 02-data

Raw data files are .gitignored, but their provenance is described in the data preparation
files in /03-analysis/00-pila_range_data.R and /03-analysis/01-import_and_process_fia.R.

This folder contains data at several stages of processing, from preliminary cleaning 
and aggregation to model-ready inputs, fitted model objects, and simulation realizations.

## 03-analysis

Contains ordered scripts used to implement every step of the data processing, 
statistical analysis, and generation of tables and figures. Not all analytical approaches 
made it into the final paper, but they are preserved here for posterity.

## 04-communication

Contains figures and tables intended for inclusion in reports, presentations, 
and the journal article, as well as the working files for publishing the article.