# tenm 0.5.1

* Improvements to documentation

# tenm 0.5.0

* Initial CRAN submission.

* Time-specific niche modeling (TENM) is a novel approach that allows 
  calibrating niche models with high temporal resolution spatial information, 
  which aims to reduce niche estimation biases. 
  What makes the tenm package stand out is its distinctive capability to 
  calibrate models by incorporating specific temporal information, 
  whether on a yearly, monthly, or even daily basis. This feature distinguishes
  it from traditional models that rely on averaged temporal data. 
  
  Some of the package functions are:
  - Time-specific spatial data thinning: data cleaning considering the 
    temporal dimension of records. 
  - Time-specific environmental data extraction: extract environmental data 
    from variables considering the temporal dimension of data.
  - Time-specific background generation: generate background points for 
    the modeling process by considering the temporal dimension of the 
    occurrence and environmental data: The number of background points for
    each year is proportional to the number of occurrences for each year
    of observation. 
  - Exporting time-specific information as Samples With Data format: 
    export the time-specific data to the Samples With Data format table. 
    This function allows users to use other modeling algorithms such 
    as MaxEnt and GLMs.
  - Time-specific model calibration: calibrate time-specific niche models 
    using minimum volume ellipsoids. It fits numerous models based on a 
    combination of user-set parameters, including different combinations of
    environmental variables.
  - Time-specific model selection: select n models from all the fitted 
    models using statistical and model performances as model 
    selection criteria.
  - Projecting time-specific niche models: project one or more of the 
    selected models onto both the environmental and/or geographical space. 
