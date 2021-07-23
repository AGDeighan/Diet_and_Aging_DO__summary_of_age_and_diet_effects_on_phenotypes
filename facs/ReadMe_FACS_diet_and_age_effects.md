# Diet and age effects on fluorescence-activated cell sorting (FACS) phenotypes

This directory contains the results of the analysis of diet and age effects on phenotypes measured by the fluorescence-activated cell sorting (FACS, percent counts of white blood cell types) assay. For each phenotype the direct effects of diet and age were estimated as well as the interaction between diet and age (i.e. the diet-specific effects of age). 

The "figures" folder contains a large figure summarizing the diet and age effects across all of the phenotypes, as well as individual plots for each phenotype that give more detailed information on the diet and age effects for each phenotype. See the "Figures" section below for more details.

The "results" folder has CSV files that contain the information on the estimated diet and age effects as well as a CSV file with general information on the phenotypes. See the "Output files" section below for more information.

The "scripts" folder contains the actual R-scripts used to run the analysis. The overall method is summarized below in the "Method" section.

<br>
<br>


## Summary of notable findings

N/A

<br>
<br>


## Method

The effects of diet and age on the FACS phenotypes were estimated using linear mixed models (R/lme4). The significance of age and diet effects were assessed with approximate F-tests using the Kenward and Roger (1997) approach (R/pbkrtest). The FACS traits are all percentages (percent of subpopultion within a parent population) and were logit-transformed before model-fitting and statistical testing.

The 28-month timepoint is imbalanced in terms of diet, but there still quite a few ad libitum mice (72), so we will use all three timepoints for the analysis. That being said, we should keep this imbalance in mind when interpreting the results of the tests for the direct diet effects and the diet-age interaction effects.

| Timepoint  | Ad Libitum | 1-day fast | 2-day fast | 20% CR | 40% CR |
| -----------| ---------- | ---------- | ---------- | ------ | ------ |  
|   5 months | 188        | 188        | 190        | 189    | 182    | 	
|  16 months | 166        | 163        | 168        | 167    | 166    |					
|  28 months |  72        |  93        |  92        | 107    | 118    |	


<br>

The direct effect of diet on each phenotype was estimated independently for each timepoint. When estimating these diet effects, collection date (test batch) was included as a random effect and diet was treated as a fixed effect. The significance of the diet effect was assessed by comparing the full and null models shown below:  
 - Full model: Phenotype ~ (1 | DateCollect) + Diet
 - Null model: Phenotype ~ (1 | DateCollect)
 
For phenotypes significantly (p < 0.05) affected by diet at a timepoint, the significance of the pairwise differences between each diet at that timepoint were assessed using a Tukey correction for multiple tests (R/emmeans).

<br>

 The direct effect of age on each phenotype was estimated using a linear mixed model with diet, collection date (batch), and mouse ID as random effects and age in months as a fixed effect. The significance of the age effect was assessed by comparing the full and null models shown below:
 - Full model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + (1 | Diet) + Age
 - Null model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + (1 | Diet)

 <br>

 The interaction effect of age and diet on each phenotype was estimated using a linear mixed model with collection date (batch) and mouse ID as random effects and diet, age in months, and the interaction between diet and age in months as fixed effects. The significance of the age-diet interaction effect was assessed by comparing the full and null models shown below:
 - Full model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + Diet + Age + Age:Diet
 - Null model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + Diet + Age
 
 For phenotypes whose age effect significantly (p < 0.05) varies by diet, the significance of the pairwise differences between each diet were assessed using a Tukey correction for multiple tests (R/emmeans).



<br>
<br>

## Figures

 - histograms_of_phenotypes_and_transformations.pdf: histograms of the untransformed and transformed (if a transform was performed) phenotype values
 - diet_and_age_effects_summarized.pdf: summarizes the age and diet effects on the phenotypes  
   - The top figures show the estimated effect (intercept effect, a shift in the mean) of each diet on the mean of each phenotype at each timepoint. All the phenotypes were Z-scaled (mean 0 and standard deviation 1), so the units of the estimated effects are in standard deviations. A positive diet coefficient for a given phenotype at a given timepoint, indicates that the specified diet group had a higher mean relative to the other diet groups. The Holm (1979) correction for a family wise error rate (FWER) for 26 tests (the number of phenotypes tested) was made separately for each timepoint. Phenotypes that were significantly affected (FWER < 0.05) by diet at a given timepoint are indicated by an asterisk. Note that the first timepoint (5 months) is actually before the start of the intervention (the diet intervention started at 6 months), so we do not expect to see any diet effects at this timepoint.
   - The bottom-left plot shows the overall/main age effects. The plot shows the age coefficient in units of standard deviations per month +/- 1 standard error. Thus, an age coefficient of -0.05 indicates that the mean value of a phenotype across diets decreases 0.05 standard deviations per month. Again, the Holm correction (26 tests) for a family wise error rate was used to determine the statistical significance of the age effects. Phenotypes that are significantly (FWER < 0.05) affected by age are indicated by an asterisk on the left edge of the plot.
   - The bottom-right plot shows the diet-specific age affects (diet-specific age coefficients +/- 1 standard error), also in units of standard deviations per month. Again, the Holm correction (26 tests) for a family wise error rate was used to determine the statistical significance. Phenotypes for which the age effect varies significantly (FWER < 0.05) by diet are indicated by an asterisk on the right edge of the plot.
- individual_phenotype_plots: this folder contains a separate PDF for each phenotype, each PDF contains two figures
   - The top figure shows the distribution of the phenotype at each timepoint with box-scatter plots grouped by diet.
   - The bottom figure shows the mean and standard error of the phenotype at each timepoint.
   - The caption below the two figures gives the detailed results of the diet and age effect analysis for the phenotype.
   

<br>
<br>

## Output files

There are four primary output files which contain the results of the analysis and information on the phenotype(s)  
 - phenotype_info_table.csv: contains general information on the phenotypes
    - Phenotype: the name of the phenotype
    - Description: description of the phenotype
    - Assay: the name of the assay to which the phenotype belongs
    - NAll_11months: The total number of non-missing samples/values from the 11 month timepoint 	
    - NAll_23months: The total number of non-missing samples/values from the 23 month timepoint 	
    - NAll_35months: The total number of non-missing samples/values from the 35 month timepoint 	
    - NAL_11months: The number of non-missing samples/values from the 11 month timepoint belonging to mice of the ad libitum group	
    - NAL_23months: The number of non-missing samples/values from the 23 month timepoint belonging to mice of the ad libitum group	
    - NAL_35months: The number of non-missing samples/values from the 35 month timepoint belonging to mice of the ad libitum group	
    - N1D_11months: The number of non-missing samples/values from the 11 month timepoint belonging to mice of the 1-day fasting group		
    - N1D_23months: The number of non-missing samples/values from the 23 month timepoint belonging to mice of the 1-day fasting group		
    - N1D_35months: The number of non-missing samples/values from the 35 month timepoint belonging to mice of the 1-day fasting group		
    - N2D_11months: The number of non-missing samples/values from the 11 month timepoint belonging to mice of the 2-day fasting group		
    - N2D_23months: The number of non-missing samples/values from the 23 month timepoint belonging to mice of the 2-day fasting group		
    - N2D_35months: The number of non-missing samples/values from the 35 month timepoint belonging to mice of the 2-day fasting group		
    - N20_11months: The number of non-missing samples/values from the 11 month timepoint belonging to mice of the 20% caloric restriction group		
    - N20_23months: The number of non-missing samples/values from the 23 month timepoint belonging to mice of the 20% caloric restriction group		
    - N20_35months: The number of non-missing samples/values from the 35 month timepoint belonging to mice of the 20% caloric restriction group		
    - N40_11months: The number of non-missing samples/values from the 11 month timepoint belonging to mice of the 40% caloric restriction group		
    - N40_23months: The number of non-missing samples/values from the 23 month timepoint belonging to mice of the 40% caloric restriction group		
    - N40_35months: The number of non-missing samples/values from the 35 month timepoint belonging to mice of the 40% caloric restriction group	
    - Skewness: The skewness of the untransformed data
    - Kurtosis: The kurtosis of the untransformed data
    - Transform: Transform performed on the data before model fitting (if any)
        - none: No transform performed
        - log: natural-log transform ln[x + 0.01], 0.01 was added to account for any zeroes
        - logit: logit transform ln[((x + 0.01)/100) / (1 - (x + 0.01)/100))], 0.01 was added to account for any zeroes and the values were divided by 100 to convert from percent to proportion
    - SD: The standard deviation of the phenotype (after transform if one was performed and excluding values from timepoints not included in the analysis)
    - Mean: The mean of the phenotype (after transform if one was performed and excluding values from timepoints not included in the analysis)
    - Median: The median of the phenotype (after transform if one was performed and excluding values from timepoints not included in the analysis)
    - Min: The minimum value of the phenotype (after transform if one was performed and excluding values from timepoints not included in the analysis)
    - Max: The maximum value of the phenotype (after transform if one was performed and excluding values from timepoints not included in the analysis)
    - REF_Cov: Covariates that were included as random effect terms in the models
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;  Note:   
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   when estimating the direct effect of age, diet was additionally included as a random effect term   
       &nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;&nbsp;   when estimating the direct effect of diet at each timepoint, mouse ID was not included as a random effect term 
    - FEF_Cov: Covariates that were included as fixed effect terms in the models  
    <br>

 - diet_effects_stratified_by_timepoint.csv: 
    - Phenotype: the name of the phenotype
    - Timepoint: the timepoint tested
    - AL_Mean: the effect (intercept effect so it is a mean) of the ad libitum diet on the phenotype at the indicated timepoint	
    - AL_SE: the standard error of the effect	
    - 1D_Mean: the effect (intercept effect so it is a mean) of the 1-day fasting diet on the phenotype at the indicated timepoint	
    - 1D_SE: the standard error of the effect	
    - 2D_Mean: the effect (intercept effect so it is a mean) of the 2-day fasting diet on the phenotype at the indicated timepoint	
    - 2D_SE: the standard error of the effect	
    - 20_Mean: the effect (intercept effect so it is a mean) of the 20% caloric restriction diet on the phenotype at the indicated timepoint	
    - 20_SE: the standard error of the effect	
    - 40_Mean: the effect (intercept effect so it is a mean) of the 40% caloric restriction diet on the phenotype at the indicated timepoint	
    - 40_SE: the standard error of the effect
    - Diet_FTestStat: calculated F-statistic
    - Diet_FTestNumDF: numerator degrees of freedom
    - Diet_FTestDenDF: denominator degrees of freedom
    - Diet_FTestScaling: factor by which the F-statistic is scaled (multiplied)	
    - Diet_FTestPValue: p-value of the main diet effect	
    - Diet_PairwiseDifferences: diet pairs that are significantly different (pTukey < 0.05) from each other (only calculated if Diet_FTestPValue < 0.05)  
    <br>

 - age_effects_conditioned_on_diet.csv: 
    - Phenotype: the name of the phenotype
    - Age_Slope: the mean change in the phenotype (transformed values if a transform was made) per month (age coefficient)
    - Age_SE: the standard error of the effect (standard error of age coefficient)
    - Age_FTestStat: calculated F-statistic
    - Age_FTestNumDF: numerator degrees of freedom
    - Age_FTestDenDF: denominator degrees of freedom
    - Age_FTestScaling: factor by which the F-statistic is scaled (multiplied)	
    - Age_FTestPValue: p-value of the main age effect	

    <br>  
 - diet_age_interaction_effects.csv: 
    - Phenotype: the name of the phenotype
    - AL_Slope: the diet-specific age coefficient for mice on the ad libitum diet	
    - AL_SE: the standard error of the diet-specific age coefficient	
    - 1D_Slope: the diet-specific age coefficient for mice on the 1-day fasting diet	
    - 1D_SE: the standard error of the diet-specific age coefficient	
    - 2D_Slope: the diet-specific age coefficient for mice on the 2-day fasting diet	
    - 2D_SE: the standard error of the diet-specific age coefficient	
    - 20_Slope: the diet-specific age coefficient for mice on the 20% caloric restriction diet	
    - 20_SE: the standard error of the diet-specific age coefficient	
    - 40_Slope: the diet-specific age coefficient for mice on the 40% caloric restriction diet	
    - 40_SE: the standard error of the diet-specific age coefficient
    - DietAge_FTestStat: calculated F-statistic
    - DietAge_FTestNumDF: numerator degrees of freedom
    - DietAge_FTestDenDF: denominator degrees of freedom
    - DietAge_FTestScaling: factor by which the F-statistic is scaled (multiplied)	
    - DietAge_FTestPValue: p-value of the diet-age interaction effect	
    - DietAge_PairwiseDifferences: diet pairs for which the age effects are significantly different (pTukey < 0.05) from each other (only calculated if DietAge_FTestPValue < 0.05)  

<br>