# Diet and age effects on frailty phenotypes

This directory contains the results of the analysis of diet and age effects on phenotypes measured by the frailty assay. For each phenotype the direct effects of diet and age were estimated as well as the interaction between diet and age (i.e. the diet-specific effects of age). 

The "figures" folder contains a large figure summarizing the diet and age effects across all of the phenotypes, as well as individual plots for each phenotype that give more detailed information on the diet and age effects for each phenotype. See the "Figures" section below for more details.

The "results" folder has CSV files that contain the information on the estimated diet and age effects as well as a CSV file with general information on the phenotypes. See the "Output files" section below for more information.

The "scripts" folder contains the actual R-scripts used to run the analysis. The overall method is summarized below in the "Method" section.

<br>
<br>


## Notes

*Keep in mind: higher scores are worse*  

While there are no diet effects on the frailty data at the 5-month timepoint (because the diets were started at 6 months), diet-group has a significant effect on mean frailty at 5 months. Fortunately, The ad libitum group has the lowest overall frailty score at 5 months, so any increase in frailty in the ad libitum mice compared to the other diet groups cannot be attributed to baseline differences. The two diet-groups that have significantly (Tukey honestly significant difference) different mean frailty values at 5 months are the mice assigned to the ad libitum group and the mice assigned to the 20% CR group (the mean frailty is higher in the 20% CR group). Additionally, when looking at the individual frailty components at 5 months, the mice assigned to the 2-day fasting group have significantly more dermatitis than the other four groups and the mice assigned to the 20% CR group have more loss of whiskers than the other four groups. The higher dermatitis in the 2-day fasting group does not persist across all the timepoints (perhaps because dermatitis is a reason for euthanasia), but the 20% CR mice consistently are among the groups with high loss of whisker scores across the timepoints.


<br>
<br>


## Summary of results  

*Keep in mind: higher scores are worse*  

The overall mean frailty score (mean of component indices) increases with age and the rate of increase is most rapid for the ad libitum mice. The rate of increase for the ad libitum mice is significantly (Tukey p-value < 0.05) greater than the rate of increase for the 2-day fasting mice, the 20% CR mice, and the 40% CR mice. In addition to a difference in the overall frailty score, several of the individual component indices are affected by diet.  

   - Dermatitis and alopecia both increase most rapidly in the ad libitum mice. Alopecia increases significantly (Tukey p-value < 0.05) more quickly in the ad libitum mice than any of the other diet groups, and significantly (Tukey p-value < 0.05) more slowly in the 40% CR group than any of the other groups. The rate of increase in dermatitis is significantly (Tukey p-value < 0.05) higher in the ad libitum mice than the 2-day fasting, 20% CR, and 40% CR mice. In addition to the baseline differences described in the **Notes** section above, another point to keep in mind with dermatitis is that it was one of our conditions for euthanisia. Therefore, mice that develop dermatitis may often be shortly thereafter euthanized and this dermatitis may not appear to increase as quickly with age as it would if the mice were not euthanized.  
   - The change in the distended abdomen and tumor scores with age is also affected by diet. The rate of increase for the tumor score is significantly (Tukey p-value < 0.05) lower in the 2-day fasting and 40% CR groups than the other groups (the difference between the 20% CR group and the ad libitum and 1-day fasting groups is not significant). Similarly, the distended abdoman score increases significantly (Tukey p-value < 0.05) faster in the ad libitum and 1-day fasting groups than the other groups.  
   - Gait disorders increase significantly (Tukey p-value < 0.05) more rapidly with age for the ad libitum mice than the 2-day, 20% CR, and 40% CR groups.
   - Abnormalities in breathing rate and/or depth increases significantly (Tukey p-value < 0.05) more rapidly for the ad libitum mice than the other groups.
   - The 40% CR mice actually perform *worse* in terms of tremors and body condition. The tremor score increases significantly (Tukey p-value < 0.05) more rapidly in the 40% CR mice than in the 1-day and 2-day fasting mice (the difference from ad libitum and 20% CR is not significant). The body condition score increases significantly (Tukey p-value < 0.05) more rapidly in the 40% CR mice than in the 1-day, 2-day, and 20% CR mice (the difference from ad libitum is not significant).

<br>
<br>


## Method

The effects of diet and age on the frailty phenotypes were estimated using linear mixed models (R/lme4). The significance of age and diet effects were assessed with approximate F-tests using the Kenward and Roger (1997) approach (R/pbkrtest). The traits were not transformed before model fitting and testing. 

Only the first six timepoints (5, 8, 10, 16, 22, and 28 months) were used for testing the direct diet effects and the diet-age interaction effects because by 34 months there is a large diet-imbalance (20 ad libitum mice and 91 40% CR mice). However, for testing the direct age effects, all six timepoints were used. The table below shows the sample counts for each timepoint (see the "phenotype_info_table.csv" file in the results folder for more information).

| Timepoint  | Ad Libitum | 1-day fast | 2-day fast | 20% CR | 40% CR |
| -----------| ---------- | ---------- | ---------- | ------ | ------ |  
|   5 months | 158        | 159        | 158        | 151    | 144    | 	
|   8 months |  30        |  28        |  32        |  30    |  30    |					
|  10 months | 185        | 177        | 184        | 184    | 179    |	
|  16 months | 169        | 164        | 169        | 166    | 166    | 	
|  22 months | 127        | 137        | 141        | 142    | 157    |					
|  28 months |  74        |  94        |  92        | 108    | 121    |				
|  34 months |  20        |  38        |  51        |  60    |  91    |			
|  37 months |   7        |  14        |  29        |  35    |  68    |	


<br>

The direct effect of diet on each phenotype was estimated independently for each timepoint. When estimating these diet effects, collection date (test batch), technician, and coat color were included as a random effects and diet was treated as a fixed effect. The significance of the diet effect was assessed by comparing the full and null models shown below:  
 - Full model: Phenotype ~ (1 | DateCollect) + (1 | Tech) + (1 | Coat) + Diet
 - Null model: Phenotype ~ (1 | DateCollect) + (1 | Tech) + (1 | Coat) 
 
For phenotypes significantly (p < 0.05) affected by diet at a timepoint, the significance of the pairwise differences between each diet at that timepoint were assessed using a Tukey correction for multiple tests (R/emmeans).

<br>

 The direct effect of age on each phenotype was estimated using a linear mixed model with diet, collection date (batch), technician, coat color, and mouse ID as random effects and age in months as a fixed effect. The significance of the age effect was assessed by comparing the full and null models shown below:
 - Full model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + (1 | Tech) + (1 | Coat) + (1 | Diet) + Age
 - Null model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + (1 | Tech) + (1 | Coat) + (1 | Diet)

 <br>

 The interaction effect of age and diet on each phenotype was estimated using a linear mixed model with collection date (batch), technician, coat color, and mouse ID as random effects and diet, age in months, and the interaction between diet and age in months as fixed effects. The significance of the age-diet interaction effect was assessed by comparing the full and null models shown below:
 - Full model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + (1 | Tech) + (1 | Coat) + Diet + Age + Age:Diet
 - Null model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + (1 | Tech) + (1 | Coat) + Diet + Age
 
 For phenotypes whose age effect significantly (p < 0.05) varies by diet, the significance of the pairwise differences between each diet were assessed using a Tukey correction for multiple tests (R/emmeans).



<br>
<br>

## Figures

 - histograms_of_phenotypes_and_transformations.pdf: histograms of the untransformed and transformed (if a transform was performed) phenotype values
 - diet_and_age_effects_summarized.pdf: summarizes the age and diet effects on the phenotypes  
   - The top figures show the estimated effect (intercept effect, a shift in the mean) of each diet on the mean of each phenotype at each timepoint. All the phenotypes were Z-scaled (mean 0 and standard deviation 1), so the units of the estimated effects are in standard deviations. A positive diet coefficient for a given phenotype at a given timepoint, indicates that the specified diet group had a higher mean relative to the other diet groups. The Holm (1979) correction for a family wise error rate (FWER) for 26 tests (the number of phenotypes tested) was made separately for each timepoint. Phenotypes that were significantly affected (FWER < 0.05) by diet at a given timepoint are indicated by an asterisk.
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
    - NAll_05months: The total number of non-missing samples/values from the 5 month timepoint 	
    - NAll_08months: The total number of non-missing samples/values from the 8 month timepoint 	
    - NAll_10months: The total number of non-missing samples/values from the 10 month timepoint 
    - NAll_16months: The total number of non-missing samples/values from the 16 month timepoint 	
    - NAll_22months: The total number of non-missing samples/values from the 22 month timepoint 	
    - NAll_28months: The total number of non-missing samples/values from the 28 month timepoint 	
    - NAll_29amonths: The total number of non-missing samples/values from the 34 month timepoint 	
    - NAll_37months: The total number of non-missing samples/values from the 37 month timepoint 
    - NAL_05months: The number of non-missing samples/values from the 5 month timepoint belonging to mice of the ad libitum group	
    - NAL_08months: The number of non-missing samples/values from the 8 month timepoint belonging to mice of the ad libitum group	
    - NAL_10months: The number of non-missing samples/values from the 10 month timepoint belonging to mice of the ad libitum group	
    - NAL_16months: The number of non-missing samples/values from the 16 month timepoint belonging to mice of the ad libitum group	
    - NAL_22months: The number of non-missing samples/values from the 22 month timepoint belonging to mice of the ad libitum group	
    - NAL_28months: The number of non-missing samples/values from the 28 month timepoint belonging to mice of the ad libitum group	
    - NAL_29amonths: The number of non-missing samples/values from the 34 month timepoint belonging to mice of the ad libitum group	
    - NAL_37months: The number of non-missing samples/values from the 37 month timepoint belonging to mice of the ad libitum group	
    - N1D_05months: The number of non-missing samples/values from the 5 month timepoint belonging to mice of the 1-day fasting group		
    - N1D_08months: The number of non-missing samples/values from the 8 month timepoint belonging to mice of the 1-day fasting group		
    - N1D_10months: The number of non-missing samples/values from the 10 month timepoint belonging to mice of the 1-day fasting group	
    - N1D_16months: The number of non-missing samples/values from the 16 month timepoint belonging to mice of the 1-day fasting group		
    - N1D_22months: The number of non-missing samples/values from the 22 month timepoint belonging to mice of the 1-day fasting group		
    - N1D_28months: The number of non-missing samples/values from the 28 month timepoint belonging to mice of the 1-day fasting group	
    - N1D_34months: The number of non-missing samples/values from the 34 month timepoint belonging to mice of the 1-day fasting group		
    - N1D_37months: The number of non-missing samples/values from the 37 month timepoint belonging to mice of the 1-day fasting group	
    - N2D_05months: The number of non-missing samples/values from the 5 month timepoint belonging to mice of the 2-day fasting group		
    - N2D_08months: The number of non-missing samples/values from the 8 month timepoint belonging to mice of the 2-day fasting group		
    - N2D_10months: The number of non-missing samples/values from the 10 month timepoint belonging to mice of the 2-day fasting group		
    - N2D_16months: The number of non-missing samples/values from the 16 month timepoint belonging to mice of the 2-day fasting group		
    - N2D_22months: The number of non-missing samples/values from the 22 month timepoint belonging to mice of the 2-day fasting group		
    - N2D_28months: The number of non-missing samples/values from the 28 month timepoint belonging to mice of the 2-day fasting group	
    - N2D_34months: The number of non-missing samples/values from the 34 month timepoint belonging to mice of the 2-day fasting group		
    - N2D_37months: The number of non-missing samples/values from the 37 month timepoint belonging to mice of the 2-day fasting group		
    - N20_05months: The number of non-missing samples/values from the 5 month timepoint belonging to mice of the 20% caloric restriction group		
    - N20_08months: The number of non-missing samples/values from the 8 month timepoint belonging to mice of the 20% caloric restriction group		
    - N20_10months: The number of non-missing samples/values from the 10 month timepoint belonging to mice of the 20% caloric restriction group		
    - N20_16months: The number of non-missing samples/values from the 16 month timepoint belonging to mice of the 20% caloric restriction group		
    - N20_22months: The number of non-missing samples/values from the 22 month timepoint belonging to mice of the 20% caloric restriction group		
    - N20_28months: The number of non-missing samples/values from the 28 month timepoint belonging to mice of the 20% caloric restriction group	
    - N20_34months: The number of non-missing samples/values from the 34 month timepoint belonging to mice of the 20% caloric restriction group		
    - N20_37months: The number of non-missing samples/values from the 37 month timepoint belonging to mice of the 20% caloric restriction group	
    - N40_05months: The number of non-missing samples/values from the 5 month timepoint belonging to mice of the 40% caloric restriction group		
    - N40_08months: The number of non-missing samples/values from the 8 month timepoint belonging to mice of the 40% caloric restriction group		
    - N40_10months: The number of non-missing samples/values from the 10 month timepoint belonging to mice of the 40% caloric restriction group	
    - N40_16months: The number of non-missing samples/values from the 16 month timepoint belonging to mice of the 40% caloric restriction group		
    - N40_22months: The number of non-missing samples/values from the 22 month timepoint belonging to mice of the 40% caloric restriction group		
    - N40_28months: The number of non-missing samples/values from the 28 month timepoint belonging to mice of the 40% caloric restriction group		
    - N40_34months: The number of non-missing samples/values from the 34 month timepoint belonging to mice of the 40% caloric restriction group		
    - N40_37months: The number of non-missing samples/values from the 37 month timepoint belonging to mice of the 40% caloric restriction group	
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
