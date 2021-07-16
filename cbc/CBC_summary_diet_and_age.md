
## Summary of notable findings

Hematocrit, hemoglobin, and RBC count all decrease with age, but the decrease is most rapid in the 2-day fasting group. Also, while the increase in the number and percent of reticulocytes with age is not significantly affected by diet, the means at each timepoint are significantly affected by diet and are highest in the 2-day fasting mice. The CR groups have the highest red blood cell counts, hematocrit, and hemoglobin at both timepoints.

Mean RDW is significantly (Tukey p < 0.05) higher in the 2-day fasting group than all the other diet groups at both 11 months and 23 months. The mean RDW level in the 40% CR mice is also significantly higher (Tukey p < 0.05) than the ad libitum, 1-day fast, and 20% CR groups, but the rate of increase with age is lowest for the 40% CR (though the difference is not statistically significant) and by 23 months the mean RDW level in the 40% CR mice is not significantly different from the ad libitum, 1-day fast, and 20% CR groups.

Mean HDW at 11 months is significantly (Tukey p < 0.05) higher in the 40% CR and 2-day fasting groups, but the rate of increase in HDW is slowest (though the difference is not statistically significant) for the 40% CR mice and by 23 months mean HDW is not higher in the 40% CR mice than the ad libitum, 1-day fast, and 20% CR groups. At 23 months, mean HDW is still highest in the 2-day fasting mice, but the diet effect is no longer significant.

Neutrophil count increases with age. This increase appears to be exacerbated by the fasting regimens but reduced by the caloric restriction regimens. On the other hand, lymphocyte count decreases with age, and the rate of decrease is significantly greater in the 40% CR group than all the other diet groups. These two dynamics combine to mitigate the increase in the NLR in the 20% CR group compared to the other groups. This can also be seen in the percent counts. For all diets, the percent lymphocytes decrease with age and the percent neutrophils increase with age but the age effects on both are least pronounced in the 20% CR group; however, the age-diet interaction effect is not statistically significant (p > 0.05) for these two phenotypes.

For all diets except the 40% CR group monocyte count increases with age (note that the pairwise difference is only significant between the 40% CR group the 2-day fasting group which showed the greatest increase in monocyte count with age).

<br>
<br>


## Method

The effects of diet and age on the CBC phenotypes were estimated using linear mixed models (R/lme4). The significance of age and diet effects were assessed with approximate F-tests using the Kenward and Roger (1997) approach (R/pbkrtest). Before fitting the models, percent cell counts were logit-transformed and phenotypes with a skewness (calculated using R/e1071) greater than 1 were log-transformed.

Only the first two timepoints, 11 and 23 months, were used for testing diet and age effects because at 35 months there is notable survivor bias for the diet groups (only 14 ad libitum mice but 79 40% CR mice).

| Timepoint  | Ad Libitum | 1-day fast | 2-day fast | 20% CR | 40% CR |
| -----------| ---------- | ---------- | ---------- | ------ | ------ |  
|  11 months | 182        | 172        | 179        | 183    | 176    | 	
|  23 months | 113        | 129        | 133        | 137    | 153    |					
|  35 months |  14        |  28        |  37        |  50    |  79    |	


<br>

The direct effect of diet on each phenotype was estimated independently for each timepoint. When estimating these diet effects, collection date (test batch) was included as a random effect and diet was treated as a fixed effect. For traits that are affected by platelet clumping (eosinophil and platelet phenotypes), the log-transformed clump-count was included as a fixed effect. The significance of the diet effect was assessed by comparing the full and null models shown below (note that for the eosinophil and platelet counts, both the full and null models also include a fixed effect term for log-transformed clump-count that is not shown below)  
 - Full model: Phenotype ~ (1 | DateCollect) + Diet
 - Null model: Phenotype ~ (1 | DateCollect)
 
For phenotypes significantly (p < 0.05) affected by diet at a timepoint, then the significance of the pairwise differences between each diet at that timepoint were assessed using a Tukey correction for multiple tests (R/emmeans).

<br>

 The direct effect of age on each phenotype was estimated using a linear mixed model with diet, collection date (batch), and mouse ID as random effects and age in months as a fixed effect. For traits that are affected by platelet clumping (eosinophil and platelet phenotypes), the log-transformed clump-count was also included as a fixed effect. The significance of the age effect was assessed by comparing the full and null models shown below (note that for the eosinophil and platelet counts, both the full and null models also include a fixed effect term for log-transformed clump-count that is not shown below)  
 - Full model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + (1 | Diet) + Age
 - Null model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + (1 | Diet)

 <br>

 The interaction effect of age and diet on each phenotype was estimated using a linear mixed model with collection date (batch) and mouse ID as random effects and diet, age in months, and the interaction between diet and age in months as fixed effects. For traits that are affected by platelet clumping (eosinophil and platelet phenotypes), the log-transformed clump-count was also included as a fixed effect. The significance of the age-diet interaction effect was assessed by comparing the full and null models shown below (note that for the eosinophil and platelet counts, both the full and null models also include a fixed effect term for log-transformed clump-count that is not shown below)  
 - Full model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + Diet + Age + Age:Diet
 - Null model: Phenotype ~ (1 | MouseID) + (1 | DateCollect) + Diet + Age
 
 For phenotypes whose age effect significantly (p < 0.05) varies by diet, the significance of the pairwise differences between each diet were assessed using a Tukey correction for multiple tests (R/emmeans).



<br>
<br>

## Figures

The plots in "histograms_of_phenotypes_and_transformations.pdf" show the histograms of the untransformed and transformed (if a transform was performed) phenotype values (only including values from the timepoints used for analysis)

The plots in "diet_and_age_effects_summarized.pdf" summarizes the age and diet effects on the phenotypes. The top figures show boxplots of each phenotype grouped by diet and timepoint. The phenotype values were normalized (mean = 0 and standard deviation = 1) before plotting so the y-axis is in units of standard deviations (of the transformed phenotype if a transform was made). Asterisks indicate phenotypes whose means at the indicated timepoint are significantly affected by diet at a family-wise error rate (F-test p-values adjusted using the Holm method for 26 tests) of 0.05. The bottom-left plot shows the overall/main age effects. The plot shows the age coefficient in units of standard deviations per month +/- 1 standard error. Asterisks indicate phenotypes that are significantly affected by age (FWER < 0.05). The bottom-right plot shows the diet-specific age affects, also in units of standard deviations per month. The asterisks indicate phenotypes for which the age-effect significantly (FWER < 0.05) varies by diet.

The folder "individual_phenotype_plots" contains a separate PDF for each phenotype. Each PDF contains two figures. The top figure shows the distribution of the phenotype at each timepoint (including the timepoint(s) not included in analysis) with box-scatter plots grouped by diet. The bottom figure shows the mean and standard error of the phenotype at each timepoint. The caption below the two figures gives more details on results of the diet and age effect analysis for that phenotype.




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