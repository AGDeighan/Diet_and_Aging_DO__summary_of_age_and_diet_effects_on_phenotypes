# 2021-07-15
################################################################################
#
#   This script creates a phenotype information table to be used for estimating
#   diet and age effects on the FACS phenotypes
#
#
#   Author:  Andrew Deighan
#   E-mail: andrew.deighan@jax.org
#
################################################################################


################################################################################
## libraries etc #########################################################
options(na.action = 'na.exclude')
options(stringsAsFactors = FALSE)
library(tidyverse)

################################################################################
## helper functions #########################################################

standardise <- function(x){
  return(
    (x - mean(x, na.rm = TRUE))/sd(x, na.rm = TRUE)
  )
}

#####


################################################################################
## load data etc #########################################################

ANIMAL_DATA <- read_csv(
  'data/AnimalData_Processed_20210120.csv'
) %>% 
  mutate(
    Diet = factor(
      Diet,
      levels = c('AL', '1D', '2D', '20', '40')
    )
  )

PHENO_DATA <- read_csv(
  'data/Frailty_Processed_20201125.csv'
) %>% 
  select(-c(
    Frailty_AdjBatchTechCoat, Frailty_AdjBatchTechCoatBW
  )) %>% 
  rename(
    FrailtyMean = Frailty_Raw
  ) %>% 
  mutate(
    DateCollect = as.character(DateCollect),
    AgeInDays = AgeInDays/30.4,
    BW = ifelse(
      is.na(BW_Weekly), BW_Test,
      BW_Weekly
    ),
    BW = standardise(BW)
  ) %>% 
  rename(
    Age = AgeInDays
  )

#####


################################################################################
## Specify assay #########################################################

ASSAY <- 'Frailty'

#####


################################################################################
##### Select timepoints to test for diet and diet-age effects #####

# We will exclude the 34 and 37 month timepoints when estimating direct diet 
# effects and diet=age interaction effects

PHENO_DATA %>% 
  left_join(
    ANIMAL_DATA,
    by = 'MouseID'
  ) %>% 
  select(
    Timepoint, Diet
  ) %>% 
  table(useNA = 'ifany')
#            Diet
# Timepoint    AL  1D  2D  20  40
#   05 months 158 159 158 151 144
#   08 months  30  28  32  30  30
#   10 months 185 177 184 184 179
#   16 months 169 164 169 166 166
#   22 months 127 137 141 142 157
#   28 months  74  94  92 108 121
#   34 months  20  38  51  60  91
#   37 months   7  14  29  35  68

PHENO_DATA %>% 
  left_join(
    ANIMAL_DATA,
    by = 'MouseID'
  ) %>% 
  group_by(
    Timepoint, Diet
  ) %>% 
  summarise(
    N = n()
  ) %>% 
  group_by(
    Timepoint
  ) %>% 
  mutate(
    N = N/sum(N)
  ) %>% 
  ungroup() %>% 
  pivot_wider(
    names_from = 'Diet',
    values_from = 'N'
  )
#   Timepoint     AL   `1D`  `2D`  `20`  `40`
# 1 05 months 0.205  0.206  0.205 0.196 0.187
# 2 08 months 0.2    0.187  0.213 0.2   0.2  
# 3 10 months 0.204  0.195  0.202 0.202 0.197
# 4 16 months 0.203  0.197  0.203 0.199 0.199
# 5 22 months 0.180  0.195  0.200 0.202 0.223
# 6 28 months 0.151  0.192  0.188 0.221 0.247
# 7 34 months 0.0769 0.146  0.196 0.231 0.35 
# 8 37 months 0.0458 0.0915 0.190 0.229 0.444

TIMEPOINTS_FOR_ANALYSIS <- paste0(c('05', '08', '10', '16', '22', '28'), ' months')

#####


################################################################################
##### Set up info table #####

# Exclude: None

PHENO_NAMES <- c(
  'FrailtyMean',
  'Alopecia',
  'LossOfFurColor',
  'Dermatitis',
  'LossOfWhiskers',
  'CoatCondition',
  'Tumours',
  'DistendedAbdomen',
  'Kyphosis',
  'TailStiffening',
  'GaitDisorders',
  'Tremor',
  'BodyCondition',
  'RightingReflex',
  'VestibDisturbance',
  'Cataracts',
  'EyeDischargeOrSwelling',
  'Microphthalmia',
  'CornealOpacity',
  'VisionLoss',
  'MenaceReflex',
  'NasalDischarge',
  'Malocclusions',
  'RectalProlapse',
  'VaginalUterineProlapse',
  'Diarrhea',
  'BreathingRateDepth',
  'Piloerection'
)

PHENO_LEVELS <- c(
  'FrailtyMean',
  'Alopecia',
  'LossOfFurColor',
  'Dermatitis',
  'LossOfWhiskers',
  'CoatCondition',
  'Tumours',
  'DistendedAbdomen',
  'Kyphosis',
  'TailStiffening',
  'GaitDisorders',
  'Tremor',
  'BodyCondition',
  'RightingReflex',
  'VestibDisturbance',
  'Cataracts',
  'EyeDischargeOrSwelling',
  'Microphthalmia',
  'CornealOpacity',
  'VisionLoss',
  'MenaceReflex',
  'NasalDischarge',
  'Malocclusions',
  'RectalProlapse',
  'VaginalUterineProlapse',
  'Diarrhea',
  'BreathingRateDepth',
  'Piloerection'
)

DESCRIPTIONS <- tibble(
  Phenotype = PHENO_NAMES,
  Description = c(
    'composite (mean of component indices) frailty score (possible values range from 0 to 1)',
    
    'Fur loss (0, 0.5, 1)',
    'Loss of fur color (0, 0.5, 1)',
    'Skin lesions (0, 0.5, 1)',
    'Reduction in number of whiskers (0, 0.5, 1)',
    'Poor grooming (0, 0.5, 1)',
    'Presence of tumors (0, 0.5, 1)',
    'Bulging abdomen (0, 0.5, 1)',
    'Curvature of spine (0, 0.5, 1)',
    'Tail stiff, does not curl around finger (0, 0.5, 1)',
    'Abnormalities in gait (0, 0.5, 1)',
    'Presence of tremor (0, 0.5, 1)',
    'Lacking flesh around hips (0, 0.5, 1)',
    'Poor righting reflex (0, 0.5, 1)',
    'Vestibular disturbance when lowered by tail towards a flat surface (0, 0.5, 1)',
    'Presence of cataracts',
    'Eye discharge or swelling (0, 0.5, 1)',
    'Small or sunken eyes (0, 0.5, 1)',
    'White spots or clouding of cornea (0, 0.5, 1)',
    'Vision loss, height at which mouse reaches towards a surface when lowered (0, 0.5, 1)',
    'Lack of menace reflex, does not blink when object moved towards face (0, 0.5, 1)',
    'Nasal discharge (0, 0.5, 1)',
    'Uneven or overgrown teeth (0, 0.5, 1)',
    'Rectal prolapse (0, 0.5, 1)',
    'Vaginal or uterine prolapse (0, 0.5, 1)',
    'Signs of diarrhea',
    'Gasping or abnormal rate/depth of breathing (0, 0.5, 1)',
    'Piloerection (0, 0.5, 1)'
  )
)


TP_COUNT <- PHENO_DATA %>% 
  select(
    MouseID, Timepoint,
    all_of(PHENO_NAMES)
  ) %>% 
  pivot_longer(
    cols = all_of(PHENO_NAMES),
    names_to = 'Phenotype',
    values_to = 'Value'
  ) %>% 
  group_by(
    Phenotype, Timepoint
  ) %>% 
  summarise(
    N = sum(!is.na(Value))
  ) %>% 
  ungroup() %>% 
  mutate(
    Timepoint = paste0('NAll_', gsub(' ', '', Timepoint))
  ) %>% 
  pivot_wider(
    names_from = 'Timepoint',
    values_from = 'N'
  )


TP_DIET_COUNT <- ANIMAL_DATA %>% 
  select(
    MouseID, Diet
  ) %>% 
  left_join(
    PHENO_DATA %>% 
      select(
        MouseID, Timepoint,
        all_of(PHENO_NAMES)
      ),
    by = 'MouseID'
  ) %>% 
  pivot_longer(
    cols = all_of(PHENO_NAMES),
    names_to = 'Phenotype',
    values_to = 'Value'
  ) %>% 
  filter(
    !is.na(Timepoint)
  ) %>% 
  group_by(
    Phenotype, Diet, Timepoint
  ) %>% 
  summarise(
    N = sum(!is.na(Value))
  ) %>% 
  ungroup() %>% 
  mutate(
    Timepoint = paste0('N', as.character(Diet), '_', gsub(' ', '', Timepoint))
  ) %>% 
  select(
    -Diet
  ) %>% 
  pivot_wider(
    names_from = 'Timepoint',
    values_from = 'N'
  )


SMRY_STATS <- PHENO_DATA %>% 
  select(
    MouseID, all_of(PHENO_NAMES)
  ) %>% 
  pivot_longer(
    cols = all_of(PHENO_NAMES),
    names_to = 'Phenotype',
    values_to = 'Value'
  ) %>% 
  group_by(
    Phenotype
  ) %>% 
  mutate(
    Skewness = e1071::skewness(Value, na.rm = TRUE),
    Kurtosis = e1071::kurtosis(Value, na.rm = TRUE),
    Transform = 'none',
    Value = ifelse(
      Transform == 'none',
      Value,
      NA
    )
  ) %>% 
  summarise(
    Skewness = unique(Skewness),
    Kurtosis = unique(Kurtosis),
    Transform = unique(Transform),
    SD = sd(Value, na.rm = TRUE),
    Mean = mean(Value, na.rm = TRUE),
    Median = median(Value, na.rm = TRUE),
    Min = min(Value, na.rm = TRUE),
    Max = max(Value, na.rm = TRUE)
  )


COVARIATES <- tibble(
  Phenotype = PHENO_NAMES,
  REF_Cov = c('MouseID:DateCollect:Tech:Coat'),
  FEF_Cov = 'BW'
)


PHENO_INFO_TABLE <- DESCRIPTIONS %>% 
  left_join(
    TP_COUNT,
    by = 'Phenotype'
  ) %>% 
  left_join(
    TP_DIET_COUNT,
    by = 'Phenotype'
  ) %>% 
  left_join(
    SMRY_STATS,
    by = 'Phenotype'
  ) %>% 
  left_join(
    COVARIATES,
    by = 'Phenotype'
  ) %>% 
  mutate(
    Assay = ASSAY,
    Phenotype = factor(Phenotype, levels = PHENO_LEVELS)
  ) %>% 
  select(
    Assay, Phenotype,
    everything()
  ) %>% 
  arrange(
    Phenotype
  )
  
rm(DESCRIPTIONS, TP_COUNT, TP_DIET_COUNT, SMRY_STATS, COVARIATES, PHENO_NAMES, PHENO_LEVELS)



#####


################################################################################
##### List of transformation function #####

TRANSFORM_FUNCTIONS <- list(
  NA
)

#####


################################################################################
##### Plot histograms of phenotypes #####

PLOT_DATA <- PHENO_DATA

pdf(
  paste0(
    'figures/', ASSAY,
    '_histograms_of_phenotypes_and_transformations.pdf'
  ),
  width = 7, height = 4 + (length(PHENO_INFO_TABLE$Phenotype) - 1)*1.5
)
par(
  mfrow = c(length(PHENO_INFO_TABLE$Phenotype), 2)
)
for(PHENO in PHENO_INFO_TABLE$Phenotype){
  X <- PLOT_DATA[[PHENO]]
  XTF <- PLOT_DATA[[PHENO]]
  TRANSFORM <- PHENO_INFO_TABLE$Transform[PHENO_INFO_TABLE$Phenotype == PHENO]

  if(TRANSFORM != 'none'){
    XTF <- TRANSFORM_FUNCTIONS[[TRANSFORM]](XTF)
  }
  
  hist(
    X,
    main = paste0(
      PHENO, 
      ', skewness = ', 
      signif(PHENO_INFO_TABLE$Skewness[PHENO_INFO_TABLE$Phenotype == PHENO], 3), 
      ', kurtosis = ', 
      signif(PHENO_INFO_TABLE$Kurtosis[PHENO_INFO_TABLE$Phenotype == PHENO], 3)
    ),
    xlab = NULL,
    breaks = 100
  )
  
  
  
  hist(
    XTF,
    main = paste0(
      PHENO, ', ', 
      ifelse(
        TRANSFORM == 'none', 'no transform',
        paste0(TRANSFORM, ' transform')
      )
    ),
    xlab = NULL,
    breaks = 100
  )
  
  rm(X, XTF, TRANSFORM)
}
dev.off()

rm(PHENO, PLOT_DATA)

#####


################################################################################
##### Prepare phenotype and covariate table #####

DATA <- PHENO_DATA %>% 
  left_join(
    ANIMAL_DATA %>% 
      select(
        MouseID, Diet, Coat
      ),
    by = 'MouseID'
  ) %>% 
  select(
    MouseID,
    all_of(PHENO_INFO_TABLE$Phenotype),
    Age, Timepoint, Diet,
    BW,
    DateCollect, Tech, Coat
  )

#####


################################################################################
##### Combine into dataset #####

DATASET <- list(
  Assay = ASSAY,
  Timepoints = TIMEPOINTS_FOR_ANALYSIS,
  Info_Table = PHENO_INFO_TABLE,
  Data = DATA,
  Transform_Functions = TRANSFORM_FUNCTIONS
)

#####


################################################################################
#####  #####

#####


################################################################################
##### Output #####

write_csv(
  PHENO_INFO_TABLE,
  paste0(
    'results/', ASSAY,
    '_phenotype_info_table.csv'
  )
)

saveRDS(
  DATASET,
  paste0(
    'data/', ASSAY,
    '_diet_age_dataset.rds'
  )
)

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####