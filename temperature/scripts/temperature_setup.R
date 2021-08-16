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
  'data/BodyTemp_Processed_20201125.csv'
) %>% 
  mutate(
    DateCollect = as.character(DateCollect),
    AgeInDays = AgeInDays/30.4,
    BW = ifelse(
      is.na(BW_Weekly), BW_Test,
      BW_Weekly
    ),
    BW = standardise(BW),
    TemperatureCondBW = Temperature_Raw
  ) %>%
  select(-c(
    Temperature_AdjBatchCoat,
    BW_Test, BW_Weekly
  )) %>%  
  rename(
    Temperature = Temperature_Raw,
    Age = AgeInDays
  )

#####


################################################################################
## Specify assay #########################################################

ASSAY <- 'Temperature'

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

# Exclude: None, we are really only interested in "All" and "Fore" (the mean
# all paw and fore paw strengths), but we may as well test the three individual
# trials also just to verify that the results are consistent.

PHENO_NAMES <- c(
  'Temperature',
  'TemperatureCondBW'
)

PHENO_LEVELS <- c(
  'Temperature',
  'TemperatureCondBW'
)

DESCRIPTIONS <- tibble(
  Phenotype = PHENO_NAMES,
  Description = c(
    'body temperature (Celsius)',
    'body temperature (Celsius)'
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
  REF_Cov = c('MouseID:DateCollect:Coat'),
  FEF_Cov = c(NA, 'BW')
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
  width = 7, height = 4 + (length(PHENO_INFO_TABLE$Phenotype[1]) - 1)*1.5
)
par(
  mfrow = c(length(PHENO_INFO_TABLE$Phenotype[1]), 2),
  cex.main = 0.75
)
for(PHENO in PHENO_INFO_TABLE$Phenotype[1]){
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
    DateCollect, Coat
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