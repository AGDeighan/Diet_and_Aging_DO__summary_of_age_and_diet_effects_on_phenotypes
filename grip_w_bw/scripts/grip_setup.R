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
  'data/Grip_Processed_Raw_20201221.csv'
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

ASSAY <- 'Grip'

#####


################################################################################
##### Select timepoints to test for diet and diet-age effects #####

# We will exclude the 35 month timepoint when estimating direct diet 
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
#   05 months 188 188 190 188 181
#   11 months 185 162 168 184 168
#   17 months 169 163 169 166 166
#   23 months 125 133 137 134 157
#   29 months  73  92  88 107 121
#   35 months  18  34  45  56  91

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
#  Timepoint    AL   `1D`  `2D`  `20`  `40`
#  05 months 0.201  0.201 0.203 0.201 0.194
#  11 months 0.213  0.187 0.194 0.212 0.194
#  17 months 0.203  0.196 0.203 0.199 0.199
#  23 months 0.182  0.194 0.200 0.195 0.229
#  29 months 0.152  0.191 0.183 0.222 0.252
#  35 months 0.0738 0.139 0.184 0.230 0.373

TIMEPOINTS_FOR_ANALYSIS <- paste0(c('05', '11', '17', '23', '29'), ' months')

#####


################################################################################
##### Set up info table #####

# Exclude: None, we are really only interested in "All" and "Fore" (the mean
# all paw and fore paw strengths), but we may as well test the three individual
# trials also just to verify that the results are consistent.

PHENO_NAMES <- c(
  'All',   
  'Fore',  
  'All_1', 
  'All_2', 
  'All_3', 
  'Fore_1',
  'Fore_2',
  'Fore_3'
)

PHENO_LEVELS <- c(
  'All',   
  'Fore',  
  'All_1', 
  'All_2', 
  'All_3', 
  'Fore_1',
  'Fore_2',
  'Fore_3'
)

DESCRIPTIONS <- tibble(
  Phenotype = PHENO_NAMES,
  Description = c(
    'The trimmed mean of all-paw grip strength (grams)',
    'The trimmed mean of fore-paw grip strength (grams)',
    'The all-paw grip strength from all-paw trial 1 (grams)',
    'The all-paw grip strength from all-paw trial 2 (grams)',
    'The all-paw grip strength from all-paw trial 3 (grams)',
    'The fore-paw grip strength from fore-paw trial 1 (grams)',
    'The fore-paw grip strength from fore-paw trial 2 (grams)',
    'The fore-paw grip strength from fore-paw trial 3 (grams)'
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
  REF_Cov = c('MouseID:DateCollect'),
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
        MouseID, Diet, Cohort
      ),
    by = 'MouseID'
  ) %>% 
  select(
    MouseID,
    all_of(PHENO_INFO_TABLE$Phenotype),
    Age, Timepoint, Diet,
    DateCollect,
    BW
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