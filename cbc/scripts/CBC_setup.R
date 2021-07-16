# 2021-07-06
################################################################################
#
#   This script creates a phenotype information table to be used for estimating
#   diet and age effects on the CBC phenotypes
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

logit_of_perc <- function(x, shift = 0.01, inverse = FALSE){
  if(inverse){
    # returns inverse log-odds of x
    ILT <- exp(x)/(1+exp(x))
    PERC <- ILT * 100
    UNSHIFT <- PERC - shift
    return(UNSHIFT)
  } else{
    # returns log-odds of x
    SHIFTED <- x + shift
    PROP <- SHIFTED/100
    LT <- log(PROP/(1-PROP))
    return(LT)
  }
}

log_tf <- function(x, shift = 0.01, inverse = FALSE){
  if(inverse){
    # returns inverse log-odds of x
    ILT <- exp(x)
    UNSHIFT <- ILT - shift
    return(UNSHIFT)
  } else{
    # returns log-odds of x
    SHIFTED <- x + shift
    LT <- log(SHIFTED)
    return(LT)
  }
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
  'data/CBC_Processed_Raw_20201006.csv'
) %>% 
  mutate(
    DateCollect = as.character(DateCollect),
    AgeInDays = AgeInDays/30.4,
    NumClumps = standardise(log_tf(NumClumps))
  ) %>% 
  rename(
    Age = AgeInDays
  )

#####


################################################################################
## Specify assay #########################################################

ASSAY <- 'CBC'

#####


################################################################################
##### Select timepoints to analyse #####

# The 35-month timepoint is imbalanced in terms of diet

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
#   11 months 182 172 179 183 176
#   23 months 113 129 133 137 153
#   35 months  14  28  37  50  79

TIMEPOINTS_FOR_ANALYSIS <- paste0(c(11, 23), ' months')

#####


################################################################################
##### Set up info table #####

# Exclude: 
#   CalcHgb: calculated/derived version of Hgb, (CHCM * NumRBC * MCV) / 1000
#   MCHC: calculated/derived version of CHCM, (Hgb / (NumRBC * MCV)) * 1000 
#   MCH: calculated/derived version of CH, (Hgb / NumRBC) * 10 
#   RDWsd: standard deviation of version of RDW, using coefficient of variation
#   HDWsd: standard deviation of version of HDW, using coefficient of variation 
#   PDWsd: standard deviation of version of PDW, using coefficient of variation

PHENO_NAMES <- c(
  'NumRBC', 'NumRetic', 'PercRetic', 
  'Hct', 'Hgb', 
  'MCV', 'CHCM',  'CH',  'CHm', 'CHr', 
  'RDWcv', 'HDWcv', 
  'NumWBC', 'NumLymph', 'NumNeut', 'NumMono', 'NumEos', 
  'NLR', 'PercLymph', 'PercNeut', 'PercMono', 'PercEos', 
  'NumPlt', 'MPV', 'MPM', 'PDWcv'
)

CLUMP_PHENOS <- c(
  'NumEos', 'PercEos',
  'NumPlt', 'MPV', 'MPM', 'PDWcv', 'PDWsd'
)

DESCRIPTIONS <- tibble(
  Phenotype = PHENO_NAMES,
  Description = c(
    'RBC count (10^6/uL)',
    'reticulocyte count (10^6/uL)',
    'percent reticulocytes of total RBCs (%) = (NumRetic / NumRBC) * 100',
    'hematocrit (%) = (NumRBC * MCV) / 10',
    'blood hemoglobin concentration (g/dL)',
    'mean corpuscular (RBC) volume (fL)',
    'mean cellular hemoglobin concentration (g/dL)',
    'mean cellular hemoglobin content (pg)',
    'mean cellular hemoglobin content of mature RBCs (pg)',
    'mean cellular hemoglobin content of reticulocytes (pg)',
    'coefficient of variation of RBC volume (%)',
    'coefficient of variation of cellular hemoglobin concentration (%) = (HDWsd / CHCM) * 100%',
    
    'total WBC count (10^3/uL) = NumLymph + NumNeut + NumMono + NumEos',
    'lymphocyte count (10^3/uL)',
    'neutrophil count (10^3/uL)',
    'monocyte count (10^3/uL)',
    'eosinophil count (10^3/uL)',
    'ratio of neutrophils to lymphocytes = NumNeut / NumLymph',
    'percent lymphocytes (%) = NumLymph / NumWBC',
    'percent neutrophils (%) = NumNeut / NumWBC',
    'percent monocytes (%) = NumMono / NumWBC',
    'percent eosinophils (%) = NumEos / NumWBC',
    
    'platelet count (10^3/uL)',
    'mean platelet volume (fL)',
    'mean platelet mass (pg)',
    'coefficient of variation of platelet volume (%)'
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

# Skewness:
#   Negative: left-skewed
#   Positive: right-skewed
# 
# Kurtosis:
#   Negative: thinner tails than normal distribution
#   Positive: thicker tails than normal distribution
SMRY_STATS <- PHENO_DATA %>% 
  filter(
    Timepoint %in% TIMEPOINTS_FOR_ANALYSIS
  ) %>% 
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
    Transform = ifelse(
      substr(Phenotype, 1, 4) == 'Perc', 'logit',
      ifelse(
        Skewness > 1, 'log',
        'none'
      )
    ),
    Value = ifelse(
      Transform == 'log', log_tf(Value),
      ifelse(
        Transform == 'logit', logit_of_perc(Value),
        Value
      )
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
  REF_Cov = c('Cohort:MouseID:DateCollect'),
  FEF_Cov = ifelse(Phenotype %in% CLUMP_PHENOS, 'NumClumps', as.character(NA))
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
    Phenotype = factor(Phenotype, levels = PHENO_NAMES)
  ) %>% 
  select(
    Assay, Phenotype,
    everything()
  ) %>% 
  arrange(
    Phenotype
  )
  
rm(DESCRIPTIONS, TP_COUNT, TP_DIET_COUNT, SMRY_STATS, COVARIATES, CLUMP_PHENOS, PHENO_NAMES)



#####


################################################################################
##### List of transformation function #####

TRANSFORM_FUNCTIONS <- list(
  logit = logit_of_perc,
  log = log_tf
)

#####


################################################################################
##### Plot histograms of phenotypes #####

PLOT_DATA <- PHENO_DATA %>% 
  filter(
    Timepoint %in% TIMEPOINTS_FOR_ANALYSIS
  )

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
    Cohort, DateCollect,
    NumClumps
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