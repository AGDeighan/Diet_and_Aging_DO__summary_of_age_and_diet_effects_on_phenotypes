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
  'data/PIXI_Processed_Raw_20201022.csv'
) %>% 
  mutate(
    DateCollect = as.character(DateCollect),
    AgeInDays = AgeInDays/30.4
  ) %>% 
  rename(
    Age = AgeInDays
  )

#####


################################################################################
## Specify assay #########################################################

ASSAY <- 'PIXI'

#####


################################################################################
##### Select timepoints to test for diet and diet-age effects #####

# The 34-month timepoint is imbalanced in terms of diet

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
#   10 months 184 173 178 184 176
#   22 months 124 133 136 140 156
#   34 months  17  35  42  56  91

TIMEPOINTS_FOR_ANALYSIS <- paste0(c('10', '22'), ' months')

#####


################################################################################
##### Set up info table #####

# Exclude: None, I don't know enough about echo data to confidently say which
# phenotypes are redundant or superfluous 

PHENO_NAMES <- c(
  'TissueArea',		
  'RST',
  
  'TTM',			      
  'FTM',			      
  'LTM',
  
  'PercFat',			  
  'PercLean',	
  
  'BoneArea',		  
  'BMC',			      
  'BMD'		      
)

PHENO_LEVELS <- c(
  'TissueArea',		
  'RST',
  
  'TTM',			      
  'FTM',			      
  'LTM',
  
  'PercFat',			  
  'PercLean',	
  
  'BoneArea',		  
  'BMC',			      
  'BMD'
)

DESCRIPTIONS <- tibble(
  Phenotype = PHENO_NAMES,
  Description = c(
    'The soft tissue area (cm^2). Directly measured',
    'The R-value of soft tissue. Directly measured',
    
    'Total soft tissue (excluding bone) mass (grams). I believe this is directly measured but it may be calculated from RST and Tissue_Area',
    'Fat tissue mass (grams). Calculated during data processing: (PercFat/100) * (TTM)',
    'Lean tissue mass (grams). Calculated during data processing: ((100 - PercFat)/100) * TTM',
    
    'Percent of total soft tissue mass that is fat tissue. Calculated from RST during initial data generation',
    'Percent of total soft tissue mass that is lean tissue. Calculated during data processing: 100 - PercFat',
    
    'Area of bone tissue (cm^2). Directly measured',
    'Bone mineral content. The total mineral content of all bone tissue (g). Directly measured',
    'Bone mineral density (~grams/cm^2). This is BMC that has been regression-adjusted (log[BMC] ~ log[BoneArea]) for bone area and then shifted to be on the same scale as simple ratio adjustment (BMC/BoneArea). See the data processing R-script for details.'
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
  FEF_Cov = NA
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
    DateCollect
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