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

# All the FACS variables are percents so we will use the logit transform

logit_of_perc_add <- function(x, shift = 0.01, inverse = FALSE){
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

logit_of_perc_sub <- function(x, shift = -0.01, inverse = FALSE){
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
  'data/FACS_Processed_Raw_20201125.csv'
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

ASSAY <- 'FACS'

#####


################################################################################
##### Select timepoints to analyse #####

# The 28-month timepoint is imbalanced in terms of diet, but there still
# quite a few ad libitum mice (72), so we will use all three timepoints for
# the analysis. That being said, we should keep this imbalance in mind when
# interpreting results

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
#   05 months 188 188 190 189 182
#   16 months 166 163 168 167 166
#   28 months  72  93  92 107 118

TIMEPOINTS_FOR_ANALYSIS <- paste0(c('05', '16', '28'), ' months')

#####


################################################################################
##### Set up info table #####

# Exclude: None, I don't know enough about FACS data to confidently say which
# phenotypes are redundant or superfluous 

PHENO_NAMES <- c(
  'Lymph_PercViable', 'Myeloid_PercViable', 'CLLlike_PercViable', 'Other_PercViable', 
  'B_PercLymph', 'B_PercViable', 'B_CD11BposPercB', 'B_CD62LposPercB', 'B_NKG2DposPercB', 
  'T_PercViable', 
  'CD4T_PercLymph', 'CD4T_PercViable', 'CD4T_CD25posPercCD4T', 'CD4T_CD62LnegCD44negPercCD4T', 'CD4T_CD62LnegCD44posPercCD4T', 'CD4T_CD62LposCD44negPercCD4T', 'CD4T_CD62LposCD44posPercCD4T', 'CD4T_NKG2DposPercCD4T', 
  'CD8T_PercLymph', 'CD8T_PercViable', 'CD8T_CD62LnegCD44negPercCD8T', 'CD8T_CD62LnegCD44posPercCD8T', 'CD8T_CD62LposCD44negPercCD8T', 'CD8T_CD62LposCD44posPercCD8T', 'CD8T_NKG2DposPercCD8T', 
  'DNT_B220posPercDNT', 'DNT_CD25posPercDNT', 'DNT_CD62LnegCD44negPercDNT', 'DNT_CD62LnegCD44posPercDNT', 'DNT_CD62LposCD44negPercDNT', 'DNT_CD62LposCD44posPercDNT', 'DNT_NKG2DposPercDNT', 
  'NKG2DposT_PercLymph', 'NKG2DposT_PercViable', 'NKG2DposT_CD4posPercNKG2DposT', 'NKG2DposT_CD8posPercNKG2DposT', 'NKG2DposT_DNPercNKG2DposT', 
  'NK_PercLymph', 'NK_PercViable', 'NK_B220posPercNK', 'NK_CD62LposPercNK', 'NK_CD11BnegCD11CnegPercNK', 'NK_CD11BnegCD11CposPercNK', 'NK_CD11BposCD11CnegPercNK', 'NK_CD11BposCD11CposPercNK', 'NK_CD11BnegCD11CnegPercCD62LposNK', 'NK_CD11BnegCD11CposPercCD62LposNK', 'NK_CD11BposCD11CnegPercCD62LposNK', 'NK_CD11BposCD11CposPercCD62LposNK', 
  'Neut_PercMyeloid', 'Neut_PercViable', 'Neut_CD62LposPercNeut', 
  'Mono_PercMyeloid', 'Mono_PercViable', 'Mono_ResidPercMono', 'Mono_InflPercMono', 'Mono_CD11CposCD62LposPercMono', 'Mono_OtherPercMono', 
  'Eos_PercMyeloid', 'Eos_PercViable'
)

PHENO_LEVELS <- c(
  # Viable cells
  'Lymph_PercViable', 'Myeloid_PercViable', 'CLLlike_PercViable', 'Other_PercViable', 
  'B_PercViable', 
  'T_PercViable', 'CD4T_PercViable', 'CD8T_PercViable', 
  'NKG2DposT_PercViable', 
  'NK_PercViable', 
  'Neut_PercViable', 'Mono_PercViable', 'Eos_PercViable',
  
  # Lymphoid cells
  'B_PercLymph', 'CD4T_PercLymph', 'CD8T_PercLymph', 
  'NKG2DposT_PercLymph', 
  'NK_PercLymph', 
  
  # Myeloid cells
  'Neut_PercMyeloid', 'Mono_PercMyeloid', 'Eos_PercMyeloid',
  
  # B cells
  'B_CD11BposPercB', 'B_CD62LposPercB', 'B_NKG2DposPercB',
  
  # CD4 T cells
  'CD4T_CD25posPercCD4T', 
  'CD4T_CD62LnegCD44negPercCD4T', 'CD4T_CD62LnegCD44posPercCD4T', 'CD4T_CD62LposCD44negPercCD4T', 'CD4T_CD62LposCD44posPercCD4T', 
  'CD4T_NKG2DposPercCD4T',
  
  # CD8 T cells
  'CD8T_CD62LnegCD44negPercCD8T', 'CD8T_CD62LnegCD44posPercCD8T', 'CD8T_CD62LposCD44negPercCD8T', 'CD8T_CD62LposCD44posPercCD8T', 
  'CD8T_NKG2DposPercCD8T',
  
  # Double negative T cells
  'DNT_B220posPercDNT', 'DNT_CD25posPercDNT', 
  'DNT_CD62LnegCD44negPercDNT', 'DNT_CD62LnegCD44posPercDNT', 'DNT_CD62LposCD44negPercDNT', 'DNT_CD62LposCD44posPercDNT', 
  'DNT_NKG2DposPercDNT',
  
  # NKG2D positive T cells
  'NKG2DposT_CD4posPercNKG2DposT', 'NKG2DposT_CD8posPercNKG2DposT', 'NKG2DposT_DNPercNKG2DposT',
  
  # Natural killer cells
  'NK_B220posPercNK', 'NK_CD62LposPercNK', 
  'NK_CD11BnegCD11CnegPercNK', 'NK_CD11BnegCD11CposPercNK', 'NK_CD11BposCD11CnegPercNK', 'NK_CD11BposCD11CposPercNK', 
  
  # CD62L positive natural killer cells
  'NK_CD11BnegCD11CnegPercCD62LposNK', 'NK_CD11BnegCD11CposPercCD62LposNK', 'NK_CD11BposCD11CnegPercCD62LposNK', 'NK_CD11BposCD11CposPercCD62LposNK',
  
  # Neutrophils
  'Neut_CD62LposPercNeut',
  
  # Monocytes
  'Mono_ResidPercMono', 'Mono_InflPercMono', 'Mono_CD11CposCD62LposPercMono', 'Mono_OtherPercMono'
)

DESCRIPTIONS <- tibble(
  Phenotype = PHENO_NAMES,
  Description = c(
    '% of viable cells that are lymphocytes', 
    '% of viable cells that are myeloid cells (CD11B+ cells)', 
    '% of viable cells that are chronic lymphocytic leukemia (CLL) like', 
    '% percent of viable cells that are not one of the other populations, in most cases RBCs that snuck through the gating on samples where there were a lot which did not lyse', 
    
    '% of lymphocytes that are B cells', 
    '% of viable cells that are B cells', 
    '% of B cells that are CD11+', 
    '% of B cells that are CD62L+', 
    '% of B cells that are NKG2D+', 
    
    '% of viable cells that are T cells', 
    
    '% of lymphocytes that are CD4+ T cells', 
    '% of viable cells that are CD4+ T cells', 
    '% of CD4+ T cells that are CD25+', 
    '% of CD4+ T cells that are CD62L- and CD44-', 
    '% of CD4+ T cells that are CD62L- and CD44+', 
    '% of CD4+ T cells that are CD62L+ and CD44-', 
    '% of CD4+ T cells that are CD62L+ and CD44+', 
    '% of CD4+ T cells that are NKG2D+', 
    
    '% of lymphocytes that are CD8+ T cells', 
    '% of viable cells that are CD8+ T cells', 
    '% of CD8+ T cells that are CD62L- and CD44-', 
    '% of CD8+ T cells that are CD62L- and CD44+', 
    '% of CD8+ T cells that are CD62L+ and CD44-', 
    '% of CD8+ T cells that are CD62L+ and CD44+', 
    '% of CD8+ T cells that are NKG2D+', 
    
    '% of double negative T cells that are B220+', 
    '% of double negative T cells that are CD25+', 
    '% of double negative T cells that are CD62L- and CD44-', 
    '% of double negative T cells that are CD62L- and CD44+', 
    '% of double negative T cells that are CD62L+ and CD44-', 
    '% of double negative T cells that are CD62L+ and CD44+', 
    '% of double negative T cells that are NKG2D+', 
    
    '% of lymphocytes that are NKG2D+ T cells', 
    '% of of viable cells that are NKG2D+ T cells', 
    '% of NKG2D+ T cells that are CD4+', 
    '% of NKG2D+ T cells that are CD8+', 
    '% of NKG2D+ T cells that are double negative (CD4- and CD8-)', 
    
    '% of lymphocytes that are natural killer cells', 
    '% of viable cells that are natural killer cells', 
    '% of natural killer cells that are B220+', 
    '% of natural killer cells that are CD62L+', 
    '% of natural killer cells that are CD11B- and CD11C-', 
    '% of natural killer cells that are CD11B- and CD11C+', 
    '% of natural killer cells that are CD11B+ and CD11C-', 
    '% of natural killer cells that are CD11B+ and CD11C+', 
    '% of CD62L+ natural killer cells that are CD11B- and CD11C-', 
    '% of CD62L+ natural killer cells that are CD11B- and CD11C+', 
    '% of CD62L+ natural killer cells that are CD11B+ and CD11C-', 
    '% of CD62L+ natural killer cells that are CD11B+ and CD11C+', 
    
    '% of myeloid cells that are neutrophils', 
    '% of viable cells that are neutrophils', 
    '% of neutrophils that are CD62L+', 
    
    '% of myeloid cells that are monocytes', 
    '% of viable cells that are monocytes', 
    '% of monocytes that are resident monocytes (CD11C+ and CD62L-)', 
    '% of monocytes that are inflammatory monocytes (CD11C- and CD62L+)', 
    '% of monocytes that are CD11C+ and CD62L+', 
    '% of monocytes that are "other" monocytes (CD11C- and CD62L-)', 
    
    '% of myeloid cells that are eosinophils', 
    '% of viable cells that are eosinophils'
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
    Transform = ifelse(
      max(Value, na.rm = TRUE) >= 99.99,
      'logit_neg_0.01',
      'logit_pos_0.01'
    ),
    Value = ifelse(
      Transform == 'logit_neg_0.01',
      logit_of_perc_sub(Value),
      logit_of_perc_add(Value)
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
  logit_pos_0.01 = logit_of_perc_add,
  logit_neg_0.01 = logit_of_perc_sub
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
        MouseID, Diet
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