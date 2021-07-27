# 2021-07-06
################################################################################
#
#   This estimates the direct effects of diet and age and the diet and age 
#   interaction (diet-specific age effects) on each phenotype
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
library(e1071)
library(lme4)
library(pbkrtest)
library(emmeans)

################################################################################
## helper functions #########################################################

#####


################################################################################
## load dataset #########################################################

DATASET <- readRDS('data/Grip_diet_age_dataset.rds')

#####


################################################################################
##### Estimate diet effects stratified by timepoint #####

DIET_EFFECTS <- tibble(
  Phenotype = vector(),
  Timepoint = vector(),
  AL_Mean = vector(),
  AL_SE = vector(),
  `1D_Mean` = vector(),
  `1D_SE` = vector(),
  `2D_Mean` = vector(),
  `2D_SE` = vector(),
  `20_Mean` = vector(),
  `20_SE` = vector(),
  `40_Mean` = vector(),
  `40_SE` = vector(),
  Diet_FTestStat = vector(),
  Diet_FTestNumDF = vector(),
  Diet_FTestDenDF = vector(),
  Diet_FTestScaling = vector(),
  Diet_FTestPValue = vector(),
  Diet_PairwiseDifferences = vector()
)

for(PHENO in DATASET$Info_Table$Phenotype){
  cat('\n\n'); cat(PHENO); cat('\n')
  
  for(TP in DATASET$Timepoints){
    cat(TP); cat('\n')
    ### Set-up data set
    VOI <- 'Diet'
    REF_TERMS <- as.character(str_split(
      DATASET$Info_Table$REF_Cov[DATASET$Info_Table$Phenotype == PHENO],
      ':',
      simplify = TRUE
    ))
    REF_TERMS <- REF_TERMS[REF_TERMS != 'MouseID']
    FEF_TERMS <- as.character(na.omit(str_split(
      DATASET$Info_Table$FEF_Cov[DATASET$Info_Table$Phenotype == PHENO],
      ':',
      simplify = TRUE
    )))
    
    MODEL_DATA <- DATASET$Data %>%
      filter(
        Timepoint == TP
      ) %>%
      select(
        MouseID,
        all_of(na.omit(c(PHENO, VOI, REF_TERMS, FEF_TERMS)))
      )
    
    ### Create formulas
    REF_TERMS <- paste0('(1|', REF_TERMS, ')', collapse = ' + ')
    FEF_TERMS <- paste0(FEF_TERMS, collapse = ' + ')
    if(FEF_TERMS == ''){FEF_TERMS <- NULL}
    
    NULL_FORMULA <- as.formula(paste0(
      PHENO,
      ' ~ ',
      paste0(
        c(REF_TERMS, FEF_TERMS),
        collapse = ' + '
      )
    ))
    
    FULL_FORMULA <- as.formula(paste0(
      PHENO,
      ' ~ -1 + ',
      paste0(
        c(REF_TERMS, FEF_TERMS, VOI),
        collapse = ' + '
      )
    ))
    
    
    ### Transform phenotype if a transform is designated
    TRANSFORM <- DATASET$Info_Table$Transform[DATASET$Info_Table$Phenotype == PHENO]
    if(TRANSFORM != 'none'){
      MODEL_DATA[[PHENO]] <- DATASET$Transform_Functions[[TRANSFORM]](MODEL_DATA[[PHENO]])
    }
    
    
    ### Standardise phenotype
    PHENO_SD <- DATASET$Info_Table$SD[DATASET$Info_Table$Phenotype == PHENO]
    PHENO_MEAN <- DATASET$Info_Table$Mean[DATASET$Info_Table$Phenotype == PHENO]
    MODEL_DATA[[PHENO]] <- (MODEL_DATA[[PHENO]] - PHENO_MEAN)/PHENO_SD
    
    
    ### Fit models
    cat('\nFull model formula\n')
    print(FULL_FORMULA)
    FULL_MODEL <- lmer(
      FULL_FORMULA,
      lmerControl(
        optimizer = 'bobyqa'
      ),
      REML = TRUE,
      data = MODEL_DATA
    )
    cat('\nFull model random effect standard deviations\n')
    print(VarCorr(FULL_MODEL))
    
    cat('\nNull model formula\n')
    print(NULL_FORMULA)
    NULL_MODEL <- lmer(
      NULL_FORMULA,
      lmerControl(
        optimizer = 'bobyqa'
      ),
      REML = TRUE,
      data = MODEL_DATA
    )
    cat('\nNull model random effect standard deviations\n')
    print(VarCorr(NULL_MODEL))
    
    
    ### Pull estimated effects and standard errors from model
    EFFECT <- c(
      as.numeric(fixef(FULL_MODEL)['DietAL']),
      as.numeric(fixef(FULL_MODEL)['Diet1D']),
      as.numeric(fixef(FULL_MODEL)['Diet2D']),
      as.numeric(fixef(FULL_MODEL)['Diet20']),
      as.numeric(fixef(FULL_MODEL)['Diet40'])
    ) 
    SE <- c(
      sqrt(as.numeric(vcov(FULL_MODEL)['DietAL', 'DietAL'])),
      sqrt(as.numeric(vcov(FULL_MODEL)['Diet1D', 'Diet1D'])),
      sqrt(as.numeric(vcov(FULL_MODEL)['Diet2D', 'Diet2D'])),
      sqrt(as.numeric(vcov(FULL_MODEL)['Diet20', 'Diet20'])),
      sqrt(as.numeric(vcov(FULL_MODEL)['Diet40', 'Diet40']))
    )
    names(EFFECT) <- names(SE) <- c('AL', '1D', '2D', '20', '40')
    
    
    ### Un-standardize
    EFFECT <- (EFFECT * PHENO_SD) + PHENO_MEAN
    SE <- (SE * PHENO_SD)
    
    
    ### Test significance of effect
    TEST_RESULTS <- KRmodcomp(
      largeModel = FULL_MODEL,
      smallModel = NULL_MODEL
    )
    cat('\nTest results\n')
    print(TEST_RESULTS$test)
    
    FTEST_STAT <- TEST_RESULTS$test['Ftest', 'stat']
    FTEST_NUMDF <- TEST_RESULTS$test['Ftest', 'ndf']
    FTEST_DENDF <- TEST_RESULTS$test['Ftest', 'ddf']
    FTEST_SCALE <- TEST_RESULTS$test['Ftest', 'F.scaling']
    FTEST_PV <- TEST_RESULTS$test['Ftest', 'p.value']
    
    
    ### If the effect is significant, perform tukey test for different groups
    PAIRWISE_DIFF <- NA
    if(FTEST_PV < 0.05){
      PCOMP <- pairs(emmeans(FULL_MODEL, 'Diet'), adjust = "tukey") %>% 
        data.frame() %>% 
        filter(p.value < 0.05)
      PAIRWISE_DIFF <- paste0(gsub(' ', '', PCOMP$contrast), collapse = ':')
      rm(PCOMP)
    }
    
    
    ### Add results to data frame
    DF_ROW <- tibble(
      Phenotype = PHENO,
      Timepoint = TP,
      AL_Mean = EFFECT['AL'],
      AL_SE = SE['AL'],
      `1D_Mean` = EFFECT['1D'],
      `1D_SE` = SE['1D'],
      `2D_Mean` = EFFECT['2D'],
      `2D_SE` = SE['2D'],
      `20_Mean` = EFFECT['20'],
      `20_SE` = SE['20'],
      `40_Mean` = EFFECT['40'],
      `40_SE` = SE['40'],
      Diet_FTestStat = FTEST_STAT,
      Diet_FTestNumDF = FTEST_NUMDF,
      Diet_FTestDenDF = FTEST_DENDF,
      Diet_FTestScaling = FTEST_SCALE,
      Diet_FTestPValue = FTEST_PV,
      Diet_PairwiseDifferences = PAIRWISE_DIFF
    )
    
    DIET_EFFECTS <- rbind(
      DIET_EFFECTS,
      DF_ROW
    )
    
    
    ### Clean up
    rm(
      MODEL_DATA,
      VOI, FEF_TERMS, REF_TERMS, NULL_FORMULA, FULL_FORMULA,
      TRANSFORM,
      FULL_MODEL, NULL_MODEL, TEST_RESULTS,
      PHENO_SD, PHENO_MEAN,
      EFFECT, SE, 
      FTEST_STAT, FTEST_NUMDF, FTEST_DENDF, FTEST_SCALE, FTEST_PV,
      PAIRWISE_DIFF,
      DF_ROW
    )
  }
  
  cat('\n# ----------------------------- #\n')
}
rm(PHENO)

DIET_EFFECTS <- DIET_EFFECTS %>% 
  mutate(
    Phenotype = factor(Phenotype, levels = levels(DATASET$Info_Table$Phenotype))
  )

DATASET$Diet_Effects <- DIET_EFFECTS

#####


################################################################################
##### Estimate age effects #####

AGE_EFFECTS <- tibble(
  Phenotype = vector(),
  Age_Slope = vector(),
  Age_SE = vector(),
  Age_FTestStat = vector(),
  Age_FTestNumDF = vector(),
  Age_FTestDenDF = vector(),
  Age_FTestScaling = vector(),
  Age_FTestPValue = vector()
)

for(PHENO in DATASET$Info_Table$Phenotype){
  cat('\n\n'); cat(PHENO); cat('\n')
  
  ### Set-up data set
  VOI <- 'Age'
  REF_TERMS <- as.character(str_split(
    DATASET$Info_Table$REF_Cov[DATASET$Info_Table$Phenotype == PHENO],
    ':',
    simplify = TRUE
  ))
  REF_TERMS <- c('Diet', REF_TERMS)
  FEF_TERMS <- as.character(na.omit(str_split(
    DATASET$Info_Table$FEF_Cov[DATASET$Info_Table$Phenotype == PHENO],
    ':',
    simplify = TRUE
  )))
  
  MODEL_DATA <- DATASET$Data %>%
    # filter(
    #   Timepoint %in% DATASET$Timepoints
    # ) %>%
    select(
      MouseID, Timepoint,
      all_of(na.omit(c(PHENO, VOI, REF_TERMS, FEF_TERMS)))
    )
  
  ### Create formulas
  REF_TERMS <- paste0('(1|', REF_TERMS, ')', collapse = ' + ')
  FEF_TERMS <- paste0(FEF_TERMS, collapse = ' + ')
  if(FEF_TERMS == ''){FEF_TERMS <- NULL}
  
  NULL_FORMULA <- as.formula(paste0(
    PHENO,
    ' ~ ',
    paste0(
      c(REF_TERMS, FEF_TERMS),
      collapse = ' + '
    )
  ))
  
  FULL_FORMULA <- as.formula(paste0(
    PHENO,
    ' ~ ',
    paste0(
      c(REF_TERMS, FEF_TERMS, VOI),
      collapse = ' + '
    )
  ))
  
  
  ### Transform phenotype if a transform is designated
  TRANSFORM <- DATASET$Info_Table$Transform[DATASET$Info_Table$Phenotype == PHENO]
  if(TRANSFORM != 'none'){
    MODEL_DATA[[PHENO]] <- DATASET$Transform_Functions[[TRANSFORM]](MODEL_DATA[[PHENO]])
  }
  
  
  ### Standardise age
  AGE_SD <- sd(MODEL_DATA[[VOI]], na.rm = TRUE)
  AGE_MEAN <- mean(MODEL_DATA[[VOI]], na.rm = TRUE)
  MODEL_DATA[[VOI]] <- (MODEL_DATA[[VOI]] - AGE_MEAN)/AGE_SD
  
  
  ### Standardise phenotype
  PHENO_SD <- DATASET$Info_Table$SD[DATASET$Info_Table$Phenotype == PHENO]
  PHENO_MEAN <- DATASET$Info_Table$Mean[DATASET$Info_Table$Phenotype == PHENO]
  MODEL_DATA[[PHENO]] <- (MODEL_DATA[[PHENO]] - PHENO_MEAN)/PHENO_SD
  
  
  # Fit models
  cat('\nFull model formula\n')
  print(FULL_FORMULA)
  FULL_MODEL <- lmer(
    FULL_FORMULA,
    lmerControl(
      optimizer = 'bobyqa'
    ),
    REML = TRUE,
    data = MODEL_DATA
  )
  cat('\nFull model random effect standard deviations\n')
  print(VarCorr(FULL_MODEL))
  
  cat('\nNull model formula\n')
  print(NULL_FORMULA)
  NULL_MODEL <- lmer(
    NULL_FORMULA,
    lmerControl(
      optimizer = 'bobyqa'
    ),
    REML = TRUE,
    data = MODEL_DATA
  )
  cat('\nNull model random effect standard deviations\n')
  print(VarCorr(NULL_MODEL))
  
  
  ### Pull estimated effect and standard error
  EFFECT <- as.numeric(fixef(FULL_MODEL)[VOI])
  SE <- sqrt(as.numeric(vcov(FULL_MODEL)[VOI, VOI]))
  
  
  ### Un-standardize
  EFFECT <- (EFFECT / AGE_SD) * PHENO_SD
  SE <- (SE/ AGE_SD) * PHENO_SD
  
  
  ### Test significance of the effect
  TEST_RESULTS <- KRmodcomp(
    largeModel = FULL_MODEL,
    smallModel = NULL_MODEL
  )
  cat('\nTest results\n')
  print(TEST_RESULTS$test)
  
  FTEST_STAT <- TEST_RESULTS$test['Ftest', 'stat']
  FTEST_NUMDF <- TEST_RESULTS$test['Ftest', 'ndf']
  FTEST_DENDF <- TEST_RESULTS$test['Ftest', 'ddf']
  FTEST_SCALE <- TEST_RESULTS$test['Ftest', 'F.scaling']
  FTEST_PV <- TEST_RESULTS$test['Ftest', 'p.value']
  
  
  ### Add results to dataframe
  DF_ROW <- tibble(
    Phenotype = PHENO,
    Age_Slope = EFFECT,
    Age_SE = SE,
    Age_FTestStat = FTEST_STAT,
    Age_FTestNumDF = FTEST_NUMDF,
    Age_FTestDenDF = FTEST_DENDF,
    Age_FTestScaling = FTEST_SCALE,
    Age_FTestPValue = FTEST_PV
  )
  
  AGE_EFFECTS <- rbind(
    AGE_EFFECTS,
    DF_ROW
  )
  
  ### Clean up
  rm(
    MODEL_DATA,
    VOI, FEF_TERMS, REF_TERMS, NULL_FORMULA, FULL_FORMULA,
    TRANSFORM,
    FULL_MODEL, NULL_MODEL, TEST_RESULTS,
    PHENO_SD, PHENO_MEAN,
    AGE_SD, AGE_MEAN, 
    EFFECT, SE, 
    FTEST_STAT, FTEST_NUMDF, FTEST_DENDF, FTEST_SCALE, FTEST_PV,
    DF_ROW
  )
  
  cat('\n# ----------------------------- #\n')
}
rm(PHENO)

AGE_EFFECTS <- AGE_EFFECTS %>% 
  mutate(
    Phenotype = factor(Phenotype, levels = levels(DATASET$Info_Table$Phenotype))
  )

DATASET$Age_Effects <- AGE_EFFECTS

#####


################################################################################
##### Estimate diet:age interaction effect #####

DIETAGE_EFFECTS <- tibble(
  Phenotype = vector(),
  AL_Slope = vector(),
  AL_SE = vector(),
  `1D_Slope` = vector(),
  `1D_SE` = vector(),
  `2D_Slope` = vector(),
  `2D_SE` = vector(),
  `20_Slope` = vector(),
  `20_SE` = vector(),
  `40_Slope` = vector(),
  `40_SE` = vector(),
  DietAge_FTestStat = vector(),
  DietAge_FTestNumDF = vector(),
  DietAge_FTestDenDF = vector(),
  DietAge_FTestScaling = vector(),
  DietAge_FTestPValue = vector(),
  DietAge_PairwiseDifferences = vector()
)

for(PHENO in DATASET$Info_Table$Phenotype){
  cat('\n\n'); cat(PHENO); cat('\n')
  
  ### Set-up data set
  VOI <- 'Diet:Age'
  REF_TERMS <- as.character(str_split(
    DATASET$Info_Table$REF_Cov[DATASET$Info_Table$Phenotype == PHENO],
    ':',
    simplify = TRUE
  ))
  FEF_TERMS <- as.character(na.omit(str_split(
    DATASET$Info_Table$FEF_Cov[DATASET$Info_Table$Phenotype == PHENO],
    ':',
    simplify = TRUE
  )))
  FEF_TERMS <- c('Diet', FEF_TERMS)
  
  MODEL_DATA <- DATASET$Data %>%
    filter(
      Timepoint %in% DATASET$Timepoints
    ) %>%
    select(
      MouseID, Timepoint,
      all_of(na.omit(c(
        PHENO, 
        unlist(str_split(VOI, pattern = ':')), 
        REF_TERMS, 
        FEF_TERMS
      )))
    )
  
  ### Create formulas
  REF_TERMS <- paste0('(1|', REF_TERMS, ')', collapse = ' + ')
  FEF_TERMS <- paste0(FEF_TERMS, collapse = ' + ')
  
  NULL_FORMULA <- as.formula(paste0(
    PHENO,
    ' ~ ',
    paste0(
      c(REF_TERMS, FEF_TERMS, 'Age'),
      collapse = ' + '
    )
  ))
  
  FULL_FORMULA <- as.formula(paste0(
    PHENO,
    ' ~ ',
    paste0(
      c(REF_TERMS, FEF_TERMS, VOI),
      collapse = ' + '
    )
  ))
  
  
  ### Transform phenotype if a transform is designated
  TRANSFORM <- DATASET$Info_Table$Transform[DATASET$Info_Table$Phenotype == PHENO]
  if(TRANSFORM != 'none'){
    MODEL_DATA[[PHENO]] <- DATASET$Transform_Functions[[TRANSFORM]](MODEL_DATA[[PHENO]])
  }
  
  
  ### Standardise age
  AGE_SD <- sd(MODEL_DATA[['Age']], na.rm = TRUE)
  AGE_MEAN <- mean(MODEL_DATA[['Age']], na.rm = TRUE)
  MODEL_DATA[['Age']] <- (MODEL_DATA[['Age']] - AGE_MEAN)/AGE_SD
  
  
  ### Standardise phenotype
  PHENO_SD <- DATASET$Info_Table$SD[DATASET$Info_Table$Phenotype == PHENO]
  PHENO_MEAN <- DATASET$Info_Table$Mean[DATASET$Info_Table$Phenotype == PHENO]
  MODEL_DATA[[PHENO]] <- (MODEL_DATA[[PHENO]] - PHENO_MEAN)/PHENO_SD
  
  
  ### Fit models
  cat('\nFull model formula\n')
  print(FULL_FORMULA)
  FULL_MODEL <- lmer(
    FULL_FORMULA,
    lmerControl(
      optimizer = 'bobyqa'
    ),
    REML = TRUE,
    data = MODEL_DATA
  )
  cat('\nFull model random effect standard deviations\n')
  print(VarCorr(FULL_MODEL))
  
  cat('\nNull model formula\n')
  print(NULL_FORMULA)
  NULL_MODEL <- lmer(
    NULL_FORMULA,
    lmerControl(
      optimizer = 'bobyqa'
    ),
    REML = TRUE,
    data = MODEL_DATA
  )
  cat('\nNull model random effect standard deviations\n')
  print(VarCorr(NULL_MODEL))
  
  
  ### Pull estimated effects and standard errors from model
  EFFECT <- c(
    as.numeric(fixef(FULL_MODEL)['DietAL:Age']),
    as.numeric(fixef(FULL_MODEL)['Diet1D:Age']),
    as.numeric(fixef(FULL_MODEL)['Diet2D:Age']),
    as.numeric(fixef(FULL_MODEL)['Diet20:Age']),
    as.numeric(fixef(FULL_MODEL)['Diet40:Age'])
  ) 
  SE <- c(
    sqrt(as.numeric(vcov(FULL_MODEL)['DietAL:Age', 'DietAL:Age'])),
    sqrt(as.numeric(vcov(FULL_MODEL)['Diet1D:Age', 'Diet1D:Age'])),
    sqrt(as.numeric(vcov(FULL_MODEL)['Diet2D:Age', 'Diet2D:Age'])),
    sqrt(as.numeric(vcov(FULL_MODEL)['Diet20:Age', 'Diet20:Age'])),
    sqrt(as.numeric(vcov(FULL_MODEL)['Diet40:Age', 'Diet40:Age']))
  )
  names(EFFECT) <- names(SE) <- c('AL', '1D', '2D', '20', '40')
  
  
  ### Un-standardize
  EFFECT <- (EFFECT / AGE_SD) * PHENO_SD
  SE <- (SE / AGE_SD) * PHENO_SD
  
  
  ### Test significance of effect
  TEST_RESULTS <- KRmodcomp(
    largeModel = FULL_MODEL,
    smallModel = NULL_MODEL
  )
  cat('\nTest results\n')
  print(TEST_RESULTS$test)
  
  FTEST_STAT <- TEST_RESULTS$test['Ftest', 'stat']
  FTEST_NUMDF <- TEST_RESULTS$test['Ftest', 'ndf']
  FTEST_DENDF <- TEST_RESULTS$test['Ftest', 'ddf']
  FTEST_SCALE <- TEST_RESULTS$test['Ftest', 'F.scaling']
  FTEST_PV <- TEST_RESULTS$test['Ftest', 'p.value']
  
  
  ### If the effect is significant, perform tukey test for different groups
  PAIRWISE_DIFF <- NA
  if(FTEST_PV < 0.05){
    PCOMP <- emtrends(FULL_MODEL, pairwise ~ Diet, var = "Age", adjust = "tukey")
    PCOMP <- PCOMP$contrasts %>% 
      data.frame() %>% 
      filter(p.value < 0.05)
    PAIRWISE_DIFF <- paste0(gsub(' ', '', PCOMP$contrast), collapse = ':')
    rm(PCOMP)
  }
  
  
  ### Add results to data frame
  DF_ROW <- tibble(
    Phenotype = PHENO,
    AL_Slope = EFFECT['AL'],
    AL_SE = SE['AL'],
    `1D_Slope` = EFFECT['1D'],
    `1D_SE` = SE['1D'],
    `2D_Slope` = EFFECT['2D'],
    `2D_SE` = SE['2D'],
    `20_Slope` = EFFECT['20'],
    `20_SE` = SE['20'],
    `40_Slope` = EFFECT['40'],
    `40_SE` = SE['40'],
    DietAge_FTestStat = FTEST_STAT,
    DietAge_FTestNumDF = FTEST_NUMDF,
    DietAge_FTestDenDF = FTEST_DENDF,
    DietAge_FTestScaling = FTEST_SCALE,
    DietAge_FTestPValue = FTEST_PV,
    DietAge_PairwiseDifferences = PAIRWISE_DIFF
  )
  
  DIETAGE_EFFECTS <- rbind(
    DIETAGE_EFFECTS,
    DF_ROW
  )
  
  
  ### Clean up
  rm(
    MODEL_DATA,
    VOI, FEF_TERMS, REF_TERMS, NULL_FORMULA, FULL_FORMULA,
    TRANSFORM,
    FULL_MODEL, NULL_MODEL, TEST_RESULTS,
    PHENO_SD, PHENO_MEAN,
    AGE_SD, AGE_MEAN,
    EFFECT, SE, 
    FTEST_STAT, FTEST_NUMDF, FTEST_DENDF, FTEST_SCALE, FTEST_PV,
    PAIRWISE_DIFF,
    DF_ROW
  )
  
  cat('\n# ----------------------------- #\n')
}
rm(PHENO)

DIETAGE_EFFECTS <- DIETAGE_EFFECTS %>% 
  mutate(
    Phenotype = factor(Phenotype, levels = levels(DATASET$Info_Table$Phenotype))
  )

DATASET$DietAge_Effects <- DIETAGE_EFFECTS

#####


################################################################################
##### Output #####

saveRDS(
  DATASET,
  paste0(
    'data/', DATASET$Assay,
    '_diet_age_dataset.rds'
  )
)

write_csv(
  DIET_EFFECTS,
  paste0(
    'results/', DATASET$Assay,
    '_diet_effects_stratified_by_timepoint.csv'
  )
)

write_csv(
  AGE_EFFECTS,
  paste0(
    'results/', DATASET$Assay,
    '_age_effects_conditioned_on_diet.csv'
  )
)

write_csv(
  DIETAGE_EFFECTS,
  paste0(
    'results/', DATASET$Assay,
    '_diet_age_interaction_effects.csv'
  )
)

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####