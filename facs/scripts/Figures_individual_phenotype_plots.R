# 2021-07-08
################################################################################
#
#   This creates a separate figure for each phenotype that gives more details 
#   on the effects of the age and diet on that phenotype
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
library(cowplot)
DIET_COLORS <- c("seashell4", "skyblue", "royalblue4", "orange", "firebrick")

################################################################################
## helper functions #########################################################

#####


################################################################################
## load dataset #########################################################

DATASET <- readRDS('data/FACS_diet_age_dataset.rds')

#####


################################################################################
##### Transform traits that need to be transformed #####

TRANSFORMED_DATA <- DATASET$Data

for(PHENO in as.character(DATASET$Info_Table$Phenotype)){
  TRANSFORM <- DATASET$Info_Table$Transform[DATASET$Info_Table$Phenotype == PHENO]
  if(TRANSFORM != 'none'){
    TRANSFORMED_DATA[[PHENO]] <- DATASET$Transform_Functions[[TRANSFORM]](TRANSFORMED_DATA[[PHENO]])
  }
  
  rm(TRANSFORM)
}


#####


################################################################################
##### Plot traits #####


for(PHENO in as.character(DATASET$Info_Table$Phenotype)){
  TITLE <- paste0(
    'Diet and age effects on ',
    DATASET$Info_Table$Description[DATASET$Info_Table$Phenotype == PHENO],
    '\n'
  )
  
  TRANSFORMED_DATA[['Pheno']] <- TRANSFORMED_DATA[[PHENO]]
  DATA <- TRANSFORMED_DATA %>% 
    select(
      MouseID, Timepoint, Age, Diet,
      Pheno
    )
  TRANSFORMED_DATA <- TRANSFORMED_DATA %>% 
    select(-Pheno)
  
  PHENO_LABEL <- ifelse(
    DATASET$Info_Table$Transform[DATASET$Info_Table$Phenotype == PHENO] == 'none', PHENO,
    paste0(
      DATASET$Info_Table$Transform[DATASET$Info_Table$Phenotype == PHENO],
      '[', PHENO, ']'
    )
  )
  
  TP_BOXPLOT <- DATA %>%
    ggplot() +
    theme_minimal(
      base_size = 9
    ) +
    geom_point(
      aes(x = Timepoint, y = Pheno, color = Diet),
      position = position_jitterdodge(jitter.width = 0.4),
      alpha = 1/3
    ) +
    geom_boxplot(
      aes(x = Timepoint, y = Pheno, color = Diet),
      outlier.shape = NA,
      alpha = 0
    ) +
    scale_color_manual(
      values = DIET_COLORS
    ) +
    theme(
      legend.position = 'none'
    ) +
    labs(
      y = PHENO_LABEL,
      x = 'Timepoint'
    )
  
  MEAN_SE_PLOT <- DATA %>% 
    group_by(
      Diet, Timepoint
    ) %>% 
    summarise(
      Mean = mean(Pheno, na.rm = TRUE),
      SE = sd(Pheno, na.rm = TRUE)/sqrt(sum(!is.na(Pheno)))
    ) %>% 
    ungroup() %>% 
    mutate(
      Position = as.numeric(substr(Timepoint, 1, 2)) +
        ((as.numeric(Diet) - 1)/(length(unique(Diet)) - 1) - 0.5)*1/2
    ) %>% 
    ggplot() +
    theme_minimal(
      base_size = 9
    ) +
    geom_errorbar(
      aes(x = Position, ymin = Mean - SE, ymax = Mean + SE, color = Diet),
      width = 0.1
    ) +
    geom_line(
      aes(x = Position, y = Mean, color = Diet, group = Diet)
    ) +
    geom_point(
      aes(x = Position, y = Mean, color = Diet),
      shape = 20
    ) +
    scale_color_manual(
      values = DIET_COLORS,
      guide = FALSE
    ) +
    scale_x_continuous(
      breaks = sort(as.numeric(gsub(' months', '', unique(DATA$Timepoint))))
    ) +
    labs(
      x = 'Age (months)',
      y = PHENO_LABEL
    )
  
  
  LEGEND <- get_legend(
    DATA %>% 
      ggplot() +
      theme_minimal() +
      geom_point(
        aes(x = Timepoint, y = Age, color = Diet),
        shape = 15, size = 4
      ) +
      scale_color_manual(
        values = DIET_COLORS,
        labels = c(
          'Ad libitum  ',
          '1-day fast  ',
          '2-day fast  ',
          '20% caloric restriction  ',
          '40% caloric restriction  '
        )
      ) +
      theme(
        legend.position = 'top'
      ) +
      labs(
        color = NULL
      )
  )
  
  
  TITLE <- ggdraw() + 
    draw_label(
      TITLE,
      size = 11,
      x = 0,
      y = 0.8,
      hjust = 0,
      vjust = 1,
      angle = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 3)
    )
  
  
  DIET_PV <- DATASET$Diet_Effects %>% 
    filter(
      Phenotype == PHENO
    ) %>% 
    select(
      Timepoint, Diet_FTestPValue, Diet_PairwiseDifferences
    )
  
  AGE_PV <- DATASET$Age_Effects %>% 
    filter(
      Phenotype == PHENO
    ) %>% 
    select(
      Age_FTestPValue
    )
  
  DIETAGE_PV <- DATASET$DietAge_Effects %>% 
    filter(
      Phenotype == PHENO
    ) %>% 
    select(
      DietAge_FTestPValue, DietAge_PairwiseDifferences
    )
  
  CAPTION_TEXT <- paste0(
    'Only the following timepoints were used when testing for direct diet and age-diet interaction effects (all timepoints were used when testing for direct age effects): ',
    paste0(
      paste0(
        DATASET$Timepoints[1:(length(DATASET$Timepoints) - 1)], 
        collapse = ', '
      ),
      ' and ',
      DATASET$Timepoints[length(DATASET$Timepoints)],
      '. '
    ),
    'The effects of age, diet, and the age-diet interaction were estimated using mixed linear models and the significance of the effects were assessed with an approximate F-test using the Kenward and Roger (1997) approach. ',
    'The p-values for the diet effect at each timepoint are: ',
    paste0(
      paste0(
        DATASET$Timepoints[1:(length(DATASET$Timepoints) - 1)],
        ' = ',
        signif(DIET_PV$Diet_FTestPValue[1:(nrow(DIET_PV) - 1)], digits = 3),
        collapse = '; '
      ),
      ' and ',
      paste0(
        DATASET$Timepoints[length(DATASET$Timepoints)],
        ' = ',
        signif(DIET_PV$Diet_FTestPValue[nrow(DIET_PV)], digits = 3),
        collapse = '; '
      ),
      '.'
    )
  )
  
  for(i in 1:nrow(DIET_PV)){
    if(!is.na(DIET_PV$Diet_PairwiseDifferences[i])){
      DIFFERENCES <- DIET_PV$Diet_PairwiseDifferences[i]
      DIFFERENCES <- unlist(str_split(DIFFERENCES, pattern = ':'))
      if(length(DIFFERENCES) > 1){
        DIFFERENCES <- paste0(
          paste0(
            DIFFERENCES[1:(length(DIFFERENCES) - 1)],
            collapse = ', '
          ),
          ' and ',
          DIFFERENCES[length(DIFFERENCES)]
        )
      }
      TEXT <- paste0(
        'The diet pairs that have significantly different (Tukey p-value < 0.05) means at ',
        DIET_PV$Timepoint[i],
        ' are ',
        DIFFERENCES,
        '.'
      )
      CAPTION_TEXT <- paste0(CAPTION_TEXT, ' ', TEXT)
      rm(DIFFERENCES, TEXT)
    }
  }
  rm(i)
  
  CAPTION_TEXT <- paste0(
    CAPTION_TEXT,
    ' ',
    'The p-value for the direct effect of age on ',
    PHENO,
    ' is ',
    signif(AGE_PV$Age_FTestPValue, digits = 3), 
    '. ',
    'The p-value for the effect of the interaction between age and diet on ',
    PHENO,
    ' is ',
    signif(DIETAGE_PV$DietAge_FTestPValue, digits = 3), 
    '.'
  )
  
  if(!is.na(DIETAGE_PV$DietAge_PairwiseDifferences)){
    DIFFERENCES <- DIETAGE_PV$DietAge_PairwiseDifferences
    DIFFERENCES <- unlist(str_split(DIFFERENCES, pattern = ':'))
    if(length(DIFFERENCES) > 1){
      DIFFERENCES <- paste0(
        paste0(
          DIFFERENCES[1:(length(DIFFERENCES) - 1)],
          collapse = ', '
        ),
        ' and ',
        DIFFERENCES[length(DIFFERENCES)]
      )
    }
    TEXT <- paste0(
      'The diet pairs that have significantly different (Tukey p-value < 0.05) rates of change with age are ',
      DIFFERENCES,
      '.'
    )
    CAPTION_TEXT <- paste0(CAPTION_TEXT, ' ', TEXT)
    rm(DIFFERENCES, TEXT)
  }
  
  
  i <- 0
  for(INDEX in seq(130, nchar(CAPTION_TEXT), by = 130)){
    POSITION <- INDEX + i
    if(substr(CAPTION_TEXT, POSITION, POSITION) != ' '){
      for(j in 1:20){
        if(substr(CAPTION_TEXT, POSITION - j, POSITION - j) == ' '){
          POSITION <- POSITION - j
          break
        }
        if(substr(CAPTION_TEXT, POSITION + j, POSITION + j) == ' '){
          POSITION <- POSITION + j
          break
        }
      }
      rm(j)
    }
    CAPTION_TEXT <- paste0(
      substr(CAPTION_TEXT, 1, POSITION),
      '\n',
      substr(CAPTION_TEXT, POSITION + 1, 999999)
    )
    i <- i + 2
    rm(POSITION)
  }
  rm(i, INDEX)
  
  CAPTION <- ggdraw() + 
    draw_label(
      CAPTION_TEXT,
      size = 9,
      x = 0,
      y = 0.8,
      hjust = 0,
      vjust = 1,
      angle = 0
    ) +
    theme(
      # add margin on the left of the drawing canvas,
      # so title is aligned with left edge of first plot
      plot.margin = margin(0, 0, 0, 3)
    )
  
  PLOT <- plot_grid(
    TITLE,
    LEGEND,
    TP_BOXPLOT,
    MEAN_SE_PLOT,
    CAPTION,
    ncol = 1,
    rel_heights = c(0.15, 0.05, 1, 1, 2/3)
  )
  
  
  WIDTH <- 6 + (length(unique(DATA$Timepoint)) - 1)
  HEIGHT <- 7
  
  pdf(
    paste0(
      'figures/individual_phenotype_plots/',
      PHENO,
      '.pdf'
    ),
    width = WIDTH, height = HEIGHT
  )
  plot(PLOT)
  dev.off()
  
  rm(
    DATA, 
    PHENO_LABEL, 
    TP_BOXPLOT, MEAN_SE_PLOT, 
    LEGEND, TITLE,
    DIET_PV, AGE_PV, DIETAGE_PV, 
    CAPTION_TEXT, CAPTION, 
    PLOT, WIDTH, HEIGHT
  )
}
rm(PHENO)


#####


################################################################################
#####  #####

#####


################################################################################
## clear workspace ###########################################

rm(list = ls())
pacman::p_unload('all')
graphics.off()

#####