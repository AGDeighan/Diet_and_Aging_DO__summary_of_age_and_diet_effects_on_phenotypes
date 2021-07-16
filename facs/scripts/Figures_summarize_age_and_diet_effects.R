# 2021-07-06
################################################################################
#
#   This script creates a figure that summarizes the diet and age effects for
#   all the phenotypes
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
##### Diet effects by TP #####

STANDARDISED_DATA <- DATASET$Data %>% 
  filter(
    Timepoint %in% DATASET$Timepoints
  )

for(PHENO in as.character(DATASET$Info_Table$Phenotype)){
  TRANSFORM <- DATASET$Info_Table$Transform[DATASET$Info_Table$Phenotype == PHENO]
  if(TRANSFORM != 'none'){
    STANDARDISED_DATA[[PHENO]] <- DATASET$Transform_Functions[[TRANSFORM]](STANDARDISED_DATA[[PHENO]])
  }
  
  MEAN <- DATASET$Info_Table$Mean[DATASET$Info_Table$Phenotype == PHENO]
  SD <- DATASET$Info_Table$SD[DATASET$Info_Table$Phenotype == PHENO]
  
  STANDARDISED_DATA[[PHENO]] <- (STANDARDISED_DATA[[PHENO]] - MEAN)/SD
  
  rm(TRANSFORM, MEAN, SD)
}
rm(PHENO)

PLOT_DATA <- STANDARDISED_DATA %>% 
  pivot_longer(
    cols = as.character(DATASET$Info_Table$Phenotype),
    values_to = 'Value',
    names_to = 'Phenotype'
  ) %>% 
  mutate(
    Phenotype = factor(Phenotype, levels = levels(DATASET$Info_Table$Phenotype))
  ) %>% 
  arrange(
    Phenotype
  )

LIMITS <- STANDARDISED_DATA %>% 
  pivot_longer(
    cols = as.character(DATASET$Info_Table$Phenotype),
    values_to = 'Value',
    names_to = 'Phenotype'
  ) %>% 
  group_by(
    Phenotype, Timepoint, Diet
  ) %>% 
  summarise(
    IQR = IQR(Value, na.rm = TRUE),
    LB = quantile(Value, 0.25, na.rm = TRUE) - 1.5 * IQR,
    UB = quantile(Value, 0.75, na.rm = TRUE) + 1.5 * IQR
  ) %>% 
  ungroup() %>% 
  summarise(
    LB = min(LB),
    UB = max(UB)
  ) %>% 
  unlist() %>% 
  as.numeric() %>% 
  abs() %>% 
  max()
LIMITS <- ceiling((LIMITS + 0.25)/0.1)*0.1
LIMITS <- c(-1*LIMITS, LIMITS)


DIET_BY_TP_PLOT <- PLOT_DATA %>% 
  ggplot() +
  theme_minimal(
    base_size = 9
  ) +
  geom_boxplot(
    aes(x = Phenotype, y = Value, color = Diet),
    outlier.shape = NA,
    lwd = 0.25
  ) +
  geom_point(
    data = DATASET$Diet_Effects %>% 
      select(
        Phenotype, Timepoint, Diet_FTestPValue
      ) %>% 
      group_by(
        Timepoint
      ) %>% 
      mutate(
        FWER = p.adjust(Diet_FTestPValue, 'holm')
      ) %>% 
      ungroup() %>% 
      filter(
        FWER < 0.05
      ) %>% 
      mutate(
        Position = LIMITS[2] - 0.25
      ),
    aes(x = Phenotype, y = Position),
    shape = '*',
    size = 5
  ) +
  scale_y_continuous(
    breaks = seq(-100, 100, 2),
    minor_breaks = seq(-100, 100, 1),
    limits = LIMITS
  ) +
  scale_color_manual(
    values = DIET_COLORS
  ) +
  theme(
    axis.text.x = element_text(
      angle = 45, vjust = 1, hjust=1
    ),
    panel.grid.major.x = element_blank(),
    plot.caption = element_text(hjust = 0),
    legend.position = 'none'
  ) +
  facet_grid(
    Timepoint ~ .
  ) +
  labs(
    title = 'Diet effect on trait values by timepoint',
    y = 'Standardized trait value',
    x = NULL,
    color = NULL,
    caption = paste0(
      'The boxplots show the distribution of the normalized (mean = 0, standard deviation = 1) phenotype at each timepoint grouped by diet',
      '\n',
      'The asterisks indicate phenotypes for which the diet effect is significant at a FWER of 0.05',
      '\n'
    )
  )

rm(
  STANDARDISED_DATA, PLOT_DATA, LIMITS
)

DIET_TP_WIDTH <- 6 + (length(DATASET$Info_Table$Phenotype) - 1)*0.2
DIET_TP_HEIGHT <- 3 + (length(DATASET$Timepoints) - 1)*2


#####


################################################################################
##### Age effects  #####

PLOT_DATA <- rbind(
  DATASET$Age_Effects %>% 
    mutate(
      Coef = Age_Slope/DATASET$Info_Table$SD[DATASET$Info_Table$Phenotype == Phenotype],
      SE = Age_SE/DATASET$Info_Table$SD[DATASET$Info_Table$Phenotype == Phenotype],
      FWER = p.adjust(Age_FTestPValue, 'holm'),
      Phenotype = fct_rev(Phenotype),
      Position = as.numeric(Phenotype),
      Diet = NA,
      Group = 'Overall'
    ) %>% 
    select(
      Phenotype, Group, Diet, Position, Coef, SE, FWER
    ),
  DATASET$DietAge_Effects %>% 
    mutate(
      FWER = p.adjust(DietAge_FTestPValue, 'holm'),
      Phenotype = fct_rev(Phenotype)
    ) %>% 
    select(
      Phenotype,
      AL_Slope:`40_SE`,
      FWER
    ) %>% 
    pivot_longer(
      cols = AL_Slope:`40_SE`,
      values_to = 'V',
      names_to = 'K'
    ) %>% 
    separate(
      K,
      into = c('Diet', 'K'),
      sep = '_'
    ) %>% 
    pivot_wider(
      values_from = 'V',
      names_from = 'K'
    ) %>% 
    rowwise() %>% 
    mutate(
      Coef = Slope/DATASET$Info_Table$SD[DATASET$Info_Table$Phenotype == Phenotype],
      SE = SE/DATASET$Info_Table$SD[DATASET$Info_Table$Phenotype == Phenotype],
      Diet = factor(Diet, levels = rev(levels(DATASET$Data$Diet))),
      Position = as.numeric(Phenotype) + ((as.numeric(Diet) - 1)/4 - 0.5)*1/3,
      Group = 'Diet specific'
    ) %>% 
    ungroup() %>% 
    select(
      Phenotype, Group, Diet, Position, Coef, SE, FWER
    )
) %>% 
  mutate(
    Diet = factor(Diet, levels = levels(DATASET$Data$Diet)),
    Group = factor(Group, levels = c('Overall', 'Diet specific'))
  )
  
  


LIMITS <- PLOT_DATA %>% 
  summarise(
    M = max(abs(c(min(Coef - SE), max(Coef + SE))))
  ) %>% 
  unlist() %>% 
  as.numeric()
LIMITS <- ceiling((LIMITS + 0.01)/0.05)*0.05
LIMITS <- c(-1*LIMITS, LIMITS)

AGE_EFFECTS_PLOT <- PLOT_DATA %>% 
  ggplot() +
  theme_minimal(
    base_size = 9
  ) +
  geom_vline(
    xintercept = 0,
    color = 'gray60',
    linetype = 3
  ) +
  geom_point(
    data = PLOT_DATA %>% 
      group_by(
        Phenotype, Group, FWER
      ) %>% 
      summarise() %>% 
      ungroup() %>% 
      filter(
        FWER < 0.05
      ) %>% 
      mutate(
        Position = as.numeric(Phenotype),
        X_Position = ifelse(
          Group == 'Overall', LIMITS[1],
          LIMITS[2]
        )
      ),
    aes(x = X_Position, y = Position),
    shape = '*',
    size = 5
  ) +
  geom_errorbarh(
    aes(xmin = Coef - SE, xmax = Coef + SE, y = Position, color = Diet),
    height = 1/10
  ) +
  geom_point(
    aes(x = Coef, y = Position, color = Diet)
  ) +
  scale_x_continuous(
    breaks = seq(-10, 10, 0.05),
    minor_breaks = seq(-10, 10, 0.01),
    limits = LIMITS
  ) +
  scale_y_continuous(
    breaks = unique(as.numeric(PLOT_DATA$Phenotype)),
    minor_breaks = seq(0, 500, 0.5),
    labels = rev(levels(PLOT_DATA$Phenotype)),
    limits = c(0.49, max(as.numeric(PLOT_DATA$Phenotype)) + 0.51),
    expand = c(0, 0)
  ) +
  scale_color_manual(
    values = DIET_COLORS,
    na.value = 'black'
  ) +
  facet_grid(
    . ~ Group
  ) +
  theme(
    panel.grid.major.y = element_blank(),
    plot.caption = element_text(hjust = 0),
    legend.position = 'none'
  ) +
  labs(
    title = 'Age effect on trait values',
    x = 'Age coefficient',
    y = NULL,
    caption = paste0(
      'The age coefficient for each phenotype indicates the estimated change with age in units of standard deviations per month',
      '\n',
      'The bars around the coefficient show +/- 1 standard error of the coefficient',
      '\n',
      'The asterisks indicate phenotypes for which the direct age effect is significant (left) or the age-diet interaction is signficant (right) at a FWER of 0.05'
    )
  )

rm(
  PLOT_DATA, LIMITS
)

#####


################################################################################
##### Final plot  #####

LEGEND <- get_legend(
  DATASET$Data %>% 
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


DIET_BY_TP_PLOT <- plot_grid(
  LEGEND,
  DIET_BY_TP_PLOT,
  ncol = 1,
  rel_heights = c(0.1, 1)
)

PLOT <- plot_grid(
  DIET_BY_TP_PLOT,
  AGE_EFFECTS_PLOT,
  ncol = 1,
  rel_heights = c(DIET_TP_HEIGHT, DIET_TP_WIDTH)
)




pdf(
  paste0(
    'figures/', DATASET$Assay,
    '_diet_and_age_effects_summarized.pdf'
  ),
  height = DIET_TP_WIDTH + DIET_TP_HEIGHT,
  width = DIET_TP_WIDTH
)
plot(PLOT)
dev.off()

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