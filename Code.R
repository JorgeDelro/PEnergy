# BAyesian Hierarchical analysis for PE manuscript

library(tidyverse)
library(rstanarm)
library(bayesplot)
library(tidybayes)
library(modelr)
library(ggpubr)

# Read database
db <- read.csv2("PEdb.csv")

# Fixing variables for analyses
db$Height <- db$Height * 100
db$Sex <- factor(db$Sex, levels = c(1,2), labels = c("Boys", "Girls"))
db$BMIGroup <- factor(db$BMIGroup, levels = c("Underweight", "Normalweight", "Overweight", "Obese"),
                      labels = c("Underweight", "Normalweight", "Overweight", "Obese"))

# Z-score of outcomes and predictors
db$Pow_Duncan_stand <- as.vector(scale(db$Pow_Duncan))
db$Pow_Gomez_stand <- as.vector(scale(db$Pow_Gomez))
db$CMJcm_stand <- as.vector(scale(db$CMJcm))
db$PE_CMJ_stand <- as.vector(scale(db$PE_CMJ))


### Table S1 - Descriptives
# Sex
db$Sex %>% table()
db$Sex %>% table() %>% prop.table()
db %>% filter(BMIGroup == "Normalweight") %>% select(Sex) %>% table()
db %>% filter(BMIGroup == "Normalweight") %>% select(Sex) %>% table() %>% prop.table()

# Age, Heaight, Weight, CMJ, PE_CMJ, CMJ_Duncan, CMJ_Gomez
db %>%
  #group_by(BMIGroup) %>%
  select(Age, Height, Weight, BMI, CMJcm, PE_CMJ, Pow_Duncan, Pow_Gomez) %>%
  summarise(M_Age = mean(Age, na.rm = T),
            SD_Age = sd(Age, na.rm = T),
            M_Height = mean(Height, na.rm = T),
            SD_Height = sd(Height, na.rm = T),
            M_Weight = mean(Weight, na.rm = T),
            SD_Weight = sd(Weight, na.rm = T),
            M_BMI = mean(BMI, na.rm = T),
            SD_BMI = sd(BMI, na.rm = T),
            M_CMJcm = mean(CMJcm, na.rm = T),
            SD_CMJcm = sd(CMJcm, na.rm = T),
            M_PE_CMJ = mean(PE_CMJ, na.rm = T),
            SD_PE_CMJ = sd(PE_CMJ, na.rm = T),
            M_Pow_Duncan = mean(Pow_Duncan, na.rm = T),
            SD_Pow_Duncan = sd(Pow_Duncan, na.rm = T),
            M_Pow_Gomez = mean(Pow_Gomez, na.rm = T),
            SD_Pow_Gomez = sd(Pow_Gomez, na.rm = T)) %>% 
  glimpse()

## Plots - Data Exploration
desCMJDuncan <- db %>%
  ggplot(aes(x = CMJcm, y = Pow_Duncan)) +
  geom_point() +
  labs(x = "CMJ (cm)", y = expression(CMJ[Duncan]*" (W)")) +
  theme(axis.text=element_text(size=12),
        axis.title.y = element_text(size = 14, face="bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size = 14, face="bold", margin = margin(t = 10, r = 20, b = 0, l = 0))) +
  theme_bw()

desPEDuncan <- db %>%
  ggplot(aes(x = CMJcm, y = Pow_Gomez)) +
  geom_point() + 
  labs(x = "CMJ (cm)", y = expression(CMJ[Gomez]*" (W)")) +
  theme(axis.text=element_text(size=12),
        axis.title.y = element_text(size = 14, face="bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size = 14, face="bold", margin = margin(t = 10, r = 20, b = 0, l = 0))) +
  theme_bw()

desCMJGomez <- db %>%
  ggplot(aes(x = PE_CMJ, y = Pow_Duncan)) +
  geom_point() +
  labs(x = expression(PE[CMJ]*" (J)"), y = expression(CMJ[Duncan]*" (W)")) +
  theme(axis.text=element_text(size=12),
        axis.title.y = element_text(size = 14, face="bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size = 14, face="bold", margin = margin(t = 10, r = 20, b = 0, l = 0))) +
  theme_bw()

desPEGomez <- db %>%
  ggplot(aes(x = PE_CMJ, y = Pow_Gomez)) +
  geom_point() +
  labs(x = expression(PE[CMJ]*" (J)"), y = expression(CMJ[Gomez]*" (W)")) +
  theme(axis.text=element_text(size=12),
        axis.title.y = element_text(size = 14, face="bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
        axis.title.x = element_text(size = 14, face="bold", margin = margin(t = 10, r = 20, b = 0, l = 0))) +
  theme_bw()

# Figure S1
ggarrange(desCMJDuncan,
          desPEDuncan,
          desCMJGomez,
          desPEGomez,
          ncol = 2,
          nrow = 2)

ggsave("FigureS1.png", height=7, width=7, units='in', dpi=600)

##
##### Fit Multilevel models
##
options(mc.cores = parallel::detectCores())

##### DUNCAN - CMJ
# Variying intercept - slope model for CJMcm
bmod1 <- stan_glmer(formula = Pow_Duncan_stand ~ CMJcm_stand + (CMJcm_stand|BMIGroup),
                     data = db,
                     family = gaussian,
                     adapt_delta = 0.99,
                     seed = 1234)

# Prior distribution
prior_summary(bmod1)
posterior_vs_prior(bmod1)

# Results
print(bmod1) # Median
summary(bmod1) # Mean

# R-Squared
hist(bayes_R2(bmod1))
median(bayes_R2(bmod1)) # 0.53
HDInterval::hdi(bayes_R2(bmod1)) # (0.49 - 0.56)

# Plot Fitted draws
#db %>%
#  group_by(BMIGroup) %>%
#  data_grid(PE_CMJ = seq_range(PE_CMJ, n = 101)) %>%
#  add_fitted_draws(bmod8, n = 100) %>%
#  ggplot(aes(x = PE_CMJ, y = Pow_Gomez)) +
#  geom_point(data = db, alpha = 0.7) +
#  geom_line(aes(y = .value, group = paste(BMIGroup, .draw)), alpha = 0.7, color = "blue") +
#  facet_grid(BMIGroup ~ ., space = "free_x", scales = "free_x") +
#  labs(x = expression(PE[CMJ]*" (J)"), y = expression(CMJ[Gomez]*" (W)")) +
#  theme_bw()

# Plot Predictive draws
DuncanCMJ_pred <- db %>%
  group_by(BMIGroup) %>%
  data_grid(CMJcm_stand = seq_range(CMJcm_stand, n = 101)) %>%
  add_predicted_draws(bmod1, n = 100) %>%
  ggplot(aes(x = CMJcm_stand, y = Pow_Duncan_stand)) +
  geom_point(data = db, alpha = 0.8) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 0.3, show.legend = F, color = "blue") +
  scale_y_continuous(limits = c(-4, 5)) +
  facet_grid(BMIGroup ~ ., space = "free_x", scales = "free_x") +
  labs(x = "CMJ (z-score)", y = expression(CMJ[Duncan]*" (z-score)")) +
  theme_bw()

# Save it
# ggsave("Duncan_CMJ_predictive.png", height=7, width=9, units='in', dpi=600)


##### DUNCAN _ PE
bmod2 <- stan_glmer(formula = Pow_Duncan_stand ~ PE_CMJ_stand + (PE_CMJ_stand|BMIGroup),
                    data = db,
                    family = gaussian,
                    adapt_delta = 0.99,
                    seed = 1234)

# Prior distribution
prior_summary(bmod2)
posterior_vs_prior(bmod2)

# Results
print(bmod2)
summary(bmod2)

# R-Squared
hist(bayes_R2(bmod2))
median(bayes_R2(bmod2)) # 0.95
HDInterval::hdi(bayes_R2(bmod2)) # (0.94 - 0.95)

# Plot Predictive draws
DuncanPE_pred <- db %>%
  group_by(BMIGroup) %>%
  data_grid(PE_CMJ_stand = seq_range(PE_CMJ_stand, n = 101)) %>%
  add_predicted_draws(bmod2, n = 100) %>%
  ggplot(aes(x = PE_CMJ_stand, y = Pow_Duncan_stand)) +
  geom_point(data = db, alpha = 0.8) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 0.3, show.legend = F, color = "blue") +
  scale_y_continuous(limits = c(-4, 5)) +
  facet_grid(BMIGroup ~ ., space = "free_x", scales = "free_x") +
  labs(x = expression(PE[CMJ]*" (z-score)"), y = expression(CMJ[Duncan]*" (z-score)")) +
  theme_bw()

# Save it
# ggsave("Duncan_PE_predictive.png", height=7, width=9, units='in', dpi=600)

##
## Compute loo to compare predicitive accuracy
##
loo1 <- loo(bmod1)
loo2 <- loo(bmod2, k_threshold = 0.7)
compare_models(loo1, loo2)



##### GOMEZ - CMJ
bmod3 <- stan_glmer(formula = Pow_Gomez_stand ~ CMJcm_stand + (CMJcm_stand|BMIGroup),
                    data = db,
                    family = gaussian,
                    adapt_delta = 0.99,
                    seed = 1234)

# Prior distribution
prior_summary(bmod3)

# Results
print(bmod3)
summary(bmod3)

# R-Squared
hist(bayes_R2(bmod3))
median(bayes_R2(bmod3)) # 0.67
HDInterval::hdi(bayes_R2(bmod3)) # (0.65 - 0.69)

# Plot Predictive draws
GomezCMJ_pred <- db %>%
  group_by(BMIGroup) %>%
  data_grid(CMJcm_stand = seq_range(CMJcm_stand, n = 101)) %>%
  add_predicted_draws(bmod3, n = 100) %>%
  ggplot(aes(x = CMJcm_stand, y = Pow_Gomez_stand)) +
  geom_point(data = db, alpha = 0.8) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 0.3, show.legend = F, color = "blue") +
  scale_y_continuous(limits = c(-4, 5)) +
  facet_grid(BMIGroup ~ ., space = "free_x", scales = "free_x") +
  labs(x = "CMJ (z-score)", y = expression(CMJ[Gomez]*" (z-score)")) +
  theme_bw()

# Save it
# ggsave("Gomez_CMJ_predictive.png", height=7, width=9, units='in', dpi=600)




##### GOMEZ - PE
bmod4 <- stan_glmer(formula = Pow_Gomez_stand ~ PE_CMJ_stand + (PE_CMJ_stand|BMIGroup),
                    data = db,
                    family = gaussian,
                    adapt_delta = 0.99,
                    seed = 1234)

# Prior distribution
prior_summary(bmod4)

# Results
print(bmod4)
summary(bmod4)

# R-Squared
hist(bayes_R2(bmod4))
median(bayes_R2(bmod4)) # 0.97
HDInterval::hdi(bayes_R2(bmod4)) # (0.97 - 0.97)

# Plot Predictive draws
GomezPE_pred <- db %>%
  group_by(BMIGroup) %>%
  data_grid(PE_CMJ_stand = seq_range(PE_CMJ_stand, n = 101)) %>%
  add_predicted_draws(bmod4, n = 100) %>%
  ggplot(aes(x = PE_CMJ_stand, y = Pow_Gomez_stand)) +
  geom_point(data = db, alpha = 0.8) +
  stat_lineribbon(aes(y = .prediction), .width = c(.95), alpha = 0.3, show.legend = F, color = "blue") +
  scale_y_continuous(limits = c(-4, 5)) +
  facet_grid(BMIGroup ~ ., space = "free_x", scales = "free_x") +
  labs(x = expression(PE[CMJ]*" (z-score)"), y = expression(CMJ[Gomez]*" (z-score)")) +
  theme_bw()

# Save it
# ggsave("Gomez_PE_predictive.png", height=7, width=9, units='in', dpi=600)


# Compute loo to compare predicitive accuracy
loo3 <- loo(bmod3)
loo4 <- loo(bmod4)
compare_models(loo3, loo4)

## Arrange predictive drwas figures
ggarrange(DuncanCMJ_pred,
          DuncanPE_pred,
          GomezCMJ_pred,
          GomezPE_pred,
          ncol = 2,
          nrow = 2)

ggsave("Figure1.png", height=9, width=9, units='in', dpi=600)


## Z-score plot
scale_var <- function(x){
  (x - mean(x, na.rm=TRUE)) / sd(x, na.rm=TRUE)
}

db %>%
  as.tibble() %>%
  mutate(Pow_Duncan_scale = scale_var(Pow_Duncan),
         Pow_Gomez_scale = scale_var(Pow_Gomez),
         PE_CMJ_scale = scale_var(PE_CMJ),
         CMJcm_scale = scale_var(CMJcm)) %>%
  select(BMIGroup, Pow_Duncan_scale, Pow_Gomez_scale, PE_CMJ_scale, CMJcm_scale) %>%
  group_by(BMIGroup) %>%
  summarize(Pow_Duncan_scale_mean = mean(Pow_Duncan_scale),
         Pow_Gomez_scale_mean = mean(Pow_Gomez_scale),
         PE_CMJ_scale_mean = mean(PE_CMJ_scale),
         CMJcm_scale_mean = mean(CMJcm_scale)) %>%
  gather(Pow_Duncan_scale_mean, Pow_Gomez_scale_mean, PE_CMJ_scale_mean, CMJcm_scale_mean, 
         key = "Variable", 
         value = "Value") %>%
  mutate(Variable = factor(Variable, levels = c("Pow_Duncan_scale_mean",
                                                "Pow_Gomez_scale_mean",
                                                "PE_CMJ_scale_mean",
                                                "CMJcm_scale_mean"))) %>%
  ggplot(aes(x = BMIGroup, y = Value, group = Variable, shape = Variable, linetype = Variable)) +
    geom_line() +
    geom_point(aes(size = 2)) +
    scale_shape_manual(values=c(0, 1, 2, 3),
                                     breaks=c("Pow_Duncan_scale_mean", 
                                              "Pow_Gomez_scale_mean",
                                              "PE_CMJ_scale_mean",
                                              "CMJcm_scale_mean"),
                                     labels=c(expression(CMJ[Duncan]), 
                                              expression(CMJ[Gomez]), 
                                              expression(PE[CMJ]),
                                              "CMJ")) +
    scale_linetype_manual(values=c(1, 2, 3, 4),
                          breaks=c("Pow_Duncan_scale_mean", 
                                   "Pow_Gomez_scale_mean",
                                   "PE_CMJ_scale_mean",
                                   "CMJcm_scale_mean"),
                          labels=c(expression(CMJ[Duncan]), 
                                   expression(CMJ[Gomez]), 
                                   expression(PE[CMJ]),
                                   "CMJ")) +
    labs(x = "BMI Group", y = "Z-Score") +
    theme(legend.position="top",
          legend.title = element_blank(),
          legend.spacing.x = unit(0.5, 'cm'),
          axis.text=element_text(size=12),
          axis.title.y = element_text(size = 14, face="bold", margin = margin(t = 0, r = 15, b = 0, l = 0)),
          axis.title.x = element_text(size = 14, face="bold", margin = margin(t = 10, r = 20, b = 0, l = 0))) +
  scale_size(guide = 'none')

ggsave("ZScore_plot.png", height=7, width=10, units='in', dpi=600)


#### Prior predictive checkings
priorpbmod1 <- posterior_vs_prior(bmod1, prob = 0.95) + 
  ggtitle(expression(CMJ[Duncan]*" ~ CMJ")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text.x  = element_text(size = 6))

priorpbmod2 <- posterior_vs_prior(bmod2, prob = 0.95) + 
  ggtitle(expression(CMJ[Duncan]*" ~ "*PE[CMJ])) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text.x  = element_text(size = 6))

priorpbmod3 <- posterior_vs_prior(bmod3, prob = 0.95) + 
  ggtitle(expression(CMJ[Gomez]*" ~ CMJ")) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text.x  = element_text(size = 6))

priorpbmod4 <- posterior_vs_prior(bmod4, prob = 0.95) + 
  ggtitle(expression(CMJ[Gomez]*" ~ "*PE[CMJ])) +
  theme(legend.position = "none",
        plot.title = element_text(size = 12),
        axis.text.x  = element_text(size = 6))

ggarrange(priorpbmod1,
          priorpbmod2,
          priorpbmod3,
          priorpbmod4,
          ncol = 2,
          nrow = 2)

ggsave("FigureS2.png", height=7, width=7, units='in', dpi=600)

#### Posterior predictive checkings
ppbmod1 <- pp_check(bmod1) +
  ggtitle(expression("Outcome: "*CMJ[Duncan])) +
  xlab("CMJ") 
  
ppbmod2 <- pp_check(bmod2) +
  ggtitle(expression("Outcome: "*CMJ[Duncan])) +
  xlab(expression(PE[CMJ])) 

ppbmod3 <- pp_check(bmod3) +
  ggtitle(expression("Outcome: "*CMJ[Gomez])) +
  xlab("CMJ") 

ppbmod4 <- pp_check(bmod4) +
  ggtitle(expression("Outcome: "*CMJ[Gomez])) +
  xlab(expression(PE[CMJ])) 

ggarrange(ppbmod1,
          ppbmod2,
          ppbmod3,
          ppbmod4,
          ncol = 2,
          nrow = 2)
  
ggsave("FigureS3.png", height=7, width=7, units='in', dpi=600)
  
