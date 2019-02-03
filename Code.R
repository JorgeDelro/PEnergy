## Bayesian Hierarchical model

library(tidyverse)
library(ggpubr)
library(brms)

# Read database
db <- read.csv2("PEnergy_db.csv")
str(db)


# Variables
db$BMIGroup <- factor(db$BMIGroup, levels = c("bajopeso", "normopeso", "sobrepeso", "obeso"),
                      labels = c("Underweight", "Normalweight", "Overweight", "Obese"))
db$Sex <- factor(db$Sex, levels = c(1, 2),
                 labels = c("chico", "chica"))
db$CMJcm <- db$CMJ*100

# Descriptives
# CMJ
CMJ_box <-  ggplot(data = db, aes(x = BMIGroup, y = CMJcm)) +
  geom_boxplot() +
  xlab("BMI") +
  ylab("CMJ (m)") +
  theme_bw()

# PE
PEcmj_box <-  ggplot(data = db, aes(x = BMIGroup, y = PEcmj)) +
  geom_boxplot() +
  xlab("BMI") +
  ylab(expression(PE[CMJ]*" "*(Joules))) +
  theme_bw()

# Power by Duncan formula
PDuncan_box <-  ggplot(data = db, aes(x = BMIGroup, y = POWCMJDuncan)) +
  geom_boxplot() +
  xlab("BMI") +
  ylab(expression(CMJ[Duncan]*" "*(W))) +
  theme_bw()

# Power by Gomez formula
PGomez_box <-  ggplot(data = db, aes(x = BMIGroup, y = POWCMJGomez)) +
  geom_boxplot() +
  xlab("BMI") +
  ylab(expression(CMJ[Gomez]*" "*(W))) +
  theme_bw()

ggarrange(CMJ_box,
          PEcmj_box,
          PDuncan_box,
          PGomez_box,
          nrow = 2,
          ncol = 2)

# Population effects model
bmod0 <- brm(formula = POWCMJDuncan ~ Age + Sex + CMJcm + BMIGroup,
            data = db,
            family = student(link = "identity"),
            chains = 4,
            warmup = 1000,
            iter = 2000,
            sample_prior = c("no"),
            control = list(adapt_delta = 0.9))

summary(bmod0)

# Random intercepts
bmod1 <- update(bmod0, formula = POWCMJDuncan ~ Age + Sex + CMJcm + (1|BMIGroup),
                control = list(adapt_delta = 0.9))
summary(bmod1)


# Full model
bmod1 <- update(bmod0, formula = POWCMJDuncan ~ Age + Sex + (CMJcm|BMIGroup),
                control = list(adapt_delta = 0.9))
summary(bmod1)



