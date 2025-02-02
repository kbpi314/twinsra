# glm for twins ra
df <- read.csv('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_glm_R.csv')

# glms
model <- lm(formula = shannon_entropy_diff ~ Abx_diff + smoking.previously_diff + Alcohol.current_diff + NSAIDS.current_diff + NSAIDS.last.3.mo_diff, data=df)
summary(model)

model <- lm(formula = Bray_Curtis ~ Abx_diff + smoking.previously_diff + Alcohol.current_diff + NSAIDS.current_diff + NSAIDS.last.3.mo_diff, data=df)
summary(model)


# 
library(lme4)
df <- read.csv('/Users/KevinBu/Desktop/clemente_lab/Projects/twinsra/inputs/df_lmer.csv')


# glms
model <- lm(formula = shannon_entropy_diff ~ MTX, data=df)
summary(model)

# glms
model <- lm(formula = Bray_Curtis ~ MTX, data=df)
summary(model)
