install.packages("aod")

library(aod)
library(ggplot2)

data <- read.csv("logistic regression.csv", head = TRUE)
head(data)

summary(data)
sapply(data, sd)

#contingency table for discrete/categorical variables
xtabs(~pCR + drug, data = data)
xtabs(~pCR + agreement, data = data)
#convert categorical variable to factor
data$drug <- factor(data$drug)
data$agreement <- factor(data$agreement)

#calculate coefficiencts and p-values
#logit by drug
mylogit1 <- glm(pCR ~ cell_viability + speed_of_cancer_spread + max_z_position + drug, data = data, family = "binomial")
summary(mylogit1)
#logit by agreement
mylogit2 <- glm(pCR ~ cell_viability + speed_of_cancer_spread + max_z_position + agreement, data = data, family = "binomial")
summary(mylogit2)
#logit by drug and agreement
mylogit3 <- glm(pCR ~ cell_viability + speed_of_cancer_spread + max_z_position + drug + agreement, data = data, family = "binomial")
summary(mylogit3)

#confidence intervals for log likelihood
confint(mylogit1)
confint(mylogit2)
confint(mylogit3)
# CIs using standard errors
confint.default(mylogit1)
confint.default(mylogit2)
confint.default(mylogit3)

#compare response categories from device with Wald Test
wald.test(b = coef(mylogit), Sigma = vcov(mylogit3), Terms = 5:6)
wald.test(b = coef(mylogit), Sigma = vcov(mylogit3), Terms = 7)
#if p<0.05, the effects of the categorical variable are significant


#compare drug 2 with drug 3
#l <- cbind(0, 0, 0, 0, 1, -1,0)
#wald.test(b = coef(mylogit3), Sigma = vcov(mylogit3), L = l)
#if p<0.05, the 2 categories are signficantly different
#similar comparisons can be made between any 2 categories


#get odds ratios
exp(coef(mylogit1))
exp(cbind(OR = coef(mylogit1), confint(mylogit1)))

exp(coef(mylogit2))
exp(cbind(OR = coef(mylogit2), confint(mylogit2)))

exp(coef(mylogit3))
exp(cbind(OR = coef(mylogit3), confint(mylogit3)))


#graphing probability of success (pCR) for each category (agreement, drug)
#can only graph 1 categorical variable at a time, so change code when comparing drugs
#drug=factor(1:3)

newdata1 <- with(data, data.frame(cell_viability = mean(cell_viability), 
		speed_of_cancer_spread = mean(speed_of_cancer_spread), 
		max_z_position = mean(max_z_position),
		agreement=factor(1:2)))
newdata1

#agreementP column will indicate the probability that each category acheives pCR
newdata1$agreementP <- predict(mylogit2, newdata = newdata1, type = "response")
newdata1

#create points for graphing predicted pCR probabilities for a categorical variable
#ex. cell viability is varied and speed of spread and max z position are kept at mean
#do this for any variable you want to plot while controlling for others at their mean

newdata2 <- with(data, data.frame(cell_viability = rep(seq(from = 0, to = 100, length.out = 100),2), 
		speed_of_cancer_spread = mean(speed_of_cancer_spread), 
		max_z_position = mean(max_z_position),
		agreement = factor(rep(1:2, each = 100))))

newdata3 <- cbind(newdata2, predict(mylogit2, newdata = newdata2, type = "link",
    se = TRUE))
newdata3 <- within(newdata3, {
    PredictedProb <- plogis(fit)
    LL <- plogis(fit - (1.96 * se.fit))
    UL <- plogis(fit + (1.96 * se.fit))
})

ggplot(newdata3, aes(x = cell_viability, y = PredictedProb)) + geom_ribbon(aes(ymin = LL,
    ymax = UL, fill = agreement), alpha = 0.2) + geom_line(aes(colour = agreement),
    size = 1)

#check model fit and compare models
#deviance (chi-square)
with(mylogit2, null.deviance - deviance)
#degrees of freedom
with(mylogit2, df.null - df.residual)
#p-value
with(mylogit2, pchisq(null.deviance - deviance, df.null - df.residual, lower.tail = FALSE))
#if p<0.05, the model is better than a null model
#check model log likelihood
logLik(mylogit2)