install.packages("pROC")
library(pROC)

data <- read.csv('ROC data.csv')
head(data)

#area under the curve
roc(data$outcome, data$cell_viability)

#smoothing (only if enough patients available for)
roc(outcome ~ cell_viability, data, smooth=TRUE)

#make curve for cell viability
roc1 <- roc(data$outcome,
            data$cell_viability, percent=TRUE, 
            # arguments for auc
            partial.auc=c(100, 80), partial.auc.correct=TRUE,
            partial.auc.focus="spec",
            # arguments for ci
            ci=TRUE, boot.n=100, ci.alpha=0.9, stratified=FALSE,
            # arguments for plot
            plot=TRUE, auc.polygon=TRUE, max.auc.polygon=TRUE, grid=TRUE,
            print.auc=TRUE, show.thres=TRUE)

# Add curves for speed of cancer spread and max z position
roc2 <- roc(data$outcome, data$speed_of_cancer_spread,
            plot=TRUE, add=TRUE, percent=roc1$percent)

roc3 <- roc(data$outcome, data$max_z_position,
            plot=TRUE, add=TRUE, percent=roc1$percent)

#specifity, negative predictive value, sensitivity and positive predictive value
coords(roc1, "best", ret=c("threshold", "specificity", "npv","sensitivity", "ppv"))

# Confidence interval of AUC
ci(roc1)

#plot confidence interval of curve
sens.ci <- ci.se(roc1, specificities=seq(0, 100, 5))
plot(sens.ci, type="shape", col="green")
plot(sens.ci, type="bars")

#compare 2 curves (ex. cell viability vs max z position as predictiors of pCR)
roc.test(roc1, roc3, reuse.auc=FALSE, partial.auc=c(100, 80),	
	partial.auc.focus="spec", partial.auc.correct=TRUE)

#population (cases and controls) needed for 0.9 power, p<0.05
power.roc.test(auc=0.8, power=0.8)

#calculate p-value against a null hypotesis of 50/50 
roc4 <- roc(data$outcome, data$control,
            plot=TRUE, add=TRUE, percent=roc1$percent)

roc.test(roc1, roc4, reuse.auc=FALSE, partial.auc=c(100, 80),	
	partial.auc.focus="spec", partial.auc.correct=TRUE)


