### This script is used to generate a death curve after infection ###

library(survival)

colors1=c("darkorchid3","gray0")
time = c(7,7,7,8,15,12,7,4,7,9,8,15,6,5,4,7,6,9,9,6,4,6,4,5) #Time of event
status = c(1,1,1,1,0,1,1,1,1,1,1,0,1,1,1,1,1,1,1,1,1,1,1,1) #Nature of event
group = c("F","F","F","F","F","F","F","F","F","F","F","F","S","S","S","S","S","S","S","S","S","S","S","S") #What group the event belongs to
sdata=data.frame(time,status,group)
formula1 <- Surv(time, status == 1) ~ group
test_result <- survdiff(formula = formula1, data = sdata, rho = 0)
test_result
plot(survfit(formula1), lty=c(1:2), col=colors1,xlab = "Time (days)", ylab = "Survival probability")
legend("topright", legend=c("Fed","Starved"),col=colors1, lty=1, cex=0.6)
