#survcurve
nm <- read.csv("http://www.sgi.com/tech/mlc/db/churn.names", 
               skip = 4, 
               colClasses = c("character", "NULL"), 
               header = FALSE, 
               sep = ":")[[1]]

dat <- read.csv("http://www.sgi.com/tech/mlc/db/churn.data", 
                header = FALSE, 
                col.names = c(nm, "Churn"),
                colClasses = c("factor",
                               "numeric",
                               "factor",
                               "character",
                               rep("factor", 2),
                               rep("numeric", 14),
                               "factor"))
# test data

test <- read.csv("http://www.sgi.com/tech/mlc/db/churn.test", 
                 header = FALSE, 
                 col.names = c(nm, "Churn"),
                 colClasses = c("factor",
                                "numeric",
                                "factor",
                                "character",
                                rep("factor", 2),
                                rep("numeric", 14),
                                "factor"))

library(survival)
dat$Churn <- as.numeric(dat$Churn) - 1

s <- with(dat, 
          Surv(account.length, Churn))

head(s, n = 100)

## Kaplan-Meier estimator. The "log-log" confidence interval is preferred.
km.as.one <- survfit(s ~ international.plan, #or area.code,
                     data = dat, 
                     conf.type = "log-log")

## Show object
km.as.one

library(survminer)

ggsurvplot(km.as.one,
           break.time.by = 12,
           color = "red",
           surv.scale = "percent",
           xlab = "Tenure- months",
           ylab = "% survived",
           ylim = c(0.2, 1),
           censor = T,
           legend = "top",
           risk.table = T,
           ggtheme = theme_gray())
