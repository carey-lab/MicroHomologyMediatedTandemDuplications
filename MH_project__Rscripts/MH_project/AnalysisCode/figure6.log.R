# load libraries
library(magrittr)
#library(count)
library(plyr)
library(dplyr)
library(seqinr)

# set working directory and load data
setwd("E:\\work\\SynologyDrive\\Projects\\2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU\\Sarah\\MH_project\\ProcessedData")
setwd('~/SynologyDrive/Projects/2019__MicroHomologyMediatedIndels__XiangweHe_ZhejiangU/Sarah/MH_project/ProcessedData') # for Lucas
MHlen_Dup <- read.table("10k.sign.count.tsv.MHlen_Dup.txt",header = F, comment.char = "#",as.is = T)
MH_seq <- read.table("10k.sign.count.tsv.mh_seq.txt",header = F, comment.char = "#",as.is = T)
MH_seq <- read.table("10k.sign.count.tsv.MHseq.txt",header = F, comment.char = "#",as.is = T) # for Lucas
MH_inter <- read.table("10k.sign.count.tsv.inter_MH.txt",header = F, comment.char = "#",as.is = T)
MH_inter <- read.table("10k.sign.count.tsv.interMH.txt",header = F, comment.char = "#",as.is = T) # for Lucas
dup <- MHlen_Dup[,2]
MHlen <- MHlen_Dup[,1]
GC <- apply(MH_seq,1,function(x){sum(unlist(strsplit(x,''))=="C")+sum(unlist(strsplit(x,''))=="G")})



re_model <- zelig(y ~ x1, data = train, model = "relogit", cite = FALSE )
z.out <- zelig(Y ~ X1 + X2, model = "relogit", tau = NULL,
               case.control = c("prior", "weighting"),
               bias.correct = TRUE,
               data = mydata, ...)
x.out <- setx(z.out)
s.out <- sim(z.out, x = x.out)
Arguments




###lu
set.seed(123456)
signal_to_noise <- .2
n <- 750
prop_ones <- .03

beta_1 <- sqrt(signal_to_noise)


x1 = rnorm(n)
intercept = qnorm(prop_ones, sd = sqrt(beta_1 ^ 2 +1))
y_signal = beta_1 * x1
y_star = intercept + y_signal + rnorm(n)
y = y_star >= 0
true_prob = pnorm(y_signal + intercept)
dat <- data.frame("x1"=x1,"intercept"=intercept,"y_signal"=y_signal,"y_star"=y_star,"y"=y,"true_prob"=true_prob)

dat %>%
  count(y) %>%
  mutate(prop = round(n/sum(n),2))

train <- dat %>%
  sample_frac(.66)

test <- dat %>%
  anti_join(train, by = c("x1", "y_star", "y", "true_prob"))

X_train <- train %>%
  select(x1)

X_test <- test %>%
  select(x1)

fit <- glm(as.factor(y) ~ x1, data = train, family = "binomial")  

coefs <- fit %>%
  coefficients() %>%     ###coef()
  as.numeric()

V <- vcov(fit)

logisticPred <- function(X, coef) {
  
  X %>%
    na.omit() %>%
    mutate(int = 1) %>%
    select(int, everything()) %>%
    as.matrix(.) %*% coef %>%
    as.vector() %>%
    (function(x) 1 / (1 + exp(-x)))
}


pred <- X_train %>%
  logisticPred(coefs)

test_pred <- X_test %>%
  logisticPred(coefs)


X_matrix <- X_train %>%
  mutate(bias = 1) %>%
  select(bias, everything()) %>%
  as.matrix()

X_test_matrix <- X_test %>%
  mutate(bias = 1) %>%
  select(bias, everything()) %>%
  as.matrix()

W <- diag(pred * (1 - pred))
Q <- X_matrix %*% solve(t(X_matrix) %*% W %*% X_matrix) %*% t(X_matrix)
e <- 0.5 * diag(Q) * (2 *pred - 1)

bias <- (solve(t(X_matrix) %*% W %*% X_matrix) %*% t(X_matrix) %*% W %*% e) %>%
  as.numeric()


unbiased_coefs <- coefs - bias
updated_var <- (nrow(X_train) / (nrow(X_train) + ncol(X_train) +1)) ^ 2 * V

library(Zelig)
re_model <- zelig(y ~ x1, data = train, model = "relogit", cite = FALSE )
unbiased_coefs
coef(re_model)
sqrt(diag(updated_var))
sqrt(diag(as.matrix(vcov(re_model)[[1]])))
zelig(y ~ x1, data = train, model = "relogit", cite = FALSE) %>%
  summary

glm(y ~ x1, data = train, family = "binomial") %>%
  summary

unbiased_coef_pred <- X_test %>%
  logisticPred(unbiased_coefs)


new_X_test <- test %>%
  as_data_frame %>%
  mutate(p = unbiased_coef_pred) %>%
  mutate(zeta = (.5-p)*p*(1-p),
         eta = (1*updated_var[1,1] + x1*updated_var[1,2])*1 +
           (1*updated_var[2,1] + x1*updated_var[2,2])*x1) %>%
  mutate(updated_p = p + zeta*eta) %>%
  mutate(event = updated_p > .5)

corrected_pred <- unlist(new_X_test$updated_p)
