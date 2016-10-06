# Showing rate of change distributions as histograms
# Best option will be to fit a distribution to each variable and 
# then generate a probability distribution from that fit
library(fitdistrplus)

# figure out the best distribution for each model/variable combo
gpp.ed  <- ecosys.change[ecosys.change$Model=="ed2" & ecosys.change$var=="GPP" & ecosys.change$Year<1850, ]
gpp.clm <- ecosys.change[ecosys.change$Model=="clm.cn" & ecosys.change$var=="GPP" & ecosys.change$Year<1850, ]


# Adding a relativized derivative -- relative to current state
ecosys.change$deriv.rel <- ecosys.change$deriv.mean / ecosys.change$Y
summary(ecosys.change)



rel.range <- c(quantile(ecosys.change[ecosys.change$var=="GPP" & ecosys.change$Year<1850, "deriv.rel"], 0.05), quantile(ecosys.change[ecosys.change$var=="GPP" & ecosys.change$Year<1850, "deriv.rel"], 0.95))
ggplot(data= ecosys.change[ecosys.change$var=="GPP" & ecosys.change$Year<1850, ]) +
  facet_wrap(~Model, scales="fixed") +
  geom_histogram(aes(x=deriv.rel, fill=Model)) +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_x_continuous(limits=rel.range)

dist.out <- NULL
v="LAI"
x.range <- seq(quantile(ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850, "deriv.mean"], 0.05), 
               quantile(ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850, "deriv.mean"], 0.95), 
               length.out=100)
x.rel   <- seq(quantile(ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850, "deriv.rel"], 0.05), 
               quantile(ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850, "deriv.rel"], 0.95), 
               length.out=100)
ggplot(data= ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850, ]) +
  facet_wrap(~Model, scales="free") +
  geom_histogram(aes(x=deriv.mean, fill=sig)) +
  geom_vline(xintercept=0, linetype="dashed")
ggplot(data= ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850, ]) +
  facet_wrap(~Model, scales="fixed") +
  geom_histogram(aes(x=deriv.rel, fill=sig)) +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_x_continuous(limits=rel.range)

distr="normal"
for(m in unique(ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850, "Model.Order"])){
  fit.rel <- fitdistr(ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850 & ecosys.change$Model.Order==m, "deriv.rel"], distr)
  prob.rel <- dnorm(x.rel, mean=fit.rel$estimate[[1]], sd=fit.rel$estimate[[2]])
  prob.rel <- prob.rel / sum(prob.rel)
  
  fit.raw <- fitdistr(ecosys.change[ecosys.change$var==v & ecosys.change$Year<1850 & ecosys.change$Model.Order==m, "deriv.mean"] , distr)
  prob.raw <- dnorm(x.range, mean=fit.raw$estimate[[1]], sd=fit.raw$estimate[[2]])
  prob.raw <- prob.raw / sum(prob.raw)
  
  dist.out <- rbind(dist.out, data.frame(Model=rep(m, length(prob.raw)), var=v, x=x.range, x.rel=x.rel, prob.raw=prob.raw, prob.rel=prob.rel))
}

ggplot(data=dist.out) +
  geom_line(aes(x=x, y=prob.raw, color=Model)) +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_x_continuous(expand=c(0,0))

ggplot(data=dist.out) +
  geom_line(aes(x=x.rel, y=prob.rel, color=Model)) +
  geom_vline(xintercept=0, linetype="dashed") +
  scale_x_continuous(expand=c(0,0))


test1 <- fitdist(gpp.ed$deriv.mean, "gamma", method = "mme")
test2 <- fitdist(gpp.ed$deriv.mean, "norm", method = "mme")
test3 <- fitdist(gpp.ed$deriv.mean, "exp", method = "mme")
test4 <- fitdist(gpp.ed$deriv.mean, "geom", method = "mme")

gofstat(list(test1, test2, test3, test4), fitnames=c("gamma", "normal", "exponential", "geometric"))
cdfcomp(list(test1, test2, test3, test4), fitnames=c("gamma", "normal", "exponential", "geometric"))


ks.test(gpp.ed$deriv.mean, "pnorm", mean=test2$estimate[[1]], sd=test2$estimate[[2]])
ks.test(gpp.ed$deriv.mean, "pgamma", shape=test1$estimate[[1]], rate=test1$estimate[[2]])



library(MASS)
fit.ed <- fitdistr(gpp.ed$deriv.mean, "normal")
mean.ed = fit.ed$estimate[[1]]
sd.ed = fit.ed$estimate[[2]]
prob.ed <- dnorm(x.range, mean=mean.ed, sd=sd.ed)
prob.ed <- prob.ed / sum(prob.ed)




fit.clm <- fitdistr(gpp.clm$deriv.mean, "normal")
mean.clm = fit.clm$estimate[[1]]
sd.clm = fit.clm$estimate[[2]]
prob.clm <- dnorm(x.range, mean=mean.clm, sd=sd.clm)
prob.clm <- prob.clm / sum(prob.clm)

deriv.dist <- data.frame(Model=c(rep("ed2", length(prob.ed)), rep("clm", length(prob.clm))), x=rep(x.range, 2), prob=c(prob.ed, prob.clm))

# add to results
results[i,] = c(gf_shape, est_mean, est_sd, ks$statistic, ks$p.value)


dummy.dat <- 
  
  res = fitData(dat.test, fit=c("logistic","normal","exponential","poisson"),
                sample=1)
ed.gpp <- fitdist(abs(dat.test*1e10), distr="beta")
