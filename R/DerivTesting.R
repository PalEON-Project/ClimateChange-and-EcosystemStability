df.tmp$deriv.fac <- as.vector(lD$fac.fit$deriv)
df.tmp$deriv.Year <- as.vector(lD$Year$deriv)
df.tmp$deriv.sd <- as.vector(lD$Year$se.deriv)
df.tmp$lwr2 <- df.tmp$deriv.Year - 2*df.tmp$deriv.sd
df.tmp$upr2 <- df.tmp$deriv.Year + 2*df.tmp$deriv.sd


plot(mean ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], ylim=range(df.tmp[df.tmp$fac.fit=="PHA", c("mean", "deriv.Year")]), type="l", lwd=2)
lines(deriv.Year ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], type="l", lwd=2, col="blue")

plot(mean ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], ylim=range(df.tmp[df.tmp$fac.fit=="PHA", c("mean", "lwr", "upr", "lwr2", "upr2")]), type="l", lwd=2)
lines(lwr ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], type="l", lty="dashed", lwd=1, col="black")
lines(upr ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], type="l", lty="dashed", lwd=1, col="black")
abline(h=0, col="red")
# lines(deriv.fac ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], type="l", lwd=2, col="red")
plot(deriv.Year ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], ylim=range(df.tmp[df.tmp$fac.fit=="PHA", c("mean", "lwr", "upr", "lwr2", "upr2")]), type="l", lwd=2, col="blue")
lines(lwr2 ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], type="l", lty="dashed", lwd=1, col="blue")
lines(upr2 ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], type="l", lty="dashed", lwd=1, col="blue")
# lines(deriv3 ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], type="l", lwd=2, col="orange")
abline(h=0, col="red")


plot(V1 ~ Year, data=df.tmp[df.tmp$fac.fit=="PHA",], ylim=range(df.tmp[df.tmp$fac.fit=="PHA", c("mean", "lwr", "upr")]), type="l", lwd=2)



plot(mean ~ Year, data=df.tmp[df.tmp$fac.fit=="PBL",], ylim=range(df.tmp[df.tmp$fac.fit=="PBL", c("mean", "lwr", "upr")]), type="l", lwd=2)
lines(lwr ~ Year, data=df.tmp[df.tmp$fac.fit=="PBL",], type="l", lty="dashed", lwd=1, col="red")
lines(upr ~ Year, data=df.tmp[df.tmp$fac.fit=="PBL",], type="l", lty="dashed", lwd=1, col="red")
abline(h=0, col="blue")


plot(mean ~ Year, data=mod.predict[mod.predict$Site=="PHA",], ylim=range(mod.predict[mod.predict$Site=="PHA", c("mean", "lwr", "upr")]), type="l", lwd=2)
lines(lwr ~ Year, data=mod.predict[mod.predict$Site=="PHA",], type="l", lty="dashed", lwd=1, col="red")
lines(upr ~ Year, data=mod.predict[mod.predict$Site=="PHA",], type="l", lty="dashed", lwd=1, col="red")
abline(h=mean(mod.predict[mod.predict$Site=="PHA","mean"], na.rm=T), col="blue")
