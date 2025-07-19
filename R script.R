library(quadprog) #Needed for solve.QP
library(reshape2) # Used to melt the data
library(tseries)
library(quantmod)
library(PerformanceAnalytics)
library(ggplot2) # Used to graph efficient frontier
library(IntroCompFinR)

tckk <- c("AAPL", "TSLA","MSFT","AMD") #Choose ticker symbols.
numtk <- length(tckk);
ustart <- "2020-01-01"; #Define start date.
uend <- "2024-1-20" #Define end date.
all_dat <- list();
for(i in 1:numtk) {
  all_dat[[i]] <- xxx <- get.hist.quote(
    instrument = tckk[i],
    start=ustart, end=uend,
    quote = c("Close"), provider = "yahoo", compression = "m")
}


APP <- as.xts(all_dat[[1]]) #Apple
TSLA <- as.xts(all_dat[[2]]) #Tesla
MSFT <- as.xts(all_dat[[3]]) #Microsoft
AMD <- as.xts(all_dat[[4]]) #AMD
plot(APP)
plot(TSLA)
plot(AMD)

ret <- cbind(APP,TSLA,MSFT,AMD) 
returns <-apply(log(ret), 2, diff) 
colnames(returns) <- c("APP", "TSLA","MSFT", "AMD") 
colnames(returns) 
#1 "APP"  "TSLA" "MSFT" "AMD"
# plot return 
Return.APP  <- apply(log(APP), 2, diff) 
Return.TSLA <- apply(log(TSLA), 2, diff) 
Return.MSFT <- apply(log(MSFT), 2, diff) 
Return.AMD <- apply(log(AMD), 2, diff) 
plot(xts(Return.APP,order.by = as.Date(row.names(Return.APP))), 
     main = "The montly return of APP")
plot(xts(Return.TSLA,order.by = as.Date(row.names(Return.TSLA))), 
     main = "The montly return of TSLA") 
plot(xts(Return.MSFT,order.by = as.Date(row.names(Return.MSFT))), 
     main = "The montly return of MSFT") 

covmat<-cov(returns) #variance-covariance matrix
covmat
mearR <- colMeans(returns)  #expected return 
mearR 

Gport <- globalMin.portfolio(er=mearR , cov.mat=covmat)
summary (Gport)
plot(Gport)


#no ness
sapply(as.data.frame(returns), summary) 
ns <- ncol(returns) 
er <- apply(returns, 2, mean) #expected return 
std <- apply(returns, 2, sd) #standard deviation 
er
std
covm <- cov(returns) #variance-covariance matrix 
covm 
#Global minimum variance (gmv) portfolio 
Am <- rbind(2*covm, rep(1, ns)) 
Am <- cbind(Am, c(rep(1, ns), 0)) 
b <- c( rep(0, ns), 1) 
w.gmv <- solve(Am) %*% b #last value is lambda constraint. Hence not relevant 
w.gmv <- w.gmv[-(ns+1), ] #sum of weights are 1 
sum(w.gmv) 

var.gmv <- t(w.gmv) %*% covm %*% w.gmv  
sigma <-sqrt(var.gmv) 
ret.gmv <- t(w.gmv) %*% er 
uv <- rep(1, ns) #unit vector 
w.gmv2 <- (solve(covm) %*% uv) * c(solve( (t(uv) %*% solve(covm) %*% uv))) 
#the results yield the same weight, expected return, and standard deviation 
ret.gmv #return of portfolio 


u0 <- 0.05 #desired expected return of the portfolio  
if(u0 < ret.gmv) { 
  message("#u0 should be greater than return on minimum variance portfolio
 ") 
} 

M <- cbind(er, uv)
B <-  t(M) %*% solve(covm) %*% M 
mu.tilde <- c(u0, 1) 
#weight of efficient portfolio 
w.ep <- solve(covm) %*% M %*% solve(B) %*% mu.tilde 
w.ep 
#portfolio expected return 
ret.ep <- t(er) %*% w.ep  
ret.ep 
#portofolio standard deviation 
var.ep <- t(w.ep) %*% covm %*% w.ep  
sigma_EP <- sqrt(var.ep) 
sigma_EP 

w1 <- .4
w <- seq(from = -5, to = 5, by = .01) #various combination of weights. 
eff <- function(w1) { 
  z <- w1 * w.gmv + (1 - w1) * w.ep 
  c(ret = t(z) %*% er  * 12 , sd = sqrt(t(z) %*% covm %*% z * 12)) 
} #annualize expected return and standard deviation. 
comb <- cbind(w, t(sapply(w, FUN = eff))) 

efficient <- comb[, "ret"] > c(ret.gmv * 12) # month expected return and standard deviation. 
xlim <- c(min(comb[, "sd"]), max(comb[, "sd"])) 
ylim <- c(min(comb[, "ret"]), max(comb[, "ret"])) 
col <- ifelse(efficient, "blue", "red") 
plot(comb[ , "ret"] ~ comb[, "sd"],  
     col = col, 
     xlim = xlim,  
     ylim = ylim,  
     xlab = "Portfolio RisK",  
     ylab = "Portfolio Return",  
     pch = 16,  
     main = "Efficient Frontier",  
     cex = .7)



#SR = (mu - rf) / sigma 
rfree <- .01 / 12  #define the risk-free interest rate and adjust to the same data frequency.  
(w.sr <- solve(covm) %*% (er - rfree) / c(t(uv) %*% solve(covm) %*% (er - rfree))) #weights of optimal portfolio 
(ret.sr <- t(w.sr) %*% er ) #expected return 
(var.sr <- t(w.sr)%*% covm %*% w.sr) #variance of portfolio 
#plotting results 
comb <- rbind(comb, c(NA, ret.sr * 12, sqrt(var.sr * 12))) #annualize exp return and standard deviation. 
col <- c(ifelse(efficient, "blue", "red"), "green") 
xlim <- c(0.15, max(comb[, "sd"])) 
ylim <- c(min(comb[,"ret"]), max(comb[, "ret"])) 
cex <- c(ifelse(efficient, .6, .6), 1.2) 
pch <- c(rep(1, nrow(comb) - 1), 16) 
plot(comb[ , "ret"] ~ comb[, "sd"],  
     col = col,  
     xlim = xlim,  
     ylim = ylim,  
     xlab = "Portfolio RisK",  
     ylab = "Portfolio Return",  
     pch = pch,  
     main = "Efficient Frontier",  
     cex = cex) 
#capital market line  
abline(a = rfree, b = ret.sr * 12 / sqrt(var.sr * 12), lty = 2) #annualize expected return and standard deviation.




