library('forecast')
library('zoo')
library('tseries')
library('fUnitRoots')
data <- read.csv2('data/valeurs_mensuelles.csv')
data <- data[4:nrow(data),]
data <- data.frame(date = as.yearmon(data$Libellé),
                   x =as.numeric(data$Indice.CVS.CJO.de.la.production.industrielle..base.100.en.2021....Commerce.d.électricité..NAF.rév..2..niveau.classe..poste.35.14.))
x <- zoo(data$x, order.by = data$date)
dx <- diff(x,1)

plot(x)
plot(dx)

Qtests <- function(series, k, fitdf=0) {
  pvals <- apply(matrix(1:k), 1, FUN=function(l) {
    pval <- if (l<=fitdf) NA else Box.test(series, lag=l, type="Ljung-Box", fitdf=fitdf)$p.value
    return(c("lag"=l,"pval"=pval))
  })
  return(t(pvals))
}

adfLag <-function(series,kmax,type){ #ADF tests until no more autocorrelated residuals
  k <- 0
  noautocorr <- 0
  while (noautocorr==0){
    
    adf <- adfTest(series,lags=k,type=type)
    pvals <- Qtests(adf@test$lm$residuals,24,fitdf=length(adf@test$lm$coefficients))[,2]
    if (sum(pvals<0.05,na.rm=T) == 0) {
      noautocorr <- 1}
    k <- k + 1
  }
  return(k-1)
}

summary(lm(dx ~ index(dx)))

#adf test on dx
k <- suppressWarnings(adfLag(dx,60,"nc"))
adf <- suppressWarnings(adfTest(dx,lags=k,type='nc'))
adf
#pp test on dx
pp.test(dx)
#kpss test on dx
kpss.test(dx)

#acf dx
acf <- acf(zoo(coredata(dx), order.by = seq_along(index(dx))), lag.max = 24, plot = TRUE, main = 'ACF of first difference')
which(abs(acf$acf)>= 1.96/sqrt(acf$n.used)) # select q_max =14

#pacf
pacf <- pacf(zoo(coredata(dx), order.by = seq_along(index(dx))), lag.max = 24, plot = TRUE,main = 'PACF of first difference')
which(abs(pacf$acf)>= 1.96/sqrt(pacf$n.used)) # select p_max =13

#Calculate AICs & BICs for p<=pmax and q<=qmax
pmax = 13
qmax = 14
mat <- matrix(NA,nrow=pmax+1,ncol=qmax+1) #empty matrix to fill
rownames(mat) <- paste0("p=",0:pmax) #renames lines
colnames(mat) <- paste0("q=",0:qmax) #renames columns
AICs <- mat #AIC matrix not filled non remplie
BICs <- mat #BIC matrix not filled non remplie
pqs <- expand.grid(0:pmax,0:qmax) #all possible combinations of p and q
for (row in 1:dim(pqs)[1]){ #loop for each (p,q)
  p <- pqs[row,1] #gets p
  q <- pqs[row,2] #gets q
  estim <- try(arima(x,c(p,1,q),include.mean = F,optim.control = list(maxit = 50000))) #tries ARIMA estimation
  AICs[p+1,q+1] <- if (class(estim)=="try-error") NA else estim$aic #assigns the AIC
  BICs[p+1,q+1] <- if (class(estim)=="try-error") NA else BIC(estim) #assigns the BIC
  print(row/dim(pqs)[1])
}
AICs
BICs
#best choice according to AICs
AICs == min(AICs, na.rm = TRUE) # best choice is ARIMA(12,1,10)
                                # AIC = 2322.348
                                # BIC = 2415.328
#best choice according to BICs
BICs == min(BICs, na.rm = TRUE) # best choice is ARIMA(1,1,1)
                                # AIC = 2330.655
                                # BIC =2342.783

#function to Analyse the significativity of the coef of ARIMA
signif <- function(estim){
  coef <- estim$coef
  se <- sqrt(diag(estim$var.coef))
  t <- coef/se
  pval <- (1-pnorm(abs(t)))^2
  return(rbind(coef,se,pval))
}
#Best choice according to AICs have parameter which aren't signif (ex : ma5)
arima_AIC <- arima(x,c(12,1,10),include.mean = F,optim.control = list(maxit = 50000))
signif(arima_AIC)

#all coef of best BICs ARIMA are signif
arima_BIC <- arima(x,c(1,1,1),include.mean = F)
signif(arima_BIC)

#Arima_bic model's residuals are not correlated 
#(using ljung-box test with null hypothesis: no autocorrelation for 1 to r; 
#all test for r=2 to 24 fail to reject null hypothesis )
Qtests(arima_BIC$residuals,24,fitdf=2)
#only ljung-box test for lag 24
Box.test(arima_BIC$residuals, lag = 24, type = "Ljung-Box", fitdf = 2)

#we choose arima_bic model
arima_BIC

#residual analysis
qqnorm(arima_BIC$residuals)
qqline(arima_BIC$residuals, col = "blue", lty = 2)
plot(arima_BIC$residuals, ylab = 'residuals')
jarque.bera.test(arima_BIC$residuals)#residuals do not follow gaussian


