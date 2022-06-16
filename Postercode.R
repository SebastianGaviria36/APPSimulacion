library(survival)
library(gamlss)
library(ggplot2)
library(dplyr)
library(plotly)

#Building auxiliar functions
censdata <- function(n, censrate, dist){
  #dist is a number generator function of 
  #a specified distribution. must be parameterized 
  #and its only argument must be n 
  nobs <- n/2
  g1 <- dist(nobs)
  g2 <- rweibull(nobs,3,200)
  groups <- rep(c(1,2),each=nobs)
  statusg1 <- ifelse(quantile(g1, 1-censrate)<g1,0,1)
  statusg2 <- ifelse(quantile(g2, 1-censrate)<g2,0,1)
  data.frame(survtime = c(g1,g2),
             status = c(statusg1,statusg2),
             group = groups)
}

mysurvdiff <- function(data,test){
  #test=0 is log rank, test=1 is peto
  chisq <- survdiff(Surv(survtime,status)~group,data=data,rho=test)$chisq
  ifelse(1-pchisq(chisq,1)<0.05,1,0)
}

#Setting distribution to simulate data
survdata <- survival::ovarian #the simulation will based in ovarian data
bestdist <- fitDist(survdata$futime,type = "realplus")$family #best dist is weibull
params <- MASS::fitdistr(survdata$futime,"weibull")$estimate #finding distribution parameters
plot(density(survdata$futime)) #cheking adjust
lines(density(rweibull(1000,params[1],params[2])), col = "blue")

#Building simulation function
logranksim <- function(n,censrate,test,iter=10000){
  f <- function(n, censrate,test){
    datoss <- censdata(n,censrate, function(n) rweibull(n,params[1],params[2]))
    mysurvdiff(datoss,test)}
  mean(replicate(iter, f(n,censrate,test)))
}

grid <- expand.grid(n = c(10,20,40,80,160), censrate = seq(0,0.8,0.2), test = 0:1)

power <- c()

for (i in 1:nrow(grid)){
  power[i] <- logranksim(grid[i,1],grid[i,2],grid[i,3])
}

simres <- cbind(grid,power)

simres <- simres %>%
  mutate(test = if_else(test==0,"LogRank","Peto"),
         n = as.factor(n),
         censrate = as.factor(censrate))

#-----------------------Finish simulation----------------------------------------

ggplot(simres, aes(x=censrate, y = power, group = n, colour = n)) + 
  facet_wrap(~test) +
  geom_line(cex=0.8) +
  labs(x = "Tasa de censura", y = "Potencia", title = "Potencia Test:") + 
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5))



plot_ly(x=simres[simres$test=="LogRank",2], 
        y=simres[simres$test=="LogRank",1], 
        z=simres[simres$test=="LogRank",4], 
        type="mesh3d",
        intensity = simres[simres$test=="LogRank",4],
        color = simres[simres$test=="LogRank",4]) %>%
  layout(title = "Potencia LogRank",
         scene = list(xaxis=list(title = "Tasa de censura"),
                      yaxis=list(title = "n"),
                      zaxis=list(title = "Potencia")))


plot_ly(x=simres[simres$test=="Peto",2], 
        y=simres[simres$test=="Peto",1], 
        z=simres[simres$test=="Peto",4], 
        type="mesh3d",
        intensity = simres[simres$test=="Peto",4],
        color = simres[simres$test=="Peto",4]) %>%
  layout(title = "Potencia Peto",
         scene = list(xaxis=list(title = "Tasa de censura"),
                      yaxis=list(title = "n"),
                      zaxis=list(title = "Potencia")))



