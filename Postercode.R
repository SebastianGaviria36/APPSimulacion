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

grid <- expand.grid(n = c(10,20,40,80,160,320,640), censrate = seq(0,0.8,0.1), test = 0:1)

power <- c()

for (i in 1:nrow(grid)){
  power[i] <- logranksim(grid[i,1],grid[i,2],grid[i,3])
}

simres <- cbind(grid,power)

saveRDS(simres, "simres2.Rds")

simres <- simres %>%
  mutate(test = if_else(test==0,"LogRank","Peto"),
         n = as.factor(n),
         censrate = as.factor(censrate))

#-----------------------Finish simulation----------------------------------------
 
simres <- simres[simres$n!=640,] #n = 640 not necessary
ggplot(simres, aes(x=censrate, y = power, group = n, colour = n)) + 
  facet_wrap(~test) +
  geom_line(cex=0.8) +
  labs(x = "Tasa de censura", y = "Potencia") + 
  theme_light() +
  theme(strip.text.x = element_text(
    size = 12, color = "black", face = "bold.italic"),
    plot.title = element_text(hjust = 0.5),
    text = element_text(size = 20),
    axis.text.x= element_text(size= 10))


plot_ly(x=simres[simres$test=="LogRank",2], 
        y=simres[simres$test=="LogRank",1], 
        z=simres[simres$test=="LogRank",4], 
        type="mesh3d",
        intensity = simres[simres$test=="LogRank",4],
        color = simres[simres$test=="LogRank",4]) %>%
  layout(title = "Potencia Test LogRank", font=list(size = 12),
         scene = list(xaxis=list(title = "Tasa de censura", 
                                 tickfont = list(size = 15),
                                 titlefont = list(size = 20)),
                      yaxis=list(title = "n", 
                                 tickfont = list(size = 15),
                                 titlefont = list(size = 20)),
                      zaxis=list(title = "Potencia", 
                                 tickfont = list(size = 15),
                                 titlefont = list(size = 20))))


plot_ly(x=simres[simres$test=="Peto",2], 
        y=simres[simres$test=="Peto",1], 
        z=simres[simres$test=="Peto",4], 
        type="mesh3d",
        intensity = simres[simres$test=="Peto",4],
        color = simres[simres$test=="Peto",4]) %>%
  layout(title = "Potencia test de Peto", font=list(size = 12),
          scene = list(xaxis=list(title = "Tasa de censura", 
                                  tickfont = list(size = 15),
                                  titlefont = list(size = 20)),
                       yaxis=list(title = "n", 
                                  tickfont = list(size = 15),
                                  titlefont = list(size = 20)),
                       zaxis=list(title = "Potencia", 
                                  tickfont = list(size = 15),
                                  titlefont = list(size = 20))))

distframe <- data.frame(dist=c(rweibull(10000,params[1],params[2]),
                               rweibull(10000,3,200)),
                        Distribution = factor(c(rep("Ovarian",10000),
                                                rep("Setup",10000))))

distframe%>%
  ggplot(aes(x=dist, fill=Distribution)) +
  geom_density(alpha=0.3, cex=0.7)+ 
  labs(x = "Tiempos de falla", y = "Densidad", title = "Densidad de ambas muestras") + 
  theme_light() +
  theme(plot.title = element_text(hjust = 0.5),
        text = element_text(size = 20),
        axis.text.x= element_text(size= 10),
        axis.text.y= element_text(size= 10))

