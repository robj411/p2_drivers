fn <- function(){}
folder <- utils::getSrcDirectory(fn)
setwd(ifelse(folder=='','.',folder))

rnms <- c(paste0('ihr',1:17),paste0('ifr',1:17),
          'ps','Tlat','Tay','Tsr','Tsh','Threc','Thd','Ti','red','R0')
pp <- data.frame(SARS=c(c(0.0578, 0.0578, 0.0578, 0.0578,	
           0.0816, 0.0816, 0.0816, 0.0816,	
           0.3026, 0.3026 ,0.3026, 0.3026,	
           0.8670, 0.8670, 0.8670, 0.8670, 0.6018),
0.867*c(0.017, 0.017, 0.017, 0.017,
                  0.024, 0.024, 0.024, 0.024,
                  0.089, 0.089, 0.089, 0.089, 
                  0.255, 0.255, 0.255, 0.255, 0.177),
c(0.867, 4.6, 2.1, 4.0, 3.75, 26.5-3.75,23.7-3.75,365,0.58,1.7500)),


Flu2009=c(0.669*c(0.00697, 0.00274, 0.00274, 0.00274,
                  0.00274, 0.00561, 0.00561, 0.00561,
                  0.00561, 0.00561, 0.01060, 0.01060,
                  0.01060, 0.01546, 0.01546, 0.01546, 0.01546),
0.669*c(0.000276, 0.00011, 0.00011, 0.00012,
                  0.00012, 0.00030, 0.00030, 0.00030,
                  0.00030, 0.00065, 0.00065, 0.00065,
                  0.00065, 0.00980, 0.00980, 0.00980, 0.00980),
c(0.669, 1.1,2.5, 2.5,2.5, 5.0,5.0,365,0.58, 1.5800)),

Flu1957=c(0.669*13.5*c(0.0001, 0.0001, 0.0001, 0.0001,
                       0.0001, 0.0001, 0.0001, 0.0001,
                       0.0001, 0.0025, 0.0025, 0.0025,
                       0.0025, 0.0200, 0.0200, 0.0200, 0.0200),   
0.669*c(0.0001, 0.0001, 0.0001, 0.0001,
                  0.0001, 0.0001, 0.0001, 0.0001,
                  0.0001, 0.0025, 0.0025, 0.0025,
                  0.0025, 0.0200, 0.0200, 0.0200, 0.0200),
c(0.669,1.1, 2.5, 2.5, 2.5, 5.0, 5.0, 365, 0.58,1.8000)),
Flu1918=c(0.669*8*c(0.02284, 0.00398, 0.00478, 0.00983,
                    0.01700, 0.02922, 0.02470, 0.02205,
                    0.01647, 0.01195, 0.01647, 0.01169,
                    0.03081, 0.04144, 0.04941, 0.04941, 0.04941),
0.669*c(0.02284, 0.00398, 0.00478, 0.00983,0.01700, 0.02922, 0.02470, 0.02205,0.01647, 0.01195, 0.01647, 0.01169,0.03081, 0.04144, 0.04941, 0.04941, 0.04941),
c(0.669, 1.1,2.5,2.5,2.5,5.0,5.0,365,0.58,2.5000)),
CovidWT = c(c(0.000016, 0.000016, 0.000408, 0.000408,0.010400, 0.010400, 0.034300, 0.034300,0.042500, 0.042500, 0.081600, 0.081600, 0.118000, 0.118000, 0.166000, 0.166000, 0.184000),
c(0.000016, 0.000016, 0.000070, 0.000070,0.000309, 0.000309, 0.000844, 0.000844,0.001610, 0.001610, 0.005950, 0.005950,0.019300, 0.019300, 0.042800, 0.042800, 0.078000),
c(0.595, 4.6,2.1,4.0,4.0,12.0,12.0,365, 0.58,2.8700)),
CovidOM = c(
1.85*c(0.000016, 0.000016, 0.000408, 0.000408,0.010400, 0.010400, 0.034300, 0.034300,0.042500, 0.042500, 0.081600, 0.081600,0.118000, 0.118000, 0.166000, 0.166000, 0.184000)*
c(1.10, 1.10, 0.78, 0.78,0.43, 0.43, 0.31, 0.31,0.20, 0.20, 0.14, 0.14,0.14, 0.14, 0.20, 0.20, 0.33),
1.85*c(0.000016, 0.000016, 0.000070, 0.000070,0.000309, 0.000309, 0.000844, 0.000844,0.001610, 0.001610, 0.005950, 0.005950,0.019300, 0.019300, 0.042800, 0.042800, 0.078000)*
  c(1.10, 1.10, 0.78, 0.78,0.43, 0.43, 0.31, 0.31,0.20, 0.20, 0.14, 0.14,0.14, 0.14, 0.20, 0.20, 0.33),
c(0.595,4.0,2.1,4.0,4.0,5.5,5.5,365,0.58,5.9436)),
CovidDE = c(
 1.85*c(0.000016, 0.000016, 0.000408, 0.000408,0.010400, 0.010400, 0.034300, 0.034300,0.042500, 0.042500, 0.081600, 0.081600,0.118000, 0.118000, 0.166000, 0.166000, 0.184000),
 1.85*c(0.000016, 0.000016, 0.000070, 0.000070,0.000309, 0.000309, 0.000844, 0.000844,0.001610, 0.001610, 0.005950, 0.005950,0.019300, 0.019300, 0.042800, 0.042800, 0.078000),
c(0.595,4.0,2.1,4.0,4.0,7.6,7.6,365,0.58,5.0800)))

start_ages <- seq(0,5*16,by=5)
end_ages <- c(start_ages[-1] - 1,'+')
agegroups <- paste(start_ages,end_ages,sep='-')
rownames(pp) <- c(paste0('IHR, ages ',agegroups),
  paste0('IFR, ages ',agegroups),
  'Probability asymptomatic',
  'Latent period','Time to recovery (asymptomatic)','Time to recovery (symptomatic)','Time to hospitalisation','Time to recovery (hospital)','Time to death','Immunity waning','Reduced infectivity asymptomatic','R0')

pairs(t(pp[35:nrow(pp),]))

pp <- rbind(pp,colMeans(pp[1:17,]),colMeans(pp[18:34,]))

tpp <- as.data.frame(t(pp))
colnames(tpp) <- c(rnms,'ihr','ifr')
pairs(tpp[,-c(1:34,42:43)])

write.csv(tpp,'../data/sevenpathogens.csv',row.names = F)


## ihr #######################################
library(brms)
library(pracma)
library(GGally)
library(gtools)
library(reshape2)
library(MASS)

hfr <- tpp[,1:17+17]/tpp[,1:17]

# logitihr <- logit(tpp[,1:17])
rrtpp <- tpp[,1:17]/t(repmat(unlist(tpp[,1]),17,1))
logihr <- log(rrtpp)
x <- (1:length(logihr))*5 - 3

dat <- melt(t(logihr))
dat$x <- rep(1:17,7)
colnames(dat)[colnames(dat)=='Var2'] <- 'variable'

fit8 <- brm(value ~ -1+ variable + s(x,by=variable), data=subset(dat,x>1), chains = 2)
me8 <- brms:::conditional_effects.brmsfit(fit8, ndraws = 200, spaghetti = TRUE)
x11(); plot(me8, ask = FALSE, points = TRUE)
png('../figures/trainingihr.png'); matplot(t(logihr),typ='l',xlab='Age group index',ylab='Log relative risk: IHR'); dev.off()

fit1 <- brm(value ~ -1 + variable + s(x,by=variable,k=4), data=subset(dat,x>1), chains = 2)
me1 <- brms:::conditional_effects.brmsfit(fit1, ndraws = 200, spaghetti = TRUE)
x11(); plot(me1, ask = FALSE, points = F) 

ps <- posterior_samples(fit1)
names(ps)
meltps <- melt(ps)
diseases <- unique(dat$variable)
meltps$disease <- NA
for(d in diseases) meltps$disease[grepl(d,meltps$variable)] <- d
meltps$parameter <- meltps$variable
for(d in diseases) meltps$parameter <- gsub(d,'',meltps$parameter)
meltps$sample <- rep(1:2000)
castps <- dcast(subset(meltps,!is.na(disease)),formula=sample+disease~parameter,value.var='value')
pairsp <- ggpairs(castps,columns=3:6,aes(colour=disease))
png("../figures/pairsplot.png", height = 900, width = 900)
print(pairsp)
dev.off()

stan_data <- standata(fit1)
cov_model_matrix <- stan_data$Xs[1:16,1]
spline_model <- stan_data$Zs_1_1[1:16,]
ps[1,1] + cov_model_matrix*ps[1,8] + spline_model %*% t(ps[1,23:24])
x <- 7
ps[1,x] + cov_model_matrix*ps[1,7+x] + spline_model %*% t(ps[1,22+(2*x-1):(2*x)])
i <- 2

allsamples <- lapply(1:7,function(x)t(sapply(1:nrow(ps),function(i)
  c(0,ps[i,x] + cov_model_matrix*ps[i,7+x] + spline_model %*% t(ps[i,22+(2*x-1):(2*x)])))
))
allsampleslong <- do.call(rbind,allsamples)
allsampleslonglong <- melt(allsampleslong)
allsampleslonglong$age <- rep(1:17,each=2000*7)
allsampleslonglong$disease <- rep(diseases,each=2000)
ggplot(allsampleslonglong) + 
  geom_line(aes(x=Var2,y=value,group=Var1,colour=disease),alpha=1,size=1.5) +
  labs(x='Age group index',y='Log relative risk: IHR',colour='Pathogen') -> pathogenplot
ggsave(pathogenplot,filename='pathogenplot.png',width=9,height=9)
png("../figures/pathogenplot.png", height = 600, width = 600)
print(pathogenplot)
dev.off()

x11(); matplot(t(allsampleslong),typ='l',
               col=rep(c('grey','navyblue','skyblue1','maroon','darkorange1','hotpink','turquoise'),each=2000))


covmat <- cov(castps[,3:7])
parametertoihr <- function(x) c(0,x[1] + cov_model_matrix*x[2] + spline_model %*% (x[3:4]))

nsamples <- 4096
newpoints <- samples <- c()
for(i in 1:nsamples){
  newpoint <- mvrnorm(1,colMeans(castps[,3:7]),covmat)
  newpoint[1] <- min(newpoint[1],1.1*max(castps[,3]))
  funeval <- parametertoihr(newpoint)
  newpoints <- rbind(newpoints,newpoint)
  samples <- rbind(samples,funeval)
}
x11(); matplot(t(samples),typ='l')
write.csv(samples,'ihrrr.csv',row.names = F)


meltnew <- melt(samples[1:100,])
meltnew$sample <- 1:100
ggplot() + 
  geom_line(data=allsampleslonglong,aes(x=Var2,y=value,group=Var1),colour='grey',alpha=1,size=1) +
  geom_line(data=meltnew,aes(x=Var2,y=value,group=sample),colour='navyblue',alpha=.5,size=1) +
  labs(x='Age group index',y='Log relative risk: IHR') -> samplelogplot
png("../figures/samplelogplot.png", height = 600, width = 600)
print(samplelogplot)
dev.off()


## hfr #####################################

hfr <- tpp[,1:17+17]/tpp[,1:17]
hfr <- hfr/repmat(t(t(hfr[,1])),1,17)
loghfr <- log(hfr)
matplot(t(loghfr),typ='l',xlab='Age group index',ylab='Log relative risk');


dat <- melt(t(loghfr))
dat$x <- rep(1:17,7)
colnames(dat)[colnames(dat)=='Var2'] <- 'variable'

fit8 <- brm(value ~ -1+ variable + s(x,by=variable), data=subset(dat,x>1), chains = 2)
me8 <- brms:::conditional_effects.brmsfit(fit8, ndraws = 200, spaghetti = TRUE)
x11(); plot(me8, ask = FALSE, points = TRUE)
png('../figures/traininghfr.png'); matplot(t(loghfr),typ='l',xlab='Age group index',ylab='Log relative risk: HFR'); dev.off()

fit1 <- brm(value ~ -1 + variable + s(x,by=variable,k=4), data=subset(dat,x>1), chains = 2)
me1 <- brms:::conditional_effects.brmsfit(fit1, ndraws = 200, spaghetti = TRUE)
x11(); plot(me1, ask = FALSE, points = F) 

ps <- posterior_samples(fit1)
names(ps)
meltps <- melt(ps)
diseases <- unique(dat$variable)
meltps$disease <- NA
for(d in diseases) meltps$disease[grepl(d,meltps$variable)] <- d
meltps$parameter <- meltps$variable
for(d in diseases) meltps$parameter <- gsub(d,'',meltps$parameter)
meltps$sample <- rep(1:2000)
castps <- dcast(subset(meltps,!is.na(disease)),formula=sample+disease~parameter,value.var='value')
pairsp <- ggpairs(castps,columns=3:6,aes(colour=disease))
png("../figures/pairsplothfr.png", height = 900, width = 900)
print(pairsp)
dev.off()

stan_data <- standata(fit1)
cov_model_matrix <- stan_data$Xs[1:16,1]
spline_model <- stan_data$Zs_1_1[1:16,]
ps[1,1] + cov_model_matrix*ps[1,8] + spline_model %*% t(ps[1,23:24])
x <- 7
ps[1,x] + cov_model_matrix*ps[1,7+x] + spline_model %*% t(ps[1,22+(2*x-1):(2*x)])
i <- 2

allsamples <- lapply(1:7,function(x)t(sapply(1:nrow(ps),function(i)
  c(0,ps[i,x] + cov_model_matrix*ps[i,7+x] + spline_model %*% t(ps[i,22+(2*x-1):(2*x)])))
))
allsampleslong <- do.call(rbind,allsamples)
allsampleslonglong <- melt(allsampleslong)
allsampleslonglong$age <- rep(1:17,each=2000*7)
allsampleslonglong$disease <- rep(diseases,each=2000)
ggplot(allsampleslonglong) + 
  geom_line(aes(x=Var2,y=value,group=Var1,colour=disease),alpha=1,size=1.5) +
  labs(x='Age group index',y='Log relative risk: HFR',colour='Pathogen') -> pathogenplot
png("../figures/pathogenplothfr.png", height = 600, width = 600)
print(pathogenplot)
dev.off()

x11(); matplot(t(allsampleslong),typ='l',
               col=rep(c('grey','navyblue','skyblue1','maroon','darkorange1','hotpink','turquoise'),each=2000))


covmat <- cov(castps[,3:7])
parametertoihr <- function(x)
  c(0,x[1] + cov_model_matrix*x[2] + spline_model %*% (x[3:4]))

newpoints <- c()
for(i in 1:nsamples){
  newpoint <- mvrnorm(1,colMeans(castps[,3:7]),covmat)
  funeval <- parametertoihr(newpoint)
  newpoints <- rbind(newpoints,funeval)
}
x11(); matplot(t(newpoints),typ='l')
write.csv(newpoints,'hfrrr.csv',row.names = F)


meltnew <- melt(newpoints[1:100,])
meltnew$sample <- 1:100
ggplot() + 
  geom_line(data=allsampleslonglong,aes(x=Var2,y=value,group=Var1),colour='grey',alpha=1,size=1) +
  geom_line(data=meltnew,aes(x=Var2,y=value,group=sample),colour='navyblue',alpha=.5,size=1) +
  labs(x='Age group index',y='Log relative risk: HFR') -> samplelogplot
png("../figures/samplelogplothfr.png", height = 600, width = 600)
print(samplelogplot)
dev.off()

