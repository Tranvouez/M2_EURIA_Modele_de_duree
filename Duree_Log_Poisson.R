# Le modèle log-Poisson  
##########################

# Human Mortality Database  (HMD)
# http://www.mortality.org

setwd("~/Documents/Prog_R")

# Décès
de=read.csv("DeathsFrance1.csv",header = TRUE, sep = ";")
# remarque : la classe d'âge "110" est en réalité "110 et plus".

# Expositions
ex=read.csv("ExposuresFrance1.csv",header = TRUE, sep = ";")

# 
##########################
# Le modèle log-Poisson  
# Il s'agit d'un modèle que l'on peut calibrer avec la fonction gnm.
# (generalized nonlinear model)
library(gnm)

ind=which((de$Age>44)&(de$Age<100)&(de$Year>1949)&(de$Year<2013))
annee=1950:2012 ; nc=length(annee)
age=45:99 ; nl=length(age)

D=de$Male[ind]
E=ex$Male[ind]
x=as.factor(ex$Age[ind])
t=as.factor(ex$Year[ind])

x=ex$Age[ind]
t=ex$Year[ind]
# x et t sont des vecteurs de taille nl*nc

regp <- gnm(D~0+as.factor(x)+Mult(as.factor(x),as.factor(t)), offset=log(E),
            family=poisson(link="log"))

nomvar=names(regp$coefficients)
nomvar # permet de localiser les paramètres

plot(45:99,regp$coefficients[1:55], main="coefficients alpha_x",xlab="Ages")
plot(45:99,regp$coefficients[56:110], main="coefficients beta_x",xlab="Ages")
plot(1950:2012,regp$coefficients[111:173], main="coefficients k_t",xlab="Années")

alpha=regp$coefficients[1:55]
k=regp$coefficients[111:173]
beta=regp$coefficients[56:110]

# On "normalise" les paramètres comme pour Lee-Carter :
sb=sum(beta)
beta=beta/sb
mk=mean(k)
k=(k-mk)*sb
alpha=alpha+beta*mk

plot(45:99,alpha, main="coefficients alpha_x",xlab="Ages")
plot(45:99,beta, main="coefficients beta_x",xlab="Ages")
plot(1950:2012,k, main="coefficients k_t",xlab="Années")

# calcul des log(mu_{x,t})
logmu=matrix(NA,nrow=55,ncol=63)
for (i  in 1:55)
{for (j in 1:63)
{logmu[i,j]=alpha[i]+beta[i]*k[j]}}

# Comparaison avec Lee Carter classique
library(demography)
muh=matrix( de$Male[ind]/ ex$Male[ind],nl,nc)
poph=matrix(ex$Male[ind],nl,nc)
Baseh=demogdata(data=muh,pop=poph,ages=age,years=annee,type="mortality",label='France',name='Hommes',lambda=1)
lch=lca(Baseh)
plot(lch)

# Comparaison paramètres Log Poisson / Lee Carter
plot(45:99,alpha, main="coefficients alpha_x",xlab="Ages",type='l',col='blue')
lines(45:99, lch$ax,col='red')
legend(50, -1, legend=c("Log Poisson", "Lee Carter"),col=c("blue","red"), lty=1, cex=0.8)

plot(45:99,beta, main="coefficients beta_x",xlab="Ages",type='l',col='blue')
lines(45:99, lch$bx,col='red')
legend(50, 0.01, legend=c("Log Poisson", "Lee Carter"),col=c("blue","red"), lty=1, cex=0.8)

plot(1950:2012,k, main="coefficients k_t",xlab="Années",type='l',col='blue')
lines(1950:2012,lch$kt,col='red')
legend(1955, -10, legend=c("Log Poisson", "Lee Carter"),col=c("blue","red"), lty=1, cex=0.8)

predh=lch$fitted$y # c'est log(mu_{x,t}) qui est prédit

# comparaison Lee-Carter / données réelles
rmseh=sqrt(sum((log(muh)-(predh))^2)/(nl-1)/nc)
rmseh

# comparaison Log Poisson / données réelles
rmseh1=sqrt(sum((logmu-(log(muh)))^2)/(nl-1)/nc)
rmseh1

# comparaison Log Poisson / Lee-Carter
rmseh2=sqrt(sum((logmu-predh)^2)/(nl-1)/nc)
rmseh2

# comparaison graphique \mu_{x,t} Lee-Carter / Log Poisson
# pour une année donnée
plot(45:99,exp(logmu[,63]),main="mu_{x,2012}",type='l',xlab="Ages",col='blue')
lines(45:99,exp(predh[,63]), col='red')
legend(50, 0.3, legend=c("Log Poisson", "Lee Carter"),col=c("blue","red"), lty=1, cex=0.5)

plot(45:99,logmu[,63],main="log(mu_{x,2012})",type='l',xlab="Ages",col='blue')
lines(45:99,predh[,63], col='red')
lines(45:99,log(muh[,63]),col='green')
legend(50, -1, legend=c("Log Poisson", "Lee Carter", "obs."),col=c("blue","red", "green"), lty=1, cex=0.5)

      