###########################################################################################################
########### Construction of the MAR(1) model 2  ###########################################################
########### Bensebaini Meriem Cyria    ####################################################################
########### coded 11/05/2023    ###########################################################################
###########################################################################################################

rm(list=ls())
RootDir<-'~/Cyria/R/MEDITS_Data/MICE/Git_Statistical_Interaction_Network/'
InputDir<-paste0(RootDir,'Input/')
OutputDir<-paste0(RootDir,'Output/')
ScriptDir<-paste0(RootDir,'Script/')

###
### Load packages 
###

library(readr)
library(MARSS)
library(questionr)
library(tidyverse)

###
### load Functions
###

source(file=paste0(ScriptDir,'MAR_Functions.R'))

###
### Function to Log transform densities
###
scale.fct<-function(x){
  log.x<-log(x)
  log.cent<-log.x-mean(log.x,na.rm=TRUE)
  return(log.cent)
}

###
### load data
###

# Biotic time series
load(paste0(InputDir,"tab_all.Rdata"))
tab_all

# Abiotic time series
Env.param<- read_csv(paste0(InputDir,"Env.param_year_MEDITS.csv"))
Env.param <-Env.param %>% column_to_rownames(v='Med_year')
Env.param
name.env<-c('SST','WeMOI','MLD','Nitrate con.','Phosphate con.','R.flow','Dem.trawlers','Pel.trawlers','Seiner')
colnames(Env.param)<-name.env
Env.param


###
### Second MAR(1): BB. angler, Red mullet, Anchovy, and Nitrate, SST, Demersal trawlers   
###

## Input matrices
data.var<-data.frame(row.names = rownames(tab_all),
  x1=tab_all$Lophius_budegassa_juv,
  x2=tab_all$Mullus_barbatus_adu,
  x3=tab_all$Engraulis_encrasicolus)


data.cov<-data.frame(row.names=rownames(Env.param), 
                     c1= Env.param$`Nitrate con.`,
                     c2= Env.param$SST,
                     c3= Env.param$Dem.trawlers
                     )
str.y<-2001
f.y<-2021

data.var<-data.var %>% filter(rownames(data.var)>=str.y & rownames(data.var)<=f.y )
data.var
data.cov<-data.cov %>% filter(rownames(data.cov)>=str.y & rownames(data.cov)<=f.y )
data.cov

## standardize data
data.var%>% mutate_each(funs(scale.fct))->newdata.var
data.cov%>% mutate_each(funs(scale))->newdata.cov 

## matrice X : transpose data
x.matrix<-t(newdata.var)
c.matrix<-t(newdata.cov) 
## set n and p
n<-3 # n: is the number of series input into the matrix y
p<-3# m: number of series input in matrix of covariates C

### list of parameters for a MAR model - the library fit state-space models,
### but we are directly fitting the state equation from the data - no need for
### the space counterpart, which is why many params (Z,A,R) are set to unity or 0
listparam<-list(Z=diag(n),                         ## conversion matrice from space to state
                A="zero", ## space eq. intercept
                R="zero", ## space eq. var-cov matrix
                c=c.matrix, ## covariate
                U="zero", ## state eq. intercept
                Q="diagonal and unequal", ## state eq. var-cov matrix
                x0="zero",## mean of the distrib from which starting values are drawn
                V0=diag(n),
                tinitx=1)      ## var-cov matrix of the distrib from which starting values are drawn

## Specify the "Full model" matrix for the MAR community matrix

B.full=matrix(list('b11','b12',0,
                  0,'b22','b23',
                  0,0,'b33'),
             nrow=n,ncol=n,byrow=T)

### Specify the "Full model" vector for the env effect in the MAR model

C.full=matrix(list('c11',0,'c13',
                  0,'c22','c23',
                  0,0,0),nrow=n,ncol=p,byrow=T)


## The corresponding "Null model" matrix and vector, if needed
Bzero=matrix(0,nrow=n,ncol=n,byrow=T)
Czero=matrix(0,nrow=n,ncol=p,byrow=T)

# model 
model.gen=list(Z=listparam$Z,
               A=listparam$A,
               R=listparam$R,
               B=B.full,
               C=C.full,
               c=listparam$c,
               U=listparam$U,
               Q=listparam$Q,
               x0=listparam$x0,
               V0=listparam$V0)
m.try<-MARSS(x.matrix,model=model.gen,control=list(conv.test.slope.tol=0.05))
m.ci<-MARSSparamCIs(m.try,method="hessian",alpha=0.05,nboot=1000)
result<-getpval(m.ci)
result

## Save results
tab.result.MAR<-data.frame(coef=names(m.try$coef)[c(1:5,9:12)],Est=result$est[c(1:5,9:12)])
write_csv(tab.result.MAR,file=paste0(OutputDir,'Coef_MAR_Model_2.csv'))


res_without_newCov<-data.frame(coef=rownames(result)[c(1:5,9:12)],Est_without_newCov=result$est[c(1:5,9:12)])

#######################################################################################
####################################################################################### see the results with the new "shift" covariate

###
### create the "shift" covariate
###
new_cov<-data.frame(new_cov=c(rep(-2,13),-1.2,-0.4,0.4,1.2,rep(2,12)))
row.names(new_cov)<-1993:2021
plot(row.names(new_cov),new_cov$new_cov,type='l')

# Add the new covariate
Env.param<- cbind(Env.param,new_cov)

###
### Second MAR(1): BB. angler, Red mullet, Anchovy, and Nitrate, SST, Demersal trawlers   
###

## Input matrices
data.var<-data.frame(row.names = rownames(tab_all),
                     x1=tab_all$Lophius_budegassa_juv,
                     x2=tab_all$Mullus_barbatus_adu,
                     x3=tab_all$Engraulis_encrasicolus)


data.cov<-data.frame(row.names=rownames(Env.param), 
                     c1= Env.param$`Nitrate con.`,
                     c2= Env.param$SST,
                     c3= Env.param$Dem.trawlers,
                     c4= Env.param$new_cov)
str.y<-2001
f.y<-2021

data.var<-data.var %>% filter(rownames(data.var)>=str.y & rownames(data.var)<=f.y )
data.var
data.cov<-data.cov %>% filter(rownames(data.cov)>=str.y & rownames(data.cov)<=f.y )
data.cov

## standardize data
data.var%>% mutate_each(funs(scale.fct))->newdata.var
data.cov%>% mutate(across(colnames(data.cov)[1:3],scale))->newdata.cov

## matrice X : transpose data
x.matrix<-t(newdata.var)
c.matrix<-t(newdata.cov) 
## set n and p
n<-3 # n: is the number of series input into the matrix y
p<-4# m: number of series input in matrix of covariates C

### list of parameters for a MAR model - the library fit state-space models,
### but we are directly fitting the state equation from the data - no need for
### the space counterpart, which is why many params (Z,A,R) are set to unity or 0
listparam<-list(Z=diag(n),                         ## conversion matrice from space to state
                A="zero", ## space eq. intercept
                R="zero", ## space eq. var-cov matrix
                c=c.matrix, ## covariate
                U="zero", ## state eq. intercept
                Q="diagonal and unequal", ## state eq. var-cov matrix
                x0="zero",## mean of the distrib from which starting values are drawn
                V0=diag(n),
                tinitx=1)      ## var-cov matrix of the distrib from which starting values are drawn

## Specify the "Full model" matrix for the MAR community matrix
B.full=matrix(list('b11','b12',0,
                       0,'b22','b23',
                       0,0,'b33'),
                  nrow=n,ncol=n,byrow=T)

## Specify the "Full model" vector for the env effect in the MAR model
C.full=matrix(list('c11',0,'c13',0,
                       0,'c22','c23','c24',
                       0,0,0,'c34'),nrow=n,ncol=p,byrow=T)
# C.no.const=matrix(list('c11',0,'c13','c14',
#                        0,'c22','c23','c24',
#                        0,0,0,'c34'),nrow=n,ncol=p,byrow=T)

## The corresponding "Null model" matrix and vector, if needed
Bzero=matrix(0,nrow=n,ncol=n,byrow=T)
Czero=matrix(0,nrow=n,ncol=p,byrow=T)

# model 
model.gen=list(Z=listparam$Z,
               A=listparam$A,
               R=listparam$R,
               B=B.full,
               C=C.full,
               c=listparam$c,
               U=listparam$U,
               Q=listparam$Q,
               x0=listparam$x0,
               V0=listparam$V0)
m.try<-MARSS(x.matrix,model=model.gen,control=list(conv.test.slope.tol=0.05))
m.ci<-MARSSparamCIs(m.try,method="hessian",alpha=0.05,nboot=1000)
result<-getpval(m.ci)
result

res_with_newCov<-data.frame(coef=rownames(result)[c(1:5,9:14)],Est_with_newCov=result$est[c(1:5,9:14)])
#res_with_newCov<-data.frame(coef=rownames(result)[c(1:5,9:15)],Est_with_newCov=result$est[c(1:5,9:15))


Table_comp<-full_join(res_without_newCov, res_with_newCov, by = "coef") 
Table_comp<- Table_comp %>% mutate(delta=round(Est_without_newCov-Est_with_newCov,3))
Table_comp

write_csv(Table_comp, paste0(OutputDir,'Tab_comp_ShiftCov_MAR2.csv'))


