###########################################################################################################
########### Pairwise Granger causality test : Abiotic drivers  ############################################
########### Bensebaini Meriem Cyria    ####################################################################
########### coded 03/11/2021    ###########################################################################
###########################################################################################################
rm(list=ls())
RootDir<-'~/Cyria/R/MEDITS_Data/MICE/Git_Statistical_Interaction_Network/'
InputDir<-paste0(RootDir,'Input/')
OutputDir<-paste0(RootDir,'Output/')

###
### Load packages 
###

library(readr)
library(vars)
library(tidyverse)


###
### Granger function
###
Gran.test.fct<-function(data,lag_order=1){
  newdf<-data
  ## Using VAR function with Schwarz criterion as a yardstick to estimate simultaneously model order 
  varpp<-VAR(y=newdf, type="none",ic="SC",p=lag_order)
  
  ## For x cause y
  
  # retrieve the p-value using grangertest function (Wald p-value) 
  gxy = grangertest(newdf$x,newdf$y,order = lag_order) # x causes y 
  Pval_xcausey=gxy$`Pr(>F)`[2]  # if p_v < alpha we reject H0, alpha = 5 % and H0: x do not Granger-cause y
  
  ## Let's compute log ratio (Effect size)
  ar_y=ar(newdf$y,order=lag_order,aic=F,method="ols") 
  log_xy_inter=log(sum((ar_y$resid)^2,na.rm=T)/sum((varpp$varresult$y$residuals)^2,na.rm=T))# For x granger cause y
  
  ## Let's extract parameters from x equation
  ar_x=ar(newdf$x,order=lag_order,aic=F,method="ols")
  
  ## retrieval and storage of results    
  result<- c(lag_order,
             round(Pval_xcausey,3),# p-value
             round(log_xy_inter,2),# G
             round(varpp$varresult$y$coefficients[1],2),# b11 : x->y
             round(varpp$varresult$y$coefficients[2],2),# b21 : y->y
             round(ar_y$ar[,,1],2)) # a1: y->y without interaction
  
  return(result)}


###
### Function to Log transform densities
###
scale.fct<-function(x){
  log.x<-log(x)
  log.cent<-log.x-mean(log.x)
  return(log.cent)
}

###
### load data
###

# Biotic time series
load(paste0(InputDir,"tab_all.Rdata"))
tab_all
dim(tab_all)

Selected_Biotic_Interactions<- read_csv(paste0(OutputDir,"Selected_Biotic_Interactions.csv"))

# Abiotic time series 
Env.param<- read_csv(paste0(InputDir,"Env.param_year_MEDITS.csv"))
Env.param <-Env.param %>% column_to_rownames(v='Med_year')
Env.param
name.env<-c('SST','WeMOI','MLD','Nitrate con.','Phosphate con.','R.flow','Dem.trawlers','Pel.trawlers','Seiner')
colnames(Env.param)<-name.env
Env.param

###
### Plot data
###
myfile<-paste(OutputDir,"Abiotic_TimeSeries",".jpg",sep='')
jpeg(filename=myfile,units="mm",res=500,width=170,height=200,bg="white")
par(mfrow=c(5,2),mar=c(2,3.1,1.5,0.5),omi=c(0,0,0,0))
for(v in 1: length(colnames(Env.param))){
  ts<-data.frame(var=Env.param[,v],
                 year=rownames(Env.param))
  ts<-ts%>% drop_na()
  plot(ts$year,ts$var,type='b',pch=19,xlab='',ylab='')
  title(main=name.env[v], cex.main=1.1, line=0.2)
}
dev.off()

###
### Test the correlation between covariates
###
Env.param.cor.test<-drop_na(Env.param) 
Env.param.cor.test<- Env.param.cor.test %>% mutate(across(colnames(Env.param.cor.test),scale))
cor_test<-data.frame(cor(Env.param.cor.test))
cor_test
write_csv(cor_test,file=paste0(OutputDir,'Correlation_Test_AbioticDrivers.csv'))

### Check only if >= 0.7 
Pos.Covariate.1<-7
Pos.Covariate.2<-8
par(mfrow=c(1,1))
plot(rownames(Env.param.cor.test), Env.param.cor.test[,Pos.Covariate.1], type='l', ylim=c(-2,3), col='red')
lines(rownames(Env.param.cor.test), Env.param.cor.test[,Pos.Covariate.2], col='green')
legend(2010, 3, legend=c(colnames(Env.param.cor.test)[Pos.Covariate.1],
                           colnames(Env.param.cor.test)[Pos.Covariate.2]),
                           col=c("red", "green"), lty=1, cex=0.8)

###
### Preparing lists for conditions 
###

# Which species are concerned by pairwise GC tests  
List.sp.ToTest<-unique(c(Selected_Biotic_Interactions$X,Selected_Biotic_Interactions$Y))

# Lists for conditions

list.sp.dem<-c('Merluccius_merluccius_juv','Merluccius_merluccius_adu','Mullus_barbatus_adu',
               'Eledone_cirrhosa_juv','Eledone_cirrhosa_adu','Lophius_budegassa_juv','Lophius_budegassa_adu',
               'Zeus_faber','Eutrigla_gurnardus_adu','Trisopterus_capelanus_juv') # List of demersal species
list.sp.sard<-c('Sardina_pilchardus')
list.sp.anch<-c('Engraulis_encrasicolus')
list.FE.pel.anch<-c('Pel.trawlers')
list.FE.pel.sar<-c('Seiner')
list.FE.dem<-c('Dem.trawlers')
list.env<-c('SST','WeMOI','MLD','R.flow','Nitrate con.','Phosphate con.')


###
### Pairwise Granger test 
###

n.cov<-dim(Env.param)[2]
n.lag<-c(0:1)

## storage object
results<-data.frame()


## Loop 

# Set the variable
for (v in 1:length(List.sp.ToTest)){
  
var<-colnames(tab_all[List.sp.ToTest[v]])

tab.1<- data.frame(year=as.numeric(rownames(tab_all[List.sp.ToTest[v]])),
                var=tab_all[List.sp.ToTest[v]][,1],
                row.names = NULL) %>% drop_na()

# Set the covariate
for (c in 1:dim(Env.param)[2]){
  cov<- colnames(Env.param[c])
  tab.2<-data.frame(year=as.numeric(rownames(Env.param[c])),
                    cov=Env.param[c][,1],
                    row.names = NULL)  %>% drop_na()

  start.cov<-min(as.numeric(tab.2$year))
  end.cov<-max(as.numeric(tab.2$year))
  start.var<-min(as.numeric(tab.1$year))
  end.var<-max(as.numeric(tab.1$year))
  
  ## Reorganize time series so that time lag is set to 0
  
  if(start.cov==start.var){
    tab<-data.frame(var=tab.1$var[tab.1$year>=start.var & tab.1$year<=(end.cov-1)],
                    cov=tab.2$cov[tab.2$year>=(start.cov+1)])
    rownames(tab)<-c((start.var+1):(end.cov))
    tab$var<-scale.fct(tab$var)
    tab$cov<-scale(tab$cov)
  } 
  
      if(start.cov<start.var){
        tab<-data.frame(var=tab.1$var[tab.1$year>=start.var & tab.1$year<=(end.cov-1)],
                        cov=tab.2$cov[tab.2$year>=(start.var+1)])
        rownames(tab)<-c((start.var+1):(end.cov))
        tab$var<-scale.fct(tab$var)
        tab$cov<-scale(tab$cov)
      } 
  
    if(start.cov>start.var){
      tab<-data.frame(var=tab.1$var[tab.1$year>=(start.cov-1) & tab.1$year<=(end.cov-1)],
                      cov=tab.2$cov[tab.2$year>=start.cov])
      rownames(tab)<-c(start.cov:end.cov)
      tab$var<-scale.fct(tab$var)
      tab$cov<-scale(tab$cov)
    }
 
  
  ## Conditions : we want to keep only the tests we're intersted in
  
  ## all environmental drivers are tested for all species 
  cond1<-var %in% list.sp.dem &  (cov %in% list.FE.dem | cov %in% list.env) # The effect of demersal trawling is only tested on demersal species
  cond2<-var %in% list.sp.sard &  (cov %in% list.FE.pel.sar | cov %in% list.env) # Testing the effect of seiners on sardine
  cond3<-var %in% list.sp.anch &  (cov %in% list.FE.pel.anch |cov %in% list.env) # Testing the effect of pelagic trawling on anchovy
  
  if (cond1 |cond2| cond3) {
    
    ## rename the variables
    tab<-tab%>% rename(x=cov,y=var) %>% dplyr::select(x,y)
    ## Test 
    res<-Gran.test.fct(data=tab,lag_order = 1)
    
    results<-rbind(results,data.frame(x=cov,
                              y=var,
                              Leng_series=dim(tab)[1],
                              P_val_WT=round(res[2],3),
                              Effect_Size=round(res[3],2),
                              Coef_b11=round(res[4],2),
                              Coef_b21=round(res[5],2),
                              Coef_a1=round(res[6],2)))
    
  }
  }}
results
dim(results)
write_csv(results,file=paste0(OutputDir,'Abiotic_drivers.csv'))

###
### Results analysis
###

### Choose the G_threshold and alpha for selection
G_th<- 0.28 # G_th ≥ln(1/0.75)≈0.28 : explains at least 25% of the variance
alpha<- 0.1 # Deduced from the power test (see the Power_test code)


covar_tab <-results %>% rename(G_obs=Effect_Size, p_value=P_val_WT) 

myfile <- paste0(OutputDir,"P_vs_G_AbioEffect.png")
png(file = myfile, width =85, height =85 ,bg = "White",
    units="mm",res=500)
ggplot(covar_tab, aes(G_obs, p_value))+
  geom_point()+
  geom_vline(xintercept=0.28)+
  geom_hline(yintercept=0.1)+
  theme_bw()
dev.off()

### 1st filter according to the G_obs and p-value : interactions with big effect and significant
dim(covar_tab)[1]
fil.1<-covar_tab %>% arrange(desc(G_obs)) %>% filter(G_obs>=G_th & p_value<=alpha )  %>%
  dplyr::select(x,y,G_obs,p_value,Coef_b11)
fil.1
dim(fil.1)[1]

write_csv(fil.1,file=paste0(OutputDir,'Selected_Abiotic_Effects.csv'))

### 2nd filter according to the G_obs and p-value : interactions with a weak effect but significant 
fil.2<-covar_tab %>% arrange(desc(G_obs)) %>% filter(G_obs<G_th & p_value<=alpha )  %>%
  dplyr::select(x,y,G_obs,p_value,Coef_b11)
dim(fil.2)[1]
fil.2
write_csv(fil.2,file=paste0(OutputDir,'Fil.2_Abiotic_Effects.csv'))

### 3d filter according to the G_obs and p-value : interactions with a weak effect and non significant 
fil.3<-covar_tab %>% arrange(desc(G_obs)) %>% filter(G_obs<G_th & p_value>alpha)  %>%
  dplyr::select(x,y,G_obs,p_value,Coef_b11)
dim(fil.3)[1]
fil.3
write_csv(fil.3,file=paste0(OutputDir,'Fil.3_Abiotic_Effects.csv'))



