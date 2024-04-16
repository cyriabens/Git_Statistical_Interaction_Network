###########################################################################################################
########### Pairwise Granger causality test : Biotic interactions  ########################################
########### Bensebaini Meriem Cyria    ####################################################################
########### coded 06/04/2021    ###########################################################################
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
load(paste0(InputDir,"tab_all.Rdata"))
Tab_densties<-tab_all
dim(Tab_densties)
Tab_densties

###
### See the time series
###
par(mfrow=c(2,2))
for(sp in 1:dim(Tab_densties)[2]){
  plot(as.numeric(rownames(Tab_densties)),Tab_densties[,sp],type='b',xlab='Year',ylab='',main=colnames(Tab_densties)[sp])
}


###
### Pairwise Granger test 
###

Bio_inter  <-data.frame(X=NULL,
                               Y=NULL,
                               Leng_T.series=NULL,
                               Effect_sens=NULL,
                               lag=NULL,
                               i=NULL,
                               j=NULL,
                               P_val_WT=NULL,
                               Effect_Size=NULL,
                               Coef_b11=NULL,
                               Coef_b21=NULL,
                               Coef_a1=NULL)


### Choose the targeted species. 
Target_sp<-c(1,4,5,11,14,15,16,22,23)
colnames(tab_all[Target_sp])

## A loop for x
for (j in 1: length(Target_sp)) {
  y<-tab_all[Target_sp[j]] %>% rownames_to_column(var='year')
  
  ## A loop for y
  
  for (i in 1:length(colnames(tab_all))) {
    
    if(i!=Target_sp[j]){
      x<-tab_all[i] %>% rownames_to_column(var='year') 
      
      
      ## standardize data
      df<-full_join(x,y,by=c("year")) 
      df %>%
        filter(!is.na(df[,2]) & !is.na(df[,3])) %>%
        column_to_rownames(var='year') %>%
        mutate(across(colnames(df)[2:3],scale.fct)) -> newdf
      names(newdf)=c("x","y")
      res<- Gran.test.fct(newdf,lag_order = 1)
      ## retrieval and storage of results    (effect of columns on rows, and only significant results)
      Bio_inter<- rbind(Bio_inter, data.frame(X=colnames(tab_all)[i],
                                                                Y=colnames(tab_all)[Target_sp[j]],
                                                                Leng_T.series=dim(newdf)[1],
                                                                Effect_sens='x_cause_y',
                                                                lag=res[1],
                                                                i=i,
                                                                j=Target_sp[j],
                                                                P_val_WT=round(res[2],3),
                                                                Effect_Size=round(res[3],2),
                                                                Coef_b11=round(res[4],2),
                                                                Coef_b21=round(res[5],2),
                                                                Coef_a1=round(res[6],2)))
      
    }
    
  }
}

# save the data frame
Bio_inter
write.csv(Bio_inter,paste0(OutputDir,"Biotic_interactions.csv"), row.names = FALSE)


###
### Results analysis
###

### Choose the G_threshold and alpha for selection
G_th<- 0.28 # G_th ≥ln(1/0.75)≈0.28 : explains at least 25% of the variance
alpha<- 0.1 # Deduced from the power test (see the Power_test code)

Biotic_interactions<-Bio_inter

Biotic_interactions %>% arrange (desc(Effect_Size), P_val_WT)

Biotic_interactions <-Biotic_interactions %>% rename(G_obs=Effect_Size, p_value=P_val_WT)

myfile <- paste0(OutputDir,"P_vs_G_BioInter.png")
png(file = myfile, width =85, height =85 ,bg = "White",
    units="mm",res=500)
ggplot(Biotic_interactions, aes(G_obs, p_value))+
  geom_point()+
  geom_vline(xintercept=G_th)+
  geom_hline(yintercept=alpha)+
  theme_bw()
dev.off()

### 1st filter according to the G_obs and p-value : interactions with big effect and significant
dim(Biotic_interactions)[1]
fil.1<-Biotic_interactions %>% filter(G_obs>=G_th & p_value<=alpha) %>% arrange(desc(G_obs)) %>%
  select(X,Y,G_obs,p_value,Coef_b11)
dim(fil.1)[1]
fil.1

### 2nd filter according to the G_obs and p-value : interactions with a weak effect but significant 
dim(Biotic_interactions)[1]
fil.2<-Biotic_interactions %>%  filter(G_obs<G_th & p_value<=alpha) %>% arrange(desc(G_obs), abs(Coef_b11))%>%
  select(X,Y,G_obs,p_value,Coef_b11)
dim(fil.2)[1]
fil.2

### 3d filter according to the G_obs and p-value : interactions with a weak effect and non significant 
dim(Biotic_interactions)[1]
fil.3<-Biotic_interactions %>% filter(G_obs<G_th & p_value>alpha) %>% arrange(desc(G_obs),abs(Coef_b11)) %>%
  select(X,Y,G_obs,p_value,Coef_b11)
dim(fil.3)[1]
fil.3

# Save selected interactions
write.csv(fil.1,paste0(OutputDir,"/Selected_Biotic_Interactions.csv"), row.names = FALSE)
write.csv(fil.2,paste0(OutputDir,"/Fil.2_Biotic_Interaction.csv"), row.names = FALSE)
write.csv(fil.3,paste0(OutputDir,"/Fil.3_Biotic_Interaction.csv"), row.names = FALSE)

