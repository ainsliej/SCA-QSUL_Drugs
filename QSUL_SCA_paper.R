# My SCA analysis for the QSUL data PAPER VERSION
# Ainslie Johnstone
# 27.05.20

# Load packages + data ------
setwd('~/OneDrive - University College London/3_QSULoutcomes/')
library(tidyverse)
library(broom)
library(cowplot)
library(MuMIn)
library(ddpcr)
datatable<- read.table("400_QSUL_all_drugs.csv", header=T, sep=",")
datatable<-datatable[-c(1,2),]
#datatable$X<-1
datatable$NACA_antiepileptic<-factor(datatable$NACA_antiepileptic)
datatable$Antidepressants<-factor(datatable$Antidepressants)
datatable$GABA_agonists<-factor(datatable$GABA_agonists)
datatable$Gender<-factor(datatable$Gender)
datatable$Affected.limb<-factor(datatable$Affected.limb)
datatable$DominantUL<-factor(datatable$DominantUL)
datatable$NFI<-as.numeric(as.character(datatable$NFI))
datatable<-datatable[!is.na(datatable$Age),]

datatable$DomAffected<-factor(ifelse(datatable$Affected.limb==datatable$DominantUL,"Y", "N") )
datatable$time_since_stroke[is.na(datatable$time_since_stroke)]<-mean(datatable$time_since_stroke, na.rm=TRUE)

datatable$FM1_tot_pc<-datatable$FM1_tot/54
datatable$ARAT1_tot_pc<-datatable$ARAT1_tot/57
datatable$CAHAI1_tot_pc<-datatable$CAHAI1_tot/91
datatable$FM2_tot_pc<-datatable$FM2_tot/54
datatable$ARAT2_tot_pc<-datatable$ARAT2_tot/57
datatable$CAHAI2_tot_pc<-datatable$CAHAI2_tot/91
datatable$FM3_tot_pc<-datatable$FM3_tot/54
datatable$ARAT3_tot_pc<-datatable$ARAT3_tot/57
datatable$CAHAI3_tot_pc<-datatable$CAHAI3_tot/91
datatable$FM4_tot_pc<-datatable$FM4_tot/54
datatable$ARAT4_tot_pc<-datatable$ARAT4_tot/57
datatable$CAHAI4_tot_pc<-datatable$CAHAI4_tot/91

datatable$FM1.FM4_pc<-datatable$FM1.FM4/54
datatable$AR1.AR4_pc<-datatable$AR1.AR4/57
datatable$CA1.CA4_pc<-datatable$CA1.CA4/91

datatable$relFM1.FM4_pc<-datatable$FM1.FM4_pc/(1-datatable$FM1_tot_pc)
datatable$relAR1.AR4_pc<-datatable$AR1.AR4_pc/(1-datatable$ARAT1_tot_pc)
datatable$relCA1.CA4_pc<-datatable$CA1.CA4_pc/(1-datatable$CAHAI1_tot_pc)


alldata<- datatable[!is.na(datatable$Benzos)&!is.na(datatable$HADS_T)&!is.na(datatable$NFI)&
                      !is.na(datatable$FM1_tot)&!is.na(datatable$FM2_tot),]

alldata$time_since_stroke[is.na(alldata$time_since_stroke)]<-mean(alldata$time_since_stroke, na.rm=TRUE)
alldata23<-alldata[!is.na(alldata$FM2.FM3)&!is.na(alldata$AR2.AR3),]
alldata4<-alldata[!is.na(alldata$FM4_tot)&!is.na(alldata$AR1.AR4)&!is.na(alldata$relAR1.AR4_pc)&is.finite(alldata$relAR1.AR4_pc),]

alldata4$HighHADS<-factor(ifelse(alldata4$HADS_T>12,1,0))
alldata4$DepressionGroups<-factor(ifelse(alldata4$Antidepressants==1,'OnAD',
                                         ifelse(alldata4$HADS_T>=12, 'OffAD_High','OffAD_Low')))


#alldata4$Antidepressants<-factor(alldata4$SSRI) Testing something for reviewers

#### ADMISSION FUNCTIONS #####   -------
# Create models through loop 
AdmissionSCA<- function(drug_var_n){
#set na.action for dredge
options(na.action = "na.fail")
  
  if (drug_var_n=="GABA_agonists"){
    drug_covA_n="Antidepressants"
    drug_covB_n="NACA_antiepileptic"
   } else if (drug_var_n=="Antidepressants"){
    drug_covA_n="GABA_agonists"
    drug_covB_n="NACA_antiepileptic"
   } else if (drug_var_n=="NACA_antiepileptic"){
     drug_covA_n="GABA_agonists"
     drug_covB_n="Antidepressants"}
    

# run full model
for (meas in 1:3){
  for (outliers in 1:2){
    if (outliers==1){this_data<-alldata4
    }else{
      this_data<- filter(alldata4, !(abs(alldata4$FM1.FM4_pc- median(alldata4$FM1.FM4_pc)) > 2.5*IQR(alldata4$FM1.FM4_pc))&
                           !(abs(alldata4$AR1.AR4_pc- median(alldata4$AR1.AR4_pc)) > 2.5*IQR(alldata4$AR1.AR4_pc))&
                           !(abs(alldata4$CA1.CA4_pc- median(alldata4$CA1.CA4_pc)) > 2.5*IQR(alldata4$CA1.CA4_pc)))}
    
    if (meas==1){ y.var<-this_data$FM1_tot_pc}
    else if (meas==2){y.var<-this_data$ARAT1_tot_pc}
    else {y.var<-this_data$CAHAI1_tot_pc}
    
    
    this_data<-cbind(this_data,y.var,drug_covA=get(drug_covA_n,this_data),
                     drug_covB=get(drug_covB_n,this_data),drug_var=get(drug_var_n,this_data))
    full.model=lm(y.var ~ Age+ time_since_stroke+ DomAffected+  Gender + NFI+ 
                    HADS_T + drug_var + drug_covA + drug_covB, data= this_data)    
    
    models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                             subset=(((Age&Gender)|(!Age&!Gender))&((NFI&HADS_T)|(!NFI&!HADS_T))&
                                       ((drug_covA&drug_covB)|(!drug_covA&!drug_covB))&
                                       ((time_since_stroke&DomAffected)|(!time_since_stroke&!DomAffected))))
    
    
    #run all possible nested models + filtering for only ones I want
    if (meas==1){
      FMons<- as.numeric(rep(1, times=nrow(models.sca)))
      ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
      CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
    }else if (meas==2){
      FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
      ARons<- as.numeric(rep(1, times=nrow(models.sca)))
      CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
    }else {   
      FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
      ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
      CAons<- as.numeric(rep(1, times=nrow(models.sca)))
    }
    
    if(outliers==2){
      OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
    }  else {
      OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
    }
    
    
    
    # extract parameter estimate
    (model_params = MuMIn::get.models(models.sca, subset = TRUE) %>%
        tibble() %>%
        dplyr::rename("model" = ".") %>%
        mutate(tidied = purrr::map(model, broom::tidy),
               model_num = row_number()) %>%
        select(model_num, tidied) %>%
        unnest() %>%
        select(model_num, term, estimate) %>%
        spread(term, estimate)) %>%
      select(-starts_with("sd"))
    
    model_params$model_num<-model_params$model_num+(nrow(models.sca)*(((meas-1)*2)+(outliers-1)))
    model_params<-cbind(model_params,FMons, ARons, CAons, OutsOut)
    
    (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
        tibble() %>%
        dplyr::rename("model" = ".") %>%
        mutate(tidied = purrr::map(model, broom::tidy),
               model_num = row_number()) %>%
        select(model_num, tidied) %>%
        unnest() %>%
        filter(term ==  paste("drug_var1", sep="")) %>%
        ungroup() %>%
        select(model_num, estimate, std.error, p.value))
    
    model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*2)+(outliers-1)))
    
    if (outliers==1&meas==1){
      model_params_all<-model_params
      model_ps_all<-model_ps
    }else{
      model_params_all<-full_join(model_params_all, model_params)
      model_ps_all<-full_join(model_ps_all, model_ps)
    }
  }}


print(test_sig<-(model_ps_all$p.value<=0.05) %>%
        sum())
print(summary(model_ps_all$estimate[model_ps_all$p.value<=0.05]))
print(summary(model_ps_all$estimate))
return(list(model_ps_all, model_params_all, test_sig))
}
# Plotting the outcomes 
plotAdmissionSCA<-function(drug_var_n){
  
  if (drug_var_n=="GABA_agonists"){
    drug_name="GABA agonist presciption,"
    thiscolour='#FF0000'
    drug_covA_n="Antidepressants"
    drug_covB_n="NACA_antiepileptic"
  } else if (drug_var_n=="Antidepressants"){
    drug_name="Antidepressant presciption,"
    thiscolour='#00A08A'
    drug_covA_n="GABA_agonists"
    drug_covB_n="NACA_antiepileptic"
  } else if (drug_var_n=="NACA_antiepileptic"){
    drug_name="Antiepileptic presciption,"
    thiscolour='#F2AD00'
    drug_covA_n="Antidepressants"
    drug_covB_n="GABA_agonists"}
  
  AdmissOutput<-AdmissionSCA(drug_var_n)
  model_ps_all<-AdmissOutput[[1]]
  model_params_all<-AdmissOutput[[2]]
  # merge and tidy for plotting
  plot.data= left_join(model_ps_all, model_params_all, by = "model_num") %>%
    arrange(estimate) %>%
    mutate(specification = row_number(),
           significant.p = ifelse(p.value < .05, "yes", "no")) %>%
    gather(variable, value, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p) %>% 
    mutate(variable = gsub("[()]", "", variable),
           variable = gsub("Intercept", "intercept", variable)) %>%
    spread(variable, value)  
  
  
  # get names of variables included in model
  variable.names = names(select(plot.data, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p,
                                -FMons, -ARons, -CAons, -drug_var1, -intercept, -Age, -drug_covA1,
                                -time_since_stroke, -NFI, -OutsOut))
  measure.names =c("FMons", "ARons", "CAons")
  model.names=c( 'OutsOut')
  
  
  # plot top panel
  top = plot.data %>%
    ggplot(aes(specification, estimate, color = significant.p, size=significant.p)) +
    geom_point(shape = "|") +scale_size_manual(values=c(2,4))+
    geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
    scale_color_manual(limits=c("no", "yes"), values = c("black", thiscolour)) +
    labs(x = "", y = paste(drug_name ,"\nregression coefficient", sep = "")) + 
    theme_minimal(base_size = 14) +
    scale_x_continuous(breaks = c(0,50,100))+
    ylim(-0.17,0.02)+
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 9),
          axis.text = element_text(color = "black"),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  
  # set plotting order for variables based on number of times it's included in better fitting models
  
  order.vars=data.frame(variable=rbind("Gender","DomAffectedY", "HADS_T", "drug_covB1"), order=rbind(4,3,2,1))
  order.meas=data.frame(variable=rbind("FMons", "ARons", "CAons"), order=rbind(3,2,1))
  
  
  
  # rename variables and plot middle panel
  middle1 = plot.data %>%
    gather(variable, value, eval(measure.names)) %>% 
    left_join(., order.meas, by = "variable") %>%
    mutate(value = ifelse(!is.na(value),"|", ""),
           variable = ifelse(variable == "FMons", "Fugl-Meyer", 
                             ifelse(variable == "ARons", "ARAT",
                                    ifelse(variable == "CAons", "CAHAI", variable)))) %>%
    ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
    scale_size_manual(values=c(2,4))+
    geom_text(aes(label = value)) +
    scale_color_manual(limits=c("no", "yes"),values = c("black", thiscolour)) +
    labs( y = "measures\n") + 
    theme_minimal(base_size = 14) +
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 9),
          axis.text = element_text(color = "black"),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  middle2 = plot.data %>%
    gather(variable, value, eval(model.names)) %>% 
    mutate(value = ifelse(!is.na(value),"|", ""),
           variable = ifelse(variable == "OutsOut", "Outliers Removed", variable)) %>%
    ggplot(aes(specification,variable,  color = significant.p, size=significant.p)) +
    geom_text(aes(label = value)) +
    scale_size_manual(values=c(2,4))+
    scale_color_manual(limits=c("no", "yes"),values = c("black", thiscolour)) +
    theme_minimal(base_size = 14) +
    labs( y = " ") + 
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 9),
          axis.text = element_text(color = "black"),
          axis.line = element_line(colour = "black"),
          axis.line.x = element_blank(),
          axis.text.x = element_blank(),
          axis.title.x = element_blank(),
          axis.ticks.x = element_blank(),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  # rename variables and plot bottom panel
  bottom = plot.data %>%
    gather(variable, value, eval(variable.names)) %>% 
    left_join(., order.vars, by = "variable") %>%
    mutate(value = ifelse(!is.na(value),"|", ""),
           variable = ifelse(variable == "GenderM", "Demographic Info", 
                             ifelse(variable == "HADS_T", "Subjective Scores", 
                                    ifelse(variable == "drug_covB1", "Other Drugs",
                                           ifelse(variable == "DomAffectedY", "Stroke Info",
                                                  variable))))) %>%
    ggplot(aes(specification, reorder(variable, order), color = significant.p, size=significant.p)) +
    geom_text(aes(label = value)) +
    scale_color_manual(limits=c("no", "yes"),values = c("black", thiscolour)) +
    scale_size_manual(values=c(2,4))+
    labs(x = "\nmodel number", y = "covariates\n") + 
    theme_minimal(base_size = 14) +
    scale_x_continuous(breaks = c(0,50,100))+
    theme(legend.title = element_text(size = 10),
          legend.text = element_text(size = 9),
          axis.text = element_text(color = "black"),
          axis.line = element_line(colour = "black"),
          legend.position = "none",
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          panel.border = element_blank(),
          panel.background = element_blank())
  
  
  # join panels
  return(list(AdSCA = cowplot::plot_grid(top, middle1, middle2,  bottom, ncol = 1, align = "v", 
                                         rel_heights = c(1,0.30, 0.20, 0.60))))
  
}
# Permutation testing 
permuteAdmissionSCA<-function(drug_var_n){
  
  if (drug_var_n=="GABA_agonists"){
    drug_covA_n="Antidepressants"
    drug_covB_n="NACA_antiepileptic"
    thiscolour='#FF0000'
  } else if (drug_var_n=="Antidepressants"){
    drug_covA_n="GABA_agonists"
    drug_covB_n="NACA_antiepileptic"
    thiscolour='#00A08A'
  } else if (drug_var_n=="NACA_antiepileptic"){
    drug_covA_n="GABA_agonists"
    drug_covB_n="Antidepressants"
    thiscolour='#F2AD00'}
# For reproducability 
set.seed(1992)

# No. of permutations
p<-500
n<-nrow(alldata4)
list<-seq(1,n,1)

# Makes a progress bar to display how were getting on 
pb = txtProgressBar(min = 1, max = p, initial = 0) 

options(na.action = "na.fail")

# Initialise a matrix to hold all the permutation data 
permlists<-matrix(0, nrow=n, ncol=p)

# take samples 
for (perms in 1:p){
  ##Permutation test a list of numbers, then reorder all the different data by the permuted list...
  permlists[,perms]<-sample(list, size=n, replace=FALSE)
}


# run full model with permutations
for (perms in 1:p){
  outcomedata<-arrange(cbind(alldata4[c("FM1_tot_pc", "ARAT1_tot_pc", "CAHAI1_tot_pc")],
                             permlists[,perms]), permlists[,perms])
  permdata<-cbind(select(alldata4, -FM1_tot_pc, -ARAT1_tot_pc, -CAHAI1_tot_pc),outcomedata)
  
  for (meas in 1:3){
    for (outliers in 1:2){
      if (outliers==1){this_data<-permdata
      
      }else{
        this_data<- filter(permdata, !(abs(permdata$FM1.FM4_pc- median(permdata$FM1.FM4_pc)) > 2.5*IQR(permdata$FM1.FM4_pc))&
                             !(abs(permdata$AR1.AR4_pc- median(permdata$AR1.AR4_pc)) > 2.5*IQR(permdata$AR1.AR4_pc))&
                             !(abs(permdata$CA1.CA4_pc- median(permdata$CA1.CA4_pc)) > 2.5*IQR(permdata$CA1.CA4_pc)))}
      
      if (meas==1){ y.var<-this_data$FM1_tot_pc
      } else if (meas==2){y.var<-this_data$ARAT1_tot_pc
      } else {y.var<-this_data$CAHAI1_tot_pc}
      
      this_data<-cbind(this_data,y.var,drug_covA=get(drug_covA_n,this_data),
                       drug_covB=get(drug_covB_n,this_data),drug_var=get(drug_var_n,this_data))
      full.model=lm(y.var ~ Age+ time_since_stroke+ DomAffected+  Gender + NFI+ 
                      HADS_T + drug_var + drug_covA + drug_covB, data= this_data)    
      
      
      quiet((
        models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                                 subset=(((Age&Gender)|(!Age&!Gender))&((NFI&HADS_T)|(!NFI&!HADS_T))&
                                           ((DomAffected&time_since_stroke)|(!DomAffected&!time_since_stroke))&
                                           ((drug_covA&drug_covB)|(!drug_covA&!drug_covB))))) 
        , all=TRUE)
      
      #run all possible nested models + filtering for only ones I want
      if (meas==1){
        FMons<- as.numeric(rep(1, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else if (meas==2){
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(1, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else {   
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(1, times=nrow(models.sca)))
      }
      
      if(outliers==2){
        OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
      }  else {
        OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
      }
      
      
      # extract parameter estimate
      (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
          tibble() %>%
          rename("model" = ".") %>%
          mutate(tidied = purrr::map(model, broom::tidy),
                 model_num = row_number()) %>%
          select(model_num, tidied) %>%
          unnest() %>%
          filter(term == paste("drug_var1", sep="")) %>%
          ungroup() %>%
          select(model_num, estimate, std.error, p.value))
      
      model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*2)+(outliers-1)))
      
      if (outliers==1&meas==1){
        model_ps_all_perms<-model_ps
      }else{
        model_ps_all_perms<-rbind(model_ps_all_perms, model_ps)
      }
      
    }} #end measure, model and outlier loops
  
  if (perms==1){
    perms_ps<-model_ps_all_perms$p.value
  }else{
    perms_ps<-cbind(perms_ps, model_ps_all_perms$p.value)
  }
  setTxtProgressBar(pb,perms)
  
}

AdmissOutput<-AdmissionSCA(drug_var_n)
test_sig<-AdmissOutput[[3]]

perm_sig<-(perms_ps<=0.05) %>%
  colSums

print('P-value for Admission as caclulated by permutation testing:')
print(drug_var_n)
print(pval<-sum(perm_sig>=test_sig)/p)

}
##### IMPROVEMENT FUNCTIONS #####   -------
# Create models through loop
ImproveSCA<- function(drug_var_n){
#set na.action for dredge
options(na.action = "na.fail")

  if (drug_var_n=="GABA_agonists"){
    drug_covA_n="Antidepressants"
    drug_covB_n="NACA_antiepileptic"
  } else if (drug_var_n=="Antidepressants"){
    drug_covA_n="GABA_agonists"
    drug_covB_n="NACA_antiepileptic"
  } else if (drug_var_n=="NACA_antiepileptic"){
    drug_covA_n="GABA_agonists"
    drug_covB_n="Antidepressants"}
  
  
  
# run full model
for (meas in 1:3){
  for (modtype in 1:3){
    for (outliers in 1:2){
      if (outliers==1){this_data<-alldata4
      }else{
        this_data<- filter(alldata4, !(abs(alldata4$FM1.FM4_pc- median(alldata4$FM1.FM4_pc)) > 2.5*IQR(alldata4$FM1.FM4_pc))&
                             !(abs(alldata4$AR1.AR4_pc- median(alldata4$AR1.AR4_pc)) > 2.5*IQR(alldata4$AR1.AR4_pc))&
                             !(abs(alldata4$CA1.CA4_pc- median(alldata4$CA1.CA4_pc)) > 2.5*IQR(alldata4$CA1.CA4_pc)))}
      
      if (meas==1){ 
        if (modtype==1){  
          cov.var<-this_data$FM1_tot_pc
          y.var<-this_data$FM4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relFM1.FM4_pc
        } else {
          cov.var<-1
          y.var<-this_data$FM1.FM4_pc}
      } else if (meas==2){
        if (modtype==1){ 
          cov.var<-this_data$ARAT1_tot_pc
          y.var<-this_data$ARAT4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relAR1.AR4_pc
        } else {
          cov.var<-1
          y.var<-this_data$AR1.AR4_pc}
      } else {
        if (modtype==1){ 
          cov.var<-this_data$CAHAI1_tot_pc
          y.var<-this_data$CAHAI4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relCA1.CA4_pc
        } else {
          cov.var<-1
          y.var<-this_data$CA1.CA4_pc}}
      
      this_data<-cbind(this_data,y.var,cov.var,drug_covA=get(drug_covA_n,this_data),
                       drug_covB=get(drug_covB_n,this_data),drug_var=get(drug_var_n,this_data))
      full.model=lm(y.var ~ cov.var + Age+ time_since_stroke+ DomAffected+ Gender + NFI+ 
                      HADS_T + drug_var + drug_covA + drug_covB, data= this_data)    
      
      models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                               subset=(((Age&Gender)|(!Age&!Gender))&cov.var&((NFI&HADS_T)|(!NFI&!HADS_T))&
                                         ((DomAffected&time_since_stroke)|(!DomAffected&!time_since_stroke))&
                                         ((drug_covA&drug_covB)|(!drug_covA&!drug_covB))))
      
      
      #run all possible nested models + filtering for only ones I want
      if (meas==1){
        FMons<- as.numeric(rep(1, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else if (meas==2){
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(1, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else {   
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(1, times=nrow(models.sca)))
      }
      
      if(outliers==2){
        OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
      }  else {
        OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
      }
      
      if(modtype==1){
        ChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(1, times=nrow(models.sca)))
      }  else if (modtype==2){
        ChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(1, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
      } else{ 
        ChangeMod<- as.numeric(rep(1, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(NA, times=nrow(models.sca)))}
      
      
      # extract parameter estimate
      (model_params = MuMIn::get.models(models.sca, subset = TRUE) %>%
          tibble() %>%
          dplyr::rename("model" = ".") %>%
          mutate(tidied = purrr::map(model, broom::tidy),
                 model_num = row_number()) %>%
          select(model_num, tidied) %>%
          unnest() %>%
          select(model_num, term, estimate) %>%
          spread(term, estimate) %>%
          select(-starts_with("sd")))
      
      model_params$model_num<-model_params$model_num+(nrow(models.sca)*(((meas-1)*3*2)+((modtype-1)*2)+(outliers-1)))
      model_params<-cbind(model_params,FMons, ARons, CAons, OutsOut, ChangeMod, RelChangeMod, OutcomeMod)
      
      (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
          tibble() %>%
          dplyr::rename("model" = ".") %>%
          mutate(tidied = purrr::map(model, broom::tidy),
                 model_num = row_number()) %>%
          select(model_num, tidied) %>%
          unnest() %>%
          filter(term == paste("drug_var1", sep="")) %>%
          ungroup() %>%
          select(model_num, estimate, std.error, p.value))
      
      model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*3*2)+((modtype-1)*2)+(outliers-1)))
      
      if (outliers==1&meas==1&modtype==1){
        model_params_all<-model_params
        model_ps_all<-model_ps
      }else{
        model_params_all<-full_join(model_params_all, model_params)
        model_ps_all<-full_join(model_ps_all, model_ps)
      }}}}

print(test_sig<-(model_ps_all$p.value<=0.05) %>%
        sum())
print(summary(model_ps_all$estimate[model_ps_all$p.value<=0.05]))
print(summary(model_ps_all$estimate))
return(list(model_ps_all, model_params_all, test_sig))

}
# Plotting the outcomes
plotImproveSCA<-function(drug_var_n){
  
  if (drug_var_n=="GABA_agonists"){
    drug_name="GABA agonist presciption,"
    thiscolour='#FF0000'
    drug_covA_n="Antidepressants"
    drug_covB_n="NACA_antiepileptic"
  } else if (drug_var_n=="Antidepressants"){
    drug_name="Antidepressant presciption,"
    thiscolour='#00A08A'
    drug_covA_n="GABA_agonists"
    drug_covB_n="NACA_antiepileptic"
  } else if (drug_var_n=="NACA_antiepileptic"){
    drug_name="Antiepileptic presciption,"
    thiscolour='#F2AD00'
    drug_covA_n="Antidepressants"
    drug_covB_n="GABA_agonists"}
  
  ImproveOutput<-ImproveSCA(drug_var_n)
  model_ps_all<-ImproveOutput[[1]]
  model_params_all<-ImproveOutput[[2]]
# merge and tidy for plotting
plot.data= left_join(model_ps_all, model_params_all, by = "model_num") %>%
  arrange(estimate) %>%
  mutate(specification = row_number(),
         significant.p = ifelse(p.value < .05, "yes", "no")) %>%
  gather(variable, value, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p) %>% 
  mutate(variable = gsub("[()]", "", variable),
         variable = gsub("Intercept", "intercept", variable)) %>%
  spread(variable, value)  


# get names of variables included in model
variable.names = names(select(plot.data, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p,
                              -FMons, -ARons, -CAons, -cov.var, -drug_var1, -intercept, -Age, -drug_covA1,
                              -time_since_stroke, -NFI, -OutsOut, -OutcomeMod, -RelChangeMod, -ChangeMod))
measure.names =c("FMons", "ARons", "CAons")
model.names=c('ChangeMod', 'RelChangeMod', 'OutcomeMod')
out.name='OutsOut'

# plot top panel
top = plot.data %>%
  ggplot(aes(specification, estimate, color = significant.p,size = significant.p)) +
  geom_point(shape = "|") +scale_size_manual(values=c(2,4))+
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(limits=c("no", "yes"), values = c("black", thiscolour)) +
  labs(x = "") + 
  theme_minimal(base_size = 14) +
  ylim(-0.17,0.02)+
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# set plotting order for variables based on number of times it's included in better fitting models

order.vars=data.frame(variable=rbind("Gender","DomAffectedY", "HADS_T", "drug_covB1"), order=rbind(4,3,2,1))
order.meas=data.frame(variable=rbind("FMons", "ARons", "CAons"), order=rbind(3,2,1))
order.mod=data.frame(variable=rbind("RelChangeMod", "ChangeMod", "OutcomeMod"), order=rbind(3,2,1))


# rename variables and plot middle panel
middle1 = plot.data %>%
  gather(variable, value, eval(measure.names)) %>% 
  left_join(., order.meas, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "FMons", "Fugl-Meyer", 
                           ifelse(variable == "ARons", "ARAT",
                                  ifelse(variable == "CAons", "CAHAI", variable)))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p,size = significant.p)) +
  geom_text(aes(label = value)) +scale_size_manual(values=c(2,4))+
  scale_color_manual(limits=c("no", "yes"),values = c("black", thiscolour)) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

middle2 = plot.data %>%
  gather(variable, value, eval(out.name)) %>% 
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "OutsOut", "Outliers Removed", variable)) %>%
  ggplot(aes(specification, variable, color = significant.p ,size = significant.p)) +
  geom_text(aes(label = value)) +scale_size_manual(values=c(2,4))+
  scale_color_manual(limits=c("no", "yes"),values = c("black", thiscolour)) +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

middle3 = plot.data %>%
  gather(variable, value, eval(model.names)) %>% 
  left_join(., order.mod, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "ChangeMod", "Abs Recovery Model",
                           ifelse(variable == "OutcomeMod", "Outcome Model",
                                  ifelse(variable == "RelChangeMod", " Rel Recovery Model", variable)))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p,size = significant.p)) +
  geom_text(aes(label = value)) +scale_size_manual(values=c(2,4))+
  scale_color_manual(limits=c("no", "yes"),values = c("black", thiscolour)) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# rename variables and plot bottom panel
bottom = plot.data %>%
  gather(variable, value, eval(variable.names)) %>% 
  left_join(., order.vars, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "GenderM", "Demographic Info", 
                                  ifelse(variable == "HADS_T", "Subjective Scores", 
                                         ifelse(variable == "drug_covB1", "Other Drugs",
                                                       ifelse(variable == "DomAffectedY", "Stroke Info",
                                                              variable))))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p,size = significant.p)) +
  geom_text(aes(label = value)) +scale_size_manual(values=c(2,4))+
  scale_color_manual(limits=c("no", "yes"),values = c("black", thiscolour)) +
  labs(x = "\nmodel number") + 
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# join panels
return(list(ImpSCA = cowplot::plot_grid(top, middle1, middle2,middle3,  bottom, ncol = 1, align = "v", 
                           rel_heights = c(1,0.245,0.125, 0.245, 0.485))))


}
# Plot Both 
plotBothSCA<-function(drug_var_n){
  AdSCA<-plotAdmissionSCA(drug_var_n)[[1]]
  ImpSCA<-plotImproveSCA(drug_var_n)[[1]]
(BothSCA=cowplot::plot_grid(AdSCA, ImpSCA, ncol=2, align = 'h', rel_widths = c(1.9,4),
                            labels="AUTO",label_size = 20))}
# Permutation testing 
permuteImproveSCA<-function(drug_var_n){
  
  if (drug_var_n=="GABA_agonists"){
    drug_covA_n="Antidepressants"
    drug_covB_n="NACA_antiepileptic"
    thiscolour='#FF0000'
  } else if (drug_var_n=="Antidepressants"){
    drug_covA_n="GABA_agonists"
    drug_covB_n="NACA_antiepileptic"
    thiscolour='#00A08A'
  } else if (drug_var_n=="NACA_antiepileptic"){
    drug_covA_n="GABA_agonists"
    drug_covB_n="Antidepressants"
    thiscolour='#F2AD00'}
  
  # For reproducability 
set.seed(1992)

# No. of permutations
p<-500
n<-nrow(alldata4)
list<-seq(1,n,1)

# Makes a progress bar to display how were getting on 
pb = txtProgressBar(min = 1, max = p, initial = 0) 

options(na.action = "na.fail")

# Initialise a matrix to hold all the permutation data 
permlists<-matrix(0, nrow=n, ncol=p)

# take samples 
for (perms in 1:p){
  ##Permutation test a list of numbers, then reorder all the different data by the permuted list...
  permlists[,perms]<-sample(list, size=n, replace=FALSE)
}


# run full model with permutations
for (perms in 1:p){
  outcomedata<-arrange(cbind(alldata4[c("FM4_tot_pc", "ARAT4_tot_pc", "CAHAI4_tot_pc", "relFM1.FM4_pc", "relAR1.AR4_pc",
                                        "relCA1.CA4_pc", "FM1.FM4_pc", "AR1.AR4_pc", "CA1.CA4_pc")],
                             permlists[,perms]), permlists[,perms])
  permdata<-cbind(select(alldata4, -FM4_tot_pc, -ARAT4_tot_pc, -CAHAI4_tot_pc, -relFM1.FM4_pc, -relAR1.AR4_pc, 
                         -relCA1.CA4_pc, -FM1.FM4_pc, -AR1.AR4_pc, -CA1.CA4_pc),outcomedata)
  
  for (meas in 1:3){
    for (modtype in 1:3){
      for (outliers in 1:2){
        if (outliers==1){this_data<-permdata
        
        }else{
          this_data<- filter(permdata, !(abs(permdata$FM1.FM4_pc- median(permdata$FM1.FM4_pc)) > 2.5*IQR(permdata$FM1.FM4_pc))&
                               !(abs(permdata$AR1.AR4_pc- median(permdata$AR1.AR4_pc)) > 2.5*IQR(permdata$AR1.AR4_pc))&
                               !(abs(permdata$CA1.CA4_pc- median(permdata$CA1.CA4_pc)) > 2.5*IQR(permdata$CA1.CA4_pc)))}
        
        if (meas==1){ 
          if (modtype==1){  
            cov.var<-this_data$FM1_tot_pc
            y.var<-this_data$FM4_tot_pc
          } else if (modtype==2) { 
            cov.var<-1
            y.var<-this_data$relFM1.FM4_pc
          } else {
            cov.var<-1
            y.var<-this_data$FM1.FM4_pc}
        } else if (meas==2){
          if (modtype==1){ 
            cov.var<-this_data$ARAT1_tot_pc
            y.var<-this_data$ARAT4_tot_pc
          } else if (modtype==2) { 
            cov.var<-1
            y.var<-this_data$relAR1.AR4_pc
          } else {
            cov.var<-1
            y.var<-this_data$AR1.AR4_pc}
        } else {
          if (modtype==1){ 
            cov.var<-this_data$CAHAI1_tot_pc
            y.var<-this_data$CAHAI4_tot_pc
          } else if (modtype==2) { 
            cov.var<-1
            y.var<-this_data$relCA1.CA4_pc
          } else {
            cov.var<-1
            y.var<-this_data$CA1.CA4_pc}}
        
        this_data<-cbind(this_data,y.var,cov.var,drug_covA=get(drug_covA_n,this_data),
                         drug_covB=get(drug_covB_n,this_data),drug_var=get(drug_var_n,this_data))
        full.model=lm(y.var ~ cov.var + Age+ time_since_stroke+ DomAffected+ Gender + NFI+ 
                        HADS_T + drug_var + drug_covA + drug_covB, data= this_data)    
        
        quiet((
          models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                                   subset=(((Age&Gender)|(!Age&!Gender))&cov.var&((NFI&HADS_T)|(!NFI&!HADS_T))&
                                             ((DomAffected&time_since_stroke)|(!DomAffected&!time_since_stroke))&
                                             ((drug_covA&drug_covB)|(!drug_covA&!drug_covB))))) , all=TRUE)
        
        #run all possible nested models + filtering for only ones I want
        if (meas==1){
          FMons<- as.numeric(rep(1, times=nrow(models.sca)))
          ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
          CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
        }else if (meas==2){
          FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
          ARons<- as.numeric(rep(1, times=nrow(models.sca)))
          CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
        }else {   
          FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
          ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
          CAons<- as.numeric(rep(1, times=nrow(models.sca)))
        }
        
        if(outliers==2){
          OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
        }  else {
          OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
        }
        
        if(modtype==1){
          ChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
          RelChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
          OutcomeMod<- as.numeric(rep(1, times=nrow(models.sca)))
        }  else if (modtype==2){
          ChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
          RelChangeMod<- as.numeric(rep(1, times=nrow(models.sca)))
          OutcomeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        } else{ 
          ChangeMod<- as.numeric(rep(1, times=nrow(models.sca)))
          RelChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
          OutcomeMod<- as.numeric(rep(NA, times=nrow(models.sca)))}
        
        
        # extract parameter estimate
        (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
            tibble() %>%
            rename("model" = ".") %>%
            mutate(tidied = purrr::map(model, broom::tidy),
                   model_num = row_number()) %>%
            select(model_num, tidied) %>%
            unnest() %>%
            filter(term == paste("drug_var1", sep="")) %>%
            ungroup() %>%
            select(model_num, estimate, std.error, p.value))
        
        model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*3*2)+((modtype-1)*2)+(outliers-1)))
        
        if (outliers==1&meas==1&modtype==1){
          model_ps_all_perms<-model_ps
        }else{
          model_ps_all_perms<-rbind(model_ps_all_perms, model_ps)
        }
        
      }}} #end measure, model and outlier loops
  
  if (perms==1){
    perms_ps<-model_ps_all_perms$p.value
  }else{
    perms_ps<-cbind(perms_ps, model_ps_all_perms$p.value)
  }
  setTxtProgressBar(pb,perms)
  
}

ImproveOutput<-ImproveSCA(drug_var_n)
test_sig<-ImproveOutput[[3]]

perm_sig<-(perms_ps<=0.05) %>%
  colSums

print('P-value for Recovery as caclulated by permutation testing:')
print(drug_var_n)
print(pval<-sum(perm_sig>=test_sig)/p)
}
##### RUN DRUG SCAs ##### -----
plotBothSCA("GABA_agonists")
plotBothSCA("NACA_antiepileptic")
plotBothSCA("Antidepressants")

permuteAdmissionSCA("GABA_agonists")
permuteAdmissionSCA("NACA_antiepileptic")
permuteAdmissionSCA("Antidepressants")

permuteImproveSCA("GABA_agonists")
permuteImproveSCA("NACA_antiepileptic")
permuteImproveSCA("Antidepressants")


##### DEMOGRAPHICS #### -----
allincluded<-datatable[!is.na(datatable$Benzos)&!is.na(datatable$relAR1.AR4)&!is.na(datatable$HADS_T)&!is.na(datatable$NFI)&is.finite(datatable$relAR1.AR4),]
withT4<- datatable[!is.na(datatable$AR1.AR4)&!is.na(datatable$Age),]
withT4HADS<- datatable[!is.na(datatable$AR1.AR4)&!is.na(datatable$HADS_T)&!is.na(datatable$NFI)&!is.na(datatable$Age),]
allexcluded<- datatable[(is.na(datatable$Benzos)|is.na(datatable$AR1.AR4)|is.na(datatable$HADS_T)|is.na(datatable$NFI))&!is.na(datatable$Age),]

t.test(allincluded$Age,allexcluded$Age)
shapiro.test(allincluded$Age)
wilcox.test(allincluded$Age,allexcluded$Age)
sd(allincluded$Age)
summary(allincluded$Age)
IQR(allincluded$Age)
sd(allexcluded$Age)
summary(allexcluded$Age)
IQR(allexcluded$Age)

gendcontab <- rbind(included=summary(allincluded$Gender),excluded=summary(allexcluded$Gender)) %>%
  as_tibble(.,rownames = NA) %>%
  select(.,-V1)%>%
  print.data.frame()
chisq.test(gendcontab)

wilcox.test(allincluded$time_since_stroke,allexcluded$time_since_stroke)
shapiro.test(allincluded$time_since_stroke)
t.test(allincluded$time_since_stroke,allexcluded$time_since_stroke)
sd(allincluded$time_since_stroke)
IQR(allincluded$time_since_stroke)
summary(allincluded$time_since_stroke)
sd(allexcluded$time_since_stroke)
IQR(allexcluded$time_since_stroke)
summary(allexcluded$time_since_stroke)

t.test(allincluded$time_since_program,allexcluded$time_since_program)
sd(allincluded$time_since_program)
summary(allincluded$time_since_program)
sd(allexcluded$time_since_program)
summary(allincluded$time_since_program)

lescontab <- rbind(included=summary(allincluded$Cause_refined),excluded=summary(allexcluded$Cause_refined)) %>%
  as_tibble(.,rownames = NA) %>%
  select(.,-V1)%>%
  print.data.frame()
chisq.test(lescontab)
AfLimbcontab <- rbind(included=summary(allincluded$Affected.limb),excluded=summary(allexcluded$Affected.limb)) %>%
  as_tibble(.,rownames = NA) %>%
  select(.,-V1)%>%
  print.data.frame()
chisq.test(AfLimbcontab)
DomAfcontab <- rbind(included=summary(allincluded$DomAffected),excluded=summary(allexcluded$DomAffected)) %>%
  as_tibble(.,rownames = NA)%>%
  print.data.frame()
chisq.test(DomAfcontab)

shapiro.test(allincluded$HADS_T)
wilcox.test(allincluded$HADS_T,allexcluded$HADS_T)
t.test(allincluded$HADS_T,allexcluded$HADS_T)
sd(allincluded$HADS_T)
IQR(allincluded$HADS_T)
summary(allincluded$HADS_T)
sd(allexcluded$HADS_T,na.rm=TRUE)
IQR(allexcluded$HADS_T, na.rm=TRUE)
summary(allexcluded$HADS_T)

shapiro.test(allincluded$NFI)
t.test(allincluded$NFI,allexcluded$NFI)
wilcox.test(allincluded$NFI,allexcluded$NFI)
summary(allincluded$NFI)
sd(allincluded$NFI)
IQR(allincluded$NFI)
summary(allexcluded$NFI)
sd(allexcluded$NFI, na.rm=TRUE)
IQR(allexcluded$NFI,na.rm=TRUE)


shapiro.test(allincluded$BARTHEL)
t.test(allincluded$BARTHEL,allexcluded$BARTHEL)
wilcox.test(allincluded$BARTHEL,allexcluded$BARTHEL)
sd(allincluded$BARTHEL, na.rm=TRUE)
summary(allincluded$BARTHEL)
IQR(allincluded$BARTHEL, na.rm=TRUE)
sd(allexcluded$BARTHEL, na.rm=TRUE)
summary(allexcluded$BARTHEL)
IQR(allexcluded$BARTHEL, na.rm=TRUE)


t.test(allincluded$MoCA,allexcluded$MoCA)
sd(allincluded$MoCA, na.rm=TRUE)
sd(allexcluded$MoCA, na.rm=TRUE)

totalEx<-allexcluded %>% filter(!is.na(Antidepressants)) %>% nrow()
AdEx<-allexcluded %>% filter(Antidepressants==1) %>% nrow()
AeEx<-allexcluded %>% filter(Antiepileptic==1) %>% nrow()
GABAEx<-allexcluded %>% filter(GABA_agonists==1) %>% nrow()
NoEx<- totalEx-(AdEx+AeEx+GABAEx)

AdIn<-allincluded %>% filter(Antidepressants==1) %>% nrow()
AeIn<-allincluded %>% filter(Antiepileptic==1) %>% nrow()
GABAIn<-allincluded %>% filter(GABA_agonists==1) %>% nrow()
NoIn<-nrow(allincluded)-(AdEx+AeEx+GABAEx)

Drugcontab <- rbind(included=cbind(GABAIn,AeIn, AdIn),
                    excluded=cbind(GABAEx,AeEx, AdEx))%>%
  as_tibble(.,rownames = NA)%>%
  print.data.frame()
chisq.test(Drugcontab)

depress.aov <- aov(HADS_T~DepressionGroups, data=alldata4)
summary(depress.aov)
TukeyHSD(depress.aov)
t.test(alldata4$HADS_T[alldata4$Antidepressants==1],alldata4$HADS_T[alldata4$Antidepressants==0])
t.test(alldata4$HADS_T[alldata4$GABA_agonists==1],alldata4$HADS_T[alldata4$GABA_agonists==0])
t.test(alldata4$HADS_T[alldata4$NACA_antiepileptic==1],alldata4$HADS_T[alldata4$NACA_antiepileptic==0])
t.test(alldata4$NFI[alldata4$Antidepressants==1],alldata4$NFI[alldata4$Antidepressants==0])

##### OTHER FIGURES ####  -----

library(reshape2)

reshapeddataFM<-melt(alldata4, id.vars= "Hospital.No.", measure.vars=c('FM1_tot', 'FM2_tot', 'FM3_tot', 'FM4_tot'),
                     variable.name = "FMtimept", value.name = "FMscore")

reshapeddataFM_full<- merge(reshapeddataFM, alldata4, by= 'Hospital.No.')



FM_Antispasticity<- ggplot(data=reshapeddataFM_full, aes(x=FMtimept, y=FMscore, color=factor(GABA_agonists), shape=factor(GABA_agonists))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  stat_summary(aes(x=FMtimept, y=FMscore, group=GABA_agonists, color=factor(GABA_agonists)),fun.y=mean, geom='line')+
  theme(legend.position='none',axis.text.x = element_blank())+scale_color_manual(values=c( 'black','#FF0000'))+ 
  scale_shape_manual(values=c(16,17))+labs( x=NULL, y='FM score')

FM_Antiepileptic<- ggplot(data=reshapeddataFM_full, aes(x=FMtimept, y=FMscore, color=factor(NACA_antiepileptic), shape=factor(NACA_antiepileptic))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(aes(x=FMtimept, y=FMscore, group=NACA_antiepileptic, color=factor(NACA_antiepileptic)),fun.y=mean, geom='line')+
  theme(legend.position='none',axis.text.x = element_blank())+scale_color_manual(values=c( 'black','#F2AD00'))+ 
  scale_shape_manual(values=c(16,17))+labs( x=NULL, y='FM score')

FM_DepressionGroups<- ggplot(data=reshapeddataFM_full, aes(x=FMtimept, y=FMscore, color=factor(DepressionGroups), shape=factor(DepressionGroups))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(fun.data = mean_se, geom = "point",size=2)+
  stat_summary(aes(x=FMtimept, y=FMscore, group=DepressionGroups, color=factor(DepressionGroups)),fun.y=mean, geom='line')+
  theme(legend.position='none', axis.text.x = element_text(size = 10))+scale_color_manual(values=c( 'grey70', 'grey50', '#00A08A'))+
  scale_shape_manual(values=c(18,15, 17))+labs( x=NULL, y='FM score')+ scale_x_discrete(labels= c('Admission','Discharge','6weeks', '6months'))


reshapeddataAR<-melt(alldata4, id.vars= "Hospital.No.", measure.vars=c('ARAT1_tot', 'ARAT2_tot', 'ARAT3_tot', 'ARAT4_tot'),
                     variable.name = "ARtimept", value.name = "ARscore")

reshapeddataAR_full<- merge(reshapeddataAR, alldata4, by= 'Hospital.No.')


AR_Antispasticity<- ggplot(data=reshapeddataAR_full, aes(x=ARtimept, y=ARscore, color=factor(GABA_agonists), shape=factor(GABA_agonists))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  stat_summary(aes(x=ARtimept, y=ARscore, group=GABA_agonists, color=factor(GABA_agonists)),fun.y=mean, geom='line')+
  theme(legend.position='none',axis.text.x = element_blank())+scale_color_manual(values=c( 'black','#FF0000'))+ 
  scale_shape_manual(values=c(16,17))+labs( x=NULL, y='ARAT score')

AR_Antiepileptic<- ggplot(data=reshapeddataAR_full, aes(x=ARtimept, y=ARscore, color=factor(NACA_antiepileptic), shape=factor(NACA_antiepileptic))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(aes(x=ARtimept, y=ARscore, group=NACA_antiepileptic, color=factor(NACA_antiepileptic)),fun.y=mean, geom='line')+
  theme(legend.position='none',axis.text.x = element_blank())+scale_color_manual(values=c( 'black','#F2AD00'))+ 
  scale_shape_manual(values=c(16,17))+labs( x=NULL, y='ARAT score')

AR_DepressionGroups<- ggplot(data=reshapeddataAR_full, aes(x=ARtimept, y=ARscore, color=factor(DepressionGroups), shape=factor(DepressionGroups))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(fun.data = mean_se, geom = "point",size=2)+
  stat_summary(aes(x=ARtimept, y=ARscore, group=DepressionGroups, color=factor(DepressionGroups)),fun.y=mean, geom='line')+
  theme(legend.position='none', axis.text.x = element_text(size = 10))+scale_color_manual(values=c('grey70', 'grey50', '#00A08A'))+
  scale_shape_manual(values=c(18,15,17))+
  labs( x=NULL, y='ARAT score')+ scale_x_discrete(labels= c('Admission','Discharge','6weeks', '6months'))


reshapeddataCA<-melt(alldata4, id.vars= "Hospital.No.", measure.vars=c('CAHAI1_tot', 'CAHAI2_tot', 'CAHAI3_tot', 'CAHAI4_tot'),
                     variable.name = "CAtimept", value.name = "CAscore")

reshapeddataCA_full<- merge(reshapeddataCA, alldata4, by= 'Hospital.No.')

CA_Antispasticity<- ggplot(data=reshapeddataCA_full, aes(x=CAtimept, y=CAscore, color=factor(GABA_agonists), shape=factor(GABA_agonists))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  stat_summary(aes(x=CAtimept, y=CAscore, group=GABA_agonists, color=factor(GABA_agonists)),fun.y=mean, geom='line')+
  theme(legend.position='bottom', axis.text.x = element_text(size = 12),legend.text = element_text(size = 12))+
  scale_color_manual(values=c( 'black','#FF0000'),name=NULL,
                                                   labels= c('Off GABA agonist', 'On GABA agonist'))+
  scale_shape_manual(values=c( 16,17),name=NULL,
                     labels= c('Off GABA agonist', 'On GABA agonist'))+
  labs( x=NULL, y='CAHAI score')+ scale_x_discrete(labels= c('Admission','Discharge','6weeks', '6months'))

CA_Antiepileptic<- ggplot(data=reshapeddataCA_full, aes(x=CAtimept, y=CAscore, color=factor(NACA_antiepileptic),shape=factor(NACA_antiepileptic))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  stat_summary(aes(x=CAtimept, y=CAscore, group=NACA_antiepileptic, color=factor(NACA_antiepileptic)),fun.y=mean, geom='line')+
  theme(legend.position='bottom', axis.text.x = element_text(size = 12),legend.text = element_text(size = 12))+
  scale_color_manual(values=c( 'black','#F2AD00'), name=NULL,
                                                     labels= c('Off Antiepileptic', 'On Antiepileptic'))+
  scale_shape_manual(values=c( 16,17),name=NULL,
                     labels= c('Off Antiepileptic', 'On Antiepileptic'))+
  labs( x=NULL, y='CAHAI score')+ scale_x_discrete(labels= c('Admission','Discharge','6weeks', '6months'))

CA_DepressionGroups<- ggplot(data=reshapeddataCA_full, aes(x=CAtimept, y=CAscore, color=factor(DepressionGroups),shape=factor(DepressionGroups))) +
  geom_violin(inherit.aes = TRUE, position='identity', linetype='blank',scale='count')+ theme_minimal()+ 
  geom_violin(inherit.aes = TRUE, position='identity', fill=NA,linetype='dotted', scale='count')+ 
  stat_summary(fun.data = mean_se, geom = "errorbar",aes(width=0.15))+
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  stat_summary(aes(x=CAtimept, y=CAscore, group=DepressionGroups, color=factor(DepressionGroups)),fun.y=mean, geom='line')+
  theme(legend.position='none', axis.text.x = element_text(size = 10),legend.text = element_text(size = 12))+
  scale_color_manual(values=c( 'grey70', 'grey50', '#00A08A'))+
  scale_shape_manual(values=c(18,15,17),name=NULL)+
  labs( x=NULL, y='CAHAI score')+ scale_x_discrete(labels= c('Admission','Discharge','6weeks', '6months'))

DepressTot<-ggplot(data=alldata4, aes(x=Antidepressants, y=HADS_T, color=factor(Antidepressants), shape=factor(Antidepressants))) +
  geom_violin(inherit.aes = TRUE, scale='count',linetype='dotted') + theme_minimal()+ 
  stat_summary( inherit.aes = TRUE,fun.data = mean_se, geom = "errorbar",
                aes(width=0.15),position=position_dodge(0.05))+
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  theme(legend.position = c(0.95, 0.91), legend.text=element_text(size=9),
        plot.margin=unit(c(5,-5,5,8), 'pt'),legend.title = element_blank(),
        axis.title.y = element_text(vjust = 3),axis.text.x=element_text(size=10) )+
  scale_color_manual(values=c( 'black','#00A08A'),labels= c('Off AD\n   All', 'On AD\n   All'))+ 
  scale_shape_manual(values=c( 16,17),labels= c('Off AD\n   All', 'On AD\n   All'))+ 
  scale_x_discrete(labels= c('Off AD\nAll', 'On AD\nAll'))+labs(x=NULL,y='HADS score')
DepressSub<-ggplot(data=alldata4[alldata4$Antidepressants==0,], 
                   aes(x=DepressionGroups, y=HADS_T, color=factor(DepressionGroups),shape=factor(DepressionGroups))) +
  geom_violin(inherit.aes = TRUE, scale='count',linetype='dotted')+theme_minimal()+
  scale_color_manual(values=c('grey70', 'grey50'),labels= c('  Off AD\nHigh HADS', '  Off AD\nLow HADS'))+ 
  scale_shape_manual(values=c(18,15),labels= c('  Off AD\nHigh HADS', '  Off AD\nLow HADS'))+ 
  stat_summary( inherit.aes = TRUE,fun.data = mean_se, geom = "errorbar",
                aes(width=0.15),position=position_dodge(0.05))+
  stat_summary(fun.data = mean_se, geom = "point", size=2)+
  theme(legend.position = c(0.60, 0.91), axis.line.x = element_blank(), legend.title =element_blank(),
        axis.text.y = element_blank(), axis.title.y = element_blank(),
        axis.ticks.y = element_blank(), plot.margin=unit(c(5,6,5,-10), 'pt'),
        axis.text.x=element_text(size=10), legend.text=element_text(size=9))+
  labs(x=NULL, y=NULL)+ scale_x_discrete(labels= c('Off AD\nHigh HADS', 'Off AD\nLow HADS'))

HADS_scores <- plot_grid(DepressTot, DepressSub, rel_widths = c(1.3,1))

plot_grid(HADS_scores,FM_DepressionGroups,AR_DepressionGroups,
          CA_DepressionGroups, nrow=2,labels="AUTO")


plot_grid( FM_Antiepileptic,AR_Antiepileptic,CA_Antiepileptic, 
          rel_heights=c(1,1,1.4), nrow=3,labels="AUTO")
plot_grid(FM_Antispasticity,AR_Antispasticity,CA_Antispasticity, 
          rel_heights=c(1,1,1.4), nrow=3,labels="AUTO")
plot_grid(FM_Antidepressant,AR_Antidepressant,CA_Antidepressant, 
          rel_heights=c(1,1,1.4), nrow=3,labels="AUTO")
plot_grid(FM_DepressionGroups,AR_DepressionGroups,CA_DepressionGroups, 
          rel_heights=c(1,1,1.3), nrow=3,labels="AUTO")


library(eulerr)
DrugsVenn<- euler(c("Antidepressants"=(length(which(alldata4$Antidepressants==1))),
                     "Antiepileptics"=(length(which(alldata4$NACA_antiepileptic==1))),
                     "GABA agonists"=(length(which(alldata4$GABA_agonists==1))),
                     "Antidepressants&Antiepileptics"=(length(which(alldata4$Antidepressants==1&alldata4$NACA_antiepileptic==1))),
                     "Antidepressants&GABA agonists"=(length(which(alldata4$Antidepressants==1&alldata4$GABA_agonists==1))),
                     "Antiepileptics&GABA agonists"=(length(which(alldata4$NACA_antiepileptic==1&alldata4$GABA_agonists==1))),
                     "Antidepressants&Antiepileptics&GABA agonists"=(length(which(alldata4$Antidepressants==1&alldata4$NACA_antiepileptic==1&alldata4$GABA_agonists==1)))))
plot(DrugsVenn, quantities = TRUE, fill=c('#47A195','#F3CF6C','#FF5555'), 
      edges = NULL)      



##### HADS SCA ##### ####-----
# Admission model-----
for (meas in 1:3){
  for (outliers in 1:2){
    if (outliers==1){this_data<-alldata4
    }else{
      this_data<- filter(alldata4, !(abs(alldata4$FM1.FM4_pc- median(alldata4$FM1.FM4_pc)) > 2.5*IQR(alldata4$FM1.FM4_pc))&
                           !(abs(alldata4$AR1.AR4_pc- median(alldata4$AR1.AR4_pc)) > 2.5*IQR(alldata4$AR1.AR4_pc))&
                           !(abs(alldata4$CA1.CA4_pc- median(alldata4$CA1.CA4_pc)) > 2.5*IQR(alldata4$CA1.CA4_pc)))}
    
    if (meas==1){ y.var<-this_data$FM1_tot_pc}
    else if (meas==2){y.var<-this_data$ARAT1_tot_pc}
    else {y.var<-this_data$CAHAI1_tot_pc}
    
    
    full.model=lm(y.var ~ Age+ time_since_stroke+ DomAffected+  Gender + NFI+ 
                    HADS_T + Antidepressants + GABA_agonists + NACA_antiepileptic, data= this_data)    
    
    models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                             subset=(((Age&Gender)|(!Age&!Gender))&
                                       ((Antidepressants&GABA_agonists&NACA_antiepileptic)|(!Antidepressants&!GABA_agonists&!NACA_antiepileptic))&
                                       ((time_since_stroke&DomAffected)|(!time_since_stroke&!DomAffected))))
    
    
    #run all possible nested models + filtering for only ones I want
    if (meas==1){
      FMons<- as.numeric(rep(1, times=nrow(models.sca)))
      ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
      CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
    }else if (meas==2){
      FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
      ARons<- as.numeric(rep(1, times=nrow(models.sca)))
      CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
    }else {   
      FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
      ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
      CAons<- as.numeric(rep(1, times=nrow(models.sca)))
    }
    
    if(outliers==2){
      OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
    }  else {
      OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
    }
    
    
    
    # extract parameter estimate
    (model_params = MuMIn::get.models(models.sca, subset = TRUE) %>%
        tibble() %>%
        dplyr::rename("model" = ".") %>%
        mutate(tidied = purrr::map(model, broom::tidy),
               model_num = row_number()) %>%
        select(model_num, tidied) %>%
        unnest() %>%
        select(model_num, term, estimate) %>%
        spread(term, estimate)) %>%
      select(-starts_with("sd"))
    
    model_params$model_num<-model_params$model_num+(nrow(models.sca)*(((meas-1)*2)+(outliers-1)))
    model_params<-cbind(model_params,FMons, ARons, CAons, OutsOut)
    
    (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
        tibble() %>%
        dplyr::rename("model" = ".") %>%
        mutate(tidied = purrr::map(model, broom::tidy),
               model_num = row_number()) %>%
        select(model_num, tidied) %>%
        unnest() %>%
        filter(term ==  paste("HADS_T", sep="")) %>%
        ungroup() %>%
        select(model_num, estimate, std.error, p.value))
    
    model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*2)+(outliers-1)))
    
    if (outliers==1&meas==1){
      model_params_all<-model_params
      model_ps_all<-model_ps
    }else{
      model_params_all<-full_join(model_params_all, model_params)
      model_ps_all<-full_join(model_ps_all, model_ps)
    }
  }}

print(test_sig<-(model_ps_all$p.value<=0.05) %>%
        sum())
summary(model_ps_all$estimate[model_ps_all$p.value<=0.05])
summary(model_ps_all$estimate)

# Plotting the outcomes ----
# merge and tidy for plotting
plot.data= left_join(model_ps_all, model_params_all, by = "model_num") %>%
  arrange(estimate) %>%
  mutate(specification = row_number(),
         significant.p = ifelse(p.value < .05, "yes", "no")) %>%
  gather(variable, value, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p) %>% 
  mutate(variable = gsub("[()]", "", variable),
         variable = gsub("Intercept", "intercept", variable)) %>%
  spread(variable, value)  


# get names of variables included in model
variable.names = names(select(plot.data, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p,
                              -FMons, -ARons, -CAons, -Antidepressants1, -intercept, -Age, -GABA_agonists1,
                              -time_since_stroke, -HADS_T, -OutsOut))
measure.names =c("FMons", "ARons", "CAons")
model.names=c( 'OutsOut')


# plot top panel
top = plot.data %>%
  ggplot(aes(specification, estimate, color = significant.p, size=significant.p)) +
  geom_point(shape = "|") +scale_size_manual(values=c(2,4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(limits=c("no", "yes"), values = c("black", "grey60")) +
  labs(x = "", y = "HADS score\n regression coefficient\n") + 
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(0,50,100))+
  ylim(-0.17,0.02)+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# set plotting order for variables based on number of times it's included in better fitting models

order.vars=data.frame(variable=rbind("Gender","DomAffectedY", "NFI", "NACA_antiepileptic1"), order=rbind(4,3,2,1))
order.meas=data.frame(variable=rbind("FMons", "ARons", "CAons"), order=rbind(3,2,1))


# rename variables and plot middle panel
middle1 = plot.data %>%
  gather(variable, value, eval(measure.names)) %>% 
  left_join(., order.meas, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "FMons", "Fugl-Meyer", 
                           ifelse(variable == "ARons", "ARAT",
                                  ifelse(variable == "CAons", "CAHAI", variable)))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  labs( y = "measures\n") + 
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

middle2 = plot.data %>%
  gather(variable, value, eval(model.names)) %>% 
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "OutsOut", "Outliers Removed", variable)) %>%
  ggplot(aes(specification, variable, color = significant.p, size= significant.p)) +
  scale_size_manual(limits=c("no", "yes"),values=c(2,4))+
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  theme_minimal(base_size = 14) +geom_text(aes(label = value)) +
  labs( y = " ") + 
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# rename variables and plot bottom panel
bottom = plot.data %>%
  gather(variable, value, eval(variable.names)) %>% 
  left_join(., order.vars, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "GenderM", "Demographic Info", 
                           ifelse(variable == "NFI", "NFI score", 
                                  ifelse(variable == "NACA_antiepileptic1", "All Drugs",
                                         ifelse(variable == "DomAffectedY", "Stroke Info",
                                                variable))))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  labs(x = "\nmodel number", y = "covariates\n") + 
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(0,50,100))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# join panels
(AdSCA = cowplot::plot_grid(top, middle1, middle2,  bottom, ncol = 1, align = "v", 
                            rel_heights = c(1,0.30, 0.20, 0.60)))


# Permutation testing ----

# For reproducability 
set.seed(1992)

# No. of permutations
p<-500
n<-nrow(alldata4)
list<-seq(1,n,1)

# Makes a progress bar to display how were getting on 
pb = txtProgressBar(min = 1, max = p, initial = 0) 

options(na.action = "na.fail")

# Initialise a matrix to hold all the permutation data 
permlists<-matrix(0, nrow=n, ncol=p)

# take samples 
for (perms in 1:p){
  ##Permutation test a list of numbers, then reorder all the different data by the permuted list...
  permlists[,perms]<-sample(list, size=n, replace=FALSE)
}


# run full model with permutations
for (perms in 1:p){
  outcomedata<-arrange(cbind(alldata4[c("FM1_tot_pc", "ARAT1_tot_pc", "CAHAI1_tot_pc")],
                             permlists[,perms]), permlists[,perms])
  permdata<-cbind(select(alldata4, -FM1_tot_pc, -ARAT1_tot_pc, -CAHAI1_tot_pc),outcomedata)
  
  for (meas in 1:3){
    for (outliers in 1:2){
      if (outliers==1){this_data<-permdata
      
      }else{
        this_data<- filter(permdata, !(abs(permdata$FM1.FM4_pc- median(permdata$FM1.FM4_pc)) > 2.5*IQR(permdata$FM1.FM4_pc))&
                             !(abs(permdata$AR1.AR4_pc- median(permdata$AR1.AR4_pc)) > 2.5*IQR(permdata$AR1.AR4_pc))&
                             !(abs(permdata$CA1.CA4_pc- median(permdata$CA1.CA4_pc)) > 2.5*IQR(permdata$CA1.CA4_pc)))}
      
      if (meas==1){ y.var<-this_data$FM1_tot_pc
      } else if (meas==2){y.var<-this_data$ARAT1_tot_pc
      } else {y.var<-this_data$CAHAI1_tot_pc}
      
      
      full.model=lm(y.var ~ Age+ time_since_stroke+ DomAffected+  Gender + NFI+ 
                      HADS_T + Antidepressants + NACA_antiepileptic + GABA_agonists, data= this_data)    
      
      
      quiet((
        models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                                 subset=(((Age&Gender)|(!Age&!Gender))&
                                           ((DomAffected&time_since_stroke)|(!DomAffected&!time_since_stroke))&
                                           ((Antidepressants&GABA_agonists&NACA_antiepileptic)|(!Antidepressants&!GABA_agonists&!NACA_antiepileptic)))))
        , all=TRUE)
      
      #run all possible nested models + filtering for only ones I want
      if (meas==1){
        FMons<- as.numeric(rep(1, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else if (meas==2){
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(1, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else {   
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(1, times=nrow(models.sca)))
      }
      
      if(outliers==2){
        OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
      }  else {
        OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
      }
      
      
      # extract parameter estimate
      (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
          tibble() %>%
          rename("model" = ".") %>%
          mutate(tidied = purrr::map(model, broom::tidy),
                 model_num = row_number()) %>%
          select(model_num, tidied) %>%
          unnest() %>%
          filter(term == paste("HADS_T", sep="")) %>%
          ungroup() %>%
          select(model_num, estimate, std.error, p.value))
      
      model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*2)+(outliers-1)))
      
      if (outliers==1&meas==1){
        model_ps_all_perms<-model_ps
      }else{
        model_ps_all_perms<-rbind(model_ps_all_perms, model_ps)
      }
      
    }} #end measure, model and outlier loops
  
  if (perms==1){
    perms_ps<-model_ps_all_perms$p.value
  }else{
    perms_ps<-cbind(perms_ps, model_ps_all_perms$p.value)
  }
  setTxtProgressBar(pb,perms)
  
}

perm_sig<-(perms_ps<=0.05) %>%
  colSums

print('P-value for Admission as caclulated by permutation testing:')
print(drug_var_n)
print(pval<-sum(perm_sig>=test_sig)/p)



# Improvement model-----
for (meas in 1:3){
  for (modtype in 1:3){
    for (outliers in 1:2){
      if (outliers==1){this_data<-alldata4
      }else{
        this_data<- filter(alldata4, !(abs(alldata4$FM1.FM4_pc- median(alldata4$FM1.FM4_pc)) > 2.5*IQR(alldata4$FM1.FM4_pc))&
                             !(abs(alldata4$AR1.AR4_pc- median(alldata4$AR1.AR4_pc)) > 2.5*IQR(alldata4$AR1.AR4_pc))&
                             !(abs(alldata4$CA1.CA4_pc- median(alldata4$CA1.CA4_pc)) > 2.5*IQR(alldata4$CA1.CA4_pc)))}
      
      if (meas==1){ 
        if (modtype==1){  
          cov.var<-this_data$FM1_tot_pc
          y.var<-this_data$FM4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relFM1.FM4_pc
        } else {
          cov.var<-1
          y.var<-this_data$FM1.FM4_pc}
      } else if (meas==2){
        if (modtype==1){ 
          cov.var<-this_data$ARAT1_tot_pc
          y.var<-this_data$ARAT4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relAR1.AR4_pc
        } else {
          cov.var<-1
          y.var<-this_data$AR1.AR4_pc}
      } else {
        if (modtype==1){ 
          cov.var<-this_data$CAHAI1_tot_pc
          y.var<-this_data$CAHAI4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relCA1.CA4_pc
        } else {
          cov.var<-1
          y.var<-this_data$CA1.CA4_pc}}
      
      this_data<-cbind(this_data,y.var,cov.var)
      full.model=lm(y.var ~ cov.var + Age+ time_since_stroke+ DomAffected+ Gender + NFI+ 
                      HADS_T + Antidepressants + NACA_antiepileptic + GABA_agonists, data= this_data)    
      
      models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                               subset=(((Age&Gender)|(!Age&!Gender))&cov.var&
                                         ((DomAffected&time_since_stroke)|(!DomAffected&!time_since_stroke))&
                                         ((Antidepressants&GABA_agonists&NACA_antiepileptic)|(!Antidepressants&!GABA_agonists&!NACA_antiepileptic))))
      
      
      #run all possible nested models + filtering for only ones I want
      if (meas==1){
        FMons<- as.numeric(rep(1, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else if (meas==2){
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(1, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else {   
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(1, times=nrow(models.sca)))
      }
      
      if(outliers==2){
        OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
      }  else {
        OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
      }
      
      if(modtype==1){
        ChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(1, times=nrow(models.sca)))
      }  else if (modtype==2){
        ChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(1, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
      } else{ 
        ChangeMod<- as.numeric(rep(1, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(NA, times=nrow(models.sca)))}
      
      
      # extract parameter estimate
      (model_params = MuMIn::get.models(models.sca, subset = TRUE) %>%
          tibble() %>%
          dplyr::rename("model" = ".") %>%
          mutate(tidied = purrr::map(model, broom::tidy),
                 model_num = row_number()) %>%
          select(model_num, tidied) %>%
          unnest() %>%
          select(model_num, term, estimate) %>%
          spread(term, estimate) %>%
          select(-starts_with("sd")))
      
      model_params$model_num<-model_params$model_num+(nrow(models.sca)*(((meas-1)*3*2)+((modtype-1)*2)+(outliers-1)))
      model_params<-cbind(model_params,FMons, ARons, CAons, OutsOut, ChangeMod, RelChangeMod, OutcomeMod)
      
      (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
          tibble() %>%
          dplyr::rename("model" = ".") %>%
          mutate(tidied = purrr::map(model, broom::tidy),
                 model_num = row_number()) %>%
          select(model_num, tidied) %>%
          unnest() %>%
          filter(term == paste("HADS_T", sep="")) %>%
          ungroup() %>%
          select(model_num, estimate, std.error, p.value))
      
      model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*3*2)+((modtype-1)*2)+(outliers-1)))
      
      if (outliers==1&meas==1&modtype==1){
        model_params_all<-model_params
        model_ps_all<-model_ps
      }else{
        model_params_all<-full_join(model_params_all, model_params)
        model_ps_all<-full_join(model_ps_all, model_ps)
      }}}}

print(test_sig<-(model_ps_all$p.value<=0.05) %>%
        sum())
summary(model_ps_all$estimate[model_ps_all$p.value<=0.05])
summary(model_ps_all$estimate)



# Plotting the outcomes ----
# merge and tidy for plotting
plot.data= left_join(model_ps_all, model_params_all, by = "model_num") %>%
  arrange(estimate) %>%
  mutate(specification = row_number(),
         significant.p = ifelse(p.value < .05, "yes", "no")) %>%
  gather(variable, value, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p) %>% 
  mutate(variable = gsub("[()]", "", variable),
         variable = gsub("Intercept", "intercept", variable)) %>%
  spread(variable, value)  


# get names of variables included in model
variable.names = names(select(plot.data, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p,
                              -FMons, -ARons, -CAons, -cov.var, -Antidepressants1, -intercept, -Age, -GABA_agonists1,
                              -time_since_stroke, -HADS_T, -OutsOut, -OutcomeMod, -RelChangeMod, -ChangeMod))
measure.names =c("FMons", "ARons", "CAons")
model.names=c('ChangeMod', 'RelChangeMod', 'OutcomeMod')
out.name='OutsOut'

# plot top panel
top = plot.data %>%
  ggplot(aes(specification, estimate, color = significant.p, size=significant.p)) +
  geom_point(shape = "|") +scale_size_manual(values=c(2,4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(limits=c("no", "yes"), values = c("black", "grey60")) +
  labs(x = "") + 
  theme_minimal(base_size = 14) +
  ylim(-0.17,0.02)+
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# set plotting order for variables based on number of times it's included in better fitting models

order.vars=data.frame(variable=rbind("Gender","DomAffectedY", "NFI", "NACA_antiepileptic1"), order=rbind(4,3,2,1))
order.meas=data.frame(variable=rbind("FMons", "ARons", "CAons"), order=rbind(3,2,1))
order.mod=data.frame(variable=rbind("OutcomeMod", "ChangeMod", "RelChangeMod"), order=rbind(1,2,3))


# rename variables and plot middle panel
middle1 = plot.data %>%
  gather(variable, value, eval(measure.names)) %>% 
  left_join(., order.meas, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "FMons", "Fugl-Meyer", 
                           ifelse(variable == "ARons", "ARAT",
                                  ifelse(variable == "CAons", "CAHAI", variable)))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

middle2 = plot.data %>%
  gather(variable, value, eval(out.name)) %>% 
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "OutsOut", "Outliers Removed", variable)) %>%
  ggplot(aes(specification, variable, color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

middle3 = plot.data %>%
  gather(variable, value, eval(model.names)) %>% 
  left_join(., order.mod, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "ChangeMod", "Abs Recovery Model",
                           ifelse(variable == "OutcomeMod", "Outcome Model",
                                  ifelse(variable == "RelChangeMod", " Rel Recovery Model", variable)))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black","grey60")) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# rename variables and plot bottom panel
bottom = plot.data %>%
  gather(variable, value, eval(variable.names)) %>% 
  left_join(., order.vars, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "GenderM", "Demographic Info", 
                           ifelse(variable == "NFI", "NFI score", 
                                  ifelse(variable == "NACA_antiepileptic1", "All Drugs",
                                         ifelse(variable == "DomAffectedY", "Stroke Info",
                                                variable))))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  labs(x = "\nmodel number") + 
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# join panels
(ImpSCA = cowplot::plot_grid(top, middle1, middle2,middle3,  bottom, ncol = 1, align = "v", 
                             rel_heights = c(1,0.245,0.125, 0.245, 0.485)))



# Plot Both -------

(BothSCA=cowplot::plot_grid(AdSCA, ImpSCA, ncol=2, align = 'h', rel_widths = c(1.9,4),
                            labels="AUTO",label_size = 20))
# END -----


##### FOR REVIEWER- TSS SCA ##### ####-----
# Admission model-----
for (meas in 1:3){
  for (outliers in 1:2){
    if (outliers==1){this_data<-alldata4
    }else{
      this_data<- filter(alldata4, !(abs(alldata4$FM1.FM4_pc- median(alldata4$FM1.FM4_pc)) > 2.5*IQR(alldata4$FM1.FM4_pc))&
                           !(abs(alldata4$AR1.AR4_pc- median(alldata4$AR1.AR4_pc)) > 2.5*IQR(alldata4$AR1.AR4_pc))&
                           !(abs(alldata4$CA1.CA4_pc- median(alldata4$CA1.CA4_pc)) > 2.5*IQR(alldata4$CA1.CA4_pc)))}
    
    if (meas==1){ y.var<-this_data$FM1_tot_pc}
    else if (meas==2){y.var<-this_data$ARAT1_tot_pc}
    else {y.var<-this_data$CAHAI1_tot_pc}
    
    
    full.model=lm(y.var ~ Age+ time_since_stroke+ DomAffected+  Gender + NFI+ 
                    HADS_T + Antidepressants + GABA_agonists + NACA_antiepileptic, data= this_data)    
    
    models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                             subset=(((Age&Gender)|(!Age&!Gender))&
                                       ((Antidepressants&GABA_agonists&NACA_antiepileptic)|(!Antidepressants&!GABA_agonists&!NACA_antiepileptic))&
                                       ((HADS_T&NFI)|(!HADS_T&!NFI))))
    
    
    #run all possible nested models + filtering for only ones I want
    if (meas==1){
      FMons<- as.numeric(rep(1, times=nrow(models.sca)))
      ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
      CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
    }else if (meas==2){
      FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
      ARons<- as.numeric(rep(1, times=nrow(models.sca)))
      CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
    }else {   
      FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
      ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
      CAons<- as.numeric(rep(1, times=nrow(models.sca)))
    }
    
    if(outliers==2){
      OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
    }  else {
      OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
    }
    
    
    
    # extract parameter estimate
    (model_params = MuMIn::get.models(models.sca, subset = TRUE) %>%
        tibble() %>%
        dplyr::rename("model" = ".") %>%
        mutate(tidied = purrr::map(model, broom::tidy),
               model_num = row_number()) %>%
        select(model_num, tidied) %>%
        unnest() %>%
        select(model_num, term, estimate) %>%
        spread(term, estimate)) %>%
      select(-starts_with("sd"))
    
    model_params$model_num<-model_params$model_num+(nrow(models.sca)*(((meas-1)*2)+(outliers-1)))
    model_params<-cbind(model_params,FMons, ARons, CAons, OutsOut)
    
    (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
        tibble() %>%
        dplyr::rename("model" = ".") %>%
        mutate(tidied = purrr::map(model, broom::tidy),
               model_num = row_number()) %>%
        select(model_num, tidied) %>%
        unnest() %>%
        filter(term ==  paste("time_since_stroke", sep="")) %>%
        ungroup() %>%
        select(model_num, estimate, std.error, p.value))
    
    model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*2)+(outliers-1)))
    
    if (outliers==1&meas==1){
      model_params_all<-model_params
      model_ps_all<-model_ps
    }else{
      model_params_all<-full_join(model_params_all, model_params)
      model_ps_all<-full_join(model_ps_all, model_ps)
    }
  }}

print(test_sig<-(model_ps_all$p.value<=0.05) %>%
        sum())
summary(model_ps_all$estimate[model_ps_all$p.value<=0.05])
summary(model_ps_all$estimate)

# Plotting the outcomes ----
# merge and tidy for plotting
plot.data= left_join(model_ps_all, model_params_all, by = "model_num") %>%
  arrange(estimate) %>%
  mutate(specification = row_number(),
         significant.p = ifelse(p.value < .05, "yes", "no")) %>%
  gather(variable, value, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p) %>% 
  mutate(variable = gsub("[()]", "", variable),
         variable = gsub("Intercept", "intercept", variable)) %>%
  spread(variable, value)  


# get names of variables included in model
variable.names = names(select(plot.data, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p,
                              -FMons, -ARons, -CAons, -Antidepressants1, -Age -intercept, -GABA_agonists1,
                              -time_since_stroke, -HADS_T, -OutsOut))
measure.names =c("FMons", "ARons", "CAons")
model.names=c( 'OutsOut')


# plot top panel
top = plot.data %>%
  ggplot(aes(specification, estimate, color = significant.p, size=significant.p)) +
  geom_point(shape = "|") +scale_size_manual(values=c(2,4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(limits=c("no", "yes"), values = c("black", "grey60")) +
  labs(x = "", y = "regression coefficient\n") + 
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(0,50,100))+
  ylim(-0.17,0.02)+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# set plotting order for variables based on number of times it's included in better fitting models

order.vars=data.frame(variable=rbind("Gender","DomAffectedY", "NFI", "NACA_antiepileptic1"), order=rbind(4,3,2,1))
order.meas=data.frame(variable=rbind("FMons", "ARons", "CAons"), order=rbind(3,2,1))


# rename variables and plot middle panel
middle1 = plot.data %>%
  gather(variable, value, eval(measure.names)) %>% 
  left_join(., order.meas, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "FMons", "Fugl-Meyer", 
                           ifelse(variable == "ARons", "ARAT",
                                  ifelse(variable == "CAons", "CAHAI", variable)))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  labs( y = "measures\n") + 
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

middle2 = plot.data %>%
  gather(variable, value, eval(model.names)) %>% 
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "OutsOut", "Outliers Removed", variable)) %>%
  ggplot(aes(specification, variable, color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  theme_minimal(base_size = 14) +
  labs( y = " ") + 
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# rename variables and plot bottom panel
bottom = plot.data %>%
  gather(variable, value, eval(variable.names)) %>% 
  left_join(., order.vars, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "GenderM", "Demographic Info", 
                           ifelse(variable == "NFI", "NFI score", 
                                  ifelse(variable == "NACA_antiepileptic1", "All Drugs",
                                         ifelse(variable == "DomAffectedY", "Arm Affected",
                                                variable))))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  labs(x = "\nmodel number") + 
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# join panels
(AdSCA = cowplot::plot_grid(top, middle1, middle2,  bottom, ncol = 1, align = "v", 
                            rel_heights = c(1,0.30, 0.20, 0.60)))


# Improvement model-----
for (meas in 1:3){
  for (modtype in 1:3){
    for (outliers in 1:2){
      if (outliers==1){this_data<-alldata4
      }else{
        this_data<- filter(alldata4, !(abs(alldata4$FM1.FM4_pc- median(alldata4$FM1.FM4_pc)) > 2.5*IQR(alldata4$FM1.FM4_pc))&
                             !(abs(alldata4$AR1.AR4_pc- median(alldata4$AR1.AR4_pc)) > 2.5*IQR(alldata4$AR1.AR4_pc))&
                             !(abs(alldata4$CA1.CA4_pc- median(alldata4$CA1.CA4_pc)) > 2.5*IQR(alldata4$CA1.CA4_pc)))}
      
      if (meas==1){ 
        if (modtype==1){  
          cov.var<-this_data$FM1_tot_pc
          y.var<-this_data$FM4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relFM1.FM4_pc
        } else {
          cov.var<-1
          y.var<-this_data$FM1.FM4_pc}
      } else if (meas==2){
        if (modtype==1){ 
          cov.var<-this_data$ARAT1_tot_pc
          y.var<-this_data$ARAT4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relAR1.AR4_pc
        } else {
          cov.var<-1
          y.var<-this_data$AR1.AR4_pc}
      } else {
        if (modtype==1){ 
          cov.var<-this_data$CAHAI1_tot_pc
          y.var<-this_data$CAHAI4_tot_pc
        } else if (modtype==2) { 
          cov.var<-1
          y.var<-this_data$relCA1.CA4_pc
        } else {
          cov.var<-1
          y.var<-this_data$CA1.CA4_pc}}
      
      this_data<-cbind(this_data,y.var,cov.var)
      full.model=lm(y.var ~ cov.var + Age+ time_since_stroke+ DomAffected+ Gender + NFI+ 
                      HADS_T + Antidepressants + NACA_antiepileptic + GABA_agonists, data= this_data)    
      
      models.sca=MuMIn::dredge(full.model, rank='AIC', extra='BIC', 
                               subset=(((Age&Gender)|(!Age&!Gender))&cov.var&
                                         ((HADS_T&NFI)|(!HADS_T&!NFI))&
                                         ((Antidepressants&GABA_agonists&NACA_antiepileptic)|(!Antidepressants&!GABA_agonists&!NACA_antiepileptic))))
      
      
      #run all possible nested models + filtering for only ones I want
      if (meas==1){
        FMons<- as.numeric(rep(1, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else if (meas==2){
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(1, times=nrow(models.sca)))
        CAons<- as.numeric(rep(NA, times=nrow(models.sca)))
      }else {   
        FMons<- as.numeric(rep(NA, times=nrow(models.sca)))
        ARons<- as.numeric(rep(NA, times=nrow(models.sca)))
        CAons<- as.numeric(rep(1, times=nrow(models.sca)))
      }
      
      if(outliers==2){
        OutsOut<- as.numeric(rep(1, times=nrow(models.sca)))
      }  else {
        OutsOut<-as.numeric(rep(NA, times=nrow(models.sca)))
      }
      
      if(modtype==1){
        ChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(1, times=nrow(models.sca)))
      }  else if (modtype==2){
        ChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(1, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
      } else{ 
        ChangeMod<- as.numeric(rep(1, times=nrow(models.sca)))
        RelChangeMod<- as.numeric(rep(NA, times=nrow(models.sca)))
        OutcomeMod<- as.numeric(rep(NA, times=nrow(models.sca)))}
      
      
      # extract parameter estimate
      (model_params = MuMIn::get.models(models.sca, subset = TRUE) %>%
          tibble() %>%
          dplyr::rename("model" = ".") %>%
          mutate(tidied = purrr::map(model, broom::tidy),
                 model_num = row_number()) %>%
          select(model_num, tidied) %>%
          unnest() %>%
          select(model_num, term, estimate) %>%
          spread(term, estimate) %>%
          select(-starts_with("sd")))
      
      model_params$model_num<-model_params$model_num+(nrow(models.sca)*(((meas-1)*3*2)+((modtype-1)*2)+(outliers-1)))
      model_params<-cbind(model_params,FMons, ARons, CAons, OutsOut, ChangeMod, RelChangeMod, OutcomeMod)
      
      (model_ps = MuMIn::get.models(models.sca, subset = TRUE) %>%
          tibble() %>%
          dplyr::rename("model" = ".") %>%
          mutate(tidied = purrr::map(model, broom::tidy),
                 model_num = row_number()) %>%
          select(model_num, tidied) %>%
          unnest() %>%
          filter(term == paste("time_since_stroke", sep="")) %>%
          ungroup() %>%
          select(model_num, estimate, std.error, p.value))
      
      model_ps$model_num<-model_ps$model_num+(nrow(models.sca)*(((meas-1)*3*2)+((modtype-1)*2)+(outliers-1)))
      
      if (outliers==1&meas==1&modtype==1){
        model_params_all<-model_params
        model_ps_all<-model_ps
      }else{
        model_params_all<-full_join(model_params_all, model_params)
        model_ps_all<-full_join(model_ps_all, model_ps)
      }}}}

print(test_sig<-(model_ps_all$p.value<=0.05) %>%
        sum())
summary(model_ps_all$estimate[model_ps_all$p.value<=0.05])
summary(model_ps_all$estimate)



# Plotting the outcomes ----
# merge and tidy for plotting
plot.data= left_join(model_ps_all, model_params_all, by = "model_num") %>%
  arrange(estimate) %>%
  mutate(specification = row_number(),
         significant.p = ifelse(p.value < .05, "yes", "no")) %>%
  gather(variable, value, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p) %>% 
  mutate(variable = gsub("[()]", "", variable),
         variable = gsub("Intercept", "intercept", variable)) %>%
  spread(variable, value)  


# get names of variables included in model
variable.names = names(select(plot.data, -estimate, -specification, -model_num, -std.error, -p.value, -significant.p,
                              -FMons, -ARons, -CAons, -cov.var, -Antidepressants1, -intercept, -Age, -GABA_agonists1,
                              -time_since_stroke, -HADS_T, -OutsOut, -OutcomeMod, -RelChangeMod, -ChangeMod))
measure.names =c("FMons", "ARons", "CAons")
model.names=c('ChangeMod', 'RelChangeMod', 'OutcomeMod')
out.name='OutsOut'

# plot top panel
top = plot.data %>%
  ggplot(aes(specification, estimate, color = significant.p, size=significant.p)) +
  geom_point(shape = "|") +scale_size_manual(values=c(2,4)) +
  geom_hline(yintercept = 0, linetype = "dashed", color = "grey") +
  scale_color_manual(limits=c("no", "yes"), values = c("black", "grey60")) +
  labs(x = "") + 
  theme_minimal(base_size = 14) +
  ylim(-0.17,0.02)+
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# set plotting order for variables based on number of times it's included in better fitting models

order.vars=data.frame(variable=rbind("Gender","DomAffectedY", "NFI", "NACA_antiepileptic1"), order=rbind(4,3,2,1))
order.meas=data.frame(variable=rbind("FMons", "ARons", "CAons"), order=rbind(3,2,1))
order.mod=data.frame(variable=rbind("OutcomeMod", "ChangeMod", "RelChangeMod"), order=rbind(1,2,3))


# rename variables and plot middle panel
middle1 = plot.data %>%
  gather(variable, value, eval(measure.names)) %>% 
  left_join(., order.meas, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "FMons", "Fugl-Meyer", 
                           ifelse(variable == "ARons", "ARAT",
                                  ifelse(variable == "CAons", "CAHAI", variable)))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

middle2 = plot.data %>%
  gather(variable, value, eval(out.name)) %>% 
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "OutsOut", "Outliers Removed", variable)) %>%
  ggplot(aes(specification, variable, color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  theme_minimal(base_size = 14) +
  theme(axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.title.y = element_blank(),
        axis.ticks.x = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

middle3 = plot.data %>%
  gather(variable, value, eval(model.names)) %>% 
  left_join(., order.mod, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "ChangeMod", "Abs Recovery Model",
                           ifelse(variable == "OutcomeMod", "Outcome Model",
                                  ifelse(variable == "RelChangeMod", " Rel Recovery Model", variable)))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black","grey60")) +
  theme_minimal(base_size = 14) +
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.line.x = element_blank(),
        axis.text.x = element_blank(),
        axis.title.x = element_blank(),
        axis.ticks.x = element_blank(),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())

# rename variables and plot bottom panel
bottom = plot.data %>%
  gather(variable, value, eval(variable.names)) %>% 
  left_join(., order.vars, by = "variable") %>%
  mutate(value = ifelse(!is.na(value),"|", ""),
         variable = ifelse(variable == "GenderM", "Demographic Info", 
                           ifelse(variable == "NFI", "NFI score", 
                                  ifelse(variable == "NACA_antiepileptic1", "All Drugs",
                                         ifelse(variable == "DomAffectedY", "Arm Affected",
                                                variable))))) %>%
  ggplot(aes(specification, reorder(variable, order), color = significant.p, size= significant.p)) +
  scale_size_manual(values=c(2,4))+  geom_text(aes(label = value)) +
  scale_color_manual(limits=c("no", "yes"),values = c("black", "grey60")) +
  labs(x = "\nmodel number") + 
  theme_minimal(base_size = 14) +
  scale_x_continuous(breaks = c(0,50,100,150,200,250,300))+
  theme(legend.title = element_text(size = 10),
        legend.text = element_text(size = 9),
        axis.text = element_text(color = "black"),
        axis.line = element_line(colour = "black"),
        axis.title.y = element_blank(),
        legend.position = "none",
        panel.grid.major = element_blank(),
        panel.grid.minor = element_blank(),
        panel.border = element_blank(),
        panel.background = element_blank())


# join panels
(ImpSCA = cowplot::plot_grid(top, middle1, middle2,middle3,  bottom, ncol = 1, align = "v", 
                             rel_heights = c(1,0.245,0.125, 0.245, 0.485)))



# Plot Both -------

(BothSCA=cowplot::plot_grid(AdSCA, ImpSCA, ncol=2, align = 'h', rel_widths = c(1.9,4),
                            labels="AUTO",label_size = 20))
