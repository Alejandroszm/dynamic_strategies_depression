library(data.table)
setDTthreads(2)
library(nnet)
library(splitstackshape)
library(sandwich)
library(lmtest)
library(foreign)
library(nnet)
library(lmtest)
library(sandwich)
library(splines)

#source("//DON/allprojs/NMIH/programs/aszmulewicz/final_dissert/2021/primary outcome/DataCleaning_24aug21.R")

bootstrap_fn <- function(data0, num.boot){
  
  
  
  set.seed(12345)
  seed <- floor(runif(num.boot)*10^8)
  
  
  # Create an empty dataset that will eventually store all of our
  # results from the bootstrapped samples.  
  boot.results <- data.frame(sequential_monotherapy=rep(NA, 1+num.boot),
                             sequential_dual=rep(NA, 1+num.boot),
                             guidelines_based=rep(NA, 1+num.boot),
                             survdiff1=rep(NA, 1+num.boot),
                             survdiff2=rep(NA, 1+num.boot),
                             rrisk1=rep(NA, 1+num.boot),
                             rrisk2=rep(NA, 1+num.boot),
                             nnt1 = rep(NA, 1+num.boot),
                             nnt2 = rep(NA, 1+num.boot)
  )
  
  
  for (i in 0:num.boot){
    if (i==0) {
      # use all data for this
      finaldata <- copy(data0)
    } else {
      # set the seed for this bootstrap sample
      set.seed(seed[i])
      ids <- as.data.table(sample(1:data_len, data_len, replace = TRUE))
      ids[, 'bid' := 1:data_len]
      colnames(ids) <- c("id", "bid")
      l <- c('id')
      setkeyv(data0, l)
      finaldata <- data0[J(ids), allow.cartesian = TRUE][,id:=bid][,bid:=NULL ]  # create the new data set names "sample"
    }
    
    
    # 1.A.) Probability of staying on the same treatment: this model should be fitted to clinical visits
    # before week 8 without intolerance (ie, scheduled clinical visits)
    # a final pre-processing first.
    finaldata<-finaldata %>% mutate(unemployed = ifelse(empl == "Unemployed", 1, 0),
                                    employed = ifelse(empl == "Employed", 1, 0)) 
    finaldata$days_study2 = finaldata$days_study*finaldata$days_study
    
    model1<-glm(formula = as.factor(changetrt==0) ~ days_study + days_study2 + unemployed + idsctot+ I(idsctot*idsctot)+
                  age + I(age*age) + female + employed + melancholic + atypical + anxiousdepr + comorbidptsd +
                  baseline.hrsd + I(baseline.hrsd*baseline.hrsd) + self.qids + I(self.qids*self.qids) + sexual +
                  qids.fup + I(qids.fup*qids.fup),
                data = subset(finaldata, days_study=!0 & 
                                ((failure==1|failure==2) & visitweek8==0 & visit==1 & intol2se==0 & intol3se==0)|
                                (visitweek8==1 & visit==1 & intol2se==0 & intol3se==0)),
                family = quasibinomial())
    
    
    finaldata$pstay1 <- predict(model1, newdata = finaldata, type = 'response') # Pr of staying. 
    
    
    #1.B) Probability of staying on the same treatment after 8 weeks: this model should be fitted to 
    # clinical visits after 8 weeks without intolerance and with inadequate benefit 
    # (ie, scheduled clinical visits)
    
    model2<-glm(formula = as.factor(changetrt==0) ~ days_study + days_study2 + unemployed + idsctot+ I(idsctot*idsctot)+
                  age + I(age*age)+ female + employed + melancholic + atypical + anxiousdepr + comorbidptsd +
                  baseline.hrsd + I(baseline.hrsd*baseline.hrsd) + self.qids + I(self.qids*self.qids) + sexual +
                  qids.fup + I(qids.fup*qids.fup),
                data = subset(finaldata, visit==1 & visitweek8==2 & 
                                intol2se==0 & intol3se==0 & after_accept_benef==0),
                family = quasibinomial())
    
    finaldata$pstay2 <- predict(model2, newdata = finaldata, type = 'response') # Pr of staying
    
    # 1.B.) Probability of receiving, after intolerance to a second antidepressant either MIRT, NTP, or Li/T4
    model3<-multinom(secondtrt ~ days_study + days_study2 + unemployed + idsctot+I(idsctot*idsctot)+
                       age + I(age*age)+ female + employed + melancholic + atypical + anxiousdepr + comorbidptsd +
                       baseline.hrsd + I(baseline.hrsd*baseline.hrsd) + self.qids + I(self.qids*self.qids) + sexual + 
                       qids.fup + I(qids.fup*qids.fup),
                     data = subset(finaldata, failure == 2 & changetrt==1 & intol2se==1 & visit==1),
                     family = quasibinomial(), maxit = 1000)
    
    df2<- data.frame(predict(model3, type="probs", newdata=finaldata)) %>%  
      select (X1, X2, X3) %>% mutate(rownumber = seq.int(nrow(.)))
    
    
    colnames(df2)<-c("Prob_MIRT_Intol", "Prob_NTP_Intol", "Prob_Li/T4_Intol", "rownumber")
    
    finaldata$rownumber <- seq.int(nrow(finaldata)) 
    
    finaldata<-merge(finaldata, df2, by ="rownumber") %>% select(-c(rownumber))
    
    # 1.C.) Probability of receiving, after non-remission to a second antidepressant either MIRT, NTP, or Li/T4
    model4<-multinom(secondtrt ~ days_study + days_study2 + unemployed + idsctot+I(idsctot*idsctot)+
                       age + I(age*age)+ female + employed + melancholic + atypical + anxiousdepr + comorbidptsd +
                       baseline.hrsd + I(baseline.hrsd*baseline.hrsd) + self.qids + I(self.qids*self.qids) + sexual + 
                       qids.fup + I(qids.fup*qids.fup),
                     data = subset(finaldata,  failure==2 & visit==1 & changetrt==1 & intol2se==0), maxit = 1000, 
                     family = quasibinomial())
    
    df3<- data.frame(predict(model4, type="probs", newdata=finaldata)) %>%  
      select (X1, X2, X3) %>% mutate(rownumber = seq.int(nrow(.)))
    
    colnames(df3)<-c("Prob_MIRT_NR", "Prob_NTP_NR", "Prob_Li/T4_NR", "rownumber")
    
    finaldata$rownumber <- seq.int(nrow(finaldata)) 
    
    finaldata<-merge(finaldata, df3, by ="rownumber") %>% select(-c(rownumber))
    
    # 2) prob of ltfu
    # define visit at exactly the day 30 since the last clinical appointment
    finaldata<-finaldata %>% group_by(id) %>% mutate(visit2 = lead(visit)) %>% ungroup()
    
    model.ltfu<-glm(formula = visit2==1 ~ days_study + days_study2 + idsctot+I(idsctot*idsctot)+
                      age + I(age*age)+ female + employed + unemployed + melancholic + atypical + anxiousdepr + comorbidptsd +
                      baseline.hrsd + I(baseline.hrsd*baseline.hrsd) + self.qids + I(self.qids*self.qids) + sexual + 
                      qids.fup + I(qids.fup*qids.fup),
                    data = subset(finaldata,  days_since_last==27 & after_accept_benef==0),
                    family = quasibinomial())
    
    finaldata$p_ltfu <- predict(model.ltfu, newdata = finaldata, type = 'response') # Pr of staying

    ####################################################################################################################    
    # next, I will create separate datasets, simulating everybody following a given strategy
    # incorporating relevant probabilities
    always_switch_data = finaldata
    
    always_switch_ids<-finaldata %>% group_by(id) %>% filter(row_number()==1,
                                                             txassign %in% c("VEN", "BUP", "SER"))
    
    always_switch_data<-always_switch_data %>% filter(id %in% always_switch_ids$id) 
    
    always_switch_data<-always_switch_data %>% mutate (censor_dev = case_when(
      failure==2 & (txassign %in% c("CIT+LI", "CIT+THY", "SER+LI", "BUP+LI","VEN+LI", "SER+THY", "BUP+THY",
                                    "VEN+THY")) ~ 1,
      failure == 3 & (txassign == "VEN+MIRT") ~ 1,
      TRUE ~ 0),
      strategy=1)
    
    # eliminate rows after the first deviation  
    always_switch_data <- always_switch_data %>%  group_by(id) %>% mutate (touseweights = lag(censor_dev))  
    always_switch_data$touseweights[is.na(always_switch_data$touseweights)]<- 0
    
    always_switch_data<-  always_switch_data %>% group_by(id) %>%
      mutate(touseweights = cumsum (touseweights)) %>%
      filter (touseweights == 0)  %>% ungroup() %>% select (-c(touseweights))
    
    always_switch_data<- always_switch_data %>% group_by(id) %>% arrange(days_study, .by_group = T)                                            
    
    #### DATASET 2: Sequential dual therapy with non-ssri
    always_dual_data = finaldata
    
    always_dual_ids<-finaldata %>% group_by(id) %>% filter(row_number()==1) %>% 
      filter((txassign=="CIT+BUP" & intol1se==0) | (txassign ==  "VEN" & intol1se==1) | 
               (txassign== "BUP"& intol1se==1))
    
    always_dual_data<-always_dual_data %>% filter(id %in% always_dual_ids$id) 
    
    always_dual_data<-always_dual_data %>% mutate (
      censor_dev = case_when(
        failure==2 & intol2se==1 & (txassign %in% c("CIT+LI","NTP", "BUP+LI", "VEN+LI", "CIT+THY", "BUP+THY", "VEN+THY")) ~ 1,
        failure==2 & intol2se==0 & (txassign %in% c("NTP")) ~ 1,
        failure==3 & intol3se==1 & txassign == "VEN+MIRT" ~ 1,
        failure==3 & intol3se==0 & txassign == "TCP" ~ 1,
        TRUE ~ 0),
      strategy = 2)
    
    # eliminate rows after the first deviation  
    always_dual_data <- always_dual_data %>%  group_by(id) %>% mutate (touseweights = lag(censor_dev))  
    always_dual_data$touseweights[is.na(always_dual_data$touseweights)]<- 0
    
    always_dual_data<-  always_dual_data %>% group_by(id) %>%
      mutate(touseweights = cumsum (touseweights)) %>%
      filter (touseweights == 0)  %>% ungroup() %>% select (-c(touseweights))
    
    always_dual_data<- always_dual_data %>% group_by(id) %>% arrange(days_study, .by_group = T)                                            
    
    #### DATASET 3: Guidelines-based sequential therapy
    guidelines_based_data = finaldata
    
    guidelines_ids<-finaldata %>% group_by(id) %>% filter(row_number()==1) %>%
      filter((txassign %in% c("SER", "CIT+BUP") & dep.improv2>= 25 & intol1se==1) | 
               (txassign == "CIT+BUP" & dep.improv2>= 25 & intol1se==0) |
               (txassign %in% c("VEN", "BUP") & dep.improv2< 25))
    
    guidelines_based_data<-guidelines_based_data %>% filter(id %in% guidelines_ids$id) 
    
    guidelines_based_data<-guidelines_based_data %>% mutate (
      censor_dev = case_when(
        failure==2 & (dep.improv2 < 25) & (txassign %in% c("MIRT","CIT+LI", "SER+LI", "BUP+LI", "SER+THY",
                                                           "VEN+LI", "CIT+THY", "BUP+THY", "VEN+THY")) ~ 1,
        failure==2 & dep.improv2 >=25 & intol2se==1 & txassign== "NTP" ~ 1,
        failure==2 & dep.improv2 >=25 & intol2se==0 & (txassign%in%c("NTP", "MIRT")) ~ 1,
        failure==3 & dep.improv2 >=25 & intol3se==0 & txassign=="TCP" ~ 1,
        failure==3 & dep.improv2 <25  & txassign=="VEN+MIRT" ~ 1,
        TRUE ~ 0),
      strategy=3)
    
    # eliminate rows after the first deviation  
    guidelines_based_data <- guidelines_based_data %>%  group_by(id) %>% mutate (touseweights = lag(censor_dev))  
    guidelines_based_data$touseweights[is.na(guidelines_based_data$touseweights)]<- 0
    
    guidelines_based_data<-  guidelines_based_data %>% group_by(id) %>%
      mutate(touseweights = cumsum (touseweights)) %>%
      filter (touseweights == 0)  %>% ungroup() %>% select (-c(touseweights))
    
    guidelines_based_data<- guidelines_based_data %>% group_by(id) %>% arrange(days_study, .by_group = T)                                            
    
    # remove datasets
    rm(always_switch_ids, always_dual_ids, guidelines_ids, df3, df2, model1, model2, model3, model4, model.ltfu)
    
    #####################################################################################################################################
    # STEP 8. Create weights for IPW (adjustment for time-varying confounding)
    #####################################################################################################################################
    
    ## start with the separated datasets (following the rows from eTable 3)
    
    ##############################
    ## 1) always switch data:
    ##############################
    
    # first, create a variable that defines adherence
    always_switch_data<-always_switch_data %>% 
      mutate(adherence = case_when(
        violation2 == 1 | censor_dev == 1  ~ 0,
        TRUE ~ 1)) 
    
    ###print(names(always_switch_data))
    ###return(0)
    # the following rows incorporate the information from eTable 3 (see paper)
    always_switch_data <- always_switch_data %>% mutate (
      denom_weight_adherence = case_when (
        visit==1 & visitweek8==0 & intol2se==0 & intol3se==0 & (failure==1|failure==2) ~ pstay1,
        visit==1 & visitweek8==0 & intol2se==1 & intol3se==0 & failure==2 ~ Prob_MIRT_Intol + Prob_NTP_Intol,
        visit==1 & visitweek8==0 & intol2se==0 & intol3se==1 & failure==3 ~ 0.5,
        visit==1 & visitweek8==1 & intol2se==1 & intol3se==0 & failure==2 ~ Prob_MIRT_Intol + Prob_NTP_Intol,
        visit==1 & visitweek8==1 & intol2se==0 & intol3se==1 & failure==3 ~ 0.5,
        visit==1 & visitweek8==1 & intol2se==0 & intol3se==0 & after_accept_benef==0 & (failure==1 | (failure==2 & changetrt==1)) ~ pstay1 + ((1-pstay1)*(Prob_MIRT_NR + Prob_NTP_NR)),
        visit==1 & visitweek8==1 & intol2se==0 & intol3se==0 & after_accept_benef==0 & (failure==3 | (failure==2 & changetrt==0)) ~ pstay1+((1-pstay1)*0.5),
        visit==1 & visitweek8==2 & intol2se==1 & intol3se==0 & failure==2 ~ Prob_MIRT_Intol + Prob_NTP_Intol,
        visit==1 & visitweek8==2 & intol2se==0 & intol3se==1 & failure==3 ~ 0.5,
        visit==1 & visitweek8==2 & intol2se==0 & intol3se==0 & after_accept_benef==0 & failure==2 & changetrt==1 ~ (1-pstay2)*(Prob_MIRT_NR+Prob_NTP_NR),
        visit==1 & visitweek8==2 & intol2se==0 & intol3se==0 & after_accept_benef==0 & failure==3 & changetrt==1 ~ (1-pstay2)*0.5,
        days_study==0 ~ 1,
        visit==0 ~ 1,
        failure==3 & changetrt==0 ~ 1,
        after_accept_benef!=0 & visitweek8!=0 & intol2se==0 & intol3se==0 ~ 1,
        TRUE ~ 1))
    
    # Numerator of unstabilized weights for adherence
    always_switch_data$k1_n <- ifelse(always_switch_data$adherence==1,1,0)
    
    # Denominator of unstabilized weights for adherence
    always_switch_data$k1_d <- with(always_switch_data, ave(denom_weight_adherence, id, FUN=cumprod))
    
    # Create  unstabilized weights for adherence
    #always_switch_data$unstabw<-with(always_switch_data, k1_n / k1_d)
    always_switch_data <- always_switch_data %>% mutate (unstabw = case_when (
      k1_n ==0 & k1_d ==0 ~ 0,
      TRUE ~ k1_n / k1_d))
    
    tau <- quantile(always_switch_data$unstabw, 0.995) # truncate weights
    always_switch_data$unstabw[always_switch_data$unstabw > tau] <- tau # replace those above 99.5% with the 99.5%
    
    summary(always_switch_data$unstabw)
    
    # weights for ltfu
    always_switch_data <- always_switch_data %>% mutate (denom_weight_ltfu = case_when (
      days_since_last==27 & after_accept_benef==0 & visit2==1  ~ p_ltfu,
      TRUE ~ 1))
    
    # Numerator of unstabilized weights for ltfu
    always_switch_data$k2_n <- ifelse(always_switch_data$censoring_event==0,1,0)
    
    # Denominator of unstabilized weights for ltfu
    always_switch_data$k2_d <- with(always_switch_data, ave(denom_weight_ltfu, id, FUN=cumprod))
    
    # Create  unstabilized weights for adherence
    always_switch_data$unstabw2<-with(always_switch_data, k2_n / k2_d)
    
    tau2 <- quantile(always_switch_data$unstabw2, 0.995) # truncate weights
    always_switch_data$unstabw2[always_switch_data$unstabw2 > tau2] <- tau2 # replace those above 99.5% with the 99.5%
    
    summary(always_switch_data$unstabw2)
    
    # multiplication of weights
    always_switch_data$final_weight = always_switch_data$unstabw * always_switch_data$unstabw2
    
    #add this line if running more than 300 iterations
    always_switch_data$final_weight[always_switch_data$final_weight > 200] <- 200 
    
    rm(tau, tau2)
    
    ##############################
    ## 2) always dual data:
    ##############################
    
    #first define adherence
    always_dual_data<-always_dual_data %>% 
      mutate(adherence = case_when(
        violation2 == 1 | censor_dev == 1  ~ 0,
        TRUE ~ 1)) 
    
    # the following incorporates information from eTable4
    always_dual_data <- always_dual_data %>% mutate (
      denom_weight_adherence = case_when (
        visit==1 & visitweek8==0 & intol2se==0 & intol3se==0 & (failure==1|failure==2) ~ pstay1,
        visit==1 & visitweek8==0 & intol2se==1 & intol3se==0 & failure==2 ~ Prob_MIRT_Intol + Prob_NTP_Intol,
        visit==1 & visitweek8==0 & intol2se==0 & intol3se==1 & failure==3 ~ 0.5,
        visit==1 & visitweek8==1 & intol2se==1 & intol3se==0 & failure==2 ~ Prob_MIRT_Intol + Prob_NTP_Intol,
        visit==1 & visitweek8==1 & intol2se==0 & intol3se==1 & failure==3 ~ 0.5,
        visit==1 & visitweek8==1 & intol2se==0 & intol3se==0 & after_accept_benef==0 & (failure==1 | (failure==2 & changetrt==1)) ~ pstay1 + ((1-pstay1)*(`Prob_Li/T4_NR`)),
        visit==1 & visitweek8==1 & intol2se==0 & intol3se==0 & after_accept_benef==0 & (failure==3 | (failure==2 & changetrt==0)) ~ pstay1+((1-pstay1)*0.5),
        visit==1 & visitweek8==2 & intol2se==1 & intol3se==0 & failure==2 ~ Prob_MIRT_Intol + Prob_NTP_Intol,
        visit==1 & visitweek8==2 & intol2se==0 & intol3se==1 & failure==3 ~ 0.5,
        visit==1 & visitweek8==2 & intol2se==0 & intol3se==0 & after_accept_benef==0 & failure==2 & changetrt==1 ~ (1-pstay2)*(`Prob_Li/T4_NR`),
        visit==1 & visitweek8==2 & intol2se==0 & intol3se==0 & after_accept_benef==0 & failure==3 & changetrt==1 ~ (1-pstay2)*0.5,
        days_study==0 ~ 1,
        visit==0 ~ 1,
        failure==3 & changetrt==0 ~ 1,
        after_accept_benef!=0 & visitweek8!=0 & intol2se==0 & intol3se==0 ~ 1,
        TRUE ~ 1))
    
    # Numerator of unstabilized weights for adherence
    always_dual_data$k1_n <- ifelse(always_dual_data$adherence==1,1,0)
    
    # Denominator of unstabilized weights for adherence
    always_dual_data$k1_d <- with(always_dual_data, ave(denom_weight_adherence, id, FUN=cumprod))
    
    # Create  unstabilized weights for adherence
    #always_dual_data$unstabw<-with(always_dual_data, k1_n/ k1_d)
    always_dual_data <- always_dual_data %>% mutate (unstabw = case_when (
      k1_n ==0 & k1_d ==0 ~ 0,
      TRUE ~ k1_n / k1_d))
    
    tau <- quantile(always_dual_data$unstabw, 0.995) # truncate weights
    always_dual_data$unstabw[always_dual_data$unstabw > tau] <- tau # replace those above 99.5% with the 99.5%
    
    summary(always_dual_data$unstabw)
    
    # weights for ltfu
    always_dual_data <- always_dual_data %>% mutate (denom_weight_ltfu = case_when (
      days_since_last==27 & after_accept_benef==0 & visit2==1  ~ p_ltfu,
      TRUE ~ 1))
    
    # Numerator of unstabilized weights for ltfu
    always_dual_data$k2_n<-ifelse(always_dual_data$censoring_event==0,1,0)
    
    # Denominator of unstabilized weights for ltfu
    always_dual_data$k2_d <- with(always_dual_data, ave(denom_weight_ltfu, id, FUN=cumprod))
    
    # Create  unstabilized weights for adherence
    always_dual_data$unstabw2<-with(always_dual_data, k2_n / k2_d)
    
    tau2 <- quantile(always_dual_data$unstabw2, 0.995) # truncate weights
    always_dual_data$unstabw2[always_dual_data$unstabw2 > tau2] <- tau2 # replace those above 99.5% with the 99.5%
    
    summary(always_dual_data$unstabw2)
    
    # multiplication of weights
    always_dual_data$final_weight = always_dual_data$unstabw * always_dual_data$unstabw2
    #add this line if running more than 300 iterations
    always_dual_data$final_weight[always_dual_data$final_weight > 200] <- 200 
    
    rm(tau, tau2)
    
    ###################################
    ## 3) guidelines-based strategy
    #########################################3
    
    #first define adherence
    guidelines_based_data<-guidelines_based_data %>% 
      mutate(adherence = case_when(
        violation2 == 1 | censor_dev == 1  ~ 0,
        TRUE ~ 1))
    
    # the following incorporates information from eTable 5
    guidelines_based_data <- guidelines_based_data %>% mutate (after56 = ifelse(day_with_atd>=56,1,0),
                                                               denom_weight_adherence = case_when (
                                                                 visit==1 & visitweek8==0 & intol2se==0 & intol3se==0 & (failure==1|failure==2) ~ pstay1,
                                                                 visit==1 & visitweek8==0 & intol2se==1 & intol3se==0 & failure==2 & dep.improv2>=25 ~ Prob_MIRT_Intol + `Prob_Li/T4_Intol`,
                                                                 visit==1 & visitweek8==0 & intol2se==1 & intol3se==0 & failure==2 & dep.improv2<25 ~ Prob_NTP_Intol,
                                                                 visit==1 & visitweek8==0 & intol2se==0 & intol3se==1 & failure==3 & dep.improv2>=25 ~ 1,
                                                                 visit==1 & visitweek8==0 & intol2se==0 & intol3se==1 & failure==3 & dep.improv2<25 ~ 0.5,
                                                                 visit==1 & visitweek8==1 & intol2se==1 & intol3se==0 & failure==2 & dep.improv2>=25 ~ Prob_MIRT_Intol + `Prob_Li/T4_Intol`,
                                                                 visit==1 & visitweek8==1 & intol2se==1 & intol3se==0 & failure==2 & dep.improv2<25 ~ Prob_NTP_Intol,
                                                                 visit==1 & visitweek8==1 & intol2se==0 & intol3se==1 & failure==3 & dep.improv2>=25 ~ 1,
                                                                 visit==1 & visitweek8==1 & intol2se==0 & intol3se==1 & failure==3 & dep.improv2<25 ~ 0.5,
                                                                 visit==1 & visitweek8==1 & intol2se==0 & intol3se==0 & after_accept_benef==0 & (failure==1 | (failure==2 & changetrt==1)) & dep.improv2>=25 ~ pstay1 + ((1-pstay1)*(`Prob_Li/T4_NR`)),
                                                                 visit==1 & visitweek8==1 & intol2se==0 & intol3se==0 & after_accept_benef==0 & (failure==1 | (failure==2 & changetrt==1)) & dep.improv2<25 ~ pstay1 + ((1-pstay1)*(Prob_NTP_NR)),
                                                                 visit==1 & visitweek8==1 & intol2se==0 & intol3se==0 & after_accept_benef==0 & (failure==3 | (failure==2 & changetrt==0)) & dep.improv2>=25 ~ pstay1+((1-pstay1)*0.5),
                                                                 visit==1 & visitweek8==1 & intol2se==0 & intol3se==0 & after_accept_benef==0 & (failure==3 | (failure==2 & changetrt==0)) & dep.improv2<25 ~ pstay1+((1-pstay1)*0.5),
                                                                 visit==1 & visitweek8==2 & intol2se==1 & intol3se==0 & failure==2 & dep.improv2>=25 ~ Prob_MIRT_Intol + `Prob_Li/T4_Intol`,
                                                                 visit==1 & visitweek8==2 & intol2se==1 & intol3se==0 & failure==2 & dep.improv2<25 ~ Prob_NTP_Intol,
                                                                 visit==1 & visitweek8==2 & intol2se==0 & intol3se==1 & failure==3 & dep.improv2>=25 ~ 1,
                                                                 visit==1 & visitweek8==2 & intol2se==0 & intol3se==1 & failure==3 & dep.improv2<25 ~ 0.5,
                                                                 visit==1 & visitweek8==2 & intol2se==0 & intol3se==0 & after_accept_benef==0 & failure==2 & changetrt==1 & dep.improv2>=25 ~ (1-pstay2)*(`Prob_Li/T4_NR`),
                                                                 visit==1 & visitweek8==2 & intol2se==0 & intol3se==0 & after_accept_benef==0 & failure==2 & changetrt==1 & dep.improv2<25 ~  (1-pstay2)* Prob_NTP_NR,
                                                                 visit==1 & visitweek8==2 & intol2se==0 & intol3se==0 & after_accept_benef==0 & failure==3 & changetrt==1 ~ (1-pstay2)*0.5,
                                                                 days_study==0 ~ 1,
                                                                 visit==0 ~ 1,
                                                                 failure==3 & changetrt==0 ~ 1,
                                                                 after_accept_benef!=0 & visitweek8!=0 & intol2se==0 & intol3se==0 ~ 1,
                                                                 TRUE ~ 1))
    
    # Numerator of unstabilized weights for adherence
    guidelines_based_data$k1_n <- ifelse(guidelines_based_data$adherence==1,1,0)
    
    # Denominator of unstabilized weights for adherence
    guidelines_based_data$k1_d <- with(guidelines_based_data, ave(denom_weight_adherence, id, FUN=cumprod))
    
    # Create  unstabilized weights for adherence
    #guidelines_based_data$unstabw<-with(guidelines_based_data, k1_n / k1_d)
    guidelines_based_data <- guidelines_based_data %>% mutate (unstabw = case_when (
      k1_n ==0 & k1_d ==0 ~ 0,
      TRUE ~ k1_n / k1_d))
    
    tau <- quantile(guidelines_based_data$unstabw, 0.995) # truncate weights
    guidelines_based_data$unstabw[guidelines_based_data$unstabw > tau] <- tau # replace those above 99.5% with the 99.5%
    
    summary(guidelines_based_data$unstabw)
    
    # weights for ltfu
    guidelines_based_data <- guidelines_based_data %>% mutate (denom_weight_ltfu = case_when (
      days_since_last==27 & after_accept_benef==0 & visit2==1  ~ p_ltfu,
      TRUE ~ 1))
    
    # Numerator of unstabilized weights for ltfu
    guidelines_based_data$k2_n <- ifelse(guidelines_based_data$censoring_event==0,1,0)
    
    # Denominator of unstabilized weights for ltfu
    guidelines_based_data$k2_d <- with(guidelines_based_data, ave(denom_weight_ltfu, id, FUN=cumprod))
    
    # Create  unstabilized weights for adherence
    guidelines_based_data$unstabw2<-with(guidelines_based_data, k2_n / k2_d)
    
    tau2 <- quantile(guidelines_based_data$unstabw2, 0.995) # truncate weights
    guidelines_based_data$unstabw2[guidelines_based_data$unstabw2 > tau2] <- tau2 # replace those above 99.5% with the 99.5%
    
    summary(guidelines_based_data$unstabw2)
    
    # multiplication of weights
    guidelines_based_data$final_weight = guidelines_based_data$unstabw * guidelines_based_data$unstabw2
   
     #add this line if running more than 300 iterations
    guidelines_based_data$final_weight[guidelines_based_data$final_weight > 200] <- 200 
    
    rm(tau, tau2)
    
    # now, merge the datasets with the weights
    
    finaldata<-rbind(always_dual_data, always_switch_data, guidelines_based_data)
    
    finaldata$days_study2 = finaldata$days_study * finaldata$days_study

    # because in the finaldata we have more than one individual with the same id, create a clone id
    finaldata<-finaldata %>% mutate(clone_id = case_when (
      strategy==1 ~ as.numeric(id),
      strategy==2 ~ as.numeric(id)+ 10000,
      strategy==3 ~ as.numeric(id) + 20000))
    
    model.surv <- glm(primary_outcome==0 ~ as.factor(strategy) + days_study + days_study2 +
                        I(strategy*days_study) + I(strategy*days_study2), 
                      data=finaldata,
                      weights = final_weight,
                      family=quasibinomial())
    
    ##################################################################################
    # Curves
    #################################################################################
    baseline <- finaldata[finaldata$days_study==0,]
    
    gf.strategy0 <- expandRows(baseline, count=270, count.is.col=F) 
    gf.strategy0$days_study <- rep(seq(0, 269), nrow(baseline))
    gf.strategy0$days_study2 <- gf.strategy0$days_study^2
    gf.strategy0$strategy <- 1
    
    gf.strategy1 <- gf.strategy0
    gf.strategy1$strategy <- 2
    
    gf.strategy2 <- gf.strategy0
    gf.strategy2$strategy <- 3
    
    gf.strategy0$p.noevent0 <- predict(model.surv, gf.strategy0, type="response")
    gf.strategy1$p.noevent1 <- predict(model.surv, gf.strategy1, type="response")
    gf.strategy2$p.noevent2 <- predict(model.surv, gf.strategy2, type="response")
    
    gf.strategy0.surv <- gf.strategy0 %>% group_by(clone_id) %>% mutate(surv0 = cumprod(p.noevent0),
                                                                        ci = 1-surv0)
    gf.strategy1.surv <- gf.strategy1 %>% group_by(clone_id) %>% mutate(surv1 = cumprod(p.noevent1),
                                                                        ci = 1-surv1)
    gf.strategy2.surv <- gf.strategy2 %>% group_by(clone_id) %>% mutate(surv2 = cumprod(p.noevent2),
                                                                        ci = 1-surv2)
    
    gf.surv0 <- aggregate(gf.strategy0.surv, by=list(gf.strategy0.surv$days_study), FUN=mean)[c("strategy", "days_study", "ci")]
    gf.surv1 <- aggregate(gf.strategy1.surv, by=list(gf.strategy1.surv$days_study), FUN=mean)[c("strategy", "days_study", "ci")]
    gf.surv2 <- aggregate(gf.strategy2.surv, by=list(gf.strategy2.surv$days_study), FUN=mean)[c("strategy", "days_study", "ci")]
    
    all <- rbind(gf.surv0, gf.surv1, gf.surv2)
    all <- all[,c('strategy', 'days_study', 'ci')]
    
    # Calculate the mean survival at each visit within each strategy
    results <- aggregate(ci ~ days_study + strategy, FUN=mean, data=all)
    
    # Add a variable that treats strategy as a factor
    results$strategy <- factor(results$strategy, labels = c("Sequential monotherapy", "Sequential dual", "Guidelines-based"))
    
    ##
    # 270-remission differences
    ##
    results2<-results %>% filter(days_study==269) # change to180?
    results2<- results2 %>% spread (strategy, ci)
    results2$survdiff1 <- results2$`Sequential dual`-results2$`Sequential monotherapy`
    results2$survdiff2 <- results2$`Guidelines-based`-results2$`Sequential monotherapy`
    results2$rrisk1 <- results2$`Sequential dual` / results2$`Sequential monotherapy`
    results2$rrisk2 <- results2$`Guidelines-based` / results2$`Sequential monotherapy`
    results2$nnt1 = 1/results2$survdiff1
    results2$nnt2 = 1/results2$survdiff2
    
    a<-results2%>%select(`Sequential monotherapy`, `Sequential dual`,`Guidelines-based`, survdiff1, 
                         survdiff2, rrisk1, rrisk2, nnt1, nnt2)
    
    if (0){
      
      
    } # end of commented out code 
    
    boot.results[i+1, ] <- a[1,]
    print(i)
    
  } # end of bootstrap loop 
  return(boot.results)
} # end of bootstrap function 


# get data set up for bootstrapping and keeping id as first variable
# play around with values and name of first column
finaldata$days_study2 = finaldata$days_study * finaldata$days_study

finaldata0<-finaldata 
rm(finaldata)
DT <- as.data.table(finaldata0)
data_len <- length(unique(DT[,id])) 
newid <- as.data.table(unique(DT[,id]))[,newid:=seq(.N)]

setnames(newid,"V1","id")
DT1 <- DT[newid,on="id"][,id:=newid][,newid:=NULL]

test<-bootstrap_fn(DT1, num.boot=100)
print(test)

quantile(test$sequential_monotherapy, probs = c(0.025, 0.975))
quantile(test$sequential_dual, probs = c(0.025, 0.975))
quantile(test$guidelines_based, probs = c(0.025, 0.975))

quantile(test$survdiff1, probs = c(0.025, 0.975))
quantile(test$survdiff2, probs = c(0.025, 0.975))

quantile(test$rrisk1, probs = c(0.025, 0.975))
quantile(test$rrisk2, probs = c(0.025, 0.975))


quantile(a$nnt1, probs = c(0.025, 0.975))
quantile(a$nnt2, probs = c(0.025, 0.975))





