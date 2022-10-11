# DATA CLEANING
# "Emulating a Target Trial of Sequential Treatment Decisions in Major Depressive Disorder 
#  Using Data from the Sequenced Treatment Alternatives to Relieve Depression (STAR*D) Trial"

library(tidyverse)
library(dplyr)
library(tidyr)

#####################################################################################################################################
# STEP 0. Identify eligible population from the original trial  (i.e., meeting inclusion criteria and consenting)
# 4,041 patients consented.  

# Here, I select IDs for the eligible individuals who started Level 2 of the trial (i.e., failed to CIT)
# 1439 patients entered level 2 
#####################################################################################################################################

# select eligible IDs 
elig <- read.delim("//DON/allproj/NMIH/data/STARD/el01.txt")
elig <- elig[-1,] 

elig <- elig %>% filter (pstat == 4) %>% select (src_subject_id) # list of 4,041 eligible IDs


# select IDs that started at Level 2

ivra<-read.delim("//DON/allproj/NMIH/data/STARD/ivra01.txt")
ivra<-ivra[-1,]
ivra <- ivra %>% filter ((level == "Level 2") & (week == "0")) %>% select (src_subject_id)
ivra<-as.data.frame(ivra[!duplicated(ivra[1]),])
names(ivra)[1]<-paste("src_subject_id")

elig<-elig %>% filter (src_subject_id %in% ivra$src_subject_id)
rm(ivra)

#####################################################################################################################################
# STEP 1. Extract baseline variables for eligible 1,439 patients
# The following baseline variables are formatted to mimic Rush et al. (2006)
# Baseline variables: age, gender, ethnicity, marital status, employment status, education years, 
#                     melancholic features, number of DSM criteria for MDD met, age at onset, 
#                     comorbid anxiety disorders (AGO-OCD-PTSD-GAD), 
#                     duration of current MDE and baseline HDRS

# Too many missing values to be included: race ,  number of prior episodes (281)
#####################################################################################################################################

dem <- read.delim("//DON/allproj/NMIH/data/STARD/dm01.txt")
dem <- dem[-1,]

# select eligible ids and carrying forward time-static confounders within IDs 
dem <- dem %>% filter (src_subject_id %in% elig$src_subject_id) %>% 
  mutate(age = as.numeric( as.character( interview_age ) ) / 12) %>%
  group_by(src_subject_id) %>% arrange(src_subject_id, days_baseline)

bas1<- dem %>% select (c(src_subject_id, days_baseline, age, gender, empl, educat, marital)) %>% 
  mutate(days_baseline = as.numeric(as.character(days_baseline)),
         female = case_when(
           gender == "M" ~ "0",
           TRUE ~ "1"),
         empl = as.numeric(as.character(empl)),
         empl = case_when(empl == 1 ~ "Unemployed",
                          empl == 2 ~ "Unemployed",
                          empl == 3 ~ "Employed",
                          empl == 4 ~ "Employed",
                          empl == 5 ~ "Employed",
                          empl == 6 ~ "Retired",
                          TRUE ~ NA_character_),
         educat = as.numeric(as.character(educat)),
         marital = as.numeric(as.character(marital)),
         marital = case_when(marital == 1 ~ "Single (never married)",
                             marital == 2 ~ "Married (cohabiting)",
                             marital == 3 ~ "Married (cohabiting)",
                             marital == 4 ~ "Divorced or separated",
                             marital == 5 ~ "Divorced or separated",
                             marital == 6 ~ "Widowed",
                             marital == 7 ~ "Married (cohabiting)",
                             marital == 8 ~ "Single (never married)",
                             TRUE ~ NA_character_)) %>%
  group_by(src_subject_id) %>% arrange (days_baseline, .by_group = TRUE) %>%
  fill(c(age, gender, empl, educat, marital), .direction="down") %>% 
  slice (c(n())) %>% ungroup() %>% select (-c(days_baseline)) # baseline info for all 1439 patients

idsc<-read.delim("//DON/allproj/NMIH/data/STARD/idsc01.txt")                                   
idsc<-idsc[-1,]
idsc1 <- idsc %>% filter (src_subject_id %in% elig$src_subject_id) %>% group_by(src_subject_id) %>% 
  arrange(days_baseline, .by_group = T) %>% slice (c(1)) %>% ungroup() %>% 
  mutate (
    idsc1 = as.numeric(as.character(isoin)),
    idsc2 = as.numeric(as.character(imnin)),
    idsc3 = as.numeric(as.character(iemin)),
    idsc4 = as.numeric(as.character(ihysm)),
    idsc5 = as.numeric(as.character(imdsd)),
    idsc6 = as.numeric(as.character(ianx)),
    idsc7 = as.numeric(as.character(ipanc)),
    idsc8 = as.numeric(as.character(iirtb)),
    idsc9 = as.numeric(as.character(irct)),
    idsc10 = as.numeric(as.character(ivrtn)),
    idsc11 = as.numeric(as.character(iqty)),
    idsc12 = as.numeric(as.character(iapdc)),
    idsc13 = as.numeric(as.character(iapin)),
    idsc14 = as.numeric(as.character(iwtdc)),
    idsc15 = as.numeric(as.character(iwtin)),
    idsc16 = as.numeric(as.character(icntr)),
    idsc17 = as.numeric(as.character(ivwsf)),
    idsc18 = as.numeric(as.character(ivwfr)),
    idsc19 = as.numeric(as.character(isuic)),
    idsc20 = as.numeric(as.character(iintr)),
    idsc21 = as.numeric(as.character(iplsr)),
    idsc22 = as.numeric(as.character(iengy)),
    idsc23 = as.numeric(as.character(isex)),
    idsc24 = as.numeric(as.character(islow)),
    idsc25 = as.numeric(as.character(iagit)),
    idsc26 = as.numeric(as.character(ismtc)),
    idsc27 = as.numeric(as.character(isymp)),
    idsc28 = as.numeric(as.character(igas)),
    idsc29 = as.numeric(as.character(iintp)),
    idsc30 = as.numeric(as.character(ildn))) %>%  replace(is.na(.), 0) %>% 
  mutate (idsctot = rowSums(.[58:87]),
          qmood = ifelse(iqty==3, 1, 0),    
          moodvar = ifelse(((ivrtn==2 | ivrtn==3) & iwrse==1), 1, 0),
          psychred = ifelse((islow==2|islow==3), 1, 0),
          earlyins = ifelse(iemin==3,1,0),
          psychag = ifelse((iagit==2 | iagit==3), 1, 0),
          appdec = ifelse((iapdc==2|iapdc==3), 1, 1),
          appinc = ifelse(iapin ==0, 0, 1),            
          ivwsf = ifelse((ivwsf==2|ivwsf==3), 1, 0),
          leaden = ifelse((ildn==2|ildn==3), 1,0),
          weight = ifelse((iwtin==2|iwtin==3|iapin==2|iapin==3), 1, 0),
          hypersom = ifelse((ihysm==2|ihysm==3),1,0),
          interpers = ifelse((iintp==3), 1, 0),
          psychomotor = ifelse((psychred==1 | psychag==1),1,0),
          weightdec = ifelse((iwtdc==3|iapdc==2|iapdc==3), 1, 0),
          atypical = case_when (
                (irct==0|irct==1|irct==2) & (leaden+weight+hypersom+interpers>=2) ~ "1",
                TRUE~"0"),
          melancholic = case_when (
                 ((irct==3 |irct==2) | (iplsr==2|iplsr==3)) & (as.numeric(qmood)+as.numeric(moodvar)+as.numeric(iemin)+as.numeric(psychomotor)+as.numeric(weightdec)+as.numeric(ivwsf)>=3) ~ 1,
                TRUE ~ 0),
          mild = case_when(
            idsctot<=25 ~ 1,
            TRUE ~ 0),
          moderate = case_when (
            idsctot>25 & idsctot < 39 ~ 1,
            TRUE ~ 0),
          severe = case_when (
            idsctot>=39 ~ 1,
            TRUE ~ 0)) %>% select (src_subject_id, melancholic, atypical, idsctot, mild, moderate, severe) # info on melancholic/atypical features for 1,434

# race has too many missing values to be used.
race<- read.delim("//DON/allproj/NMIH/data/STARD/el01.txt")
bas2<- race  %>% select (c(src_subject_id, race, ethnicity, dsmdm, dsmdi, dsmcw, dsmih, dsmpa, dsmle, dsmfw, dsmdt, dsmtd)) %>%
  filter (src_subject_id %in% elig$src_subject_id)  %>%
  mutate (ethnicity = as.factor(case_when(ethnicity == "Hispanic or Latino" ~ "Hispanic",
                                          ethnicity == "Not Hispanic or Latino" ~ "Not hispanic",
                                          TRUE ~ NA_character_))) %>%
  mutate (numberdsm = as.numeric(as.character(dsmdm))+
            as.numeric(as.character(dsmdi))+
            as.numeric(as.character(dsmcw))+
            as.numeric(as.character(dsmih))+
            as.numeric(as.character(dsmpa))+
            as.numeric(as.character(dsmle))+
            as.numeric(as.character(dsmfw))+
            as.numeric(as.character(dsmdt))+
            as.numeric(as.character(dsmtd))) %>%
  select(c(src_subject_id, ethnicity, numberdsm)) # no missing values


psych.hist <- read.delim("//DON/allproj/NMIH/data/STARD/phx01.txt")

bas3<-psych.hist %>% select (c(src_subject_id, alcoh, days_baseline, dage, epino, dep, episode_date, pd_ag, specphob, 
                               soc_phob, ocd_phx, psd, gad_phx)) %>%
  filter (src_subject_id %in% elig$src_subject_id) %>% mutate (
    ageonset = case_when (
      as.numeric(as.character(dage)) == -9 ~ NA_character_,
      as.numeric(as.character(dage)) == -3 ~ NA_character_,
      as.numeric(as.character(dage)) < 18 ~ "1",
      as.numeric(as.character(dage)) >= 18 ~ "0",
      TRUE ~ NA_character_),
    durationmde = ifelse(abs(as.numeric(as.character(episode_date)))>=730, 1, 0),
    comorbidago = as.factor( case_when (
      pd_ag == "1" ~ "1",
      pd_ag == "0" ~ "0",
      TRUE ~ "0")),
    comorbidocd = as.factor( case_when (
      ocd_phx == "1" ~ "1",
      ocd_phx == "0" ~ "0",
      TRUE ~ "0")),
    comorbidptsd = as.factor( case_when (
      psd == "1" ~ "1",
      psd == "0" ~ "0",
      TRUE ~ "0")),
    comorbidgad = as.factor( case_when (
      gad_phx == "1" ~ "1",
      gad_phx == "0" ~ "0",
      TRUE ~ "0")),
    baseline_alcoh = case_when (
      alcoh == "1" | alcoh == "2" ~ "1",
      TRUE ~ "0")) %>% group_by(src_subject_id) %>% 
  arrange(days_baseline, .by_group = T) %>% slice (c(1)) %>% ungroup() %>%
  select (c(src_subject_id, baseline_alcoh, ageonset, durationmde, comorbidago, comorbidocd, comorbidptsd, comorbidgad))

bas3$ageonset[is.na(bas3$ageonset)] <-0

# cumulative illness rating score. 
cirs<-read.delim("//DON/allproj/NMIH/data/STARD/crs01.txt")
cirs<-cirs[-1,]
bas4<- cirs %>% filter (src_subject_id %in% elig$src_subject_id) %>% group_by(src_subject_id) %>% slice (c(1)) %>%
  mutate (cirstot = 
            as.numeric(as.character(heart))+
            as.numeric(as.character(vsclr))+
            as.numeric(as.character(hema))+
            as.numeric(as.character(eyes))+
            as.numeric(as.character(ugi))+
            as.numeric(as.character(lgi))+
            as.numeric(as.character(renal))+
            as.numeric(as.character(genur))+
            as.numeric(as.character(mskl))+
            as.numeric(as.character(neuro))+
            as.numeric(as.character(psych))+
            as.numeric(as.character(respiratory))+
            as.numeric(as.character(liverd))+
            as.numeric(as.character(endod)), na.rm = T) %>% select (src_subject_id, cirstot) # all 1439 patients

## Baseline HDRS scores at study entry
hrsd<-read.delim("//DON/allproj/NMIH/data/STARD/hrsd01.txt")
hrsd<-hrsd[-1,]   

# extract all HRDS per individual.Select last HRDS at level 1
bas5<-hrsd %>% mutate(days_baseline=as.numeric(as.character(days_baseline))) %>%
  filter (src_subject_id %in% elig$src_subject_id) %>% filter (level %in% c("Level 1", "Enrollment")) %>%
  group_by(src_subject_id) %>% arrange(days_baseline, .by_group=T)%>% slice(c(n())) %>%                                                         
  replace(is.na(.), 0)  %>%
  mutate (
    hamd1 = as.numeric(as.character(hsoin)),
    hamd2 = as.numeric(as.character(hmnin)),
    hamd3 = as.numeric(as.character(hemin)),
    hamd4 = as.numeric(as.character(hmdsd)),
    hamd5 = as.numeric(as.character(hpanx)),
    hamd6 = as.numeric(as.character(hinsg)),
    hamd7 = as.numeric(as.character(happt)),
    hamd8 = as.numeric(as.character(hwl)),
    hamd9 = as.numeric(as.character(hsanx)),
    hamd10 = as.numeric(as.character(hhypc)),
    hamd11 = as.numeric(as.character(hvwsf)),
    hamd12 = as.numeric(as.character(hsuic)),
    hamd13 = as.numeric(as.character(hintr)),
    hamd14 = as.numeric(as.character(hengy)),
    hamd15 = as.numeric(as.character(hslow)),
    hamd16 = as.numeric(as.character(hagit)),
    hamd17 = as.numeric(as.character(hsex))) %>%  mutate (hrsdtot = hamd1+hamd2+hamd3+hamd4+hamd5+hamd6+hamd7+hamd8+hamd9+hamd10+hamd11+hamd12+hamd13+hamd14+
                                                            hamd15+hamd16+hamd17,
                                                          baseline.hrsd = as.numeric(as.character(hrsdtot)),
                                                          factoranx =  hamd5+hamd9+hamd7+hamd14+hamd10+hamd6,
                                                          anxiousdepr = ifelse(factoranx>=7,1,0))%>%
  ungroup() %>% 
  select(src_subject_id, baseline.hrsd, anxiousdepr)

baseline<- merge(bas1, bas2, by = "src_subject_id", all = T)
baseline<- merge (baseline, bas3, by = "src_subject_id", all = T)
baseline<-merge(baseline, bas4, by = "src_subject_id")
baseline<-merge (baseline, idsc1, by = "src_subject_id", all = T)
baseline<-merge(baseline, bas5, by = "src_subject_id", all = T)
baseline<-baseline %>% mutate(id= src_subject_id) %>% select(-c(src_subject_id))

rm(bas1, bas2, bas3, bas4, bas5, hrsd, cirs, dem, psych.hist, race, idsc, idsc1)


# primary vs. specialty care
# setting = 1 (Primary care)
# setting = 2 (Specialty care)
clinics1 <- read.csv("//don/allproj/NMIH/programs/aszmulewicz/clinics1.csv")
clinics1 <- clinics1%>%filter(ID %in% baseline$id) %>% mutate (id = ID) %>% select(id, Setting)

baseline<-merge(baseline, clinics1, by = "id", all = T)

table(baseline$Setting)

#####################################################################################################################################
# STEP 2. Incorporate random treatment assignment at each level entry
# Plus, this will show when a patient is moved to FU
#####################################################################################################################################
ivra<-read.delim("//DON/allproj/NMIH/data/STARD/ivra01.txt")
ivra<-ivra[-1,]

trt1 <- ivra %>% mutate(txassign = ifelse (txassign=="", NA, txassign)) %>% filter (src_subject_id %in% elig$src_subject_id) %>%
  filter (! level == "Level 1") %>%
  mutate (onlyct = case_when (
    (cogaugl2 == 1 | cogswl2==1) & !medswl2==1 & !medaugl2==1  ~ 1,
    TRUE ~ 0)) %>% 
  select (src_subject_id, days_baseline, txassign, onlyct) 
trt1<-na.omit(trt1)  # will keep the rows were the treatment is assigned

trt1<-trt1[!duplicated(trt1[1:2]),]

trt1<-trt1 %>% mutate (id = src_subject_id,
                       changelevel = ifelse(txassign == "FU", NA, 1)) %>% select (-c(src_subject_id)) 

# merge with baseline data and replace FU by prior medication
merged.data1<-merge (baseline, trt1, by = "id", all=T) %>% group_by(id) %>% 
  mutate(days_baseline = as.numeric(as.character(days_baseline)),
         start_time = min(days_baseline),
         txassign = ifelse(txassign=="FU", lag(as.character(txassign), order_by = days_baseline), as.character(txassign))) %>%
  dplyr::arrange (days_baseline, .by_group  = T) %>% ungroup()

rm(ivra, trt1)

#####################################################################################################################################
# STEP 3. Expand baseline data with all clinical visits for all included individuals
# Plus, include the days_baseline that all visits occurred.
# Extract indicator of week of clinical visit as well.
#####################################################################################################################################

# select the days that all clinic visits occurred for all eligible subjects
ccv<- read.delim("//DON/allproj/NMIH/data/STARD/ccv01.txt")
ccv<-ccv[-1,]

ccv1 <- ccv %>% mutate(id = src_subject_id) %>% filter(id %in% baseline$id) %>% filter (! level == "Level 1") %>%
  mutate(days_baseline = as.numeric(as.character(days_baseline)),
         week91214 = case_when (
           week == "9" | week == "12" | week == "14" ~ 1,
           TRUE ~ 0),
         week9 = ifelse(week=="9",1,0),
         week12 = ifelse(week=="12", 1,0),
         week14 = ifelse(week=="14",1,0),
         visit = 1) %>% select (id, days_baseline, visit, remsn, week91214, week9, week12, week14, level) %>% 
  arrange(id) 

ccv1<-ccv1[!duplicated(ccv1[1:2]),] # last day of a level serves as week 0 of next level (info duplicated)

#merge and create an indicator that a clinical visit took place at that day
merged.data2 <- merge (merged.data1, ccv1, by = c("id", "days_baseline"),  all=T) %>% group_by(id) %>% 
  fill (start_time, .direction="updown") %>% 
  
  # remove rows before start time. Also, remove those with days_baseline = 0. These 3 individuals entered
  # STAR*D at level 2 (hence days=0 in level 2). These individuals do not meet eligibility criteria 
  # for our target trial emulation
  
  filter (as.numeric(as.character(days_baseline)) >= as.numeric(as.character(start_time))) %>%
  ungroup()

todelete<- merged.data2 %>% filter(days_baseline==0) %>% select(id)

merged.data2 <- merged.data2 %>% filter(!id %in% todelete$id)

rm(ccv, ccv1, merged.data1, todelete)

#####################################################################################################################################
# STEP 4. Incorporation of time-varying depressive scores (QIDS-C and QIDS-SR) to adjust for confounding
####################################################################################################################################
qids<-read.delim("//DON/allproj/NMIH/data/STARD/qids01.txt")
qids<- qids[-1,]

# first, extract depressive symptoms at every clinical visit and construct QIDS-C total score by hand
qids1<-qids %>% mutate(days_baseline = as.numeric(as.character(days_baseline)),
                       id  = src_subject_id) %>% filter(version_form == "Clinician") %>% 
  filter(id %in% merged.data2$id) %>%
  mutate  (
    qids1soin= as.numeric(as.character(vsoin)),
    qids2mnin=as.numeric(as.character(vmnin)),
    qids3emin=as.numeric(as.character(vemin)),
    qids4hysm=as.numeric(as.character(vhysm)),
    qids5mdsd=as.numeric(as.character(vmdsd)),
    qids6apdc=as.numeric(as.character(vapdc)),
    qids7apin=as.numeric(as.character(vapin)),
    qids8wtdc=as.numeric(as.character(vwtdc)),
    qids9wtin=as.numeric(as.character(vwtin)),
    qids10cntr=as.numeric(as.character(vcntr)),
    qids11vwsf=as.numeric(as.character(vvwsf)),
    qids12suic=as.numeric(as.character(vsuic)),
    qids13intr=as.numeric(as.character(vintr)),
    qids14engy=as.numeric(as.character(vengy)),
    qids15slow=as.numeric(as.character(vslow)),
    qids16agit=as.numeric(as.character(vagit))) %>% select (id, days_baseline, qids1soin:qids16agit, level)

# remove rows with missing in all
qids1<-qids1[rowSums(is.na(qids1[, 3:18])) < 16, ]
qids1<-qids1 %>% replace(is.na(.), 0)

# get the maximum score from 1-4 and 6-9 and 15-16 to get the QIDS-total score
qids1$max1 <- do.call(pmax, qids1[3:6])
qids1$max2 <- do.call(pmax, qids1[8:11])
qids1$max3 <- do.call(pmax, qids1[17:18])

qids1<-qids1 %>% mutate(qstot=max1 + qids5mdsd+ max2 + qids10cntr+qids11vwsf+qids12suic+qids13intr+
                          qids14engy+max3)

# baseline QIDS will be defined as the last QIDS measurement in Level 1 (to make sure it happened before exposure)
qids.bas <- qids1 %>% filter (level == "Level 1") %>% 
  group_by(id) %>% arrange(days_baseline, .by_group = TRUE) %>% slice (n()) %>% ungroup() %>% 
  mutate (
    baseline.qids = qstot) %>% select (id, baseline.qids) # all patients with baseline QIDS

merged.data3<-merge (merged.data2, qids.bas, by = "id", all=T )

# second, incorporate all depressive measurement as follow-up measurements.
qids.fup <- qids1  %>% filter (!level == "Level 1") %>% 
  mutate (qids.fup = qstot) %>% 
  select (id, days_baseline, qids.fup)

# remove the duplicate visits (they act as visit0 of the next level)
qids.fup<-qids.fup[!duplicated(qids.fup[1:2]),]

# Second, extract % improvement at each visit using as reference QIDS-C entering STAR*D
qids.imp <- qids %>% filter(version_form=="Clinician",
                            !is.na(as.numeric(as.character(qcimp_r)))) %>% 
  mutate(days_baseline = as.numeric(as.character(days_baseline)),
         id = src_subject_id) %>% filter(id %in% baseline$id) %>%  
  filter (!level == "Level 1") %>% group_by(id) %>% arrange(days_baseline, .by_group=T) %>%
  mutate (dep.improv = as.numeric(as.character(qcimp_r))) %>%
  select (id, days_baseline, dep.improv)

# extract the QIDS at study entry to get the first improvement
qids.bas.imp<-qids1 %>% filter(level == "Level 1") %>% group_by(id) %>% 
  arrange(days_baseline, .by_group = TRUE) %>% slice (c(1)) %>% ungroup() %>% mutate(first.qids = qstot) %>%
  select(id, first.qids)

# remove the duplicate visits (they act as visit0 of the next level)
qids.imp<-qids.imp[!duplicated(qids.imp[1:2]),]

# create a second improvement variable, using the baseline (level 2) QIDS-C as the reference
qids.follow <- merge(qids.fup, qids.imp, by = c("id", "days_baseline"), all = T)

merged.data3<-merge(merged.data3, qids.follow, by = c("id", "days_baseline"), all = T)
merged.data3<-merge(merged.data3, qids.bas.imp, by = "id", all = T)

#baseline self-reported depression (latest measurement in level 1 of the star*d)
qids.self.bas <- qids %>% mutate(days_baseline = as.numeric(as.character(days_baseline)),
                                 id  = src_subject_id) %>% filter(id %in% baseline$id) %>%
  filter (level == "Level 1" & version_form == "Self Rating") %>% 
  group_by(id) %>% arrange(days_baseline, .by_group = TRUE) %>% slice (n()) %>% ungroup() %>% 
  mutate (
    baseline.qids.self = qstot) %>% select (id, baseline.qids.self) # 1438 patients with baseline self reported QIDS

merged.data4<-merge (merged.data3, qids.self.bas, by = "id", all=T )


# follow-up self-reported depression
qids.self <- qids %>% mutate(days_baseline = as.numeric(as.character(days_baseline)),
                             id  = src_subject_id) %>% filter((!level == "Level 1") & version_form == "Self Rating") %>% 
  filter(id %in% baseline$id) %>% mutate (self.qids = as.numeric(as.character(qstot))) %>%
  select(id, days_baseline, self.qids) 

qids.self<-qids.self[!duplicated(qids.self[1:2]),]

merged.data5<-merge(merged.data4, qids.self, by = c("id", "days_baseline"), all = T) %>% mutate( 
  baseline.qids.self = as.numeric(as.character(baseline.qids.self)))
####################################

#replace first row with baseline value if missing and remove rows before time zero
merged.data5<-merged.data5 %>%
  group_by(id) %>% arrange(days_baseline, .by_group = T) %>%
  fill(start_time, .direction = "updown") %>% fill(c(qids.fup, self.qids), .direction = "down")%>%
  mutate(self.qids = case_when(
    row_number() == 1 & is.na(as.numeric(as.character(self.qids))) ~ as.numeric(as.character(baseline.qids.self)), 
    TRUE ~ as.numeric(as.character(self.qids))),
    qids.fup = case_when(
      row_number() == 1 & is.na(as.numeric(as.character(qids.fup))) ~ as.numeric(as.character(baseline.qids)),
      TRUE ~ as.numeric(as.character(qids.fup)))) %>%
  fill(c(dep.improv), .direction = "down") %>% 
  filter (as.numeric(as.character(days_baseline)) >= as.numeric(as.character(start_time))) %>% 
  mutate(qids.fup = case_when(
    row_number() == 1 & is.na(as.numeric(as.character(qids.fup))) ~ as.numeric(as.character(baseline.qids)),
    TRUE ~ as.numeric(as.character(qids.fup)))) %>%
  mutate(improv = (as.numeric(baseline.qids) - as.numeric(qids.fup) )/(as.numeric(baseline.qids)) * 100,
         improv2 = case_when(
           as.numeric(improv) < 0 ~ 0,
           as.numeric(improv) > 100 ~ 100,
           TRUE ~ as.numeric(improv)),
         improv.first = ((as.numeric(first.qids) - as.numeric(qids.fup))*100)/(as.numeric(first.qids)),
         dep.improv2 = case_when (
           row_number() == 1 & !is.na(dep.improv) ~ as.numeric(dep.improv),
           row_number() == 1 & is.na(dep.improv) ~ as.numeric(improv.first),
           TRUE ~ as.numeric(improv2)),
         visit=1) %>% ungroup()

rm(qids, qids.bas, qids.fup, qids.follow, qids1, qids.imp, qids.self, qids.self.bas)

#####################################################################################################################################
# STEP 5. Side-effects
#####################################################################################################################################
prise<-read.delim("//DON/allproj/NMIH/data/STARD/prise01.txt")
prise<- prise[-1,]

prise <- prise %>% mutate (days_baseline = as.numeric(as.character(days_baseline))) %>% 
  filter(src_subject_id %in% merged.data5$id) %>% 
  filter (! level == "Level 1") %>%
  mutate(id = src_subject_id) %>% mutate (sedation = case_when (
    (slmch ==1 | oftge==1 | odegy==1) ~ 1,
    TRUE ~ 0)) %>% mutate ( sexual = case_when (
      (sxls ==1 | sxorg==1 | sxerc ==1) ~ 1,
      TRUE ~ 0)) %>% mutate (extrapiramidal = case_when (
        (nvtrm ==1 | nvcrd==1 | orsls==1) ~ 1,
        TRUE ~0))  %>% select (c(id, days_baseline, sedation, sexual, extrapiramidal)) 

prise<- prise[!duplicated(prise[1:2]),] # info at the last visit of a given level is duplicated.

merged.data6 <- merge (merged.data5, prise, by = c("id", "days_baseline"), all = T) %>%
  # replace first row with 0 if missing
  group_by(id) %>%  fill(start_time, .direction = "updown") %>%
  filter (as.numeric(as.character(days_baseline)) >= as.numeric(as.character(start_time))) %>%
  mutate(sedation = ifelse(row_number() == 1 & (is.na(sedation)), 0, sedation),
         sexual = ifelse(row_number() == 1 & (is.na(sexual)), 0, sexual),
         extrapiramidal = ifelse(row_number() == 1 & (is.na(extrapiramidal)), 0, extrapiramidal),
         visit=1)%>% 
  ungroup() 

#####################################################################################################################################
# STEP 6. Add indicators for protocol decisions at each clinical visit
# This step is needed to detect deviations from protocol
# 1) Indicator of first clinical visit after week 8 (only clinical visits)
# 2) Indicator of acceptable benefit (QIDS < 11)
# 3) Primary outcome (using QIDS)
# 4) Treatment failure
#####################################################################################################################################

#1) Indicator of the first clinical visit after week 8
merged.data6 <- merged.data6 %>% group_by(id) %>% 
  mutate(day_diff = as.numeric(days_baseline) - lag(as.numeric(days_baseline)),
         study_days = as.numeric(days_baseline) - start_time,
         change = ifelse(is.na(changelevel), 0, changelevel),
         change2 = lag(change),
         change2 = ifelse(is.na(change2), 0, change2))

merged.data6<-merged.data6 %>% group_by(id) %>% fill(txassign, .direction = "down")%>%
  mutate(level2 = cumsum(change2),
         day_diff2 = case_when(
           level2 == 2 ~ day_diff,
           TRUE ~ 0),
         day_diff3 = cumsum(day_diff2),
         day_diff4 = case_when(
           level2 == 3 ~ day_diff,
           TRUE ~ 0),
         day_diff5 = cumsum(day_diff4),
         day_diff6 = case_when(
           level2 == 4 ~ day_diff,
           TRUE ~ 0),
         day_diff7 = cumsum (day_diff6),
         day_with_atd = case_when(
           level2 == 0 ~ 0,
           level2 == 1 ~ study_days,
           level2 == 2 ~ day_diff3,
           level2 == 3 ~ day_diff5,
           level2 == 4 ~ day_diff7)) %>% select(-c(level2, day_diff, day_diff2, day_diff3, day_diff4, day_diff5, day_diff6, day_diff7))

# indicator for the visits after 56 days with the same agent
merged.data6<- merged.data6 %>% group_by(id) %>% 
  mutate(week8 = case_when(
    day_with_atd >= 56  ~ 1,
    TRUE ~ 0))

merged.data6$visitweek8 <- with(merged.data6, ave(week8, cumsum(week8 == 0), FUN = cumsum))

#2) Indicator of acceptable benefit
merged.data6<-merged.data6 %>% group_by(id) %>% fill(qids.fup) %>% 
  mutate(accept_benef = case_when (
    qids.fup < 11 & day_with_atd>=56 ~ 1,
    TRUE ~ 0),
    after_accept_benef = cumsum(accept_benef))

# 3) indicator of primary outcome (using QIDS)
merged.data6<-merged.data6 %>% group_by(id) %>%
  mutate(prim_out_qids = case_when (
    qids.fup<=5 & lag(qids.fup<=5) ~ 1,
    TRUE ~ 0)) %>% ungroup()

# 4) Indicator of treatment failure (see definition in main paper)
merged.data6 <- merged.data6 %>% mutate(treatment_failure1 = case_when(
  visitweek8==1 & is.na(changelevel) & accept_benef==0  ~ 1,
  TRUE ~ 0)) %>%
  mutate(treatment_failure2 = case_when (
    lag(treatment_failure1)==1 & visitweek8==2 & is.na(changelevel) & accept_benef==0 ~ 1,
    TRUE ~ 0))

#####################################################################################################################################
# STEP 7. Add indicators for variables in between clinical visits (by phone)
# 1) Indicator of intolerable side effects
# 2) Indicator of use of prohibited medications
# 3) Work productivity (measured over the phone)
# 4) Outcome (measured with HDRS) 
#####################################################################################################################################
# 1) Indicator for: Intolerable side-effects at each level.
le<-read.delim("//DON/allproj/NMIH/data/STARD/le01.txt")
le<-le[-1,]

c<-le %>% filter (lestatus == "1" | lestatus == "2") %>% filter (src_subject_id %in% merged.data6$id) %>%
  filter (! level == "Level 1") %>%
  mutate(
    intol2se = case_when (
      l2se==1 | l2ctse==1 | l2ase==1  ~ 1,
      l2se=="" | l2se == "0" ~ 0,
      TRUE ~ 0),
    intol3se = case_when(
      l3se==1 ~ 1,
      TRUE ~ 0),
    id = src_subject_id) %>% select(id, days_baseline,intol2se, intol3se)
# side effects at baseline (ie, level 1)

c.bas<- le %>% filter (src_subject_id %in% merged.data6$id) %>% filter(l1se==1) %>% mutate(id=src_subject_id,
                                                                                           intol1se=l1se) %>%
  select(id, intol1se)
c.bas<-c.bas[!duplicated(c.bas$id),]

merged.data6<- merge(merged.data6, c.bas, by =c("id"), all=T)
merged.data6<- merge(merged.data6, c, by =c("id", "days_baseline"), all=T) %>% 
  group_by(id) %>%
  arrange(as.numeric(days_baseline), .by_group = TRUE) %>% 
  fill(start_time, .direction = "updown") %>%
  filter (as.numeric(as.character(days_baseline)) >= as.numeric(as.character(start_time))) %>% 
  mutate(intol1se = case_when(
    is.na(intol1se) ~ "0",
    TRUE ~ intol1se),
    intol2se = case_when(
      is.na(intol2se) ~ "0",
      TRUE ~ as.character(intol2se)),
    intol3se = case_when (
      is.na(intol3se) ~ "0",
      TRUE ~ as.character(intol3se)),
    intol = case_when (
      intol2se==1 | intol3se==1 ~ 1,
      TRUE ~ 0)) %>% ungroup()

rm(c, c.bas, notmeds, le)

# 2) Receiving a prohibited medication
roa<-read.delim("//DON/allproj/NMIH/data/STARD/roa01.txt")
notmeds<-roa %>% filter(src_subject_id %in% merged.data6$id) %>% 
  filter (medication1_name %in% c("CYMBALTA",
                                  "CYMBOLTA",
                                  "CYNBALTA",
                                  "DEPAKOTE",
                                  "ESCITALOPRAM",
                                  "FLUOXETINE",
                                  "LAMICTAL",
                                  "LEXAPRO",                                      
                                  "LEXAPRO 10MG QD",
                                  "LEXIPRO",                                              
                                  "NEURONTIN",
                                  "NEUROTIN",                                         
                                  "PAXIL",
                                  "PROZAC",
                                  "STRATERA",
                                  "TEGRETOL") |
            medication2_name %in% c("DEPAKOTE",
                                    "ESCITALOPRAM",                                       
                                    "LAMICTAL",                                               
                                    "LAMITAL",
                                    "LAWICTAL") |                                         
            medication3_name %in% c("CYMBALTA",
                                    "DEPAKOTE",
                                    "ESCITALOPRAM",
                                    "TEGRATOL")) %>%
  mutate (id = src_subject_id,
          proh.meds = 1) %>% select(id, days_baseline, proh.meds)

merged.data7<-merge(merged.data6, notmeds, by = c("id", "days_baseline"), all = T) %>% group_by(id) %>%
  arrange(as.numeric(days_baseline), .by_group = T) %>% ungroup()

# 3) Productivity scores (measured phone)
wpai<-read.delim("//DON/allproj/NMIH/data/STARD/wpai01.txt")
wpai<-wpai[-1,]

wpai<-wpai %>% mutate (id = src_subject_id) %>% mutate(days_baseline=as.numeric(as.character(days_baseline))) %>%
  filter (id %in% merged.data7$id) %>%  group_by(id) %>% arrange(days_baseline, .by_group=T) %>%
  mutate (WORK = as.numeric(as.character(wpai_pctactimp)))

#latest measurement in level 1 will serve as baseline measurement
wpai.bas<- wpai %>% filter (level == "Level 1") %>%
  group_by(id) %>% slice(c(n())) %>% ungroup() %>%  mutate(bas.WORK = WORK) %>% select(id, bas.WORK) # 1412 with baseline measurement

merged.data7<-merge(merged.data7, wpai.bas, by = "id", all=T)

#extract follow-up measurements of work productivity
wpai.fup<- wpai %>% filter (! level == "Level 1",
                            id %in% merged.data6$id) %>%
  select (id, days_baseline, WORK)

merged.data8<-merge(merged.data7, wpai.fup, by = c("id", "days_baseline"), all = T)

# replace first row with 0 if missing and eliminate rows before time zero
merged.data8<-merged.data8 %>%
  group_by(id) %>%  fill(start_time, .direction = "updown") %>%  
  filter (as.numeric(as.character(days_baseline)) >= as.numeric(as.character(start_time))) %>% group_by(id) %>%
  arrange(as.numeric(days_baseline), .by_group = T) %>%
  mutate(WORK = ifelse(row_number() == 1 & (is.na(WORK)), as.numeric(as.character(bas.WORK)), 
                       as.numeric(as.character(WORK))))%>% 
  ungroup() 

rm(baseline, 
   elig, merged.data2, merged.data3, merged.data4, merged.data5,
   merged.data6, merged.data7, notmeds, prise, qids.bas.imp,
   roa, wpai, wpai.bas, wpai.fup)


#4) Add outcome (measured with HRDS)
hrsd<-read.delim("//DON/allproj/NMIH/data/STARD/hrsd01.txt")
hrsd<-hrsd[-1,]   

# extract all HRDS per individual.
hrsd<-hrsd %>% mutate (id = src_subject_id) %>% mutate(days_baseline=as.numeric(as.character(days_baseline))) %>%
  filter (id %in% merged.data8$id) %>% filter (!level %in% c("Enrollment", "Level 1")) %>% 
  replace(is.na(.), 0)  %>%
  mutate (
    hamd1 = as.numeric(as.character(hsoin)),
    hamd2 = as.numeric(as.character(hmnin)),
    hamd3 = as.numeric(as.character(hemin)),
    hamd4 = as.numeric(as.character(hmdsd)),
    hamd5 = as.numeric(as.character(hpanx)),
    hamd6 = as.numeric(as.character(hinsg)),
    hamd7 = as.numeric(as.character(happt)),
    hamd8 = as.numeric(as.character(hwl)),
    hamd9 = as.numeric(as.character(hsanx)),
    hamd10 = as.numeric(as.character(hhypc)),
    hamd11 = as.numeric(as.character(hvwsf)),
    hamd12 = as.numeric(as.character(hsuic)),
    hamd13 = as.numeric(as.character(hintr)),
    hamd14 = as.numeric(as.character(hengy)),
    hamd15 = as.numeric(as.character(hslow)),
    hamd16 = as.numeric(as.character(hagit)),
    hamd17 = as.numeric(as.character(hsex))) %>%  mutate (hrsdtot = hamd1+hamd2+hamd3+hamd4+hamd5+hamd6+hamd7+hamd8+hamd9+hamd10+hamd11+hamd12+hamd13+hamd14+
                                                            hamd15+hamd16+hamd17)%>% mutate (hrsdtot = as.numeric(as.character(hrsdtot)))  %>%
  mutate(remission_hrsd = ifelse(hrsdtot<=7,1,0),
         measured_hrsd = 1) %>%
  select(id, days_baseline, measured_hrsd, remission_hrsd)
hrsd<- hrsd[!duplicated(hrsd[1:2]),]

merged.data9<-merge (merged.data8, hrsd, by = c("id", "days_baseline"), all=T) %>% 
  group_by(id) %>% 
  fill(start_time, .direction = "updown") %>%
  filter (as.numeric(as.character(days_baseline)) >= as.numeric(as.character(start_time))) %>%
  arrange(as.numeric(days_baseline), .by_group = T) %>%
  mutate(measured_hrsd= ifelse(is.na(measured_hrsd), 0, measured_hrsd)) %>%
  ungroup()

rm(hrsd,merged.data8)

#####################################################################################################################################
# STEP 10. Use the indicators above to define protocol violations:    
#                        a) Receiving ATD, AP or MS 
#                        b) Switching medications too early (without side effects)
#                        c) No treatment change despite treatment failure              
####################################################################################################################################

merged.data9 <- merged.data9 %>% group_by(id) %>% fill(start_time, .direction = "updown") %>%
  mutate (days_study = as.numeric(days_baseline) - as.numeric(as.character(start_time))) %>%
  group_by(id) %>% mutate(cumhdrs = cumsum(measured_hrsd),
                          censorhdrs = case_when (
                            days_study >=105 & cumhdrs==0 ~ 1,
                            TRUE ~ 0))

merged.data9 <- merged.data9 %>% group_by(id) %>% fill(start_time, .direction = "updown") %>%
  mutate (days_study = as.numeric(days_baseline) - as.numeric(as.character(start_time)),
          violation  = case_when (
            is.na(visit) ~ 0, # assessments made outside of a clinical visit have NA: those cannot be protocol deviations
            proh.meds == 1 ~ 1,# received a not allowed medication
            intol==0 & changelevel==1 & visitweek8==0 & days_study != 0 ~ 2, # changed atd before week 8
            treatment_failure2 ==1 ~ 3, # did not change by the end of the wait and see period.
            intol==1 & is.na(changelevel) ~ 4, #intolerable side effects and did not switch
            censorhdrs ==1 ~ 5, # no HDRS within 90 days
            TRUE ~ 0))

table (merged.data9$violation)

#####################################################################################################################################
# STEP 11. Expand 180 days rows per individual              
####################################################################################################################################

# truncate follow-up at 270 days
merged.data10<-merged.data9 %>% filter(days_study<=270)

temp.grid <- data.frame(id=rep(unique(merged.data9$id), each=271),
                        days_study=rep(0:270, times=length(unique(merged.data9$id)))) # grid of all weeks for all ids

merged.data10 <- merge(merged.data10, temp.grid, by=c("id", "days_study"), all = T)

finaldata = merged.data10
#####################################################################################################################################
# STEP 12. Censoring due to LTFU
# Defined as a gap between visits greater than or equal 28 days, before acceptable benefit.
# Need to recreate day with atd and after_accept_ben because they were created only in clinical visits
# while LTFU can occur at any day between 0 and 270
####################################################################################################################################
# create days since last visit
finaldata$visit[is.na(finaldata$visit)]<-0

#consecutive days without a clinical visit:
finaldata<-finaldata %>% mutate(flag = visit==1) %>% 
  group_by(id, grp = with(rle(flag), rep(seq_along(lengths), lengths))) %>% 
  mutate (days_since_last = seq_along(grp)) %>%
  ungroup() %>%
  select(-c(grp))

# days with atd
finaldata<-finaldata %>% group_by(id) %>% fill(txassign, .direction = "down") %>%
  mutate(change = ifelse(is.na(changelevel), 0, changelevel),
         change2 = lag(change),
         change2 = ifelse(is.na(change2), 0, change2),
         day_diff = as.numeric(days_study) - lag(as.numeric(days_study)),
         level2 = cumsum(change2),
         day_diff2 = case_when(
           level2 == 2 ~ day_diff,
           TRUE ~ 0),
         day_diff3 = cumsum(day_diff2),
         day_diff4 = case_when(
           level2 == 3 ~ day_diff,
           TRUE ~ 0),
         day_diff5 = cumsum(day_diff4),
         day_diff6 = case_when(
           level2 == 4 ~ day_diff,
           TRUE ~ 0),
         day_diff7 = cumsum (day_diff6),
         day_with_atd = case_when(
           level2 == 0 ~ 0,
           level2 == 1 ~ days_study,
           level2 == 2 ~ day_diff3,
           level2 == 3 ~ day_diff5,
           level2 == 4 ~ day_diff7)) %>% ungroup() %>%
  select(-c(level2, day_diff, day_diff2, day_diff3, day_diff4, day_diff5, day_diff6, day_diff7))

#censoring after 28 days without a visit (in the absence of adequate benefit, see paper)
#do not censor those with data at EXACTLY 270 days

finaldata<-finaldata %>% group_by(id) %>% fill(c(start_time), .direction = "updown") %>% fill(qids.fup, .direction= "down")%>%
  mutate(lastfollowup = max(as.numeric(days_study)),
         accept_benef = case_when (
           qids.fup < 11 & day_with_atd>=56 ~ 1,
           TRUE ~ 0),
         after_accept_benef = cumsum(accept_benef),
         censoring_event = case_when(
           days_since_last == 28 & after_accept_benef==0 & days_study<270  ~ 1,
           lastfollowup == days_study & days_study<270  ~ 1,
           TRUE ~ 0)) %>% ungroup()

finaldata <- finaldata %>% group_by(id) %>% arrange(days_study, .by_group=T) %>%
  fill (c(age, gender, empl, educat, female, ethnicity, marital, melancholic, numberdsm, ageonset, bas.WORK, self.qids, onlyct,
          idsctot, visitweek8, after_accept_benef, censoring_event, Setting,
          start_time, durationmde, comorbidago, comorbidgad, comorbidocd, comorbidptsd, baseline.qids, anxiousdepr,atypical,
          cirstot, txassign, baseline.hrsd, txassign, sedation, sexual, extrapiramidal, qids.fup, dep.improv, dep.improv2,   
          self.qids, WORK, intol1se, intol2se, intol3se, baseline.qids.self, intol, baseline_alcoh), .direction="down") 

#####################################################################################################################################
# STEP 13. Remove rows after outcome and censoring (first of gap/ltfu)     
####################################################################################################################################
finaldata<-finaldata %>% mutate (violation2 = case_when(
  violation %in% c("1", "2", "3", "4", "5") ~ 1,
  TRUE ~ 0),
  primary_outcome = case_when (
    prim_out_qids==1 | remission_hrsd==1 ~ 1,
    TRUE ~ 0))%>% group_by(id) %>%
  mutate(
    todelete1 = lag(primary_outcome),
    todelete2 = lag(violation2),
    todelete3 = lag(censoring_event)) %>% ungroup()
finaldata$todelete1[is.na(finaldata$todelete1)]<- 0
finaldata$todelete2[is.na(finaldata$todelete2)]<- 0
finaldata$todelete3[is.na(finaldata$todelete3)]<- 0


finaldata<-finaldata%>%  group_by(id) %>% mutate(todeletediscont = cumsum (todelete1),
                                                 todeleteviolation = cumsum(todelete2),
                                                 todeletecensor = cumsum(todelete3)) %>%
  filter(todeletediscont == 0 & todeletecensor==0 & todeleteviolation==0) %>% # include here: & todeleteviolation==0
  ungroup() %>% 
  select(-c(todelete1, todelete2, todelete3, todeletediscont, todeletecensor, todeleteviolation)) %>% ungroup() # include here todelete2, todeleteviolation

rm(merged.data10, merged.data9, temp.grid)

#####################################################################################################################################
# STEP 14. Carry forward covariates and exclude those with missing baseline covariates 
####################################################################################################################################

#  transformation to numerical
finaldata<-finaldata %>% 
  mutate_at(c('days_study', 'age', 'numberdsm', 'ageonset', 'cirstot', 'educat', 'baseline.hrsd', 'idsctot',
              'baseline.qids.self', 'dep.improv', 'qids.fup', 'self.qids', 'WORK', 'bas.WORK'), as.numeric) 

ids_with_missing_baseline <- 
  finaldata %>% 
  filter(days_study==0) %>% 
  {.[which(
    is.na(.$ageonset) |
      is.na(.$baseline.qids.self)|
      is.na(.$Setting)|
      is.na(.$bas.WORK) |
      is.na(.$melancholic) |
      is.na(.$ethnicity) |
      is.na(.$age) |
      is.na(.$idsctot) |
      is.na(.$female) |
      is.na(.$empl) |
      is.na(.$educat) |
      is.na(.$atypical) |
      is.na(.$numberdsm) |
      is.na(.$anxiousdepr) |
      is.na(.$durationmde)|
      is.na(.$comorbidago) |
      is.na(.$comorbidocd) |
      is.na(.$comorbidptsd) |
      is.na(.$comorbidgad) |
      is.na(.$cirstot) |
      is.na(.$baseline.hrsd)
  ),]$id}

finaldata <- 
  finaldata %>% filter(!id %in% ids_with_missing_baseline)

finaldata %>% summarise(n_distinct(id)) # 1388 included patients

# remove those with incompatble agents for target trial emulation
a<-finaldata %>% group_by(id) %>% slice(c(1)) %>% ungroup() %>% filter(txassign %in% c("CIT+BUS", "CT", "CIT+CT"))

finaldata<-finaldata %>% filter(!id %in% a$id)
finaldata %>% summarise(n_distinct(id)) # 971 included patients

rm(a, ids_with_missing_baseline, clinics1)

########################
# finaldata for analysis
########################
finaldata<-finaldata %>% 
  select(id, days_study,txassign, age, gender, empl, educat, marital, female, ethnicity, onlyct,idsctot,
         numberdsm, ageonset,durationmde, comorbidago, comorbidocd, comorbidptsd, comorbidgad, Setting,mild, moderate, severe,
         cirstot,melancholic,atypical, anxiousdepr, baseline.hrsd,baseline.qids, changelevel, baseline_alcoh,
         baseline.qids.self, bas.WORK, intol1se, intol2se,intol3se,qids.fup, self.qids, dep.improv2, 
         sedation, sexual, extrapiramidal, WORK, primary_outcome, violation, violation2, censoring_event, day_with_atd,
         visitweek8, visit, lastfollowup, intol, prim_out_qids, after_accept_benef, days_since_last)

####################################################################################################
# One final step: create variables required for the data analysis:
# 1) treatment received at the first failure "firsttrt"
# 2) treatment received at the second failure "secondtrt"
# 3) treatment received at the third failure "thirdtrt"
# 4) change of agents "changetrt"
# 6) days using atd in the final dataset
###################################################################################################

finaldata <- finaldata %>% mutate(failure = case_when(
  txassign %in% c("CIT+CT", "CT", "BUP", "CIT+BUP", "CIT+BUS", "VEN",
                  "SER") ~ "1",
  txassign %in% c("CIT+LI", "CIT+THY", "MIRT","NTP", "SER+LI", "BUP+LI",
                  "VEN+LI", "SER+THY", "BUP+THY", "VEN+THY") ~ "2",
  txassign %in% c("VEN+MIRT", "TCP") ~ "3"),
  firsttrt = case_when (
    failure == 1 & txassign == "CIT+CT" ~ "1",
    failure == 1 & txassign == "CT" ~ "2",
    failure == 1 & txassign == "BUP" ~ "3",
    failure == 1 & txassign == "CIT+BUP" ~ "4",
    failure == 1 & txassign == "CIT+BUS" ~ "5",
    failure == 1 & txassign == "VEN" ~ "6",
    failure == 1 & txassign == "SER" ~ "7",
    TRUE ~ NA_character_),
  secondtrt = case_when (
    failure == 2 & txassign == "MIRT" ~ "1",
    failure == 2 & txassign == "NTP" ~ "2",
    failure == 2 & txassign %in% c("CIT+LI", "SER+LI", "BUP+LI", "VEN+LI") ~ "3",
    failure == 2 & txassign %in% c("CIT+THY","SER+THY", "BUP+THY", "VEN+THY") ~ "3",
    TRUE ~ NA_character_),
  thirdtrt = case_when (
    failure == 3 & txassign == "VEN+MIRT" ~ "1",
    failure == 3 & txassign == "TCP" ~ "0"),
  changetrt = case_when (
    changelevel==1 & days_study!=0 ~ 1,
    TRUE ~ 0))

finaldata <- finaldata %>% group_by(id) %>% fill (c(firsttrt, secondtrt, thirdtrt), .direction = "down") %>% ungroup()

# days taking antidepressant
# first, take the day of second treatment change
a<-finaldata %>% filter(changetrt==1 & failure==2) %>% mutate(day_second_tx = days_study) %>% select(id, day_second_tx)
b<-finaldata %>% filter(changetrt==1 & failure==3) %>% mutate(day_third_tx = days_study) %>% select(id, day_third_tx)
b<-b[!duplicated(b[1:2]),]

finaldata<-merge(finaldata, a, by = "id", all = T)
finaldata<-merge(finaldata, b, by = "id", all = T)
finaldata<-finaldata%>%group_by(id)%>%arrange(days_study, .by_group = T)

finaldata <- finaldata %>% group_by(id) %>% mutate(
  day_with_atd = case_when(
    failure==1 ~ days_study,
    failure==2 & changetrt==1 ~ days_study,
    failure==2 & changetrt==0 ~ days_study - day_second_tx,
    failure==3 & changetrt==1 ~ days_study - day_second_tx,
    failure==3 & changetrt==0 ~ days_study - day_third_tx),
  change.qids = qids.fup - baseline.qids) %>% select(-c(day_second_tx, day_third_tx))

rm(a,b)
