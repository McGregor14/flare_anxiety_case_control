### cleaning script for case-control analyses

# clear workspace 
rm(list = ls())

# load packages
library('tidyverse')

#  Functions
# Not in operator
'%!in%' <- function(x,y)!('%in%'(x,y))

# set seed to be the same as the one used in analyses
set.seed(31290)

# read in the main TEDS dataset
TEDS <- read_csv(
  file = paste0("raw-data/master_dataset_teds_tm_2020-04-22.csv"), 
  locale = locale("uk"), 
  trim_ws = T
  )

# adjust TEDS columns 
TEDS <- TEDS %>% 
  mutate(
  Sex = factor(as.character(TwinSex), labels=c("female","male")),
  did_consent_factor = factor(did_consent, levels = c("yes", "no"))
  )

# starting with 5914 (received email invitation)
participation_check <- TEDS %>%
  # consent
  filter(did_consent == "yes") %>% # 2987 consented
  # screening
  filter(is.na(screened_out)|screened_out != "yes") %>% # 2925 passed screening
  # completion
  filter(is.na(DidNotCompleteAcq)|DidNotCompleteAcq != 1) %>% # 2826 completed acquisition
  filter(is.na(DidNotCompleteExt)|DidNotCompleteExt != 1) %>% # 2661 completed extinction
  filter(experiment_complete == "yes") %>% # 2543 completed
  filter(is.na(USunp)|USunp > 5) %>% # 2363 found US aversive (should = TEDS_complete)
  
  filter(is.na(Restd1_wk1)|Restd1_wk1 != 1) %>% # 2245 didn't restart
  filter(is.na(volmean)|volmean > .5) %>% # 2102 didn't lower volume
  filter(is.na(Headphone_remove)|Headphone_remove != "yes") %>%  # 1662 didn't remove headphones
  filter(is.na(ExpInstructions)|ExpInstructions != "no") # 1625 followed the instructions
  # 1625 (95%.same 50%.miss needs including, should = number of NOs for exclusion variable)

TEDS_complete <- TEDS %>%
  # consent
  filter(did_consent == "yes") %>% # 2987 consented
  # screening
  filter(is.na(screened_out)|screened_out != "yes") %>% # 2925 passed screening
  # completion
  filter(is.na(DidNotCompleteAcq)|DidNotCompleteAcq != 1) %>% # 2826 completed acquisition
  filter(is.na(DidNotCompleteExt)|DidNotCompleteExt != 1) %>% # 2661 completed extinction
  filter(experiment_complete == "yes") %>% # 2543 completed
  filter(is.na(USunp)|USunp > 5) # 2363 found US aversive (should = TEDS_complete)

# create new variables needed for analyses
df1 <- TEDS_complete %>%
  mutate(
    exclude = 
      ifelse(
        Subject_ID %!in% participation_check$Subject_ID,
        "yes",
        "no"
        ),
    
    # create affective rating differentials
    # baseline to post-acquisition
    FCQ0_1aroCSp = FCQ1AroCSp - FCQ0AroCSp,
    FCQ0_1aroCSm = FCQ1AroCSm - FCQ0AroCSm,
    FCQ0_1valCSp = FCQ1ValCSp - FCQ0ValCSp,
    FCQ0_1valCSm = FCQ1ValCSm - FCQ0ValCSm,
    FCQ0_1feaCSp = FCQ1FeaCSp - FCQ0FeaCSp,
    FCQ0_1feaCSm = FCQ1FeaCSm - FCQ0FeaCSm,
    
    # post-acquisition to post-extinction
    # 'FCQ2ValCSp' missing
    
    # FCQ1_2aroCSp = FCQ2AroCSp - FCQ1AroCSp,
    # FCQ1_2aroCSm = FCQ2AroCSm - FCQ1AroCSm,
    # FCQ1_2valCSp = FCQ2ValCSp - FCQ1ValCSp,
    # FCQ1_2valCSm = FCQ2ValCSm - FCQ1ValCSm,
    # FCQ1_2feaCSp = FCQ2FeaCSp - FCQ1FeaCSp,
    # FCQ1_2feaCSm = FCQ2FeaCSm - FCQ1FeaCSm,
    
    # rename affective rating composites to match data wrangling later on
    FCQ1AffCSp = FCQ1CSp,
    FCQ1AffCSm = FCQ1CSm,
    FCQ2AffCSp = FCQ2CSp,
    FCQ2AffCSm = FCQ2CSm,
    FCQ1AffDif = FCQ1AffCSp - FCQ1AffCSm,
    FCQ2AffDif = FCQ2AffCSp - FCQ2AffCSm,
    
    # create phenotypic variables for anxiety and depression
    # MHQ.Anx: any anxiety diagnosis based on DSM-IV and Puck's paper
    MHQ.Anx = 
      if_else(
        condition =
          MHQ.AnxGen == 1 |
          MHQ.AnxSoc == 1 | 
          MHQ.AnxPho == 1 | 
          MHQ.AnxAgo == 1 | 
          MHQ.AnxPA == 1 | 
          MHQ.PTSD == 1 | 
          MHQ.OCD == 1 | 
          MHQ.OCDoth == 1,
        true = 1,
        false = 0,
        missing = 0),
    
    # PHQ.dia: diagnosis based on PHQ severity scores of more than 2
    PHQ.dia = 
      ifelse(
        PHQ.sev >= 2,
        1,
        0
        ),
    
    # anyanx: either a self-report diagnosis via MHQ or a GAD diagnosis
    anyanx = 
      ifelse(
        MHQ.Anx == 1 | GAD.dia == 1,
        1, 
        0
        ),
    
    # anydep: either a self-report diagnosis via MHQ or a PHQ diagnosis
    anydep = 
      ifelse(
        MHQ.Dep == 1 | PHQ.dia == 1, 
        1, 
        0
        ), 
    
    # control: no MHQ, PHQ, or GAD diagnosis
    control = 
      ifelse(
        GAD.sev == 0 &
          PHQ.sev == 0 &
          MHQ.Dep == 0 &
          MHQ.pnDep == 0 &
          MHQ.Mania == 0 &
          MHQ.AnxGen == 0 & 
          MHQ.AnxSoc == 0 & 
          MHQ.AnxPho == 0 & 
          MHQ.AnxAgo == 0 & 
          MHQ.AnxPA == 0 & 
          MHQ.PTSD == 0 & 
          MHQ.OCD == 0 & 
          MHQ.OCDoth == 0 & 
          MHQ.BDD == 0 & 
          MHQ.ED.AN == 0 & 
          MHQ.ED.BN == 0 & 
          MHQ.ED.bng == 0 & 
          MHQ.SCZ == 0 & 
          MHQ.Psychosis == 0 & 
          MHQ.ASD == 0 & 
          MHQ.AD == 0 & 
          MHQ.PD == 0 & 
          MHQ.Oth == 0, 
        1, 
        0
    ),
    
    # gad.control: no MHQ, PHQ, and GAD < 10
    gad.life.control = 
      ifelse(
        GAD < 10 & 
          PHQ.sev == 0 & 
          MHQ.Dep == 0 & 
          MHQ.pnDep == 0 & 
          MHQ.Mania == 0 & 
          MHQ.AnxGen == 0 & 
          MHQ.AnxSoc == 0 & 
          MHQ.AnxPho == 0 & 
          MHQ.AnxAgo == 0 & 
          MHQ.AnxPA == 0 & 
          MHQ.PTSD == 0 & 
          MHQ.OCD == 0 & 
          MHQ.OCDoth == 0 & 
          MHQ.BDD == 0 & 
          MHQ.ED.AN == 0 & 
          MHQ.ED.BN == 0 & 
          MHQ.ED.bng == 0 & 
          MHQ.SCZ == 0 & 
          MHQ.Psychosis == 0 & 
          MHQ.ASD == 0 & 
          MHQ.AD == 0 & 
          MHQ.PD == 0 & 
          MHQ.Oth == 0, 
        1, 
        0
        ),
    
    # anx_com_cont: anxiety/comorbid/control
    anx_com_cont = 
      factor(
        case_when(
          anyanx == 1 & anydep == 0 ~ "case",
          anyanx == 1 & anydep == 1 ~ "comorbid",
          control == 1 ~ "control",
        )
      ),
    
    # anx_cont: anxiety/control
    anx_cont = 
      factor(
        case_when(
          anyanx == 1 ~ "case",
          control == 1 ~ "control"
        )
      ),
    
    # anx.current_control: current anxiety (GAD)/control
    anx.current_cont = 
      factor(
        case_when(
          GAD.dia == 1 ~ "case",
          control == 1 ~ "control"
        )
      ),
    
    # anx.life_control: anxiety lifetime not current (MHQ not GAD)/control
    anx.life_cont = 
      factor(
        case_when(
          GAD.dia == 0 & MHQ.Anx == 1 ~ "case",
          control == 1 ~ "control"
        )
      ),
    
    # anx.nodep_control: anxiety but no depression/control
    anx.nodep_cont = 
      factor(
        case_when(
          anyanx == 1 & anydep == 0 ~ "case",
          control == 1 ~ "control"
        )
      ),
    
    # general.anx_control: general anxiety/control
    general.anx_control = 
      factor(
        case_when(
          MHQ.AnxGen == 1 ~ "case",
          control == 1 ~ "control"
        )
      ),
    
    # social.anx_control: social anxiety/control
    social.anx_control = 
      factor(
        case_when(
          MHQ.AnxSoc == 1 ~ "case",
          control == 1 ~ "control"
        )
      ),
    
    # panic.anx_control: panic anxiety/control
    panic.anx_control = 
      factor(
        case_when(
          MHQ.AnxPA == 1 ~ "case",
          control == 1 ~ "control"
        )
      ),
    
    # gad.anx_gad.control: gad > 10 / gad < 10
    gad.anx_gad.control = 
      factor(
        case_when(
          GAD.dia == 1 ~ "case",
          GAD.dia == 0 ~ "control"
        )
      ),
    
    # gad.life_gad.life.control: gad > 10 / gad < 10
    gad.life_gad.life.control = 
      factor(
        case_when(
          GAD.dia == 1 | MHQ.Anx == 1 ~ "case",
          gad.life.control == 1 ~ "control"
        )
      )
  )

# create a vector of the names of all the variables of interest
vars_of_interest <- c("FCT1CSpMea", "FCT1CSmMea", "FCT1CSdifMea", 
                      "FCT3CSpMea", "FCT3CSmMea", "FCT3CSdifMea")

# create a vector of the names of all the variables of interest for thirds
vars_of_interest_thirds <- c("FCT1CSpMea.01", "FCT1CSpMea.02","FCT1CSpMea.03", 
                             "FCT1CSmMea.01", "FCT1CSmMea.02", "FCT1CSmMea.03", 
                             "FCT1CSdifMea.01", "FCT1CSdifMea.02", "FCT1CSdifMea.03", 
                             
                             "FCT3CSpMea.01", "FCT3CSpMea.02", "FCT3CSpMea.03", 
                             "FCT3CSmMea.01", "FCT3CSmMea.02", "FCT3CSmMea.03", 
                             "FCT3CSdifMea.01", "FCT3CSdifMea.02", "FCT3CSdifMea.03")

# create a vector of all the trial names for each stimulus to be used for elephant plots
fc_trials <- c("FCT1CSm_1", "FCT1CSm_2", "FCT1CSm_3", "FCT1CSm_4", 
               "FCT1CSm_5", "FCT1CSm_6", "FCT1CSm_7", "FCT1CSm_8", 
               "FCT1CSm_9", "FCT1CSm_10", "FCT1CSm_11", "FCT1CSm_12",
               
               "FCT1CSp_1", "FCT1CSp_2", "FCT1CSp_3", "FCT1CSp_4", 
               "FCT1CSp_5", "FCT1CSp_6", "FCT1CSp_7", "FCT1CSp_8", 
               "FCT1CSp_9", "FCT1CSp_10", "FCT1CSp_11", "FCT1CSp_12",
               
               "FCT3CSm_1", "FCT3CSm_2", "FCT3CSm_3", "FCT3CSm_4", 
               "FCT3CSm_5", "FCT3CSm_6", "FCT3CSm_7", "FCT3CSm_8", 
               "FCT3CSm_9", "FCT3CSm_10", "FCT3CSm_11", "FCT3CSm_12", 
               "FCT3CSm_13", "FCT3CSm_14", "FCT3CSm_15", "FCT3CSm_16", 
               "FCT3CSm_17", "FCT3CSm_18", 
               
               "FCT3CSp_1", "FCT3CSp_2", "FCT3CSp_3", "FCT3CSp_4", 
               "FCT3CSp_5", "FCT3CSp_6", "FCT3CSp_7", "FCT3CSp_8", 
               "FCT3CSp_9", "FCT3CSp_10", "FCT3CSp_11", "FCT3CSp_12", 
               "FCT3CSp_13", "FCT3CSp_14", "FCT3CSp_15", "FCT3CSp_16", 
               "FCT3CSp_17", "FCT3CSp_18")

emailed_n <- nrow(TEDS) # starting with 5914 (received email invitation)

consented_n <- TEDS %>%
  # consent
  filter(did_consent == "yes") %>% nrow() # 2987 consented
  
passed_screening_n <- TEDS %>%
  # consent
  filter(did_consent == "yes") %>% # 2987 consented
  # screening
  filter(is.na(screened_out)|screened_out != "yes") %>% nrow() # 2925 passed screening

complete_n <- nrow(TEDS_complete) # 2363 complete

# save objects to use in markdown scripts
save(df1, file = "data/df1.RData")
save(vars_of_interest, file = "data/vars_of_interest.RData")
save(vars_of_interest_thirds, file = "data/vars_of_interest_thirds.RData")
save(fc_trials, file = "data/fc_trials.RData")
save(emailed_n, file = "data/emailed_n.RData")
save(consented_n, file = "data/consented_n.RData")
save(passed_screening_n, file = "data/passed_screening_n.RData")
save(complete_n, file = "data/complete_n.RData")
save(TEDS, file = "data/TEDS.RData")

# clear workspace
rm(list = ls())
