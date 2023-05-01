library(tidyverse)
library(ggsci)
setwd('~/Documents/INTRuST_Data')

# script by Luisa Berger, April, 7th, 2023

# In this script we will first exclude participants with missing data and
# order the clinical variables for this project

# load original df provided by Lara Pankatz (November 2022, excel file)
df <- read.csv('INTRuST_Flatfile.csv', na.strings = '')
df[df == -99] <- NA

# first we get a broad summary of the dataframe
df %>% head()
df %>% names()
summary(df)

# deselect all variables that contain BIO or recent (non-MRI variables)
df <- df %>% 
  dplyr::select(-c(contains('bio'))) %>%
  dplyr::select(-c(contains('recent')))

# variable names: clean ID name and only small letters
df <- df %>% 
  mutate(ID_NAA = gsub('-', '_', ID_NAA)) %>%
  dplyr::select(ID_NAA, everything())
names(df) <- tolower((names(df)))
# Manually replace each value step by step
df$id_naa <- gsub("007_TBI_0004", "007_TBI_004", df$id_naa)
df$id_naa <- gsub("007_TBI_0007", "007_TBI_007", df$id_naa)
df$id_naa <- gsub("007_TBI_0008", "007_TBI_008", df$id_naa)
df$id_naa <- gsub("007_TBI_0012", "007_TBI_012", df$id_naa)
df$id_naa <- gsub("007_TBI_0016", "007_TBI_016", df$id_naa)
df$id_naa <- gsub("007_TBI_0019", "007_TBI_019", df$id_naa)
df$id_naa <- gsub("007_TBI_0020", "007_TBI_020", df$id_naa)
df$id_naa <- gsub("007_TBI_0022", "007_TBI_022", df$id_naa)
df$id_naa <- gsub("007_TBI_0024", "007_TBI_024", df$id_naa)
df$id_naa <- gsub("007_TBI_0030", "007_TBI_030", df$id_naa)
df$id_naa <- gsub("007_TBI_0033", "007_TBI_033", df$id_naa)
df$id_naa <- gsub("007_TBI_0034", "007_TBI_034", df$id_naa)
df$id_naa <- gsub("007_TBI_0041", "007_TBI_041", df$id_naa)

### create new pcl variable, based on the summed scores of pcl_m and pcl_c
df <- df %>%
  rowwise() %>% 
  mutate(pcl_c_luisa = ifelse(all(!is.na(c_across(starts_with('pcl_cs')))),
                              sum(c_across(starts_with('pcl_cs')), na.rm = TRUE),
                              NA)) %>%
  mutate(pcl_m_luisa = ifelse(all(!is.na(c_across(starts_with('pcl_m')))),
                              sum(c_across(starts_with('pcl_m')), na.rm = TRUE),
                              NA)) %>%
  mutate(pcl_total_luisa = if (all(is.na(c(pcl_c_luisa, pcl_m_luisa)))) {NA}
         else {sum(pcl_c_luisa, pcl_m_luisa, na.rm = TRUE)})
df <- df %>% dplyr::select(id_naa, 'ptsd.diagnosis',
                    pcl_total_luisa, pcl.c_score_naa, pcl_c_luisa, pcl_m_luisa,
                    contains('pcl_'),everything())

# Exclusion criterias - we have 4 to kick them off :-)

# Criteria 1: MRI data available
df_1 <- df %>%
      # individuals without NAA do not have MRI data (therefore delete)
    filter(!is.na(id_naa)) %>%
      # only participants with 'scandate_complete_naa' == 1 do have MRI data (therefore delete)
    filter(scandate_complete_naa == 1)

# Criteria 2: Exclude participants with missing demographic values (gender, age)
df_2a <- df_1 %>%
  filter(!is.na(gender_naa))

df_2b <- df_2a %>%
  filter(!is.na(age_naa))

# Factor 3: Exclude partcipants with missing pcl-c and/ or ctq
df_3a <- df_2b %>%
  filter(if_all(contains('ctq'), ~ !is.na(.)))

### missing pcl-c
df_3b <- df_3a %>%
  filter(!is.na(pcl_total_luisa))
# identify cases that dont have the same value in these variables
df_pcl_different <- df_3b[!(is.na(df_3b$pcl.c_score_naa)
                               & is.na(df_3b$pcl_total_luisa)) & (is.na(df_3b$pcl.c_score_naa) |
                                                                       is.na(df_3b$pcl_total_luisa) |
                                                                       (df_3b$pcl.c_score_naa != df_3b$pcl_total_luisa)), ]
### missing cdrisc
# df_3c <- df_3b %>%
# filter(!is.na(cdrisc_score_naa)) %>%

# Factor 4: Exclude participants for who harmonization was not possible/
# MRI quality was bad -> only participants with DWI available
df_4 <- df_3b %>%
  # we deselect all subjects that come from sites that were not harmonized.
  filter(!site_naa %in% c(1, 3, 6, 18))
  # we deselect all subjects with bad MRI quality. Therefore compare caselists.
# Read in the Excel file
caselist <- read.csv("~/Desktop/INTRuST_final_caselist.csv")
# Filter the original data frame to include only the IDs from the Excel data
df_clinical <- df_4 %>%
  filter(id_naa %in% caselist$ID_NAA)

# load wma measures provided by Twishi (April 12 2023, csv file)
setwd('~/Documents/INTRuST_Data')
df_mri <- read.csv('new_appended_diffusion_measures_anatomical_tracts.csv',
                   na.strings = '')
df_mri[df_mri == -99] <- NA

# variable names: clean ID names and small letters
df_mri <- df_mri %>% 
  mutate(subjectkey = gsub('harmonized_', '', subjectkey)) %>%
  mutate(subjectkey = gsub('_desc-XcUnEd_dwi_b900', '', subjectkey)) %>%
  mutate(subjectkey = gsub('reconstructed_', '', subjectkey)) %>%
  mutate(subjectkey = gsub('_desc-XcUnEd_dwi', '', subjectkey)) %>%
  dplyr::select(subjectkey, everything())
names(df_mri) <- tolower((names(df_mri)))
# deselect all variables that contain other tracts
df_mri <- df_mri %>% 
  dplyr::select(c(contains('subjectkey'),contains('t_af'),contains('t_cb'),
           contains('t_cc'),contains('t_ilf'),
           contains('t_ioff'),contains('t_uf'),contains('t_slf')))
# deselect tensor 2 because we only need tensor 1 for analysis
df_mri <- df_mri %>% 
  dplyr::select(-c(contains('tensor2')))
# rename variable 'subjectkey' to variable 'id_naa'
df_mri <- df_mri %>% 
  rename("id_naa" = subjectkey)
# merge the two data sets
df_final <- merge(df_clinical, df_mri, by = c("id_naa"))

# exclude problematic case 023_NA3_033 as num.fibers is 0 at CC5
df_final <- df_final %>% filter(!id_naa == "023_NA3_033")

# update the value of ptsd.diagnosis for id_naa = "007_TBI_041"
df_final$ptsd.diagnosis <- ifelse(df_final$id_naa == "007_TBI_041", 1,
                                  df_final$ptsd.diagnosis)

# Create ctq total
df_final <- df_final %>%
  rowwise() %>% 
  mutate(ctq_totalscore_naa = sum(c_across(starts_with(c('ctq_'))),
                                  na.rm = T))

# Create ctq subscale diagnosis (jessica)
df_final$ctq_diagnosis_emo_score_naa <- 
  ifelse(df_final$ctq_emo_score_naa >= 13, 1, 0)
df_final$ctq_diagnosis_phys_score_naa <- 
  ifelse(df_final$ctq_phys_score_naa >= 10, 1, 0)
df_final$ctq_diagnosis_sex_score_naa <- 
  ifelse(df_final$ctq_sex_score_naa >= 8, 1, 0)
df_final$ctq_diagnosis_emoneg_score_naa <- 
  ifelse(df_final$ctq_emoneg_score_naa >= 15, 1, 0)
df_final$ctq_diagnosis_physneg_score_naa <- 
  ifelse(df_final$ctq_physneg_score_naa >= 10, 1, 0)

### Create ctq diagnosis (jessica)
df_final <- df_final %>%
  mutate(ctq_diagnosis_jessica_naa =
  ifelse((ctq_diagnosis_emo_score_naa +
            ctq_diagnosis_phys_score_naa +
            ctq_diagnosis_sex_score_naa +
            ctq_diagnosis_emoneg_score_naa +
            ctq_diagnosis_physneg_score_naa)>0, 2, 0))

### Create ctq diagnosis (soraya)
df_final$ctq_diagnosis_soraya_naa <-
  ifelse(df_final$ctq_totalscore_naa >= 41, 2, 0)
df_final <- df_final %>%
  dplyr::select(id_naa, ctq_totalscore_naa,
                ctq_diagnosis_soraya_naa,
                ctq_diagnosis_jessica_naa,
                starts_with('ctq'), everything())

# create ptsd diagnosis variable on pclc measurements
df_final$ptsd_diagnosis_luisa_naa <-
  ifelse(df_final$pcl_total_luisa > 31, 1, 0)
df_final <- df_final %>%
  dplyr::select(id_naa, ptsd_diagnosis_luisa_naa, ptsd.diagnosis, pcl_total_luisa,
                ctq_diagnosis_jessica_naa,
                ctq_diagnosis_soraya_naa,
                everything())

# Create group variable (soraya)
df_final$ptsd_cm_comorbid_soraya_naa <-
  ifelse(df_final$ptsd.diagnosis == 0 & df_final$ctq_diagnosis_soraya_naa == 0, 0,  # healthy controls
              ifelse(df_final$ptsd.diagnosis == 1 & df_final$ctq_diagnosis_soraya_naa == 0, 1,  # PTSD-only
                   ifelse(df_final$ptsd.diagnosis == 0 & df_final$ctq_diagnosis_soraya_naa == 2, 2,  # CM-only
                    ifelse(df_final$ptsd.diagnosis == 1 & df_final$ctq_diagnosis_soraya_naa == 2, 3,  # CM+PTSD
                          NA))))  # if neither of the conditions above is met, set to NA

# Create group variable (soraya) with new PTSD (based on pclc)
df_final$ptsd_cm_comorbid_newptsd_soraya_naa <-
  ifelse(df_final$ptsd_diagnosis_luisa_naa == 0 & df_final$ctq_diagnosis_soraya_naa == 0, 0,  # healthy controls
         ifelse(df_final$ptsd_diagnosis_luisa_naa == 1 & df_final$ctq_diagnosis_soraya_naa == 0, 1,  # PTSD-only
                ifelse(df_final$ptsd_diagnosis_luisa_naa == 0 & df_final$ctq_diagnosis_soraya_naa== 2, 2,  # CM-only
                       ifelse(df_final$ptsd_diagnosis_luisa_naa == 1 & df_final$ctq_diagnosis_soraya_naa == 2, 3,  # CM+PTSD
                              NA))))  # if neither of the conditions above is met, set to NA

# Create group variable (jessica)
df_final$ptsd_cm_comorbid_jessica_naa <-
  ifelse(df_final$ptsd.diagnosis == 0 & df_final$ctq_diagnosis_jessica_naa == 0, 0,  # healthy controls
         ifelse(df_final$ptsd.diagnosis == 1 & df_final$ctq_diagnosis_jessica_naa == 0, 1,  # PTSD-only
                ifelse(df_final$ptsd.diagnosis == 0 & df_final$ctq_diagnosis_jessica_naa== 2, 2,  # CM-only
                       ifelse(df_final$ptsd.diagnosis == 1 & df_final$ctq_diagnosis_jessica_naa == 2, 3,  # CM+PTSD
                              NA))))  # if neither of the conditions above is met, set to NA

# Create group variable (inga)
df_final$hc_tbi_ptsd_comorbid_naa <-
  ifelse(df_final$ptsd.diagnosis == 0 & df_final$tbi.diagnosis == 0, 0,  # healthy controls
         ifelse(df_final$ptsd.diagnosis == 1 & df_final$tbi.diagnosis == 0, 1,  # PTSD-only
                ifelse(df_final$ptsd.diagnosis == 0 & df_final$tbi.diagnosis == 1, 2,  # tbi-only
                       ifelse(df_final$ptsd.diagnosis == 1 & df_final$tbi.diagnosis == 1, 3,  # CM+PTSD
                              NA))))  # if neither of the conditions above is met, set to NA

# Add labels
df_final$ptsd_cm_comorbid_soraya_naa <- factor(df_final$ptsd_cm_comorbid_soraya_naa,
          levels = c(0, 2, 1, 3),
          labels = c("Healthy controls","CM-only","PTSD-only","CM+PTSD"))
df_final$ctq_diagnosis_soraya_naa <- factor(df_final$ctq_diagnosis_soraya_naa,
          levels = c(0, 2),
          labels = c("No childhood maltreatment","Childhood maltreatment"))
df_final$ptsd_cm_comorbid_jessica_naa <- factor(df_final$ptsd_cm_comorbid_jessica_naa,
                                               levels = c(0, 2, 1, 3),
                                               labels = c("Healthy controls","CM-only","PTSD-only","CM+PTSD"))
df_final$ctq_diagnosis_jessica_naa <- factor(df_final$ctq_diagnosis_jessica_naa,
                                            levels = c(0, 2),
                                            labels = c("No childhood maltreatment","Childhood maltreatment"))
df_final$ptsd.diagnosis <- factor(df_final$ptsd.diagnosis,
                                            levels = c(0, 1),
                                            labels = c("No PTSD",
                                                       "PTSD"))
df_final$hc_tbi_ptsd_comorbid_naa <- factor(df_final$hc_tbi_ptsd_comorbid_naa,
                                               levels = c(0, 2, 1, 3),
                                               labels = c("Healthy controls","mTBI-only","PTSD-only","mTBI+PTSD"))

df_final <- df_final %>%
  dplyr::select(id_naa,
                hc_tbi_ptsd_comorbid_naa,tbi.diagnosis,ptsd.diagnosis,
                ptsd_cm_comorbid_newptsd_soraya_naa,
                ptsd_cm_comorbid_soraya_naa,
                ptsd_cm_comorbid_jessica_naa, ptsd.diagnosis,
                ctq_diagnosis_jessica_naa,
                ctq_diagnosis_soraya_naa,
                everything())

# Add labels for new ptsd
df_final$ptsd_cm_comorbid_newptsd_soraya_naa <- factor(df_final$ptsd_cm_comorbid_newptsd_soraya_naa,
                                               levels = c(0, 2, 1, 3),
                                               labels = c("Healthy controls","CM-only","PTSD-only","CM+PTSD"))
df_final$ptsd_diagnosis_luisa_naa <- factor(df_final$ptsd_diagnosis_luisa_naa,
                                  levels = c(0, 1),
                                  labels = c("No PTSD",
                                             "PTSD"))


summary(df_final$ptsd_cm_comorbid_soraya_naa)
summary(df_final$ptsd_cm_comorbid_jessica_naa)
summary(df_final$ptsd_cm_comorbid_newptsd_soraya_naa)

# Create probable military variable
df_final$probable_military_naa <-
  ifelse(!is.na(df_final$pcl_m_luisa), 1, 0)
df_final <- df_final %>%
  dplyr::select(id_naa, probable_military_naa, ptsd.diagnosis,
                ptsd_cm_comorbid_soraya_naa,
                everything())
df_final$probable_military_naa <- factor(df_final$probable_military_naa,
                                            levels = c(0, 1),
                                            labels = c("No military",
                                                       "Military"))

# check on your variable class
sapply(df_final, class)
df_final$ptsd_cm_comorbid_soraya_naa <- as.factor(df_final$ptsd_cm_comorbid_soraya_naa)
df_final$ptsd_cm_comorbid_newptsd_soraya_naa <- as.factor(df_final$ptsd_cm_comorbid_newptsd_soraya_naa)
df_final$race_naa <- as.factor(df_final$race_naa)
df_final$gender_naa<- as.factor(df_final$gender_naa)
df_final$tbi.diagnosis <- as.factor(df_final$tbi.diagnosis)
df_final$ptsd.diagnosis.source <- as.factor(df_final$ptsd.diagnosis.source)

write.csv(df_final, file = "~/Documents/INTRuST_Data/df_final.csv")

# sort your variables applicable to demographics
df_demographics <- df_final %>%
  dplyr::select(id_naa, ptsd_cm_comorbid_soraya_naa,cdrisc_score_naa,
                ptsd_cm_comorbid_newptsd_soraya_naa,
                ctq_diagnosis_soraya_naa, ptsd.diagnosis,
                race_naa, age_naa, gender_naa,
                pcl_total_luisa, ctq_totalscore_naa,
                ctq_emoneg_score_naa, ctq_physneg_score_naa,
                ctq_emo_score_naa, ctq_phys_score_naa, ctq_sex_score_naa,
                phq9_score_naa,
                tbi.diagnosis, audit_10_score_naa)
sapply(df_demographics, class)

# sort your variables applicable to MRI
df_wma <- df_final %>%
  dplyr::select(id_naa,
                ptsd_cm_comorbid_soraya_naa,age_naa, gender_naa,
                pcl_total_luisa, ctq_totalscore_naa,
                contains("fa1.mean"),
                contains("num_points"),
                contains("mean_length"),
                contains("num_fibers"),
                contains("freewater.mean"),
                contains("maxeigenvalue.mean"),
                contains("meandiffusivity.mean"),
                contains("mineigenvalue.mean"),
                everything())
sapply(df_wma, class)
names(df_wma)

### df_healthycontrols <- df_final %>%
  # filter(tbi.diagnosis == 0) %>%
  # filter(audit_10_score_naa < 8) %>%
  # filter(phq9_score_naa < 10)

### df_noHC <- df_wma %>%
  # filter(!ptsd_cm_comorbid_soraya_naa == "Healthy controls")


### try matching
library(MatchIt)
library(mosaic)
# let's match
# create a balanced data set from an ubalanced data set
df_final$matching_hc_clinical <- ifelse(
  df_final$ptsd_cm_comorbid_soraya_naa == "Healthy controls",0,1)
df_final <- df_final %>%
  dplyr::select(id_naa, matching_hc_clinical,
                ptsd_cm_comorbid_soraya_naa,
                ptsd_cm_comorbid_newptsd_soraya_naa,
                everything())
df_final$matching_hc_clinical <- as.factor(df_final$matching_hc_clinical)
summary(df_final$matching_hc_clinical)

# find matches
df_matching <- df_final %>%
  filter(!is.na(tbi.diagnosis))
  
turnout_match = matchit(matching_hc_clinical ~ age_naa + gender_naa + tbi.diagnosis,
                          data = df_matching, ratio = 1)
summary(turnout_match)
df_matching = match.data(turnout_match)

# check your age
options(contrasts = c("contr.helmert", "contr.poly"))
crf.lm.age <- lm(age_naa~ptsd_cm_comorbid_soraya_naa, data = df_matching)
Anova(crf.lm.age, type=3)

# check on gender
cont_table_gender <- table(df_matching$gender_naa,
                    df_matching$ptsd_cm_comorbid_soraya_naa)
chisq.test(cont_table_gender)

# check on mtbi
cont_table_tbi <- table(df_matching$tbi.diagnosis,
                    df_matching$ptsd_cm_comorbid_soraya_naa)
chisq.test(cont_table_tbi)

# Compute pairwise comparisons

# table_demographics <- df_matching %>%
  # select(ptsd_cm_comorbid_soraya_naa, age_naa) %>%
  # tbl_summary(by = "ptsd_cm_comorbid_soraya_naa",
                                  # percent = "row",
                                  # missing = "no",
                                  # statistic = list(
                                  # all_continuous() ~ "{mean} ± {sd}",
                                  # all_categorical() ~ "{n} ({p}%)"),
                                  # digits = all_continuous() ~ 2,) %>%
  # add_n() %>% add_overall() %>% add_p() %>%
  # add_stat(all_continuous() ~ {add_stat_pairwise}) %>%
  # as_kable()