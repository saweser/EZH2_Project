##### Co occuring mutations with EZH2

rm(list=ls());  # empty workspace
setwd("/data/htp/A07/EZH2_Project/")

library(readxl)
library(dplyr)
library(tibble)

datKlaus<-read_xlsx("/data/htp/A07/EZH2_Project/EZH2_molecular_associations.xlsx", sheet = 1)
  

dat<-read_xlsx("/data/htp/A07/EZH2_Project/EZH2_molecular_associations.xlsx", sheet = 1) %>%
  select(1:3) %>%
  column_to_rownames("Gene")

dat2<-dat %>% t() %>%
  as.data.frame() %>% 
  rownames_to_column("EZH2") %>%
  mutate(total=c(639, 25)) %>% 
  column_to_rownames("EZH2")

# interesting genes DNMT3A, TET2, RUNX1, ASXL1
# DNMT3A
datDNMT3A<-dat2 %>% 
  select(DNMT3A, total) %>%
  rownames_to_column("EZH2") %>%
  mutate(DNMT3A_WT=total-DNMT3A) %>%
  column_to_rownames("EZH2") %>%
  select(-total)

chisq.test(datDNMT3A) # p:0.2984
fisher.test(datDNMT3A) # p: 0.2737

# TET2
datTET2<-dat2 %>% 
  select(TET2, total) %>%
  rownames_to_column("EZH2") %>%
  mutate(TET2_WT=total-TET2) %>%
  column_to_rownames("EZH2") %>%
  select(-total)

chisq.test(datTET2) # p: 0.6753
fisher.test(datTET2) # p: 0.5653

# RUNX1
datRUNX1<-dat2 %>% 
  select(RUNX1, total) %>%
  rownames_to_column("EZH2") %>%
  mutate(RUNX1_WT=total-RUNX1) %>%
  column_to_rownames("EZH2") %>%
  select(-total)

chisq.test(datRUNX1, simulate.p.value = TRUE) # p:0.0001663
fisher.test(datRUNX1) # p: 0.0004579

# ASXL1
datASXL1<-dat2 %>% 
  select(ASXL1, total) %>%
  rownames_to_column("EZH2") %>%
  mutate(ASXL1_WT=total-ASXL1) %>%
  column_to_rownames("EZH2") %>%
  select(-total)

chisq.test(datASXL1, simulate.p.value = TRUE) # p:0.0004998
fisher.test(datASXL1) # p: 9.265e-05

# NPM1 (most often mutated in cohort)
datNPM1<-dat2 %>% 
  select(NPM1, total) %>%
  rownames_to_column("EZH2") %>%
  mutate(NPM1_WT=total-NPM1) %>%
  column_to_rownames("EZH2") %>%
  select(-total)

fisher.test(datNPM1) # p: 0.02816


# when there is one expected value below 5 in the table it will throw the error
# Chi-squared approximation may be incorrect
# then use Fishers exact test!


#### KDM6A is mutually exclusive?
dat646<-read_xlsx("/data/htp/A07/EZH2_Project/686-sup-document2 (1).xlsx", sheet = 2)

kdm<-dat646 %>%
  filter(Gene=="KDM6A")

ezh<-dat646 %>%
  filter(Gene=="EZH2")

nrow(kdm)
nrow(ezh)

intersect(kdm$Sample, ezh$Sample)
setdiff(kdm$Sample, ezh$Sample)

kdm$Sample
ezh$Sample

