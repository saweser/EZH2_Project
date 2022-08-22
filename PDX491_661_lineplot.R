
##### PDX
# 

rm(list=ls());  # empty workspace

library(readxl)
library(tibble)
library(ggbeeswarm)


pdx491<-read.table("/data/htp/A07/EZH2_Project/revision/PDX491.txt", header=T, sep = "", dec = ",", row.names = 1) %>%
  rownames_to_column("sample")
pdx661<-read.table("/data/htp/A07/EZH2_Project/revision/PDX661.txt", header=T, sep = "", dec = ",", row.names = 1) %>%
  rownames_to_column("sample")

pdx<-rbind(pdx491, pdx661) %>%
  gather("mutation","value",ETV6_p.P214L:JAK1_p.V658F) %>%
  mutate(clonal=c(rep("clonal",21), rep("subclonal",14), rep("clonal",21), rep("subclonal",7))) %>%
  mutate(change=c(rep("stable",21), rep("change",21), rep("stable",21))) %>%
  mutate(organism=ifelse(sample=="R1","Patient", ifelse(sample=="R2","Patient","PDX")))%>%
  mutate(ezh2=ifelse(mutation=="EZH2_p.A692G","yes","no")) %>%
  mutate(sample=gsub("R1","Relapse 1", sample)) %>%
  mutate(sample=gsub("R2","Relapse 2", sample)) %>%
  mutate(sample=gsub("AML-","", sample)) %>%
  mutate(sample=gsub("P0","Primograft", sample)) %>%
  mutate(sample=gsub("P1","1st re-Tx", sample)) %>%
  mutate(sample=gsub("P2","2nd re-Tx", sample)) %>%
  mutate(sample=gsub("P3","3rd re-Tx", sample)) %>%
  mutate(sample=factor(sample, levels=c("Relapse 1", "491_Primograft", "491_2nd re-Tx", "491_3rd re-Tx",
                                        "Relapse 2","661_Primograft","661_1st re-Tx")))
  

# ezh2 separately
ezh2<-pdx %>% filter(mutation=="EZH2_p.A692G")
# all others
others<-pdx %>% filter(!mutation=="EZH2_p.A692G")

# color palette
safe_colorblind_palette <- c("#88CCEE", "#CC6677", "#DDCC77", "#332288", "#AA4499", 
                             "#44AA99", "#999933", "#882255", "#661100", "#6699CC", "#888888")
scales::show_col(safe_colorblind_palette)
#color for ezh2: "#117733"


# ezh2 plot
p=ggplot(data=ezh2, aes(x=sample, y=value, color=mutation)) +
  geom_line(aes(group=mutation), alpha=0.8, size=1) +
  geom_point()+
  scale_color_manual(values=c("#117733"))+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black",angle = 45, hjust = 1))+
  theme(axis.title.y = element_text(size = 12, face = "bold", color = "black"))+
  theme(axis.text.y = element_text(size = 10, face = "bold", color = "black"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  labs(y="VAF [%]")+
  theme(legend.title = element_blank())+
  facet_grid(~organism,scales = "free")+
  scale_y_continuous(expand = c(0.1,0))
p
ggsave("PDX_VAF_ezh2.pdf", plot = p, width = 250, height = 70, units = "mm", device = "pdf", 
       path = "/data/htp/A07/EZH2_Project/revision")


# the others plot
p=ggplot(data=others, aes(x=sample, y=value, color=mutation)) +
  geom_line(aes(group=mutation), alpha=1, size=1) +
  geom_point()+
  scale_color_manual(values=safe_colorblind_palette )+
  theme_bw()+
  theme(axis.title.x = element_blank())+
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black",angle = 45, hjust = 1))+
  theme(axis.title.y = element_text(size = 12, face = "bold", color = "black"))+
  theme(axis.text.y = element_text(size = 10, face = "bold", color = "black"))+
  theme(panel.grid.major = element_blank())+
  theme(panel.grid.minor = element_blank())+
  labs(y="VAF [%]")+
  theme(legend.title = element_blank())+
  facet_grid(~organism,scales = "free")+
  scale_y_continuous(expand = c(0.1,0))
p
ggsave("PDX_VAF_others.pdf", plot = p, width = 250, height = 70, units = "mm", device = "pdf", 
       path = "/data/htp/A07/EZH2_Project/revision")


