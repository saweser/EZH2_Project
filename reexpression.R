
##### reexpression of EZH2 mutants in HEK293T EZH2 negative cells
# 

rm(list=ls());  # empty workspace

library(readxl)


lof<-read.table("/data/htp/A07/EZH2_Project/revision/Loss_of_Function.txt", header=T, sep = "")
gof<-read.table("/data/htp/A07/EZH2_Project/revision/Gain_of_Function.txt", header=F, sep = "")

reexp<-read_xlsx("/data/htp/A07/EZH2_Project/revision/Re-Expression Zusammenfassung  in 1g2 KO 293T.xlsx", 
                 sheet = 2) %>%
  filter(!state=="WT")

##### our mutations #############################################################################

# Y646N and Y731F are the  mutations from the other paper
# leave out WT becayse that is anyway set to 1


reEZH2<-reexp %>% filter(Protein=="EZH2") %>%
  filter(!state=="KO") %>%
  filter(!mutation %in% c("A692G"))%>%
  mutate(mutation=gsub("Y730_delins","Y730_delinsX",mutation)) %>%
  mutate(mutation=gsub("c2195_1", "Y733LfsX6", mutation))%>%
  mutate(mutation=factor(mutation, levels = c("K574E","Y730_delinsX","Y733LfsX6","Q612X","G743fs","I744fs","D293G",
                                              "Y646N", "Y731F")))
# take out the KO we anyway know that is 0 from the westernblot
# take out the mutations we dont have in the paper

reEZH2mut<-reEZH2 %>% filter(state=="MUT")
# plot the MUT and MUT+WT separately because not all mutations of MUT+WT are in the paper

##### plot points 
p=ggplot(data=reEZH2mut, aes(x=mutation, y=value )) +
  geom_point(color="black")+
  theme_classic()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major=element_blank())+
  geom_hline(yintercept = 1, size=0.1)+
  theme(legend.position = NULL)+
  labs(y="Relative EZH2 expression")+
  theme(axis.title.x = element_blank())+
  #facet_grid(~mutation, scales="free")+
  stat_summary(
    geom = "crossbar",
    fun = "mean",
    col = "black",
    width = 0.5, size=.1)+
  geom_text(aes(label = Sig, y = 1.2), size = 6,
            data = t_tests)+
  theme(axis.text.x = element_text(size = 10, face = "bold", color = "black", angle = 45, vjust = 0.97, hjust = 1))+
  theme(axis.title.y = element_text(size = 10, face = "bold", color = "black"))+
  theme(axis.text.y = element_text(size = 8, face = "bold", color = "black"))
p

ggsave("EZH2_exp_reexp_point.pdf", plot = p, width = 70, height = 80, units = "mm", device = "pdf", 
       path = "/data/htp/A07/EZH2_Project/revision")

#### plot bar
# t.test, one sample t.test 
t_tests = reEZH2mut %>%
  group_by(mutation) %>%
  summarise(P = t.test(value, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", ""),
            MaxWidth = max(value))

a<-reEZH2mut%>%filter(mutation=="G743fs")
t.test(a$value, mu=1) # padj: 0.03797

p=ggplot(data=reEZH2mut, aes(x=mutation, y=value )) +
  theme_classic()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major=element_blank())+
  scale_y_continuous(breaks = c(0.2,0.6,1),expand = c(0,0))+ 
  geom_hline(yintercept = 1, size=0.1)+
  theme(legend.position = NULL)+
  labs(y="Relative EZH2 expression")+
  theme(axis.title.x = element_blank())+
  #facet_grid(~mutation, scales="free")+
  stat_summary(
    geom = "bar",
    fun = "mean",
    color = "black",
    fill="black",
    width = 0.5, size=1)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.5, size=0.4)+
  geom_text(aes(label = Sig, y = 0.75), size = 6,
            data = t_tests)+
  theme(axis.text.x = element_text(size = 11, face = "bold", color = "black", angle = 45, vjust = 0.97, hjust = 1))+
  theme(axis.title.y = element_text(size = 11, face = "bold", color = "black"))+
  theme(axis.text.y = element_text(size = 9, face = "bold", color = "black"))

p

ggsave("EZH2_exp_reexp.pdf", plot = p, width = 80, height = 70, units = "mm", device = "pdf", 
       path = "/data/htp/A07/EZH2_Project/revision")



##### mutations of the other paper ########################################################################################
### lof
lof<-read.table("/data/htp/A07/EZH2_Project/revision/Loss_of_Function.txt", header=T, sep = "", dec=",")

lof<-lof %>% t() %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  mutate(state=c(rep("KO",3), rep("WT",3), rep("MUT",3), rep("MUT_WT",3)))%>%
  mutate(mutation=c(rep("LOF",12)))%>%
  select(2:4) %>%
  filter(!state=="WT")

colnames(lof)[1]<-"value"

#### gof
gof<-read.table("/data/htp/A07/EZH2_Project/revision/Gain_of_Function.txt", header=T, sep = "", dec = ",")
gof<-gof %>% t() %>%
  as.data.frame() %>%
  rownames_to_column("variable") %>%
  mutate(state=c(rep("KO",4), rep("WT",4), rep("MUT",4), rep("MUT_WT",4)))%>%
  mutate(mutation=c(rep("GOF",16)))%>%
  select(2:4) %>%
  filter(!state=="WT")

colnames(gof)[1]<-"value"


df<-rbind(gof,lof)

t_tests = df %>%
  group_by(mutation,state) %>%
  summarise(P = t.test(value, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", ""),
            MaxWidth = max(value),
            sd=sd(value),
            mean=mean(value),
            dist=mean+sd/2)

p=ggplot(data=df, aes(x=state, y=value )) +
  #geom_point(color="black", size=2)+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major=element_blank())+
  #geom_hline(yintercept = 1, size=0.1)+
  theme(legend.position = NULL)+
  labs(y="Relative H3K27me3")+
  theme(axis.title.x = element_blank())+
  facet_wrap(~mutation, scales="free")+
  stat_summary(
    geom = "point",
    fun = "mean",
    col = "black",
    width = 0.3, size=2)+
  stat_summary(fun.data = mean_se, geom = "errorbar", width=0.3, size=.4)+
  geom_text(aes(label = Sig, y =dist+0.05), size = 7,
          data = t_tests, color="grey40")+
  theme(axis.text.x = element_text(size = 12, face = "bold", color = "black", angle = 45, vjust = 0.97, hjust = 1))+
  theme(axis.title.y = element_text(size = 12, face = "bold", color = "black"))+
  theme(axis.text.y = element_text(size = 10, face = "bold", color = "black"))+
  scale_y_continuous(expand = c(0.1,0))
  

p

ggsave("H3K27me3_PaperMut_reexp.pdf", plot = p, width = 90, height = 80, units = "mm", device = "pdf", 
       path = "/data/htp/A07/EZH2_Project/revision")



#crossbar

dplyr::group_b


