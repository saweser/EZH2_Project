
##### reexpression of EZH2 mutants in HEK293T EZH2 negative cells
# 

rm(list=ls());  # empty workspace

library(readxl)


aa<-read.table("/data/htp/A07/EZH2_Project/revision/A692G.txt", header=T, sep = "", dec = ",")


##### without WT #############################################################################

a692<-aa %>% t() %>%
  as.data.frame()%>%
  rownames_to_column("name")%>%
  mutate(state=c(rep("KO",5),rep("WT",5),rep("MUT",5)))%>%
  filter(!state=="WT")
colnames(a692)<-c("name","value","state")

# take out the KO we anyway know that is 0 from the westernblot
# take out the mutations we dont have in the paper

t_tests = a692 %>%
  group_by(state) %>%
  summarise(P = t.test(value, mu = 1)$p.value,
            Sig = ifelse(P < 0.05 & P> 0.01, "*", ifelse(P<0.01 & P>0.001,"**", ifelse(P<0.001,"***",""))),
            MaxWidth = max(value),
            sd=sd(value),
            mean=mean(value),
            dist=mean+sd/2)

p=ggplot(data=a692, aes(x=state, y=value )) +
  geom_dotplot(binaxis = "y", stackdir = "center",position = "dodge", size=5)+
  theme_classic()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major=element_blank())+
  #geom_hline(yintercept = 1, size=0.1)+
  theme(legend.position = NULL)+
  labs(y="Relative H3K27me3")+
  theme(axis.title.x = element_blank())+
  stat_summary(
    geom = "crossbar",
    fun = "mean",
    col = "black",
    width = 0.2, size=.1)+
  stat_summary(fun.data = mean_sd, geom = "errorbar", width=0.2, size=.4)+
  geom_text(aes(label = Sig, y =1.01), size = 7,
            data = t_tests, color="grey40")+
  theme(axis.text.x = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.title.y = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.text.y = element_text(size = 12, face = "bold", color = "black"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0.1,0))

p

ggsave("A962G.pdf", plot = p, width = 70, height = 80, units = "mm", device = "pdf", 
       path = "/data/htp/A07/EZH2_Project/revision")

##### with WT ############################
a692<-aa %>% t() %>%
  as.data.frame()%>%
  rownames_to_column("name")%>%
  mutate(state=c(rep("KO",5),rep("WT",5),rep("MUT",5)))
colnames(a692)<-c("name","value","state")

# take out the KO we anyway know that is 0 from the westernblot
# take out the mutations we dont have in the paper

t_tests = a692 %>%
  group_by(state) %>%
  summarise(P = t.test(value, mu = 1)$p.value,
            Sig = ifelse(P < 0.05, "*", ""),
            MaxWidth = max(value),
            sd=sd(value),
            mean=mean(value),
            dist=mean+sd/2)

my_comparisons <- list( c("KO", "WT"), c("MUT", "WT") )


p=ggplot(data=a692, aes(x=state, y=value )) +
  geom_dotplot(binaxis = "y", stackdir = "center",position = "dodge", size=5)+
  theme_bw()+
  theme(panel.grid.minor = element_blank())+
  theme(panel.grid.major=element_blank())+
  #geom_hline(yintercept = 1, size=0.1)+
  theme(legend.position = NULL)+
  labs(y="Relative H3K27me3")+
  theme(axis.title.x = element_blank())+
  stat_compare_means(label = "p.format", method = "t.test", size=5, label.y = c(1.1, 1.2),
                     method.args = list(var.equal = TRUE), comparisons=my_comparisons)+
  stat_summary(
    geom = "crossbar",
    fun = "mean",
    col = "black",
    width = 0.2, size=.1)+
  stat_summary(fun.data = mean_sd, geom = "errorbar", width=0.2, size=.4)+
  theme(axis.text.x = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.title.y = element_text(size = 14, face = "bold", color = "black"))+
  theme(axis.text.y = element_text(size = 12, face = "bold", color = "black"))+
  theme(axis.title.x = element_blank())+
  scale_y_continuous(expand = c(0.1,0))

p

ggsave("A962G.pdf", plot = p, width = 70, height = 80, units = "mm", device = "pdf", 
       path = "/data/htp/A07/EZH2_Project/revision")











