
##### PDX
# 

rm(list=ls());  # empty workspace

library(readxl)
library(tibble)
library(ggbeeswarm)


pdx<-read.table("/data/htp/A07/EZH2_Project/revision/IC50_ComparisonPDX.txt", header=T, sep = "", dec = ",")

pdx2<-pdx %>% t() %>%
  as.data.frame()%>%
  rownames_to_column("name")%>%
  mutate(state=gsub("\\..*","",name)) %>%
  mutate(V1=as.numeric(V1))

colnames(pdx2)[2]<-"value"

 second table had one replicate more
a<-c("PDX_491.7",	79.93,"PDX_491")
pdx2<-rbind(pdx2, a) %>%mutate(value=as.numeric(value))

p491<-pdx2 %>% filter(state=="PDX_491")
p661<-pdx2 %>% filter(state=="PDX_661")

t.test(p491$value, p661$value, var.equal=TRUE)
var(p491$value)
var(p661$value)

var.test(p491$value, p661$value, var.equal=TRUE)
# here the variance is unequal so welch t-test should be used (r default)
# julia used students t test equal variance everywhere so need to do that now as well 
# 

##### plot
p=ggplot(data=pdx2, aes(x=state, y=value )) +
  geom_dotplot(binaxis = "y", stackdir = "center",position = "dodge", size=5)+
  theme_classic()+
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
  labs(y="IC50 AraC [nM]")+
  #stat_compare_means(label = "p.format", label.y=200, method = "t.test", method.args=list(alternative = "greater"))+
  stat_compare_means(label = "p.signif", label.y=250, method = "t.test", size=7,color="grey40",
                     method.args = list(var.equal = TRUE), comparisons=list(c("PDX_491","PDX_661")))+
  scale_y_continuous(expand = c(0.1,0))

p
  
ggsave("PDX_491_661_IC50.pdf", plot = p, width = 70, height = 80, units = "mm", device = "pdf", 
       path = "/data/htp/A07/EZH2_Project/revision")


