# Yuchen figures #################

library(readxl)
library(ggplot2)
library(Seurat)
library(dplyr)
library(tidyr)
#library(ggsignif)

setwd("C:/Users/annab/The University of Sydney (Staff)/Thomas O'Neil - TRM")
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/projects/TRM")
list.files()

#Data input, clean up to rename the tissues, calculate percentages we need.

{data <- read_excel("thesis resilts flow.xlsx") 
data$tissue[data$tissue == "AC"] = "ULB"
data$tissue[data$tissue == "Sig"] = "LLB"
data$tissue[data$tissue == "R"] = "LLB"

df <- data.frame(samples = paste0(data$tissue, "_", data$layer),
                  tissue= data$tissue,layer=data$layer,
                 experiment = data$exp)
  
df$CD3per <- 100*data$`CD3 #`/data$`CD45 #`
df$MNPper <- 100*data$`HLADR #`/data$`CD45 #`
df$pDCper <- 100*data$`pDC #`/data$`CD45 #`
df$CD4per <- 100*data$`CD4 #`/data$`CD3 #`
df$CD8per <- 100*data$`CD8 #`/data$`CD3 #`
df$CCR5per <- 100*data$`CCR5 #`/data$`CD4 #`
df$CCR5MFI<-data$`CCR5 MFI`
df$CD69per<-100*data$`CD69 #`/data$`CD4 #`
df$HLADRperCD3<-100*data$`CD3 DR #`/data$`CD3 #`
df$HLADRmfiCD3 <- data$`CD3DR MFI`
df$HLADRperCD4<-100*data$`CD4 DR #`/data$`CD4 #`
df$HLADRmfiCD4<-data$`CD4 DR MFI`

df$CD69SP <- 100*data$CD69pCCR5n/data$`CD4 #`
df$CCR5SP <- 100*data$CD69nCCR5p/data$`CD4 #`
df$DP <- 100*data$CD69pCCR5p/data$`CD4 #`
df$DN <- 100*data$CD69nCCR5n/data$`CD4 #`


df$donor <- paste0(data$exp,"_",data$donor)

df <- na.omit(df)
}

df$samples <- factor(df$samples, levels=c("L_UM", "L_E", 
                                          "V_UM", "V_E",
                                          "LLB_UM", "LLB_LA", "LLB_E", 
                                          "ULB_UM", "ULB_LA", "ULB_E", 
                                          "A_UM", "A_E", "B_B"
))


# final figures ############

# CD4/CD8 ratio in tissues #########
pdf("C:/Users/annab/The University of Sydney (Staff)/Thomas O'Neil - TRM/thesis_plots/CD4CD8_tissue.pdf", width=6,height=4)

## CD4/CD8 ratio in a stack plot of CD4 and CD8 percentages format, scaled up to 100%. 
#We dropped it because it was a bit misleading that CD8 and CD4 didn't actually make up all T cells in the data. 
df %>% 
  filter(tissue != "ULB")%>%
  group_by(samples) %>%summarise(CD8= 100*mean(CD8per)/(mean(CD4per)+mean(CD8per)),CD4 = 100*mean(CD4per)/(mean(CD4per)+mean(CD8per))) %>%
  pivot_longer(cols = c("CD4", "CD8"), names_to = "Subset")%>% 
  ggplot(aes(samples, value, fill=Subset))+geom_col(color = "black", size = 1, width = 0.8)+
  scale_fill_manual(values=c("#2f2f2f", "white"))+theme_classic()+
  geom_hline(yintercept = 50, size=0.7)+
  ggtitle("CD4+ T cells")+RotatedAxis()+
  ylab("% of CD3+ T cells")+
  xlab("")+
  theme(axis.line = element_line(size = 1), axis.text.x = element_text(vjust=0.7), axis.text.y = element_text(face = "bold"))+
  scale_y_continuous(expand = c(0, 0), limits = c(0, 101))
  #ylim(0,101)
dev.off()

## Changed to a log ratio format.
df %>% 
  filter(tissue != "ULB")%>%
  group_by(samples) %>%reframe(MRatio = mean(log(CD4per/CD8per))) %>%
 # pivot_longer(cols = c("CD4", "CD8"), names_to = "Subset")%>% 
  ggplot(aes(samples, MRatio))+geom_col(color = "black", fill = "#2f2f2f", width = 0.8)+
  #scale_fill_manual(values=c("#2f2f2f", "white"))+
  theme_classic()+
  geom_hline(yintercept = 0, size=1)+
  ggtitle("Log(Ratio of CD4:CD8)")+RotatedAxis()+
  ylab("Log(CD4/CD8)")+
  xlab("")+
  theme(axis.line = element_line(size = 1), axis.text.x = element_text(vjust=0.7), axis.text.y = element_text(face = "bold"))
  #scale_y_continuous(expand = c(0, 0))

## To plot the samples for each tissues as dots, and to colour-code by CD4>CD8 or CD4<CD8.

# Step 1: Create the grouped summary
grouped_df <- df %>% 
  filter(tissue != "ULB") %>%
  group_by(samples) %>%
  reframe(Ratio = mean(log(CD4per/CD8per)), SD = sd(log(CD4per/CD8per)))

# Step 2: Plot using the grouped data for bars and the original data for points
ggplot() +
  # Use grouped data for the bars and error bars
  geom_col(data = grouped_df, aes(x = samples, y = Ratio, fill = Ratio > 0), color = 'black', linewidth = 1) +
  geom_errorbar(data = grouped_df, aes(x = samples, ymin = Ratio - SD, ymax = Ratio + SD), width = 0.2, size = 0.8) +
  
  # Add individual points using original df
  geom_point(data = df %>% filter(tissue != "ULB"), 
             aes(x = samples, y = log(CD4per/CD8per)), 
             position = position_jitter(width = 0.2), 
             color = "black", size = 2) +
  
  # Additional elements
  geom_vline(xintercept = 0.5 + c(2, 4, 7, 9), alpha = 0.5) +
  scale_fill_manual(values = c("TRUE" = "red", "FALSE" = "blue")) +
  geom_hline(yintercept = 0, size = 1) +
  theme_classic() +
  ggtitle("Log(Ratio of CD4:CD8)") +
  RotatedAxis() +
  ylab("Log(CD4/CD8)") +
  xlab("") +
  theme(axis.line = element_line(size = 1), 
        axis.text.x = element_text(vjust = 0.7), 
        axis.text.y = element_text(face = "bold"))

# A little summary plot with colour gradient and dot size representing the % of CD4 and CD8 T cells.
pdf("C:/Users/annab/The University of Sydney (Staff)/Thomas O'Neil - TRM/thesis_plots/CD4CD8_tissue_dot.pdf", width = 6, height = 4)

df %>% 
  filter(tissue != "ULB")%>%
  group_by(samples) %>%summarise(CD4 = 100*mean(CD4per)/(mean(CD4per)+mean(CD8per)), CD8= 100*mean(CD8per)/(mean(CD4per)+mean(CD8per))) %>%
  #pivot_longer(cols = c("CD4", "CD8"), names_to = "Subset")%>%
  ggplot(aes(samples, 1, size=CD4, color=CD4))+geom_point()+
  scale_color_gradientn(colors= c("blue","white",'red'),
                        limits = c(0, 100))+
  theme_classic()+RotatedAxis()+ylim(0,2)
dev.off()


#CD69 and CCR5 expression ############

#CD69
pdf("C:/Users/annab/The University of Sydney (Staff)/Thomas O'Neil - TRM/thesis_plots/CD69Percentage.pdf", width=6.5, height=4)
ggplot(df[df$tissue!= "ULB",], aes(y=CD69per, samples, fill=layer))+
  geom_boxplot(outliers = F, show.legend = F, size = 0.8)+geom_jitter(alpha=0.2, show.legend = F)+
  geom_vline(xintercept = 0.5+c(2,4,7,9), alpha=0.5)+geom_hline(yintercept = c(25,50,75), alpha=0.1)+
  scale_fill_manual(values=c("grey", "#FC6A03","black","#4f4f4f"))+
  ggtitle("CD4+ T cells")+
  theme_classic()+
  RotatedAxis()+xlab("")+ylab("% CD69")+ylim(0,100)+
  theme(axis.text.x = element_text(size=12, vjust=1), 
        axis.line.x = element_line(linewidth = 1), 
        axis.line.y  = element_line(linewidth = 1),
        axis.text.y = element_text(face = 'bold'))
dev.off()

#Pre-stats
bartlett.test(CD69per ~samples, data = df %>% filter(tissue == "L"))$p.value #labia E and UM sample number not equal
bartlett.test(CD69per ~samples, data = df %>% filter(tissue == "V"))$p.value #0.0079 have difference in variances
bartlett.test(CD69per ~samples, data = df %>% filter(tissue == "LLB" & layer != "LA"))$p.value #E and UM sample numbner not equal

#Stats
#LUM and LE used unpaired test because one UM does not have matched E. And that UM was not representative for dataset mean, cannot remove it.
t.test(CD69per ~samples, data = df %>% filter(tissue == "L"))$p.value#sig
#t.test(CD69per ~samples, data = df %>% filter(tissue == "V"))$p.value

#filtering VUM and VE samples
vum <- df %>% filter(samples=="V_UM")
ve<- df%>% filter(samples=="V_E")
 #paired t-test comparing vum and ve matched donors
t.test(vum$CD69per, ve$CD69per, paired= TRUE)$p.value #sig 0.0489

#LLB UM and E due to UM and E numbers not the same, and cannot remove the data point.
t.test(CD69per ~samples, data = df %>% filter(tissue == "LLB" & layer != "LA"))$p.value

#ANOVAs across tissues for E or UM
summary(aov(CD69per~samples, data = df %>% filter(layer=="E" & tissue != "A" & tissue!= "ULB")))
TukeyHSD(aov(CD69per~samples, data = df %>% filter(layer=="E" & tissue != "A" & tissue!= "ULB")))

summary(aov(CD69per~samples, data = df %>% filter(layer=="UM" & tissue != "A" & tissue!= "ULB")))#sig
TukeyHSD(aov(CD69per~samples, data = df %>% filter(layer=="UM" & tissue != "A" & tissue!= "ULB")))#sig


#CCR5
pdf("C:/Users/annab/The University of Sydney (Staff)/Thomas O'Neil - TRM/thesis_plots/CCR5Percentage.pdf", width=6.5, height=4)
ggplot(df[df$tissue!= "ULB",], aes(y=CCR5per, samples, fill=layer))+
  geom_boxplot(outliers = F, show.legend = F, size = 0.8)+geom_jitter(alpha=0.2, show.legend = F)+
  geom_vline(xintercept = 0.5+c(2,4,7,9), alpha=0.5)+geom_hline(yintercept = c(25,50,75), alpha=0.2)+
  scale_fill_manual(values=c("grey", "#FC6A03","black","#4f4f4f"))+
  ggtitle("CD4+ T cells")+
  theme_classic()+
  RotatedAxis()+xlab("")+ylab("% CCR5")+ylim(0,100)+
  theme(axis.text.x = element_text(size=12, vjust=1), 
        axis.line.x = element_line(linewidth = 1), 
        axis.line.y  = element_line(linewidth = 1),
        axis.text.y = element_text(face = 'bold'))+ #bold text numbers
        geom_signif(comparisons = list(c("L_E","L_UM")),map_signif_level=TRUE)
                   #add significance
dev.off()

#pre-stats
with(df, shapiro.test(CCR5per[samples == "L_UM"])) #normality test example

# variance tests
bartlett.test(CCR5per ~samples, data = df %>% filter(tissue == "L"))$p.value #labia sample number not equal
bartlett.test(CCR5per ~samples, data = df %>% filter(tissue == "V"))$p.value #0.02 have difference in data distribution
bartlett.test(CCR5per ~samples, data = df %>% filter(tissue == "LLB" & layer != "LA"))$p.value #0.08 sample numbner not equal

#Stats
t.test(CCR5per ~samples, data = df %>% filter(tissue == "L"))$p.value#sig
t.test(CCR5per ~samples, data = df %>% filter(tissue == "V"))$p.value#sig
t.test(vum$CCR5per, ve$CCR5per, paired= TRUE)$p.value #sig
t.test(CCR5per ~samples, data = df %>% filter(tissue == "LLB" & layer != "LA"))$p.value

summary(aov(CCR5per~samples, data = df %>% filter(layer=="E" & tissue != "A" & tissue!= "ULB")))#sig
TukeyHSD(aov(CCR5per~samples, data = df %>% filter(layer=="E" & tissue != "A" & tissue!= "ULB")))#sig

summary(aov(CCR5per~samples, data = df %>% filter(layer=="UM" & tissue != "A" & tissue!= "ULB")))
TukeyHSD(aov(CCR5per~samples, data = df %>% filter(layer=="UM" & tissue != "A" & tissue!= "ULB")))


#CCR5 gMFI
pdf("C:/Users/annab/The University of Sydney (Staff)/Thomas O'Neil - TRM/thesis_plots/CCR5MFI.pdf", width = 6.5, height = 4)
ggplot(df[df$tissue!= "ULB",], aes(y=log(CCR5MFI), samples, fill=layer))+
  geom_boxplot(outliers = F, show.legend = F, size = 0.8)+geom_jitter(alpha=0.3, show.legend = F)+ #added jitter points
  geom_vline(xintercept = 0.5+c(2,4,7,9), alpha=0.5)+
  scale_fill_manual(values=c("grey", "#FC6A03","black","#4f4f4f"))+
  ggtitle("CD4+ T cells")+
  theme_classic()+
  RotatedAxis()+xlab("")+ylab("log(CCR5 MFI)")+
  theme(axis.text.x = element_text(size=12, vjust=1), 
        axis.line.x = element_line(linewidth = 1), 
        axis.line.y  = element_line(linewidth = 1),
        axis.text.y = element_text(face = 'bold'))

dev.off()
#variance test
bartlett.test(CCR5MFI ~samples, data = df %>% filter(tissue == "L"))$p.value
bartlett.test(CCR5MFI ~samples, data = df %>% filter(tissue == "V"))$p.value
bartlett.test(CCR5MFI ~samples, data = df %>% filter(tissue == "LLB" & layer != "LA"))$p.value

# Stats
t.test(log(CCR5MFI) ~samples, data = df %>% filter(tissue == "L"))$p.value
t.test(log(CCR5MFI) ~samples, data = df %>% filter(tissue == "V"))$p.value
t.test(vum$CCR5MFI, ve$CCR5MFI, paired= TRUE)$p.value #sig
t.test(log(CCR5MFI) ~samples, data = df %>% filter(tissue == "LLB" & layer != "LA"))$p.value
summary(aov(log(CCR5MFI)~samples, data = df %>% filter(layer== "E" & tissue!= "ULB" & tissue != "A")))
summary(aov(log(CCR5MFI)~samples, data = df %>% filter(layer== "UM" & tissue!= "ULB" & tissue != "A")))



#CCR5CD69 CD4 T cell subsets #######
pdf("CCR5CD69subsets.pdf", width= 6, height = 4)
df %>% 
  filter(tissue != "ULB")%>%
  group_by(samples) %>%summarise(CCR5SP = mean(CCR5SP),CD69SP = mean(CD69SP),DN = mean(DN),DP = mean(DP)) %>%
  pivot_longer(cols = c("DP","CD69SP","CCR5SP", "DN" ), names_to = "Subset")%>%
  mutate(Subset = factor(Subset, levels = c("DP", "CD69SP", "CCR5SP", "DN"))) %>%  # Refactor the 'Subset' column
  ggplot(aes(samples, value, fill=Subset))+geom_col(color='black', width=0.8, size = 0.8)+
  scale_fill_manual(values=c('red', "orange", "blue", 'grey'))+geom_hline(yintercept = c(25,50,75), size=0.1, alpha = 0.5)+
  ggtitle("CCR5/CD69")+
  theme_classic()+
  RotatedAxis()+xlab("")+ylab("Percent of CD4+ T cells")+
  theme(axis.text.x = element_text(size=12, vjust=1),
        axis.line.x = element_line(linewidth = 0.8),
        axis.line.y = element_line(linewidth = 0.8))
dev.off()

#### calculations########
#mean pDC percentage
df%>%
  filter(tissue=="L") %>% 
  summarise(pdcL = mean(pDCper), na.rm = T, pdcLSD = sd(pDCper)) #calculate means and sd
df%>%
  filter(tissue=="V") %>% 
  summarise(pdcV = mean(pDCper), na.rm = T, pdcVSD = sd(pDCper))
df%>%
  filter(tissue=="ULB") %>% 
  summarise(pdcULB = mean(pDCper), na.rm = T, pdcUSD = sd(pDCper))
df%>%
  filter(tissue=="LLB") %>% 
  summarise(pdcLLB = mean(pDCper), na.rm = T, pdcLBSD = sd(pDCper))
df%>%
  filter(tissue=="A") %>% 
  summarise(pdcA = mean(pDCper), na.rm = T, pdcASD = sd(pDCper))
df%>%
  filter(tissue=="B") %>% 
  summarise(pdcB = mean(pDCper), na.rm = T, pdcBSD = sd(pDCper))

#mean CD69 percentage
df%>%
  filter(tissue=="B") %>% 
  summarise(cd69perB = mean(CD69per), na.rm = T)
df%>%
  filter(samples=="LLB_UM") %>% 
  summarise(cd69perLLBum = mean(CD69per), na.rm = T)
df%>%
  filter(samples=="LLB_E") %>% 
  summarise(cd69perLLBe = mean(CD69per), na.rm = T)

#mean CCR5 percentage
df%>%
  filter(samples=="LLB_LA") %>% 
  summarise(ccr5perLA = mean(CCR5per), na.rm = T)
df%>%
  filter(tissue=="B") %>% 
  summarise(ccr5perB = mean(CCR5per), na.rm = T)