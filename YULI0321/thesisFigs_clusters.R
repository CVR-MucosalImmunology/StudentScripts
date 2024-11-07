setwd("C:/Users/annab/The University of Sydney (Staff)/Thomas O'Neil - TRM")
setwd("~/Library/CloudStorage/OneDrive-TheUniversityofSydney(Staff)/projects/TRM/")

library(readxl)
library(ggplot2)
library(tidyr)
library(dplyr)
library(Seurat)
suppressPackageStartupMessages(library(ComplexHeatmap))
library(ComplexHeatmap)
list.files()

data <- readxl::read_excel("thesis resilts flow.xlsx") 
data$tissue[data$tissue == "AC"] = "ULB"
data$tissue[data$tissue == "Sig"] = "LLB"
data$tissue[data$tissue == "R"] = "LLB"

data <- data[1:47,]

{df <- data.frame(samples = paste0(data$tissue, "_", data$layer),
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
  
  df$C1 <- 100*data$C1/data$`CD4 #`
  df$C2 <- 100*data$C2/data$`CD4 #`
  df$C3 <- 100*data$C3/data$`CD4 #`
  df$C4 <- 100*data$C4/data$`CD4 #`
  df$C5 <- 100*data$C5/data$`CD4 #`
  df$C6 <- 100*data$C6/data$`CD4 #`
  df$C7 <- 100*data$C7/data$`CD4 #`
  df$C8 <- 100*data$C8/data$`CD4 #`
  
  df$C1_CCR5 <- data$`C1_CCR5 MFI`
  df$C1_DR <- data$C1_DR
  df$C1_CD69 <- data$C1_CD69
  
  df$C2_CCR5 <- data$`C2_CCR5 MFI`
  df$C2_DR <- data$C2_DR
  df$C2_CD69 <- data$C2_CD69
  
  df$C3_CCR5 <- data$`C3_CCR5 MFI`
  df$C3_DR <- data$C3_DR
  df$C3_CD69 <- data$C3_CD69
  
  df$C4_CCR5 <- data$`C4_CCR5 MFI`
  df$C4_DR <- data$C4_DR
  df$C4_CD69 <- data$C4_CD69
  
  df$C5_CCR5 <- data$`C5_CCR5 MFI`
  df$C5_DR <- data$C5_DR
  df$C5_CD69 <- data$C5_CD69
  
  df$C6_CCR5 <-data$`C6_CCR5 MFI`
  df$C6_DR <- data$C6_DR
  df$C6_CD69 <- data$C6_CD69
  
  df$C7_CCR5 <- data$`C7_CCR5 MFI`
  df$C7_DR <- data$C7_DR
  df$C7_CD69 <- data$C7_CD69
  
  df$C8_CCR5 <- data$`C8_CCR5 MFI`
  df$C8_DR <- data$C8_DR
  df$C8_CD69 <- data$C8_CD69
  
  df$C1_PD1 <- data$C1_PD1
  df$C2_PD1 <- data$C2_PD1
  df$C3_PD1 <- data$C3_PD1
  df$C4_PD1 <- data$C4_PD1
  df$C5_PD1 <- data$C5_PD1
  df$C6_PD1 <- data$C6_PD1
  df$C7_PD1 <- data$C7_PD1
  df$C8_PD1 <- data$C8_PD1
  
  df$C1_CD28 <- data$C1_CD28
  df$C2_CD28 <- data$C2_CD28
  df$C3_CD28 <- data$C3_CD28
  df$C4_CD28 <- data$C4_CD28
  df$C5_CD28 <- data$C5_CD28
  df$C6_CD28 <- data$C6_CD28
  df$C7_CD28 <- data$C7_CD28
  df$C8_CD28 <- data$C8_CD28
  
  
  df$donor <- paste0(data$exp,"_",data$donor)
  
  df <- na.omit(df)
}
colnames(df)

df$samples <- factor(df$samples, levels=c("L_UM", "L_E", 
                                          "V_UM", "V_E",
                                          "LLB_UM", "LLB_LA"
))

rownames(df) <- 1:nrow(df)
df0 <- df
df <- df0[-27,]


### Cluster % across tissues ######
df %>% mutate(Remaining = 100 - rowSums(df[,21:28]))%>% pivot_longer(cols = c("C1", "C2", "C3", "C4", "C5", "C6", "C7", "C8"), names_to = "cluster") %>%
  group_by(samples, cluster) %>%
  summarise(mean=mean(value)) %>%
  ggplot(aes(samples, mean, fill=cluster)) + geom_col(width = 0.8)+theme_classic()+geom_col(color='black', linewidth = 0.6)+
  theme(axis.line = element_line(linewidth = 0.6),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 0.7))+
  geom_hline(yintercept = 50, alpha = 0.3)+
  ylab("% of CD4+ T cells") +
  xlab("")

##stats of clusters% across tissues
# compare LE VE
t.test(df$C1[df$samples == "L_E"], df$C1[df$samples=="V_E"])
t.test(df$C2[df$samples == "L_E"], df$C2[df$samples=="V_E"])
t.test(df$C3[df$samples == "L_E"], df$C3[df$samples=="V_E"])
t.test(df$C4[df$samples == "L_E"], df$C4[df$samples=="V_E"])
t.test(df$C5[df$samples == "L_E"], df$C5[df$samples=="V_E"])
t.test(df$C6[df$samples == "L_E"], df$C6[df$samples=="V_E"])
t.test(df$C7[df$samples == "L_E"], df$C7[df$samples=="V_E"])
t.test(df$C8[df$samples == "L_E"], df$C8[df$samples=="V_E"])

# compare L, V, LLB UMs
summary(aov(C1~samples, data = df %>% filter(layer=='UM')))
summary(aov(C2~samples, data = df %>% filter(layer=='UM')))#0.00161 **
  TukeyHSD(aov(C2~samples, data = df %>% filter(layer=='UM'))) #VUM_LUM 0.0133 *, LLBUM_LUM 0.0022**
summary(aov(C3~samples, data = df %>% filter(layer=='UM')))
summary(aov(C4~samples, data = df %>% filter(layer=='UM')))
summary(aov(C5~samples, data = df %>% filter(layer=='UM')))
summary(aov(C6~samples, data = df %>% filter(layer=='UM')))#0.00534 **
   TukeyHSD(aov(C6~samples, data = df %>% filter(layer=='UM'))) #LLBUM_LUM 0.0065**, LLBUM_VUM 0.0127*
summary(aov(C7~samples, data = df %>% filter(layer=='UM')))
summary(aov(C8~samples, data = df %>% filter(layer=='UM')))

# compare epi and underlying mucosa from same tissue
# LE vs LUM
t.test(df$C1[df$samples == "L_E"], df$C1[df$samples=="L_UM"], paired = F)$p.value
t.test(df$C2[df$samples == "L_E"], df$C2[df$samples=="L_UM"], paired = F)$p.value
t.test(df$C3[df$samples == "L_E"], df$C3[df$samples=="L_UM"], paired =F)$p.value
#0.0002***
t.test(df$C4[df$samples == "L_E"], df$C4[df$samples=="L_UM"], paired = F)$p.value
#0.0334*
t.test(df$C5[df$samples == "L_E"], df$C5[df$samples=="L_UM"], paired = F)$p.value
#0.0199*
t.test(df$C6[df$samples == "L_E"], df$C6[df$samples=="L_UM"], paired = F)$p.value
#0.0282*
t.test(df$C7[df$samples == "L_E"], df$C7[df$samples=="L_UM"], paired = F)$p.value
#0.0073**
t.test(df$C8[df$samples == "L_E"], df$C8[df$samples=="L_UM"], paired = F)$p.value
#0.0002***

#paired test for VE vs VUM
t.test(df$C1[df$samples == "V_E"], df$C1[df$samples=="V_UM"], paired = T)$p.value
t.test(df$C2[df$samples == "V_E"], df$C2[df$samples=="V_UM"], paired = T)$p.value
t.test(df$C3[df$samples == "V_E"], df$C3[df$samples=="V_UM"], paired = T)$p.value
#0.005**
t.test(df$C4[df$samples == "V_E"], df$C4[df$samples=="V_UM"], paired = T)$p.value
t.test(df$C5[df$samples == "V_E"], df$C5[df$samples=="V_UM"], paired = T)$p.value
#0.039*
t.test(df$C6[df$samples == "V_E"], df$C6[df$samples=="V_UM"], paired = T)$p.value
#0.02*
t.test(df$C7[df$samples == "V_E"], df$C7[df$samples=="V_UM"], paired = T)$p.value
t.test(df$C8[df$samples == "V_E"], df$C8[df$samples=="V_UM"], paired = T)$p.value



# DR% on clusters across tissues heatmap 
df %>% pivot_longer(grep("_DR", colnames(df), value = TRUE), names_to = "DR_Clusters", values_to = "DR_Perct") %>%
  group_by(samples, DR_Clusters) %>%
  summarise(mean = mean(DR_Perct)) %>% 
  ungroup() %>%
  pivot_wider(names_from  = DR_Clusters, values_from = mean) ->df_dr #%>% #make a wider form with DR% across tissues per clustr for heatmap

rownames(df_dr) <- df_dr$samples
df_dr %>% select(-samples)%>%
  Heatmap(cluster_columns = F, cluster_rows = F, row_labels = rownames(df_dr),
          col = colorRamp2(c(0,50,100), c("blue","white","red")))

#heatmap: % CD69 on clusters across tissues
tmp3 <- df%>% pivot_longer(grep("_CD69", colnames(df), value = TRUE), names_to = "CD69_Clusters", values_to = "CD69_Perct")
tmp3 %>% 
  group_by(samples, CD69_Clusters) %>%
summarise(mean = mean(CD69_Perct)) %>% 
ungroup() %>%
  pivot_wider(names_from  = CD69_Clusters, values_from = mean) ->df2 #%>%

rownames(df2) <- df2$samples
df2 %>% select(-samples)%>%
  Heatmap(cluster_columns = F, cluster_rows = F, row_labels = rownames(df2))



#### Characteristics of clusters: complete #######
# DR on 8 clusters
df %>%
  pivot_longer(cols = grep("_DR", colnames(df), value = TRUE), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
  ggplot(aes(CCR5MFI_Clusters, CCR5MFI_Values, fill = CCR5MFI_Clusters)) +
  geom_boxplot(outliers = FALSE, linewidth = 0.6) +  # Boxplot with full data
  stat_summary(fun = mean, geom = "point", aes(group = samples, shape=samples), size = 2, position = position_jitter(width = 0.2)) + 
  # summarise the DR in each cluster per tissue, then add samples as points on top of it
  theme_classic() + 
  theme(axis.line = element_line(linewidth = 0.6),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 0.7)) +
  RotatedAxis() +
  labs(x = "", y = "DR %")+
  geom_hline(yintercept = seq(0,80,10), alpha=0.1)

### CD69% on 8 clusters
tmp3 <- df%>% pivot_longer(grep("_CD69", colnames(df), value = TRUE), names_to = "CD69_Clusters", values_to = "CD69_Perct")

tmp3 %>%
  ggplot(aes(CD69_Clusters, CD69_Perct, fill=CD69_Clusters))+
  geom_boxplot(outliers=F, size = 0.6)+
  stat_summary(geom = "point", fun = "mean", aes(group = samples, shape = samples), position = position_jitter(width = 0.2)) +
  #facet_wrap(~samples)+
  theme_classic()+
  RotatedAxis()+labs(x="", y="CD69 %") + 
  theme(axis.line = element_line(linewidth = 0.6),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 0.7))+
  geom_hline(yintercept = seq(0,80,10), alpha=0.1)

### CCR7% on clusters
tmp3 <- df%>% pivot_longer(grep("_CCR7", colnames(df), value = TRUE), names_to = "CD69_Clusters", values_to = "CD69_Perct")

tmp3 %>%
  ggplot(aes(CD69_Clusters, CD69_Perct, fill=CD69_Clusters))+
  geom_boxplot(outliers=F, size = 0.6)+
  stat_summary(geom = "point", fun = "mean", aes(group = samples, shape = samples), position = position_jitter(width = 0.2)) +
  #facet_wrap(~samples)+
  theme_classic()+
  RotatedAxis()+labs(x="", y="CCR7 %") + 
  theme(axis.line = element_line(linewidth = 0.6),
        axis.text.y = element_text(face = "bold"),
        axis.text.x = element_text(vjust = 0.7))+
  geom_hline(yintercept = seq(0,80,10), alpha=0.1)

# PD1 & CD28 gMFI on 8 clusters 
df %>% pivot_longer(cols = grep("_PD1", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
  ggplot(aes(CCR5MFI_Clusters,log(CCR5MFI_Values), fill=CCR5MFI_Clusters))+
  geom_boxplot(outliers=F, size = 0.6, show.legend = F)+ 
  #facet_wrap(~samples)+
  theme_classic()+RotatedAxis()+labs(x="", y="PD1 gMFI")+
  stat_summary(fun=mean, geom = "point", aes(group=samples, shape=samples),  size = 2, position = position_jitter(width = 0.2), show.legend = F)+
  theme_classic() + 
  theme(axis.line = element_line(linewidth = 0.6),
        axis.text.y = element_text(face = "bold", hjust = 0.8),
        axis.text.x = element_text(vjust = 0.7)) +
  RotatedAxis()
df %>% pivot_longer(cols = grep("_CD28", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
  ggplot(aes(CCR5MFI_Clusters,log(CCR5MFI_Values), fill=CCR5MFI_Clusters))+
  geom_boxplot(outliers=F, size = 0.6, show.legend = F)+ 
  #facet_wrap(~samples)+
  theme_classic()+RotatedAxis()+labs(x="", y="CD28 gMFI")+
  stat_summary(show.legend = F, fun=mean, geom = "point", aes(group=samples, shape=samples),  size = 2, position = position_jitter(width = 0.2))+
  theme_classic() + 
  theme(axis.line = element_line(linewidth = 0.6),
        axis.text.y = element_text(face = "bold", hjust = 0.8),
        axis.text.x = element_text(vjust = 0.7)) +
  RotatedAxis()

#Streamline the boxplots
cowplot::plot_grid(ncol = 3,
                   df %>%
                     pivot_longer(cols = grep("_DR", colnames(df), value = TRUE), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
                     ggplot(aes(CCR5MFI_Clusters, CCR5MFI_Values, fill = CCR5MFI_Clusters)) +
                     geom_boxplot(outliers = FALSE, linewidth = 0.6,show.legend = F) +  # Boxplot with full data
                     stat_summary(show.legend = F, fun = mean, geom = "point", aes(group = samples, shape=samples), size = 2, position = position_jitter(width = 0.2)) + 
                     # summarise the DR in each cluster per tissue, then add samples as points on top of it
                     theme_classic() + 
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold"),
                           axis.text.x = element_text(vjust = 0.7)) +
                     RotatedAxis() +
                     labs(x = "", y = "HLADR %")+
                     geom_hline(yintercept = seq(0,80,10), alpha=0.1),
                   
                   tmp3 %>%
                     ggplot(aes(CD69_Clusters, CD69_Perct, fill=CD69_Clusters))+
                     geom_boxplot(outliers=F, size = 0.6, show.legend = F)+
                     stat_summary(show.legend = F, geom = "point", fun = "mean", aes(group = samples, shape = samples), position = position_jitter(width = 0.2)) +
                     #facet_wrap(~samples)+
                     theme_classic()+
                     RotatedAxis()+labs(x="", y="CD69 %") + 
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold"),
                           axis.text.x = element_text(vjust = 0.7))+
                     geom_hline(yintercept = seq(0,80,10), alpha=0.1),
                   
                   tmp3 %>%
                     ggplot(aes(CD69_Clusters, CD69_Perct, fill=CD69_Clusters))+
                     geom_boxplot(outliers=F, size = 0.6, show.legend = F)+
                     stat_summary(show.legend = F, geom = "point", fun = "mean", aes(group = samples, shape = samples), position = position_jitter(width = 0.2)) +
                     #facet_wrap(~samples)+
                     theme_classic()+
                     RotatedAxis()+labs(x="", y="CCR7 %") + 
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold"),
                           axis.text.x = element_text(vjust = 0.7))+
                     geom_hline(yintercept = seq(0,80,10), alpha=0.1),
                   
                   df %>% pivot_longer(cols = grep("_PD1", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
                     ggplot(aes(CCR5MFI_Clusters,log(CCR5MFI_Values), fill=CCR5MFI_Clusters))+
                     geom_boxplot(outliers=F, size = 0.6, show.legend = F)+ 
                     #facet_wrap(~samples)+
                     theme_classic()+RotatedAxis()+labs(x="", y="Log(PD1 gMFI)")+
                     stat_summary(fun=mean, geom = "point", aes(group=samples, shape=samples),  size = 2, position = position_jitter(width = 0.2), show.legend = F)+
                     theme_classic() + 
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold", hjust = 0.8),
                           axis.text.x = element_text(vjust = 0.7)) +
                     RotatedAxis(),
                   
                   df %>% pivot_longer(cols = grep("_CD28", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
                     ggplot(aes(CCR5MFI_Clusters,log(CCR5MFI_Values), fill=CCR5MFI_Clusters))+
                     geom_boxplot(outliers=F, size = 0.6, show.legend = F)+ 
                     #facet_wrap(~samples)+
                     theme_classic()+RotatedAxis()+labs(x="", y="Log(CD28 gMFI)")+
                     stat_summary(show.legend = F, fun=mean, geom = "point", aes(group=samples, shape=samples),  size = 2, position = position_jitter(width = 0.2))+
                     theme_classic() + 
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold", hjust = 0.8),
                           axis.text.x = element_text(vjust = 0.7)) +
                     RotatedAxis(),
                   align = c("hv"))



### Clusters with HLADR% on CD4 T cells, and with pDC% Correlations ########

cor(df$C1,df$HLADRperCD4)# R 0.48
summary(lm(df$C1~df$HLADRperCD4)) #R-squared:  0.2019  p-value: 0.008404 **
cor(df$C2,df$HLADRperCD4) # 0.24
summary(lm(df$C2~df$HLADRperCD4)) #R-squared:  0.02042 p 0.219
cor(df$C3,df$HLADRperCD4) #-0.20
summary(lm(df$C3~df$HLADRperCD4)) #R-squared:  0.004137  p 0.3001
cor(df$C4,df$HLADRperCD4)#-0.19
summary(lm(df$C4~df$HLADRperCD4)) #R-squared:  -0.001023, p-value: 0.3331
cor(df$C5,df$HLADRperCD4) #0.19
summary(lm(df$C5~df$HLADRperCD4)) #R^2 -0.0009131, p-value: 0.3323
cor(df$C6,df$HLADRperCD4) #-0.19
summary(lm(df$C6~df$HLADRperCD4)) #Adjusted R-squared:  0.001776 p-value: 0.3146
cor(df$C7,df$HLADRperCD4) #0.015
summary(lm(df$C7~df$HLADRperCD4)) #Adjusted R-squared:  -0.03681,p-value: 0.9398
cor(df$C8,df$HLADRperCD4) #-0.40
summary(lm(df$C8~df$HLADRperCD4)) #Adjusted R-squared:  0.1255, p-value: 0.03353 *

#Correlations with DR% CD4+ T cells
y='HLADRperCD4';cowplot::plot_grid(ncol=4, 
                                     
                   df2 %>%  
                     ggplot(aes_string('C1', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#f8766d")+ylim(0,40)+
                     ylab('HLA-DR %')+theme(panel.border = element_rect(linewidth = 1),
                                            axis.text = element_text(face = 'bold')),
                   df2 %>%  
                     ggplot(aes_string('C2', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#d39200")+ylim(0,40)+
                     ylab("")+theme(panel.border = element_rect(linewidth = 1),
                                    axis.text = element_text(face = 'bold')),
                   df2 %>%  
                     ggplot(aes_string('C3', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#7cae00")+ylim(0,40)+
                     ylab("")+theme(panel.border = element_rect(linewidth = 1),
                                    axis.text = element_text(face = 'bold')),
                   df2 %>%  
                     ggplot(aes_string('C4', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#00be67")+ylim(0,40)+
                     ylab("")+theme(panel.border = element_rect(linewidth = 1),
                                    axis.text = element_text(face = 'bold')),
                   df2 %>%  
                     ggplot(aes_string('C5', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#00bfc4")+ylim(0,40)+
                     ylab('HLA-DR %')+theme(panel.border = element_rect(linewidth = 1),
                                            axis.text = element_text(face = 'bold')),
                   df2 %>%  
                     ggplot(aes_string('C6', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#00a9ff")+ylim(0,40)+
                     ylab("")+theme(panel.border = element_rect(linewidth = 1),
                                    axis.text = element_text(face = 'bold')),
                   df2 %>%  
                     ggplot(aes_string('C7', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#c77cff")+ylim(0,40)+xlim(0,15)+
                     ylab("")+theme(panel.border = element_rect(linewidth = 1),
                                    axis.text = element_text(face = 'bold')),
                   df2 %>%  
                     ggplot(aes_string('C8', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#ff61cc")+ylim(0,40)+
                     ylab("")+theme(panel.border = element_rect(linewidth = 1),
                                    axis.text = element_text(face = 'bold'))
                   )

#correlations with pDC%
y='pDCper';cowplot::plot_grid(ncol=4, 
                   df2 %>%  
                     ggplot(aes_string('C1', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#f8766d")+ylim(0,1)+
                     ylab("% pDC"),
                   df2 %>%  
                     ggplot(aes_string('C2', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#d39200")+ylim(0,1)+
                     ylab(""),
                   df2 %>%  
                     ggplot(aes_string('C3', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#7cae00")+ylim(0,1)+
                     ylab(""),
                   df2 %>%  
                     ggplot(aes_string('C4', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#00be67")+ylim(0,1)+
                     ylab(""),
                   df2 %>%  
                     ggplot(aes_string('C5', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#00bfc4")+ylim(0,1)+
                     ylab("% pDC"),
                   df2 %>%  
                     ggplot(aes_string('C6', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#00b9e3")+ylim(0,1)+
                     ylab(""),
                   df2 %>%  
                     ggplot(aes_string('C7', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#c77cff")+ylim(0,1)+xlim(0,15)+
                     ylab(""),
                   df2 %>%  
                     ggplot(aes_string('C8', y))+
                     geom_point(size=1)+theme_bw()+
                     geom_smooth(method = 'lm', colour = "#ff61cc")+ylim(0,1)+
                     ylab("")
)


#Correlation of C1/C8 with activation #######
df2 %>%  
  ggplot(aes(log(C1/C8), HLADRperCD4))+
  geom_point(size=1)+theme_bw()+
  geom_smooth(method = 'lm', colour = "blue")+ylim(0,40)+
  geom_vline(xintercept = 0)+
  theme(panel.border = element_rect(linewidth = 1),
        axis.text = element_text(face = 'bold'))+
  xlab('Log (C1 %/ C8 %)')+
  ylab('HLA-DR %')

cor(log(df$C1/df$C8), df$HLADRperCD4)#0.52
summary(lm(log(df$C1/df$C8)~df$HLADRperCD4)) #p-value: 0.003751**

#Correlation of C1 and C8
cor(df$C1, df$C8)#-0.51
summary(lm(df$C1~df$C8)) #pval = 0.004

### CD69 DR CCR5MFI on cluster1 and 8 #########

# CCR5gMFI on clusters (we collected geometric mean, but annotated as MFI in the file) 
tmp <-df %>% pivot_longer(cols = grep("_CCR5", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
  group_by(samples, CCR5MFI_Clusters) %>%
  summarise(mean = mean(CCR5MFI_Values)) #make a template which takes the column names matched _CCR5 into a new column, 
#and take values to make another one. Calculated the means across tissue samples.

cowplot::plot_grid(ncol=4,
                   df2 %>% pivot_longer(cols = grep("_CD69", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
                     filter(CCR5MFI_Clusters%in% c("C1_CD69","C8_CD69")) %>%
                     ggplot(aes(CCR5MFI_Clusters,CCR5MFI_Values, fill=CCR5MFI_Clusters))+
                     geom_boxplot(outliers=F, size = 0.6, show.legend = F)+ 
                     theme_classic()+RotatedAxis()+labs(x="", y="CD69 %")+
                     stat_summary(fun=mean, geom = "point", aes(group=samples, shape=samples),  size = 2, position = position_jitter(width = 0.2), show.legend = F)+
                     theme_classic() + 
                     scale_fill_manual(values=c("#f8766d","#ff61cc"))+
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold", hjust = 0.8),
                           axis.text.x = element_text(vjust = 0.7)) +
                     RotatedAxis(),
                   df2 %>% pivot_longer(cols = grep("_DR", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
                     filter(CCR5MFI_Clusters %in%c("C1_DR","C8_DR")) %>%
                     ggplot(aes(CCR5MFI_Clusters,CCR5MFI_Values, fill=CCR5MFI_Clusters))+
                     geom_boxplot(outliers=F, size = 0.6, show.legend = F)+ 
                     theme_classic()+RotatedAxis()+labs(x="", y="HLADR %")+
                     stat_summary(show.legend = F, fun=mean, geom = "point", aes(group=samples, shape=samples),  size = 2, position = position_jitter(width = 0.2))+
                     theme_classic() + 
                     scale_fill_manual(values=c("#f8766d","#ff61cc"))+
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold", hjust = 0.8),
                           axis.text.x = element_text(vjust = 0.7)) +
                     RotatedAxis(),
                   df2 %>% pivot_longer(cols = grep("_CCR5", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
                     filter(CCR5MFI_Clusters %in%c("C1_CCR5","C8_CCR5")) %>%
                     ggplot(aes(CCR5MFI_Clusters,log(CCR5MFI_Values), fill=CCR5MFI_Clusters))+
                     geom_boxplot(outliers=F, size = 0.6, show.legend=F)+ 
                     theme_classic()+RotatedAxis()+labs(x="", y="log(CCR5 gMFI)")+
                     stat_summary(show.legend=F, fun=mean, geom = "point", aes(group=samples, shape=samples),  size = 2, position = position_jitter(width = 0.2))+
                     theme_classic() + 
                     scale_fill_manual(values=c("#f8766d","#ff61cc"))+
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold", hjust = 0.8),
                           axis.text.x = element_text(vjust = 0.7)) +
                     RotatedAxis(),    
                   df2 %>% pivot_longer(cols = grep("_PD1", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
                     filter(CCR5MFI_Clusters %in%c("C1_PD1","C8_PD1")) %>%
                     ggplot(aes(CCR5MFI_Clusters,log(CCR5MFI_Values), fill=CCR5MFI_Clusters))+
                     geom_boxplot(outliers=F, size = 0.6, show.legend=F)+ 
                     theme_classic()+RotatedAxis()+labs(x="", y="log(PD-1 gMFI)")+
                     stat_summary(show.legend=F, fun=mean, geom = "point", aes(group=samples, shape=samples),  size = 2, position = position_jitter(width = 0.2))+
                     theme_classic() + 
                     scale_fill_manual(values=c("#f8766d","#ff61cc"))+
                     theme(axis.line = element_line(linewidth = 0.6),
                           axis.text.y = element_text(face = "bold", hjust = 0.8),
                           axis.text.x = element_text(vjust = 0.7)) +
                     RotatedAxis(), 
                   
                   rel_widths = c(1,1,1,1),
                  align = 'h',axis = 'tb')



#stats between C1 and C8 difference ###########

df2 %>% pivot_longer(cols = grep("_CD69", colnames(df), value=T), names_to = "CCR5MFI_Clusters", values_to = "CCR5MFI_Values") %>%
  filter(CCR5MFI_Clusters%in% c("C1_CD69","C8_CD69")) -> df_c1_c8stat

t.test(df2$C1_CD69, df2$C8_CD69, paired = F)$p.value #0.004
t.test(df2$C1_DR, df2$C8_DR, paired = F)$p.value%>% 
  format(scientific = F) #p< 0.0001
t.test(df2$C1_CCR5, df2$C8_CCR5, paired = F)$p.value #0.72
t.test(df2$C1_PD1, df2$C8_PD1, paired = F)$p.value #0.054



### Looking at the correlations between clusters that are ubiqutous #######
#create names
clusters <- c("C1","C2","C5","C8")
#create an empty list for plots
plot_list<-list()
#generate Cn xCn comparisons
# Nested loop to go through all combinations of Cn and Cm
for (i in seq_along(clusters)) {
  for (j in seq_along(clusters)) {
    # Skip cases where Cn and Cm are the same
    if (i != j) {
      # Define column names for Cn and Cm
      Cn <- clusters[i]
      Cm <- clusters[j]
      # Calculate the log ratio in a temporary column
      df2$temp_ratio <- log(df2[[Cn]] / df2[[Cm]])
      # Generate the plot
      plot <- df2 %>%
        ggplot(aes(x = temp_ratio, y = HLADRperCD4)) +
        geom_point(size = 1) + 
        theme_bw() +
        geom_smooth(method = 'lm', colour = "blue") + 
        ylim(0, 40) +
        geom_vline(xintercept = 0) +
        theme(panel.border = element_rect(linewidth = 1),
              axis.text = element_text(face = 'bold')) +
        xlab(paste("Log (", Cn, " % / ", Cm, " %)", sep = "")) +
        ylab("HLA-DR %")
      
      # Add the plot to the list
      plot_list[[length(plot_list) + 1]] <- plot
    }
  }
}
#plotting
cowplot::plot_grid(plotlist = plot_list, ncol = 3)

#correlation of Cn/Cn with activation
cor(log(df$C1/df$C8), df$HLADRperCD4)#0.52
summary(lm(log(df$C1/df$C8)~df$HLADRperCD4)) #p-value: 0.003751**
cor(log(df$C1/df$C2), df$HLADRperCD4)#0.44 
summary(lm(log(df$C1/df$C2)~df$HLADRperCD4)) #p:0.02 
cor(log(df$C1/df$C5), df$HLADRperCD4)#0.16
summary(lm(log(df$C1/df$C5)~df$HLADRperCD4)) #p:0.39 

cor(log(df$C2/df$C1), df$HLADRperCD4)#-0.44
summary(lm(log(df$C2/df$C1)~df$HLADRperCD4))#p 0.02
cor(log(df$C2/df$C5), df$HLADRperCD4)#-0.06
summary(lm(log(df$C2/df$C5)~df$HLADRperCD4))#p 0.74
cor(log(df$C2/df$C8), df$HLADRperCD4)#0.39
summary(lm(log(df$C2/df$C8)~df$HLADRperCD4))# p 0.04

cor(log(df$C5/df$C1), df$HLADRperCD4)#-0.16
summary(lm(log(df$C5/df$C1)~df$HLADRperCD4))# p 0.39
cor(log(df$C5/df$C2), df$HLADRperCD4)#0.06
summary(lm(log(df$C5/df$C2)~df$HLADRperCD4))# p 0.74
cor(log(df$C5/df$C8), df$HLADRperCD4)#0.42
summary(lm(log(df$C5/df$C8)~df$HLADRperCD4)) #p 0.02

cor(log(df$C8/df$C1), df$HLADRperCD4)#-0.52
summary(lm(log(df$C8/df$C1)~df$HLADRperCD4))# p 0.004
cor(log(df$C8/df$C2), df$HLADRperCD4)#-0.39
summary(lm(log(df$C8/df$C2)~df$HLADRperCD4))# p 0.03
cor(log(df$C8/df$C5), df$HLADRperCD4)#-0.42
summary(lm(log(df$C8/df$C5)~df$HLADRperCD4)) #p 0.03




# new stats check the one extra unmatched LUM datapoint ------------------------------

#df = paired

df %>% filter(samples == "L_UM") %>% 
  dplyr::select(.,-c("tissue", 'layer', "experiment", "donor")) -> df_stat
donor2 <- df0[27,5:68]
v <- c()

list= list()
for(i in 1:length(donor2)){
  list[[i]]=ggplot(df_stat, aes_string(colnames(df_stat)[i])) +geom_hist()+geom_vline(xintercept = unlist(donor2[i]))
}

cowplot::plot_grid(plotlist=list[31:50], ncol=3)

# 2,6,7,8

donor2[14]




