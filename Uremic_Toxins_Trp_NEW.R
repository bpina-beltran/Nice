library(readxl)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(tidyverse)
library(FSA)
library(ggsignif)
library(reshape)
#set directory
setwd("C:/Users/pina-beltr.b/Desktop/Trp Models/Trp Models")

#upload plasma data
Plasma_UTs<- read_excel("Data/2023-03-10_Mice_plasma_BlancaBeltran_result.xlsx")

#STATS

# group_by(Plasma_UTs,Group) %>%
#   summarise(
#     count = n(),
#     median = median(PhenylGlucuronide, na.rm = TRUE),
#     IQR = IQR(PhenylGlucuronide, na.rm = TRUE))

Plasma_UTs$Group <- as.factor(Plasma_UTs$Group)
Plasma_UTs$Group <- factor(Plasma_UTs$Group,levels= c("Trp", "5/6_Nx_Control", "5/6_Nx_Trp", 
                                                   "Ade_Control", "Ade_Trp"))

NX <- Plasma_UTs %>%
  filter(Group == '5/6_Nx_Control' | Group == '5/6_Nx_Trp')

Ade <- Plasma_UTs %>%
  filter(Group == 'Ade_Control' | Group == 'Ade_Trp')

# ##Kruskal-Wallis (Trp comparisons) 
# THIS WAS DONE WHEN I THOUGHT THAT ALL THE GROUPS HAD THE SAME SEX, 
# BUT I HAVE BOTH FEMALES AND MALES SO I CAN ONLY USE IT TO COMPARE FEMALE GROUPS
# 
# kw_results <-  list()
# for(i in names(Trp_comparisons[,5:17])){  
#   kw_results[[i]] <- kruskal.test(formula(paste(i, "~ Group")), data = Trp_comparisons)
# }
# 
# 
# kw_results <- as.data.frame(do.call(cbind, kw_results)) %>%
#               filter(!row_number() %in% c(5)) %>%
#               t()%>% as.data.frame()
# 
# kw_p.values <- as.data.frame(unlist(kw_results$p.value))
# colnames(kw_p.values) <- "p.value"
# 
# kw_p.values%>%filter(p.value<0.05)
# 
# dunns <-  list()
# dunns <- lapply(Trp_comparisons[,c("PhenylGlucuronide","IndoxylSulphate", "Hippuric_Acid", "PhenylSulphate", "Kynurenine","Kynurenic_Acid", "Indole3AceticAcid", "TMAO")], function(x) dunnTest(x ~ Group, data =Trp_comparisons,  method ="bonferroni"))
# 
# for(i in dunns[[i]]){  
#   df_dunns[[i]] <- as.data.frame(dunns[[i]][["res"]]) %>%
#     add_column(add_column = "i")
# }
# 
# 
# df_dunns <- as.data.frame(dunns[["PhenylGlucuronide"]][["res"]]) %>% 
#             add_column(Uremic_toxin = "PhenylGlucuronide") 
#   
# df_dunns <- rbind( df_dunns, (as.data.frame(dunns[["IndoxylSulphate"]][["res"]]) %>% 
#                                 add_column(Uremic_toxin = "IndoxylSulphate")),
#             (as.data.frame(dunns[["Hippuric_Acid"]][["res"]]) %>% 
#             add_column(Uremic_toxin = "Hippuric_Acid")),   
#             (as.data.frame(dunns[["PhenylSulphate"]][["res"]]) %>% 
#             add_column(Uremic_toxin = "PhenylSulphate")),
#             (as.data.frame(dunns[["Kynurenine"]][["res"]]) %>% 
#             add_column(Uremic_toxin = "Kynurenine")), 
#             (as.data.frame(dunns[["Kynurenic_Acid"]][["res"]]) %>% 
#             add_column(Uremic_toxin = "Kynurenic_Acid")),
#             (as.data.frame(dunns[["Indole3AceticAcid"]][["res"]]) %>% 
#             add_column(Uremic_toxin = "Indole3AceticAcid")),
#             (as.data.frame(dunns[["TMAO"]][["res"]]) %>% 
#             add_column(Uremic_toxin = "TMAO")))
# 
# df_dunns_signif <- df_dunns%>%filter(P.adj<0.05)
#   

##Mann Whitney U test ("Trp", "5/6_Nx_Trp") --> Only Females
Test_Trp <- list() 
Test_Trp <- lapply(Trp_comparisons[,5:17],function(x)  wilcox.test(x ~ Group, data=Trp_comparisons, exact = FALSE)$p.value)
Test_Trp <- as.data.frame(t(do.call(cbind, Test_Trp)))
colnames(Test_Trp) <- "p.value"

Test_Trp%>%filter(p.value<0.05)

##Mann Whitney U test ("5/6_Nx_Control", "5/6_Nx_Trp") --> Only Females

Test_NX <- list() 
Test_NX <- lapply(NX[,5:17],function(x)  wilcox.test(x ~ Group, data=NX, exact = FALSE)$p.value)
Test_NX <- as.data.frame(t(do.call(cbind, Test_NX)))
colnames(Test_NX) <- "p.value"

Test_NX%>%filter(p.value<0.05)

##Mann Whitney U test ("Ade_Control", "Ade_Trp") --> Only Males
Test_Ade <- list() 
Test_Ade <- lapply(Ade[,5:17],function(x)  wilcox.test(x ~ Group, data=Ade, exact = FALSE)$p.value)
Test_Ade <- as.data.frame(t(do.call(cbind, Test_Ade)))
colnames(Test_Ade) <- "p.value"

Test_Ade%>%filter(p.value<0.05)

# write.csv(df, file="results")

# boxplot
# 
# plot_list <- list()
# 
# for (i in names(Plasma_UTs[, 5:17])){
#   plot <- ggplot(Plasma_UTs, aes(y=i, x=Group)) +
#     geom_boxplot(shape=15) +
#     geom_jitter(width = 0.2) +
#     scale_colour_manual(values = c("#A7B800", "#E55451", "#990012",  "#00AFFF", "#0000A5"))
#   plot_list [[i]] <- plot
# }

quick <- function(df, toxin){
  ggplot(df, aes(x = Group, y = toxin, fill = Group)) +
       geom_boxplot(show.legend = FALSE, outlier.shape = NA, na.rm = TRUE, alpha=0.1, width = 0.6) +
       scale_fill_manual(values = c("#A7B800", "#E55451", "#990012",  "#00AFFF", "#0000A5"))+ 
       geom_jitter(position=position_jitter(0.2), shape=21, size=3) +
       theme_classic()+
       theme(axis.title.x = element_blank(),axis.title.y = element_text(size=12,color = "black", face="bold") ) +
       geom_vline(xintercept = 3.5,  color = "black", linetype = "dashed") +
       theme(axis.text.x=element_text(size=10,color = "black", face="bold"), axis.text.y=element_text(size=10,color = "black")) 
}



quick(Plasma_UTs,Plasma_UTs$p_CresylGlucuronide)+
  labs(y="p_CresylGlucuronide (µM)") +
  ylim(0, 1.2) +
  geom_label(geom="text", x=2, y=1.2, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=1.2, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  geom_signif(comparisons = list(c("Ade_Control","Ade_Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 0.3, color = "#0000E9")
ggsave(filename = "p_CresylGlucuronide.tiff")

quick(Plasma_UTs,Plasma_UTs$PhenylGlucuronide)+
  labs(y= "PhenylGlucuronide (µM)")+
  ylim(0,2) + geom_label(geom="text", x=2, y=2, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=2, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  geom_signif(comparisons = list(c("5/6_Nx_Control","5/6_Nx_Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 1.4, color = "darkred")+
  geom_signif(comparisons = list(c("5/6_Nx_Trp","Trp")), size = 1, textsize = 10, annotation = c("*"), tip_length = 0,  y_position = 1.6, color = "black")
ggsave(filename = "PhenylGlucuronide.tiff")

quick(Plasma_UTs,Plasma_UTs$IndoxylSulphate)+
  labs(y= "IndoxylSulphate (µM)", x = "Group") + 
  geom_label(geom="text", x=2, y=112, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=112, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  ylim(0,112) +
  geom_signif(comparisons = list(c("5/6_Nx_Trp","Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 95, color = "black")
ggsave(filename = "IndoxylSulphate.tiff")

quick(Plasma_UTs,Plasma_UTs$p_CresylSulphate)+
  labs(y= "p_CresylSulphate (µM)", x = "Group") + 
  ylim(0,10) +
  geom_label(geom="text", x=2, y=9, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=9, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) 
ggsave(filename = "p_CresylSulphate.tiff")

quick(Plasma_UTs,Plasma_UTs$Hippuric_Acid)+
  labs(y= "Hippuric_Acid (µM)", x = "Group")+
  ylim(0,25) + 
  geom_label(geom="text", x=2, y=25, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=25, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  geom_signif(comparisons = list(c("Ade_Control","Ade_Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 7, color = "#0000E9") +
  geom_signif(comparisons = list(c("5/6_Nx_Trp","Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 20, color = "black") 
ggsave(filename = "Hippuric_Acid.tiff")

quick(Plasma_UTs,Plasma_UTs$PhenylSulphate)+
  labs(y= "PhenylSulphate (µM)", x = "Group")+
  ylim(0,35) +
  geom_label(geom="text", x=2, y=35, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=35, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  geom_signif(comparisons = list(c("5/6_Nx_Trp","Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 29, color = "black") 
ggsave(filename = "PhenylSulphate.tiff")

quick(Plasma_UTs,Plasma_UTs$Kynurenine)+
  labs(y= "Kynurenine (µM)", x = "Group")+
  ylim(0,10) + 
  geom_label(geom="text", x=2, y=10, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=10, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  geom_signif(comparisons = list(c("5/6_Nx_Control","5/6_Nx_Trp")), size = 1, textsize = 10, tip_length = 0,  y_position = 8,  annotation = c("*"), color = "darkred")
ggsave(filename = "Kynurenine.tiff")

quick(Plasma_UTs,Plasma_UTs$Tryptophan) +
  ylim(0,110) +
  geom_label(geom="text", x=2, y=105, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=105, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  labs(y= "Tryptophan (µM)", x = "Group")
ggsave(filename = "Tryptophan.tiff")

quick(Plasma_UTs,Plasma_UTs$Kynurenic_Acid)+
  labs(y= "Kynurenic_Acid (µM)", x = "Group")+
  ylim(0,8) +
  geom_label(geom="text", x=2, y=8, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=8, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  geom_signif(comparisons = list(c("Ade_Control","Ade_Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 0.7, color = "#0000E9") +
  geom_signif(comparisons = list(c("5/6_Nx_Control","5/6_Nx_Trp")), size = 1, textsize = 10,  annotation = c("**"),tip_length = 0,  y_position = 5, color = "darkred")+
  geom_signif(comparisons = list(c("5/6_Nx_Trp","Trp")),  annotation = c("*"), size = 1, textsize = 10, tip_length = 0,  y_position = 6, color = "black") 
ggsave(filename = "Kynurenic_Acid.tiff")

quick(Plasma_UTs,Plasma_UTs$Tyrosine)+
  ylim(0,150) + 
  geom_label(geom="text", x=2, y=150, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=150, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  labs(y= "Tyrosine (µM)", x = "Group")+
  geom_signif(comparisons = list(c("5/6_Nx_Trp","Trp")),  annotation = c("*"), size = 1, textsize = 10, tip_length = 0,  y_position = 115, color = "black") 
ggsave(filename = "Tyrosine).tiff")

quick(Plasma_UTs,Plasma_UTs$Indole3AceticAcid)+
  labs(y= "Indole3AceticAcid (µM)", x = "Group")+
  ylim(0,15) + 
  geom_label(geom="text", x=2, y=15, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=15, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  geom_signif(comparisons = list(c("5/6_Nx_Control","5/6_Nx_Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 11, color = "darkred")+
  geom_signif(comparisons = list(c("5/6_Nx_Trp","Trp")),  annotation = c("*"), size = 1, textsize = 10, tip_length = 0,  y_position = 13, color = "black") 
ggsave(filename = "Indole3AceticAcid.tiff")

quick(Plasma_UTs,Plasma_UTs$PhenylAlanine)+
  ylim(0,100) + 
  geom_label(geom="text", x=2, y=100, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=100, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  labs(y= "PhenylAlanine (µM)", x = "Group")
ggsave(filename = "PhenylAlanine.tiff")

quick(Plasma_UTs,Plasma_UTs$TMAO)+
  labs(y= "TMAO (µM)", x = "Group")+
  ylim(0,80) + 
  geom_label(geom="text", x=2, y=80, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=7, fontface=2)+
  geom_label(geom="text", x=4.7, y=80, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=7, fontface=2) +
  geom_signif(comparisons = list(c("5/6_Nx_Trp","Trp")), size = 1, textsize = 10,  annotation = c("*"), tip_length = 0,  y_position = 65, color = "black") 
ggsave(filename = "TMAO.tiff")


#Bleeding time
Bleeding_time <- read_excel("Bleeding_time.xlsx")
Bleeding_time$Group <- as.factor(Bleeding_time$Group)
Bleeding_time$Group <- factor(Bleeding_time$Group,levels= c("Trp", "5/6_Nx_Control", "5/6_Nx_Trp", 
                                                      "Ade_Control", "Ade_Trp"))
Ade_B <- Bleeding_time %>%
  filter(Group == 'Ade_Control' | Group == 'Ade_Trp')
wilcox.test(`Bleeding_time(s)`~ Group, data=Ade_B )

NX_B <- Bleeding_time %>%
  filter(Group == '5/6_Nx_Control' | Group == '5/6_Nx_Trp')
wilcox.test(`Bleeding_time(s)`~ Group, data=NX_B )

Trp_comparisons_B <- Bleeding_time %>%
  filter(Group == '5/6_Nx_Trp' | Group == 'Trp')
wilcox.test(`Bleeding_time(s)`~ Group, data = Trp_comparisons_B )

quick(Bleeding_time, Bleeding_time$`Bleeding_time(s)`) +
  ylim(0,300) + 
  geom_label(geom="text", x=2, y=2, show.legend = FALSE, label ="FEMALES", fill="white", color="Black", size=10, fontface=2)+
  geom_label(geom="text", x=4.7, y=2, show.legend = FALSE, label ="MALES", fill="white", color="Black", size=10, fontface=2) 
ggsave(filename = "Bleeding_time.tiff")


#MICE WEIGHT
Mice_Weight <- read_excel("Mice_Weight.xlsx")

Mice_Weight$Group <- as.factor(Mice_Weight$Group)
Mice_Weight$Group <- factor(Mice_Weight$Group,levels= c("Trp", "5/6_Nx_Control", "5/6_Nx_Trp", 
                                                      "Ade_Control", "Ade_Trp"))


Females_Weight <- Mice_Weight%>%filter(Sex == 'Female') 
Males_Weight <- Mice_Weight%>%filter(Sex == 'Male')%>%
  filter(!is.na(Weight)) 
Males_Weight$Weight<- as.numeric(Males_Weight$Weight)
Males_Weight$Weeks <- factor(Males_Weight$Weeks,levels= c("Day_0", "Day_5", "Day_12","Day_18",
                                                          "Day_25", "Day_30", "Day_45"))
Females_Weight$Weight<- as.numeric(Females_Weight$Weight)

ggplot(Females_Weight, aes(x=Weeks, y=Weight, fill = factor(`Group`))) +
  scale_fill_manual(values = c("#E55451", "#990012")) +
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, na.rm = TRUE, alpha=0.1, width = 0.6) +
  ylim(15,25) +
  geom_point(shape=21, size=3, position=position_jitterdodge(), aes(fill=`Group`))+
  theme_classic()+
  theme(axis.title.x = element_blank())
ggsave(filename = "Trp_Females_Weight.tiff")

ggplot(Males_Weight, aes(x=Weeks, y=Weight, fill = factor(`Group`))) +
  scale_fill_manual(values = c("#A7B800", "#00AFFF", "#0000A5"))+ 
  geom_boxplot(show.legend = FALSE, outlier.shape = NA, na.rm = TRUE, alpha=0.1,width=0.6) +
  ylim(18,34) +
  geom_point(shape=21, size=3, position=position_jitterdodge(), aes(fill=`Group`))+
  theme_classic()+
  theme(axis.title.x = element_blank())
ggsave(filename = "Trp_Males_Weight.tiff")
ggsave(filename = "Trp_Males_Weight.tiff")
