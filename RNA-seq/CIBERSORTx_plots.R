#Code is used for plotting statistical analysis of CIBERSORTx results for COVID-19 2020 project
#Code developed by Alexander Capraro, August 2020

library(tidyverse)
library(dplyr)
library(ggplot2)
library(ggpubr)
library(rstatix)

#### Stacked bar average ####
celltype_average <- read.table(file = "/Users/z3416833/Documents/miCF/SARS-CoV2_Study/CIBERSORTx/GSE102580/CIBERSORTx_output_average.txt", sep = "\t", header = T)
celltype_average <- gather(celltype_average, celltype, percent, Intermediate.Secretory...Ciliated:FOXN4., factor_key=TRUE)

#Mock All
celltype_average_mock <- celltype_average[celltype_average$Condition == "Mock", ]
#celltype_average_mockN <- celltype_average_mock[celltype_average_mock$Tissue == "Nasal", ]
celltype_average_mock$Sample <- factor(celltype_average_mock$Sample, levels = c("Mock_OA","Mock_YA","Mock_CH","Mock_B"))

ggplot(celltype_average_mock, aes(x = Sample, y = percent, fill = celltype)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() #6x7.5

#Virus All
celltype_average_virus <- celltype_average[celltype_average$Condition == "Virus", ]
#celltype_average_virusN <- celltype_average_virus[celltype_average_virus$Tissue == "Nasal", ]
celltype_average_virus$Sample <- factor(celltype_average_virus$Sample, levels = c("Virus_OA","Virus_YA","Virus_CH","Virus_B"))

ggplot(celltype_average_virus, aes(x = Sample, y = percent, fill = celltype)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() #export 6x7.5

#Mock & Virus
celltype_average$Sample <- factor(celltype_average_virus$Sample, levels = c("Mock_B","Mock_CH","Mock_YA","Mock_OA","Virus_B","Virus_CH","Virus_YA","Virus_OA"))

ggplot(celltype_average, aes(x = Sample, y = percent, fill = celltype)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  scale_fill_brewer(palette = "Paired") +
  theme_minimal() +
  facet_wrap(~Condition)#export 6x7.5


#Nasal Mock vs. Virus

celltype_average_N <- celltype_average[celltype_average$Tissue == "Nasal", ]
celltype_average_N$AgeGroup <- factor(celltype_average_N$AgeGroup, levels = c("CH","YA","OA"))

ggplot(celltype_average_N, aes(x = Condition, y = percent, fill = celltype)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~AgeGroup, ncol = 3, scales = "fixed")

# Bronchial Mock vs Virus
celltype_average_B <- celltype_average[celltype_average$Tissue == "Bronchial", ]
celltype_average_B$AgeGroup <- factor(celltype_average_B$AgeGroup, levels = c("CH","YA","OA"))

ggplot(celltype_average_B, aes(x = Condition, y = percent, fill = celltype)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() #export 7.5x6

# Bronchial vs Nasal
celltype_average_B <- celltype_average[celltype_average$AgeGroup == "CH", ]

ggplot(celltype_average_B, aes(x = Tissue, y = percent, fill = celltype)) + 
  geom_bar(position = "stack", stat = "identity", alpha = 0.9) +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() #export 7.5x6

#### Box Plots ####
celltype_wide <- read.table(file = "/Users/z3416833/Documents/miCF/SARS-CoV2_Study/CIBERSORTx/GSE102580/CIBERSORTx_output.txt", sep = "\t", header = T)
celltype <- gather(celltype_wide, celltype, percent, Intermediate.Secretory...Ciliated:FOXN4., factor_key=TRUE)
celltype$AgeGroup <- factor(celltype$AgeGroup , levels = c("OA","YA","CH","B"))
celltype_n <- celltype[celltype$Tissue == "Nasal", ]
celltype_b <- celltype[celltype$Tissue == "Bronchial", ]
#celltype_bvn <- celltype[celltype$AgeGroup == "CH",]

# Mock
celltype_mock <- celltype[celltype$Condition == "Mock", ]
stat_test <- celltype_mock %>%
  group_by(celltype) %>%
  wilcox_test(percent ~ AgeGroup, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_mock, aes(x = AgeGroup, y = percent, fill = AgeGroup)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.75), dotsize = 1) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D")) +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "fixed")

# Virus
celltype_virus <- celltype[celltype$Condition == "Virus", ]
stat_test <- celltype_virus %>%
  group_by(celltype) %>%
  wilcox_test(percent ~ AgeGroup, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_virus, aes(x = AgeGroup, y = percent, fill = AgeGroup)) + 
  geom_boxplot() +
  geom_dotplot(binaxis='y', stackdir='center', 
               position=position_dodge(0.75), dotsize = 1) +
  scale_fill_manual(values = c("#BEBADA","#FFFFB3","#8DD3C7","#F7941D")) +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "fixed")

# Mock Nasal
celltype_mockN <- celltype_n[celltype_n$Condition == "Mock", ]
my_comparisons <- list( c("CH", "OA"), c("CH", "YA"), c("YA", "OA") )

stat_test <- celltype_mockN %>%
  group_by(celltype) %>%
  t_test(percent ~ AgeGroup, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_mockN, aes(x = AgeGroup, y = percent, fill = AgeGroup)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free") +
  stat_compare_means(method = "t.test", comparisons = my_comparisons) #export pdf 6x9


#Virus Nasal
celltype_virusN <- celltype_n[celltype_n$Condition == "Virus", ]
my_comparisons <- list( c("CH", "OA"), c("CH", "YA"), c("YA", "OA") )

ggplot(celltype_virusN, aes(x = AgeGroup, y = percent, fill = AgeGroup)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free")+
  stat_compare_means(method = "t.test", comparisons = my_comparisons)

#Nasal mock vs. virus

#all

ggplot(celltype_n, aes(x = AgeGroup, y = percent, fill = Condition)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free")


#Child
celltype_childN <- celltype_n[celltype_n$AgeGroup == "CH", ]

stat_test <- celltype_childN %>%
  group_by(celltype) %>%
  wilcox_test(percent ~ Condition, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_childN, aes(x = Condition, y = percent, fill = Condition)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free")

#YA
celltype_youngadultN <- celltype_n[celltype_n$AgeGroup == "YA", ]

stat_test <- celltype_youngadultN %>%
  group_by(celltype) %>%
  wilcox_test(percent ~ Condition, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_youngadultN, aes(x = Condition, y = percent, fill = Condition)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free")

#OA
celltype_oldadultN <- celltype_n[celltype_n$AgeGroup == "OA", ]

stat_test <- celltype_oldadultN %>%
  group_by(celltype) %>%
  wilcox_test(percent ~ Condition, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_oldadultN, aes(x = Condition, y = percent, fill = Condition)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free")

#Bronchial mock vs virus

stat_test <- celltype_b %>%
  group_by(celltype) %>%
  wilcox_test(percent ~ Condition, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_b, aes(x = Condition, y = percent, fill = Condition)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free")

#Bronchial vs Nasal

celltype_bvn_mock <- celltype_bvn[celltype_bvn$Condition == "Mock", ]

stat_test <- celltype_bvn_mock %>%
  group_by(celltype) %>%
  wilcox_test(percent ~ Condition, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_bvn_mock, aes(x = Tissue, y = percent, fill = Tissue)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free")

stat_test <- celltype_bvn_virus %>%
  group_by(celltype) %>%
  wilcox_test(percent ~ Condition, p.adjust.method = "bonferroni", paired = F)
stat_test

ggplot(celltype_bvn_mock, aes(x = Tissue, y = percent, fill = Tissue)) + 
  geom_boxplot() +
  scale_fill_brewer(palette = "Set3") +
  theme_minimal() +
  facet_wrap(~celltype, ncol = 3, scales = "free")

