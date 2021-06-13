
R2_control_Aiotic  <- read.csv("Aiotic GWAS Experiment control R^2.csv")
R2_stress_Aiotic <- read.csv("Aiotic GWAS Experiment stress R^2.csv")



All_H2_Aiotic  <- (0.532814077-0.125940274)/0.532814077#The total heritability was calculated according to the analysis of variance
sig_H2_Aiotic_control <- sum(R2_control_Aiotic[-1,15] ) #The heritability of significant QTLs in control environment was calculated according to the linear mixed model of tassel software
sig_H2_Aiotic_stress <- sum(R2_stress_Aiotic[-1,15] ) #The heritability of significant QTLs in stress environment was calculated according to the linear mixed model of tassel software





R2_control_Biotic  <- read.csv("Biotic GWAS Experiment monoculture R^2.csv")
R2_stress_Biotic <- read.csv("Biotic GWAS Experiment co-culture R^2.csv")



All_H2_Biotic  <- (0.186124929-0.074096221)/0.186124929
sig_H2_Biotic_control <- sum(R2_control_Biotic[-1,15] ) 
sig_H2_Biotic_stress <- sum(R2_stress_Biotic[-1,15] ) 
