---
title: "Untitled"
author: "Zach"
date: "December 9, 2019"
output: html_document
---

```{r setup, include=FALSE}
##For Importing data into R 
total_proteins_arabidopsis <- read.csv("C:/Users/HP15/Desktop/BCB546X/total_proteins_arabidopsis.txt", sep="")
total_proteins_rice <- read.csv("C:/Users/HP15/Desktop/BCB546X/total_proteins_rice.txt", sep="")


View(total_proteins_arabidopsis) ##For viewing data
View(total_proteins_rice)

library(ggplot2) ## For proper data Visualisation in R

head(total_proteins_arabidopsis[, c(2:3)]) #To be sure that the data has the expected colums
head(total_proteins_rice[, c(2:3)])

#Prepare data for analysis by making tables into matrix
table(total_proteins_arabidopsis$`Chromosome`, total_proteins_arabidopsis$Protein_Type)
table(total_proteins_rice$`Chromosome`, total_proteins_rice$Protein_type)
tableARAB3 <- table(total_proteins_arabidopsis$`Chromosome`, total_proteins_arabidopsis$Protein_Type)
tableRICE3 <- table(total_proteins_rice$`Chromosome`, total_proteins_rice$Protein_type)
tableARAB3
tableRICE3
table_percentARAB <- tableARAB3/ rowSums(table3) * 100    ## Calculate percentage of Proten types per chromosome
table_percentRICE <- tableRICE3/ rowSums(tableRICE3) * 100
table_percentARAB
table_percentRICE

library(tidyverse) ##Import packages for data manipulation to prepare it for analysis
library(reshape2)
dat_4A <- melt(table_percentARAB, id.vars = c("Chromosome"))#Melt function to stack data into three individaul columns
dat_4R <- melt(table_percentRICE, id.vars = c("Chromosome"))
view(dat_4A)
view(dat_4R)

my_dataA <- as_tibble(dat_4A)  
my_dataR <- as_tibble(dat_4R)
view(my_dataR)

ARAB<- dat_4A %>% rename(Chromosome = Var1,Protein_Type = Var2, Percent = value)
RICE<- dat_4R %>% rename(Chromosome = Var1,Protein_type = Var2, Percent = value)
View(ARAB)
view(RICE)

ARAB$Chromosome <- as.character(ARAB$Chromosome)##To convert Chromosome column values to character
RICE$Chromosome <- as.character(RICE$Chromosome)

#Change values in chromosome colum to roman numbers as they appear in the original paper
Final <- ARAB$Chromosome[ARAB$Chromosome == "1"] <- "I"    
Final2 <- ARAB$Chromosome[ARAB$Chromosome == "2"] <- "II"
Final3 <- ARAB$Chromosome[ARAB$Chromosome == "3"] <- "III"
Final4 <- ARAB$Chromosome[ARAB$Chromosome == "4"] <- "IV"
Final5 <- ARAB$Chromosome[ARAB$Chromosome == "5"] <- "V"
FinalC <- ARAB$Chromosome[ARAB$Chromosome == "C"] <- "Chloroplast"
FinalM <- ARAB$Chromosome[ARAB$Chromosome == "M"] <- "Mitochondria"
FinalR1 <- RICE$Chromosome[RICE$Chromosome == "1"] <- "I"
FinalR2 <- RICE$Chromosome[RICE$Chromosome == "2"] <- "II"
FinalR3 <- RICE$Chromosome[RICE$Chromosome == "3"] <- "III"
FinalR4 <- RICE$Chromosome[RICE$Chromosome == "4"] <- "IV"
FinalR5 <- RICE$Chromosome[RICE$Chromosome == "5"] <- "V"
FinaLR6 <- RICE$Chromosome[RICE$Chromosome == "6"] <- "VI"
FinalR7 <- RICE$Chromosome[RICE$Chromosome == "7"] <- "VII"
FinalR8 <- RICE$Chromosome[RICE$Chromosome == "8"] <- "VIII"
FinalR9 <- RICE$Chromosome[RICE$Chromosome == "9"] <- "IX"
FinalR10 <- RICE$Chromosome[RICE$Chromosome == "10"] <- "X"
FinalR11 <- RICE$Chromosome[RICE$Chromosome == "11"] <- "XI"
FinalR12 <- RICE$Chromosome[RICE$Chromosome == "12"] <- "XII"

##Specify variables for plotting data using ggplot

p2A <- ggplot(data = ARAB, aes(x = Chromosome, y = Percent, group = Protein_Type, fill = Protein_Type))
p2R <- ggplot(data = RICE, aes(x = Chromosome, y = Percent, group = Protein_type, fill = Protein_type))

p2A <- p2A + geom_bar(stat = "identity", width = 0.5, position = "dodge") #Specify graph dimmensions 
p2R <- p2R + geom_bar(stat = "identity", width = 0.5, position = "dodge")

p2A <- p2A + theme_bw() # To specify colors
p2R <- p2R + theme_bw()
p2A <- p2A + theme(axis.text.x = element_text(angle = 90))
p2R <- p2R + theme(axis.text.x = element_text(angle = 90))

##Print graphs
p2A
p2R
```

## R Markdown

This is an R Markdown document. Markdown is a simple formatting syntax for authoring HTML, PDF, and MS Word documents. For more details on using R Markdown see <http://rmarkdown.rstudio.com>.

When you click the **Knit** button a document will be generated that includes both content as well as the output of any embedded R code chunks within the document. You can embed an R code chunk like this:

```{}



```

## Including Plots

You can also embed plots, for example:

```{r pressure, echo=FALSE}

```

Note that the `echo = FALSE` parameter was added to the code chunk to prevent printing of the R code that generated the plot.
