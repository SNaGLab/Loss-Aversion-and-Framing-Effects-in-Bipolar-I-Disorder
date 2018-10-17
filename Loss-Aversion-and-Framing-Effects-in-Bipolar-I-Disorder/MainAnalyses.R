

##libraries
library(ggplot2)
library(reshape2)
library(reshape)
library(lme4)
library(nlme)
library(plyr)
library(multcomp)
library(psych)
library(car)
library(lmSupport)
library(pbkrtest)
library(sm)
library(MASS)
library(ICC)
library(vioplot)
library(Hmisc)
library(tidyr)
library(lsr)
library(rms)
library(gdata)
library(readxl)
library(stringr)
library(taRifx)
library(interplot)
library(nnet)

# -----------------------------------------------------------------------------------------------------

## Data: upload the csv files!

cd <- "/Users/Kimmie/Dropbox/Papers/BipolarPatients/Data/OnGit/" ### Put your own current directory here.
s1 <- "Demographics_full.csv" #Name file, leave intact
s2 <- "Gigantor.csv" #Name file, leave intact
s3 <- "Pabst.csv" #Name file, leave intact
s4 <- "Nerd.csv" #Name file, leave intact
s5 <- "SymptomScores.csv" #Name file, leave intact

file1=paste(cd, s1, sep = "")
file2=paste(cd, s2, sep = "")
file3=paste(cd, s3, sep = "")
file4=paste(cd, s4, sep = "")
file5=paste(cd, s5, sep = "")

Participants=read.csv(file1,stringsAsFactors=FALSE)
colnames(Participants)= c("ID","Group","Progress","Age","Caucasian","Gender","Education","income","employment","living","marital","children")

correctID=c('1_B','2_Y','3_B','7_Y','18_B','21_B','22_Y','33_Y','36_Y','49_Y','55_Y','61_Y','65_Y','67_Y','69_Y','70_Y','76_Y','80_Y','88_B','99_Y')
y=length(c(correctID))

for (i in 1:y) {
  Participants[i,"ID"]=correctID[i]
}

Gigantor=read.csv(file2)
Pabst=read.csv(file3)
Nerd=read.csv(file4)

DepressManic=read.csv(file5,stringsAsFactors=FALSE)
DepressManic=subset(DepressManic, DepressManic$Progress==100)

# -----------------------------------------------------------------------------------------------------

## In order to merge the datafiles, some manipulations are needed as ID's do not have same format or variables names differ.

# Education variable
Participants$Education=gsub("\\mba|\\+", "", Participants$Education)
destring(Participants$Education)

# Fix some ID's: should all be in XX_Y or XX_B format

Gigantor$ID=gsub("\\-", "_", Gigantor$ID)
Gigantor$ID=gsub("\\:|\\#", "", Gigantor$ID)
Pabst$ID=gsub("\\ ", "", Pabst$ID) #remove an additional space
Nerd$ID=gsub("\\_y", "_Y", Nerd$ID)
Nerd$ID=gsub("\\y", "_Y", Nerd$ID)
Nerd$ID=gsub("\\ ", "", Nerd$ID) #remove an additional space
Gigantor$ID=gsub("\\ ", "", Gigantor$ID) #remove an additional space

#remove the additional 00 for lower number ID's to be able to combine with Bipolar file
correctID=c('1_B','2_Y','3_B','7_Y','18_B','21_B','22_Y','33_Y','36_Y','49_Y','55_Y','61_Y','65_Y','67_Y','69_Y','70_Y','76_Y','80_Y','88_B','99_Y')
y=length(c(correctID))

for (i in 1:y) {
  DepressManic[i+1,'ID']=correctID[i]
}

Gigantor[5,'ID']='61_Y' #single fix
Pabst[5,'ID']='136_Y' #single fix
Pabst[22,'ID']='141_Y' #single fix: 141_D --> 141_Y

Nerd <- subset(Nerd, Progress==100)

# Fix colnames Gigantor/Nerd to match the others
colnames(Gigantor)=gsub("\\_Q3_", "_", colnames(Gigantor))
colnames(Nerd)[12]="Q5"

# Merge three data frames
Bipolar=rbind(Gigantor, Nerd, Pabst)
Bipolar=merge(by = "ID", Bipolar,Participants) 
Bipolar=merge(by="ID",Bipolar, DepressManic)

Bipolar=subset(Bipolar, Bipolar$Group!="NA") 

# -----------------------------------------------------------------------------------------------------

## Prepare some variables: create factor variables for grouping categories

Bipolar$Group=ifelse(Bipolar$Group==1,"BD","CTL")
Bipolar$Group=factor(Bipolar$Group)

Bipolar$Gender=factor(Bipolar$Gender)
Bipolar$Caucasian=factor(Bipolar$Caucasian)
Bipolar$employment=factor(Bipolar$employment)
Bipolar$marital=factor(Bipolar$marital)

# -----------------------------------------------------------------------------------------------------

### Preliminary analyses

## Demographics and clinical characteristics (Table 1)

#Age --> quite comparable, 36.1 vs 34.9 in BD vs CTL
aggregate(Bipolar$Age, by=list(group=Bipolar$Group),mean) 
aggregate(Bipolar$Age, by=list(group=Bipolar$Group),sd) 
t.test(Bipolar$Age ~ Bipolar$Group) #p=0.629

# Gender --> Gender proportion quite equal between groups (1=M, 2=F)
table(Bipolar$Group,Bipolar$Gender)
prop.test(table(Bipolar$Group,Bipolar$Gender),correct=FALSE) 

# Ethnicity --> White proportion quite equal between groups
table(Bipolar$Group,Bipolar$Caucasian)
prop.test(table(Bipolar$Group,Bipolar$Caucasian),correct=FALSE) 

#Education --> sig difference: CTL more education years
Bipolar$Education=destring(Bipolar$Education)
Bipolar$Education[is.na(Bipolar$Education)] <- 0
BipolarEduc=subset(Bipolar,Bipolar$Education>0)
aggregate(BipolarEduc$Education, by=list(group=BipolarEduc$Group),mean) 
aggregate(BipolarEduc$Education, by=list(group=BipolarEduc$Group),sd) 
t.test(BipolarEduc$Education ~ BipolarEduc$Group) #p=0.06

#Employment
Bipolar$employment
table(Bipolar$Group,Bipolar$employment)

test <- multinom(employment ~ Group, data = Bipolar)
summary(test)
z <- summary(test)$coefficients/summary(test)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

# Income (coding: 1= Less than $10,000,2=$10,000 - $25,000,3=$26,000 - $50,000,4=$51,000 - $75,000,5=$76,000 - $100,000,6=More than $100,000)
Bipolar$income
table(Bipolar$Group,Bipolar$income)
income=lm(income ~ Group, data=Bipolar)
summary(income)

# Nr children
Bipolar$children=gsub("\\P", "", Bipolar$children)
Bipolar$children=destring(Bipolar$children)
Bipolar$children[is.na(Bipolar$children)] <- 0
aggregate(Bipolar$children, by=list(group=Bipolar$Group),mean) 
aggregate(Bipolar$children, by=list(group=Bipolar$Group),sd) 
t.test(Bipolar$children ~ Bipolar$Group)

# Marital status
Bipolar$marital
table(Bipolar$Group,Bipolar$marital)
test <- multinom(marital ~ Group, data = Bipolar)
summary(test)
z <- summary(test)$coefficients/summary(test)$standard.errors
p <- (1 - pnorm(abs(z), 0, 1)) * 2

## Clinical symptoms in BD vs control group

#Create group samples
Bipolar.sample=subset(Bipolar, Group=="BD")
Control.sample=subset(Bipolar, Group=="CTL")

#Manic scores
aggregate(Bipolar$asrm_tot, by=list(group=Bipolar$Group),mean) #mania is higher for BD
aggregate(Bipolar$asrm_tot, by=list(group=Bipolar$Group),sd)
aggregate(Bipolar$asrm_tot, by=list(group=Bipolar$Gender),mean)
aggregate(Bipolar.sample$asrm_tot, by=list(group=Bipolar.sample$Gender),mean) #Females slightly more manic in BD
aggregate(Control.sample$asrm_tot, by=list(group=Control.sample$Gender),mean) #Females quite more manic in control

t.test(Bipolar$asrm_tot ~ Bipolar$Group) #p=0.12
t.test(Bipolar$asrm_tot ~ Bipolar$Gender) #p=0.08
t.test(Bipolar.sample$asrm_tot ~ Bipolar.sample$Gender) #not sig
t.test(Control.sample$asrm_tot ~ Control.sample$Gender) #sig difference, p=0.03

# Depression scores
aggregate(Bipolar$bdisf_tot, by=list(group=Bipolar$Group),mean) #Depression score is quite higher for Bipolar
aggregate(Bipolar$bdisf_tot, by=list(group=Bipolar$Group),sd)
aggregate(Bipolar$bdisf_tot, by=list(group=Bipolar$Gender),mean) 
aggregate(Bipolar.sample$bdisf_tot, by=list(group=Bipolar.sample$Gender),mean) #similar scores between males and females in Bipolar
aggregate(Control.sample$bdisf_tot, by=list(group=Control.sample$Gender),mean) #similar scores between males and females in control

t.test(Bipolar$bdisf_tot ~ Bipolar$Group) #sig difference: BD more depressed
wilcox.test(Bipolar$bdisf_tot ~ Bipolar$Group) #sig difference: BD more depressed
t.test(Bipolar$bdisf_tot ~ Bipolar$Gender)

# Anxiety
aggregate(Bipolar$stai_state_brief_tot, by=list(group=Bipolar$Group),mean) #BD more anxious than controls
aggregate(Bipolar$stai_state_brief_tot, by=list(group=Bipolar$Group),sd)
aggregate(Bipolar.sample$stai_state_brief_tot, by=list(group=Bipolar.sample$Gender),mean) #similar scores between males and females in Bipolar
aggregate(Control.sample$stai_state_brief_tot, by=list(group=Control.sample$Gender),mean) #males quite more anxious than females in control group

t.test(Bipolar$stai_state_brief_tot ~ Bipolar$Group) #p=0.11, not sig different.
t.test(Bipolar$stai_state_brief_tot ~ Bipolar$Gender) #p=0.3962
t.test(Bipolar.sample$stai_state_brief_tot ~ Bipolar.sample$Gender) #not sig
t.test(Control.sample$stai_state_brief_tot ~ Control.sample$Gender) #p=0.17

# Lastly, do symptoms correlate? --> Anxiety and depression sig pos correlate (r=0.59), anxiety and mania sig neg correlate (r=-0.23)
x=Bipolar[,c("bdisf_tot","asrm_tot","stai_state_brief_tot")]
rc=rcorr(as.matrix(x),type="pearson") 
rc

x=Bipolar.sample[,c("bdisf_tot","asrm_tot","stai_state_brief_tot")]
rc=rcorr(as.matrix(x),type="pearson") 
rc

x=Control.sample[,c("bdisf_tot","asrm_tot","stai_state_brief_tot")]
rc=rcorr(as.matrix(x),type="pearson") 
rc

# -----------------------------------------------------------------------------------------------------

# Loss aversion task --> calculate the reject gambles so a higher score indicates more loss aversion

# Some 'reject' gambles have been coded as 2 by Qualtrics, change to 0 for consistency.

Bipolar$Q3_1=ifelse(Bipolar$Q3_1==1,1,0)
Bipolar$Q3_2=ifelse(Bipolar$Q3_2==1,1,0)
Bipolar$Q3_3=ifelse(Bipolar$Q3_3==1,1,0)
Bipolar$Q3_4=ifelse(Bipolar$Q3_4==1,1,0)
Bipolar$Q3_5=ifelse(Bipolar$Q3_5==1,1,0)
Bipolar$Q3_6=ifelse(Bipolar$Q3_6==1,1,0)

Bipolar$LA=6-(rowSums(Bipolar[,c('Q3_1','Q3_2','Q3_3','Q3_4','Q3_5','Q3_6')])) #or LA=rowSums(Bipolar[,19:24])

# Figure 1
ggplot(Bipolar, aes(x=Group, y=LA)) +
  xlab("\nGroups") + ylab("Number of reject choices") +
  scale_x_discrete(labels=c("BD" = "Bipolar", "CTL" = "Control")) +
  theme(axis.text.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.x = element_text(colour="#000000", size=10, face="bold")) +
  theme(axis.title.y = element_text(colour="#000000", size=10, face="bold")) +
  theme(title = element_text(colour="#000000", size=12, face="bold")) +
  ggtitle("Loss aversion in bipolar vs controls") +
  theme(plot.title = element_text(hjust = 0.5)) +
  geom_violin() +    
  #geom_boxplot(width=0.2) + 
  stat_summary(fun.y=mean, geom="point", shape=23, size=4, color="red") 

# quite low loss avesion in overall sample: 2.6 reject choices
summary(Bipolar$LA)
sd(Bipolar$LA)
aggregate(Bipolar$LA, by=list(group=Bipolar$Group),mean)
aggregate(Bipolar$LA, by=list(group=Bipolar$Group),sd)

# Overall sample: is sig from 3 reject gambles, which implies mean levels of loss aversion.
t.test (Bipolar$LA, mu=3) #non-sig

# Difference in loss aversion between our groups --> no
t.test(Bipolar$LA ~ Bipolar$Group) #non-sig
wilcox.test(Bipolar$LA ~ Bipolar$Group) #non-sig

#Main OLS model to control for multiple variables at once: reported in main text
GenderLA=lm(LA ~ Group*Gender + Age, data=Bipolar)
summary(GenderLA)

GenderLA=lm(LA ~ Group*Gender + Age + bdisf_tot, data=Bipolar)
summary(GenderLA)

GenderLA=lm(LA ~ Group*Gender + Age + bdisf_tot + asrm_tot + stai_state_brief_tot, data=Bipolar)
summary(GenderLA)

# Participants who choose 6 gamble options is rather strange, subset these and rerun the analyses

Loss_seekers=subset(Bipolar, Bipolar$LA==0) #15 participants who are loss seekers
Bipolar$LossSeeker=ifelse(Bipolar$LA==0,1,0)
Bipolar.nls=subset(Bipolar, Bipolar$LA>0)

# who are those loss seekers?
table(Loss_seekers$Group,Loss_seekers$Gender) # both groups and gender represented
# normal symptom scores
summary(Bipolar.nls$asrm_tot) 
wilcox.test(Bipolar$asrm_tot ~ Bipolar$LossSeeker) #non-sig
summary(Bipolar.nls$stai_state_brief_tot)
wilcox.test(Bipolar$stai_state_brief_tot ~ Bipolar$LossSeeker) #p=0.1042
summary(Bipolar.nls$bdisf_tot)
wilcox.test(Bipolar$bdisf_tot ~ Bipolar$LossSeeker) #p=0.1793

# Did loss seekers not pay attention? How did they score on the catch item?
table(Bipolar$Catch_item_exclude,Bipolar$LossSeeker) #only 1 out 15 did not pay attention.

# After removing loss seekers, more in line with standard loss averse behavior
summary(Bipolar.nls$LA)
sd(Bipolar.nls$LA)
aggregate(Bipolar.nls$LA, by=list(group=Bipolar.nls$Group),mean)
aggregate(Bipolar.nls$LA, by=list(group=Bipolar.nls$Group),sd)

#OLS model without loss seekers --> Group not sig either
GenderLA=lm(LA ~ Group*Gender + Age + bdisf_tot + asrm_tot + stai_state_brief_tot, data=Bipolar.nls)
summary(GenderLA)

# -------------------------------------------------------------------------------------------------------

# Framing effects 
# typical choice pattern is B (2) Q4 then switch to A (1) Q5
# Sig group, gender and interaction differences for smaller sample

AD1=Bipolar[,c(1,11,13,15,17,18,25,26,27,30)]
AD1$APQ='Q4'
colnames(AD1)[2]="Choice"

AD2=Bipolar[,c(1,12,13,15,17,18,25,26,27,30)]
AD2$APQ='Q5'
colnames(AD2)[2]="Choice"
AD=rbind(AD1,AD2)

AD$Choice=factor(AD$Choice)
AD$Group=factor(AD$Group)
AD$Gender=factor(AD$Gender)
AD$APQ=factor(AD$APQ)
AD$LossSeeker=factor(AD$LossSeeker)

# Smaller sample when loss seeker are removed

AD1.nls=Bipolar.nls[,c(1,11,13,15,17,18,25,26,27,30)]
AD1.nls$APQ='Q4'
colnames(AD1.nls)[2]="Choice"

AD2.nls=Bipolar.nls[,c(1,12,13,15,17,18,25,26,27,30)]
AD2.nls$APQ='Q5'
colnames(AD2.nls)[2]="Choice"
AD.nls=rbind(AD1.nls,AD2.nls)

AD.nls$Choice=factor(AD.nls$Choice)
AD.nls$Group=factor(AD.nls$Group)
AD.nls$Gender=factor(AD.nls$Gender)
AD.nls$APQ=factor(AD.nls$APQ)

# Models
baselineQ <- glm(Choice ~ APQ,family=binomial(link='logit'),data=AD) #choice super sig
summary(baselineQ)

  # Sample without loss seekers similar results
  baselineQ <- glm(Choice ~ APQ,family=binomial(link='logit'),data=AD.nls)
  summary(baselineQ)

# Group not sig
baselineGroup <- glm(Choice ~ APQ*Group,family=binomial(link='logit'),data=AD)
summary(baselineGroup)

  # Sample without loss seekers similar results
  baselineGroup <- glm(Choice ~ APQ*Group,family=binomial(link='logit'),data=AD.nls) 
  summary(baselineGroup)

# GroupxGender sig: disregard main effects of group and gender, they are influenced by significance of interaction.
baselineGroup <- glm(Choice ~ APQ*Group*Gender,family=binomial(link='logit'),data=AD)
summary(baselineGroup)

  # Sample without loss seekers similar results
  baselineGroup <- glm(Choice ~ APQ*Group*Gender,family=binomial(link='logit'),data=AD.nls) #group x gender sig
  summary(baselineGroup)

#Symptom scores included: interaction remains sig
baselineGroup <- glm(Choice ~ APQ*Group*Gender + asrm_tot + bdisf_tot + stai_state_brief_tot,family=binomial(link='logit'),data=AD)
summary(baselineGroup)

  # # Sample without loss seekers similar results
  baselineGroup <- glm(Choice ~ APQ*Group*Gender + asrm_tot + bdisf_tot + stai_state_brief_tot,family=binomial(link='logit'),data=AD.nls)
  summary(baselineGroup)

# Visualise the effect: Fig 2

# Main effect

Q1=with(AD, c(sum(APQ == "Q4"))) #80
Q2=with(AD, c(sum(APQ == "Q5")))

Q1_A=with(AD, c(sum(APQ == "Q4" & Choice==2))) #43
Q2_A=with(AD, c(sum(APQ == "Q5" & Choice==2))) #29
PercQ1_A=Q1_A/Q1
PercQ2_A=Q2_A/Q2

Choice=c('die', 'save')
PercCEchoice=c(PercQ1_A,PercQ2_A)
Barplot=data.frame(Choice,PercCEchoice)

ggplot(data=Barplot, aes(x=Choice, y=PercCEchoice, fill=Choice)) +
  scale_fill_manual(values=c("lightgrey", "darkgrey")) +               
  geom_bar(colour="black", stat="identity") +
  guides(fill=FALSE) + 
  xlab("\nFrames") + ylab("Percentage choosing risky option") +
  theme(axis.text.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.5)) +
  ggtitle("Risk taking under different frames") +
  theme(title = element_text(face="bold", colour="#000000", size=11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank()) 

# GroupXGender differences --> Fig. 3

BP.FE=with(AD, c(sum(Group=="BD" & Gender==2 & APQ == "Q4"))) #33
BP.FE_countCE=with(AD, c(sum(Group=="BD" & Gender==2 & APQ == "Q4" & Choice == 2))) #21
BP.FE_countPE=with(AD, c(sum(Group=="BD" & Gender==2 & APQ == "Q5" & Choice == 2))) #13
PercBP.FE.CE=BP.FE_countCE/BP.FE #0.69
PercBP.FE.PE=BP.FE_countPE/BP.FE #0.39

BP.MA=with(AD, c(sum(Group=="BD" & Gender==1 & APQ == "Q4"))) #16
BP.MA_countCE=with(AD, c(sum(Group=="BD" & Gender==1 & APQ == "Q4" & Choice == 2))) #6
BP.MA_countPE=with(AD, c(sum(Group=="BD" & Gender==1 & APQ == "Q5" & Choice == 2))) #5
PercBP.MA.CE=BP.MA_countCE/BP.MA #0.375
PercBP.MA.PE=BP.MA_countPE/BP.MA #0.31

CTL.FE=with(AD, c(sum(Group=="CTL" & Gender==2 & APQ == "Q4"))) #23
CTL.FE_countCE=with(AD, c(sum(Group=="CTL" & Gender==2 & APQ == "Q4" & Choice == 2))) #10
CTL.FE_countPE=with(AD, c(sum(Group=="CTL" & Gender==2 & APQ == "Q5" & Choice == 2))) #9
PercCTL.FE.CE=CTL.FE_countCE/CTL.FE #0.69
PercCTL.FE.PE=CTL.FE_countPE/CTL.FE #0.39

CTL.MA=with(AD, c(sum(Group=="CTL" & Gender==1 & APQ == "Q4"))) #8
CTL.MA_countCE=with(AD, c(sum(Group=="CTL" & Gender==1 & APQ == "Q4" & Choice == 2))) #6
CTL.MA_countPE=with(AD, c(sum(Group=="CTL" & Gender==1 & APQ == "Q5" & Choice == 2))) #2
PercCTL.MA.CE=CTL.MA_countCE/CTL.MA #0.75
PercCTL.MA.PE=CTL.MA_countPE/CTL.MA #0.25

gendergroup=c('Bipolar_female_die', 'Bipolar_female_save','Bipolar_male_die', 'Bipolar_male_save','Control_female_die', 'Control_female_save','Control_male_die', 'Control_male_save')
PercCEchoice=c(PercBP.FE.CE,PercBP.FE.PE,PercBP.MA.CE,PercBP.MA.PE,PercCTL.FE.CE,PercCTL.FE.PE,PercCTL.MA.CE,PercCTL.MA.PE)
Barplot.CE.Gender=data.frame(gendergroup,PercCEchoice)

ggplot(data=Barplot.CE.Gender, aes(x=gendergroup, y=PercCEchoice, fill=gendergroup)) +
  scale_fill_manual(values=c("lightgrey", "darkgrey","lightgrey", "darkgrey","lightgrey", "darkgrey","lightgrey", "darkgrey")) +               
  geom_bar(colour="black", stat="identity") +
  guides(fill=FALSE) + 
  xlab("\nGroups") + ylab("Percentage within each group chosing risky option") +
  theme(axis.text.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.5)) +
  ggtitle("Interaction effect of group and gender under different frames") +
  theme(title = element_text(face="bold", colour="#000000", size=11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank()) 

