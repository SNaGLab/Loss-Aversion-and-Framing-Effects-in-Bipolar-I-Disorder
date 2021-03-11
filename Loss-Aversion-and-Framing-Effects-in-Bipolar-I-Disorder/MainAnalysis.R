

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

Bipolar=subset(Bipolar,Catch_item_exclude==0)

#Remove #201Y because ineligibility
Bipolar=Bipolar[-c(26),]

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

#Age --> quite comparable, 36.18 vs 34.28 in BD vs CTL
aggregate(Bipolar$Age, by=list(group=Bipolar$Group),mean) 
aggregate(Bipolar$Age, by=list(group=Bipolar$Group),sd) 
t.test(Bipolar$Age ~ Bipolar$Group) #p=0.473

# Gender --> Gender proportion equal between groups (1=M, 2=F)
table(Bipolar$Group,Bipolar$Gender)
prop.test(table(Bipolar$Group,Bipolar$Gender),correct=FALSE) #p=0.18

# Ethnicity --> White proportion quite equal between groups
table(Bipolar$Group,Bipolar$Caucasian)
prop.test(table(Bipolar$Group,Bipolar$Caucasian),correct=FALSE) #p=0.932

#Education --> sig difference: CTL more education years (16.00 vs 17.41)
Bipolar$Education=destring(Bipolar$Education)
Bipolar$Education[is.na(Bipolar$Education)] <- 0
BipolarEduc=subset(Bipolar,Bipolar$Education>0)
aggregate(BipolarEduc$Education, by=list(group=BipolarEduc$Group),mean) 
aggregate(BipolarEduc$Education, by=list(group=BipolarEduc$Group),sd) 
t.test(BipolarEduc$Education ~ BipolarEduc$Group) #p=0.02
wilcox.test(BipolarEduc$Education ~ BipolarEduc$Group) #p=0.03

#Employment
Bipolar$employment
table(Bipolar$Group,Bipolar$employment)
chisq.test(Bipolar$Group,Bipolar$employment, correct=FALSE)

# Income (coding: 1= Less than $10,000,2=$10,000 - $25,000,3=$26,000 - $50,000,4=$51,000 - $75,000,5=$76,000 - $100,000,6=More than $100,000)
Bipolar$income
table(Bipolar$Group,Bipolar$income)
chisq.test(Bipolar$Group,Bipolar$income, correct=FALSE)

# Nr children
Bipolar$children=gsub("\\P", "", Bipolar$children)
Bipolar$children=destring(Bipolar$children)
Bipolar$children[is.na(Bipolar$children)] <- 0
aggregate(Bipolar$children, by=list(group=Bipolar$Group),mean) 
aggregate(Bipolar$children, by=list(group=Bipolar$Group),sd) 
t.test(Bipolar$children ~ Bipolar$Group) #p=0.92

# Marital status
Bipolar$marital
table(Bipolar$Group,Bipolar$marital)
chisq.test(Bipolar$Group,Bipolar$marital, correct=FALSE)

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

t.test(Bipolar$asrm_tot ~ Bipolar$Group) #p=0.17
shapiro.test(Bipolar$asrm_tot) #non-normal, so perform wilcox test to verify t-test results.
wilcox.test(Bipolar$asrm_tot ~ Bipolar$Group) #p=0.11
t.test(Bipolar$asrm_tot ~ Bipolar$Gender) #p=0.15
t.test(Bipolar.sample$asrm_tot ~ Bipolar.sample$Gender) #not sig
t.test(Control.sample$asrm_tot ~ Control.sample$Gender) #not sig

# Depression scores
aggregate(Bipolar$bdisf_tot, by=list(group=Bipolar$Group),mean) #Depression score is quite higher for Bipolar
aggregate(Bipolar$bdisf_tot, by=list(group=Bipolar$Group),sd)
aggregate(Bipolar$bdisf_tot, by=list(group=Bipolar$Gender),mean) 
aggregate(Bipolar.sample$bdisf_tot, by=list(group=Bipolar.sample$Gender),mean) #similar scores between males and females in Bipolar
aggregate(Control.sample$bdisf_tot, by=list(group=Control.sample$Gender),mean) #similar scores between males and females in control

t.test(Bipolar$bdisf_tot ~ Bipolar$Group) #sig difference: BD more depressed
shapiro.test(Bipolar$bdisf_tot) #non-normal, so perform wilcox test to verify t-test results.
wilcox.test(Bipolar$bdisf_tot ~ Bipolar$Group) #sig difference: BD more depressed
t.test(Bipolar$bdisf_tot ~ Bipolar$Gender) #not sig

# Anxiety
aggregate(Bipolar$stai_state_brief_tot, by=list(group=Bipolar$Group),mean) #BD more anxious than controls
aggregate(Bipolar$stai_state_brief_tot, by=list(group=Bipolar$Group),sd)
aggregate(Bipolar.sample$stai_state_brief_tot, by=list(group=Bipolar.sample$Gender),mean) #similar scores between males and females in Bipolar
aggregate(Control.sample$stai_state_brief_tot, by=list(group=Control.sample$Gender),mean) #males quite more anxious than females in control group

t.test(Bipolar$stai_state_brief_tot ~ Bipolar$Group) #p=0.14, not sig different.
shapiro.test(Bipolar$stai_state_brief_tot) 
wilcox.test(Bipolar$stai_state_brief_tot ~ Bipolar$Group) #non-normal, so perform wilcox test to verify t-test results.
t.test(Bipolar$stai_state_brief_tot ~ Bipolar$Gender) #p=0.354
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

# quite low loss avesion in overall sample: 2.569 reject choices
summary(Bipolar$LA)
sd(Bipolar$LA)
aggregate(Bipolar$LA, by=list(group=Bipolar$Group),mean)
aggregate(Bipolar$LA, by=list(group=Bipolar$Group),sd)

# Overall sample: is sig from 3 reject gambles, which implies lower levels of loss aversion.
t.test (Bipolar$LA, mu=3) #p<.05
wilcox.test(Bipolar$LA, mu=3) #p<0.05.

# Difference in loss aversion between our groups --> no
t.test(Bipolar$LA ~ Bipolar$Group) #non-sig
wilcox.test(Bipolar$LA ~ Bipolar$Group) #non-sig

#Main OLS model to control for multiple variables at once: reported in main text
LA1=lm(LA ~ Group + Gender, data=Bipolar)
summary(LA1)

LA2=lm(LA ~ Group*Gender + Age, data=Bipolar)
summary(LA2)

LA3=lm(LA ~ Group*Gender + Age + bdisf_tot, data=Bipolar)
summary(LA3)

LA4=lm(LA ~ Group*Gender + Age + bdisf_tot + Education, data=Bipolar)
summary(LA4)

LA5=lm(LA ~ Group*Gender + Age + bdisf_tot + asrm_tot + stai_state_brief_tot + Education, data=Bipolar)
summary(LA5)

LA6=lm(LA ~ Group*Gender + Age + bdisf_tot + asrm_tot + stai_state_brief_tot + Education + Catch_item_exclude, data=Bipolar)
summary(LA6)

# Participants who choose 6 gamble options is rather strange, subset these and rerun the analyses

Loss_seekers=subset(Bipolar, Bipolar$LA==0) #14 participants who are loss seekers
Bipolar$LossSeeker=ifelse(Bipolar$LA==0,1,0)
Bipolar.nls=subset(Bipolar, Bipolar$LA>0)

# who are those loss seekers?
table(Loss_seekers$Group,Loss_seekers$Gender) # both groups and gender represented
wilcox.test(Bipolar$LossSeeker ~ Bipolar$Group) #non-sig
wilcox.test(Bipolar$LossSeeker ~ Bipolar$Gender) #non-sig
# normal symptom scores
summary(Bipolar.nls$asrm_tot) 
wilcox.test(Bipolar$asrm_tot ~ Bipolar$LossSeeker) #non-sig
summary(Bipolar.nls$stai_state_brief_tot)
wilcox.test(Bipolar$stai_state_brief_tot ~ Bipolar$LossSeeker) #p=0.303
summary(Bipolar.nls$bdisf_tot)
wilcox.test(Bipolar$bdisf_tot ~ Bipolar$LossSeeker) #p=0.276

# Did loss seekers not pay attention? How did they score on the catch item?
table(Bipolar$Catch_item_exclude,Bipolar$LossSeeker) #only 1 out 15 did not pay attention.

# After removing loss seekers, more in line with standard loss averse behavior
summary(Bipolar.nls$LA)
sd(Bipolar.nls$LA)
aggregate(Bipolar.nls$LA, by=list(group=Bipolar.nls$Group),mean)
aggregate(Bipolar.nls$LA, by=list(group=Bipolar.nls$Group),sd)

#OLS model without loss seekers --> Group not sig either
GenderLA=lm(LA ~ Group*Gender + Age + bdisf_tot + asrm_tot + stai_state_brief_tot + Education, data=Bipolar.nls)
summary(GenderLA)

# Difference in loss aversion between our groups --> no
t.test(Bipolar.nls$LA ~ Bipolar.nls$Group) #non-sig
wilcox.test(Bipolar.nls$LA ~ Bipolar.nls$Group) #non-sig

# Overall sample: is non sig from 3 reject gambles, which implies mean levels of loss aversion.
t.test (Bipolar.nls$LA, mu=3) #non-sig
wilcox.test(Bipolar.nls$LA, mu=3) #non-sig

# -------------------------------------------------------------------------------------------------------

## Framing effects in outcomes 
## typical choice pattern is B (2) Q4 then switch to A (1) Q5

# Prepare data 

AD1=Bipolar[,c(1,11,13,15,17,18,25,26,27,28,30)]
AD1$APQ='Q4'
colnames(AD1)[2]="Choice"

AD2=Bipolar[,c(1,12,13,15,17,18,25,26,27,28,30)]
AD2$APQ='Q5'
colnames(AD2)[2]="Choice"
AD=rbind(AD1,AD2)

AD$Choice=factor(AD$Choice)
AD$Group=factor(AD$Group)
AD$Gender=factor(AD$Gender)
AD$APQ=factor(AD$APQ)
AD$LossSeeker=factor(AD$LossSeeker)

# First check for group effect on switching behavior between die/save frame 

table(AD$APQ, AD$Choice)

#Existence of framing effect: yes!
prop.test(table(AD$APQ,AD$Choice),correct=FALSE) 
chisq.test(table(AD$APQ,AD$Choice))

# actual switch

ADdie=factor(Bipolar$Q4)
ADsave=factor(Bipolar$Q5)
ADswitch=ifelse(ADdie==ADsave,0,1)

# Framing effect (switch) across groups: no sig difference
prop.test(table(ADswitch,Bipolar$Group),correct=FALSE) 
chisq.test(table(ADswitch,Bipolar$Group))

# Smaller sample when loss seeker are removed

AD1.nls=Bipolar.nls[,c(1,11,13,15,17,18,25,26,27,28,30)]
AD1.nls$APQ='Q4'
colnames(AD1.nls)[2]="Choice"

AD2.nls=Bipolar.nls[,c(1,12,13,15,17,18,25,26,27,28,30)]
AD2.nls$APQ='Q5'
colnames(AD2.nls)[2]="Choice"
AD.nls=rbind(AD1.nls,AD2.nls)

AD.nls$Choice=factor(AD.nls$Choice)
AD.nls$Group=factor(AD.nls$Group)
AD.nls$Gender=factor(AD.nls$Gender)
AD.nls$APQ=factor(AD.nls$APQ)

# Logistical models
baselineQ <- glm(Choice ~ APQ,family=binomial(link='logit'),data=AD) #frame super sig
summary(baselineQ)

# Sample without loss seekers similar results
baselineQ <- glm(Choice ~ APQ,family=binomial(link='logit'),data=AD.nls)
summary(baselineQ)

# Group sig, but its interaction with frame not sig
baselineGroup <- glm(Choice ~ Group,family=binomial(link='logit'),data=AD) #group sig
summary(baselineGroup)   

baselineGroupInt <- glm(Choice ~ APQ*Group,family=binomial(link='logit'),data=AD)
summary(baselineGroupInt)

# Sample without loss seekers similar results
baselineGroup <- glm(Choice ~ Group,family=binomial(link='logit'),data=AD.nls) #group here sig, but when introducing frame, disappears, so not strong.
summary(baselineGroup)

baselineGroupInt <- glm(Choice ~ APQ*Group,family=binomial(link='logit'),data=AD.nls) 
summary(baselineGroupInt)

# Gender
baselineGender1 <- glm(Choice ~ Gender,family=binomial(link='logit'),data=AD)
summary(baselineGender1)   

baselineGender2 <- glm(Choice ~ APQ*Group+Gender,family=binomial(link='logit'),data=AD)
summary(baselineGender2)

# this model is reported in the paper
baselineGender3 <- glm(Choice ~ APQ*Group+Gender+ Age + Education + asrm_tot + bdisf_tot + stai_state_brief_tot,family=binomial(link='logit'),data=AD)
summary(baselineGender3)

# Sample without loss seekers similar results
baselineGender4 <- glm(Choice ~ Gender,family=binomial(link='logit'),data=AD.nls) 
summary(baselineGender4)

baselineGender5 <- glm(Choice ~ APQ*Group+Gender,family=binomial(link='logit'),data=AD.nls) 
summary(baselineGender5)

baselineGender6 <- glm(Choice ~ APQ*Group+Gender + Age + Education + asrm_tot + bdisf_tot + stai_state_brief_tot,family=binomial(link='logit'),data=AD.nls) 
summary(baselineGender6)

# symptom scores do not play a role
baselinesympt1 <- glm(Choice ~ APQ*Group+Gender+ Age + Education + asrm_tot*APQ + bdisf_tot*APQ + stai_state_brief_tot*APQ,family=binomial(link='logit'),data=AD)
summary(baselinesympt1)

# Visualise the effects

# Main effect: Figure 2

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
  ggtitle("Framing effects in outcomes") +
  theme(title = element_text(face="bold", colour="#000000", size=11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank()) 

# Group differences --> general effect of group on risk taking

BP=with(AD, c(sum(Group=="BD"))) #66
BP_countrisk=with(AD, c(sum(Group=="BD" & Choice == 2))) #38
PercBP.risk=BP_countrisk/BP #0.43

HC=with(AD, c(sum(Group=="CTL"))) #46
HC_countrisk=with(AD, c(sum(Group=="CTL" & Choice == 2))) #28
PercHC.risk=HC_countrisk/HC #0.64

Bar_Group=c('Bipolar', 'Control')
PercCEchoice=c(PercBP.risk,PercHC.risk)
Barplot=data.frame(Bar_Group,PercCEchoice)

ggplot(data=Barplot, aes(x=Bar_Group, y=PercCEchoice, fill=Choice)) +
  scale_fill_manual(values=c("lightgrey", "darkgrey")) +               
  geom_bar(colour="black", stat="identity") +
  guides(fill=FALSE) + 
  xlab("\nFrames") + ylab("Percentage choosing risky option") +
  theme(axis.text.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.5)) +
  ggtitle("Risk taking between groups") +
  theme(title = element_text(face="bold", colour="#000000", size=11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank()) 

# Frame x Group effect: Figure 3

BP=with(AD, c(sum(Group=="BD" & APQ == "Q4"))) 
BP_countdie=with(AD, c(sum(Group=="BD" & APQ == "Q4" & Choice==2))) 
BP_countsave=with(AD, c(sum(Group=="BD" & APQ == "Q5" & Choice==2))) 
PercBP.die=BP_countdie/BP
PercBP.save=BP_countsave/BP

CTL=with(AD, c(sum(Group=="CTL" & APQ == "Q4"))) 
CTL_countdie=with(AD, c(sum(Group=="CTL" & APQ == "Q4" & Choice==2))) 
CTL_countsave=with(AD, c(sum(Group=="CTL" & APQ == "Q5" & Choice==2))) 
PercCTL.die=CTL_countdie/CTL
PercCTL.save=CTL_countsave/CTL

gendergroup=c('Bipolar_die', 'Bipolar_save', 'Control_die', 'Control_save')
PercCEchoice=c(PercBP.die,PercBP.save,PercCTL.die,PercCTL.save)
Barplot.AD.Group=data.frame(gendergroup,PercCEchoice)

ggplot(data=Barplot.AD.Group, aes(x=gendergroup, y=PercCEchoice, fill=gendergroup)) +
  scale_fill_manual(values=c("lightgrey", "darkgrey", "lightgrey", "darkgrey")) +               
  geom_bar(colour="black", stat="identity") +
  guides(fill=FALSE) + 
  xlab("\nGroups") + ylab("Percentage within each group chosing risky option") +
  theme(axis.text.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.5)) +
  ggtitle("Group vs Choice interaction in framing effects in outcomes") +
  theme(title = element_text(face="bold", colour="#000000", size=11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank()) 

# -------------------------------------------------------------------------------------------------------

## Framing effects in probabilities
# Allais paradox
# Typical choice pattern is A - D (certainty effect).
# Typically, we observe a switch, which is in violation with Expected utility.
# Test for choice pattern A, B and occurence of switch between BD and controls.
# Q*gender*group sig in overall sample, not in smaller sample without loss seekers, 
# but remains sig when adding loss seeker as a grouping variable without reducing sample size


# Prepare data 

# overall sample
AP1=Bipolar[,c(1,3,13,15,17,18,25,26,27,28,30)]
AP1$APQ='Q1'
colnames(AP1)[2]="Choice"

AP2=Bipolar[,c(1,4,13,15,17,18,25,26,27,28,30)]
AP2$APQ='Q2'
colnames(AP2)[2]="Choice"
AP=rbind(AP1,AP2)

# Make grouping variables factors 
AP$Choice=factor(AP$Choice)
AP$Group=factor(AP$Group)
AP$Gender=factor(AP$Gender)
AP$APQ=factor(AP$APQ)
AP$LossSeeker=factor(AP$LossSeeker)

# loss seekers removed as separate dataset, but can also add as a grouping variable to the whole data set
AP1.nls=Bipolar.nls[,c(1,3,13,15,17,18,25,26,27,28)]
AP1.nls$APQ='Q1'
colnames(AP1.nls)[2]="Choice"

AP2.nls=Bipolar.nls[,c(1,4,13,15,17,18,25,26,27,28)]
AP2.nls$APQ='Q2'
colnames(AP2.nls)[2]="Choice"
AP.nls=rbind(AP1.nls,AP2.nls)

AP.nls$Choice=factor(AP.nls$Choice)
AP.nls$Group=factor(AP.nls$Group)
AP.nls$Gender=factor(AP.nls$Gender)
AP.nls$APQ=factor(AP.nls$APQ)

# First check for group effect on switching behavior 

table(AP$APQ, AP$Choice)

#Existence of framing effect: marginally sig with 0.06
prop.test(table(AP$APQ,AP$Choice),correct=FALSE) 
chisq.test(table(AP$APQ,AP$Choice))

# actual switch

APQ1=factor(Bipolar$Q1)
APQ2=factor(Bipolar$Q2)
APswitch=ifelse(APQ1==APQ2,0,1)

# Framing effect (switch) across groups: no sig difference
prop.test(table(APswitch,Bipolar$Group),correct=FALSE) 
chisq.test(table(APswitch,Bipolar$Group))

# Logistical models
baselineAPQ <- glm(Choice ~ APQ,family=binomial(link='logit'),data=AP) #frame sig
summary(baselineAPQ)

# Sample without loss seekers similar results, albeit less sig.
baselineAPQ <- glm(Choice ~ APQ,family=binomial(link='logit'),data=AP.nls)
summary(baselineAPQ)

# Group and its interaction with frame not sig
baselineAPGroup <- glm(Choice ~ Group,family=binomial(link='logit'),data=AP) 
summary(baselineAPGroup)   

baselineAPGroupInt <- glm(Choice ~ APQ*Group,family=binomial(link='logit'),data=AP)
summary(baselineAPGroupInt)

# Sample without loss seekers similar results
baselineAPGroup <- glm(Choice ~ Group,family=binomial(link='logit'),data=AP.nls) 
summary(baselineAPGroup)

baselineAPGroupInt <- glm(Choice ~ APQ*Group,family=binomial(link='logit'),data=AP.nls) 
summary(baselineAPGroupInt)

# Gender
baselineAPGender1 <- glm(Choice ~ Gender,family=binomial(link='logit'),data=AP)
summary(baselineAPGender1)   

baselineAPGender2 <- glm(Choice ~ APQ*Group+Gender,family=binomial(link='logit'),data=AP)
summary(baselineAPGender2)

# this model is reported in the paper
baselineAPGender3 <- glm(Choice ~ APQ*Group+Gender+ Age + Education + asrm_tot + bdisf_tot + stai_state_brief_tot,family=binomial(link='logit'),data=AP)
summary(baselineAPGender3)

# symptom scores
baselineAPSymp <- glm(Choice ~ APQ*Group+Gender+ Age + Education + APQ*asrm_tot + APQ*bdisf_tot + APQ*stai_state_brief_tot,family=binomial(link='logit'),data=AP)
summary(baselineAPSymp)

# Sample without loss seekers similar results
baselineAPGender4 <- glm(Choice ~ Gender,family=binomial(link='logit'),data=AP.nls) 
summary(baselineAPGender4)

baselineAPGender5 <- glm(Choice ~ APQ*Group+Gender,family=binomial(link='logit'),data=AP.nls) 
summary(baselineAPGender5)

baselineAPGender6 <- glm(Choice ~ APQ*Group+Gender + Age + Education + asrm_tot + bdisf_tot + stai_state_brief_tot,family=binomial(link='logit'),data=AP.nls) 
summary(baselineAPGender6)

# Figures

## Frame effect: Fig 4
Q1=with(AP, c(sum(APQ == "Q1"))) #80
Q2=with(AP, c(sum(APQ == "Q2")))

Q1_A=with(AP, c(sum(APQ == "Q1" & Choice==1))) #43
Q2_A=with(AP, c(sum(APQ == "Q2" & Choice==1))) #29
PercQ1_A=Q1_A/Q1
PercQ2_A=Q2_A/Q2

Choice=c('Lottery 1', 'Lottery 2')
PercCEchoice=c(PercQ1_A,PercQ2_A)
Barplot=data.frame(Choice,PercCEchoice)

ggplot(data=Barplot, aes(x=Choice, y=PercCEchoice, fill=Choice)) +
  scale_fill_manual(values=c("lightgrey", "darkgrey")) +               
  geom_bar(colour="black", stat="identity") +
  guides(fill=FALSE) + 
  xlab("\nChoices") + ylab("Overall percentage chosing option 'A'") +
  theme(axis.text.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.5)) +
  ggtitle("Framing effects in probabilities") +
  theme(title = element_text(face="bold", colour="#000000", size=11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank()) 

# Group effect

BP=with(AP, c(sum(Group=="BD"))) 
BP_countCE=with(AP, c(sum(Group=="BD" & Choice == 1))) 
PercBP.CE=BP_countCE/BP 

HC=with(AP, c(sum(Group=="CTL")))
HC_countCE=with(AP, c(sum(Group=="CTL" & Choice == 1))) 
PercHC.CE=HC_countCE/HC 

Bar_Group=c('Bipolar', 'Control')
PercCEchoice=c(PercBP.CE,PercHC.CE)
Barplot=data.frame(Bar_Group,PercCEchoice)

ggplot(data=Barplot, aes(x=Bar_Group, y=PercCEchoice, fill=Choice)) +
  scale_fill_manual(values=c("lightgrey", "darkgrey")) +               
  geom_bar(colour="black", stat="identity") +
  guides(fill=FALSE) + 
  xlab("\nFrames") + ylab("Percentage choosing option A") +
  theme(axis.text.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.5)) +
  ggtitle("Framing effects in probabilities") +
  theme(title = element_text(face="bold", colour="#000000", size=11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank()) 

# Frame x Group effect

BP=with(AP, c(sum(Group=="BD" & APQ == "Q1"))) #49
BP_countCE=with(AP, c(sum(Group=="BD" & APQ == "Q1" & Choice==1))) #27
BP_countPE=with(AP, c(sum(Group=="BD" & APQ == "Q2" & Choice==1))) #18
PercBP.CE=BP_countCE/BP
PercBP.PE=BP_countPE/BP

CTL=with(AP, c(sum(Group=="CTL" & APQ == "Q1"))) #31
CTL_countCE=with(AP, c(sum(Group=="CTL" & APQ == "Q1" & Choice==1))) #16
CTL_countPE=with(AP, c(sum(Group=="CTL" & APQ == "Q2" & Choice==1))) #11
PercCTL.CE=CTL_countCE/CTL
PercCTL.PE=CTL_countPE/CTL

gendergroup=c('Bipolar_lottery 1', 'Bipolar_lottery 2', 'Control_lottery 1', 'Control_lottery 2')
PercCEchoice=c(PercBP.CE,PercBP.PE,PercCTL.CE,PercCTL.PE)
Barplot.AP.Group=data.frame(gendergroup,PercCEchoice)

ggplot(data=Barplot.AP.Group, aes(x=gendergroup, y=PercCEchoice, fill=gendergroup)) +
  scale_fill_manual(values=c("lightgrey", "darkgrey", "lightgrey", "darkgrey")) +               
  geom_bar(colour="black", stat="identity") +
  guides(fill=FALSE) + 
  xlab("\nGroups") + ylab("Percentage within each group chosing option 'A'") +
  theme(axis.text.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.x = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.title.y = element_text(face="bold", colour="#000000", size=10)) +
  theme(axis.line.x = element_line(colour = "black", size = 0.5)) +
  ggtitle("Group vs Choice interaction in framing effects in probabilities") +
  theme(title = element_text(face="bold", colour="#000000", size=11)) +
  theme(plot.title = element_text(hjust = 0.5)) +
  theme(axis.ticks = element_blank()) 
