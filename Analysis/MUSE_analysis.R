##*****************************************##
## Analysis of MUSE-SAP study            ####
## James Liley                             ##
## 26 October 2022                         ##
##*****************************************##
## 
## Please note that this pipeline assumes that the raw data file is saved without a password.
## 

##*****************************************##
## Packages and scripts                  ####
##*****************************************##

# Packags
library(lme4) # Linear mixed models
library(readxl) # Read excel files

# Main functions
source("Analysis/MUSE_functions.R")


# Directories
raw_data_dir="./Data/Raw/" # Look for raw data CSV in this folder
working_data_dir="./Data/Working/" # Save working data to this folder
save_dir="./Analysis/Results/Figures/" # Directory to which figures are saved
text_dir="./Analysis/Results/output.txt" # Write to this text file; if NULL, write to console

# Switches and options: general
recheck_conf_int=FALSE # Set to TRUE to redo a simulation to verify confidence interval cover. If FALSE, will check if a simulation has already been performed and saved.
redo_power_calculation=FALSE # Set to TRUE to do a power calculation for a given estimated effect size. If FALSE, will check if one has already been performed and saved.
alpha=0.05 # Type 1 error rate for confidence intervals/sample size calculations
beta=0.1 # Type 2 error rate for sample size calculations
ndig=4 # Number of significant figures for outputs
npower=1000 # Number of simulations at each sample size for power calculations

# Switches and options - specific
sf36_use_uk=TRUE # Set to TRUE to use UK data to generate MCS and PCS for SF36; otherwise use US
sf36_use_z=FALSE # Set to TRUE to express SF36 MCS and PCS as Z scores; FALSE to use T scores, where T=10Z + 50.

# Setup
if (!is.null(text_dir)) sink(text_dir)
cat("Starting...\n\n")



# Set width: should be viewed in a monospace font without line breaks
options(width=100)



# Print time and date

cat("Date:\n")
print(Sys.Date())

cat("\n\nTime:\n")
print(Sys.time())

cat("\n\n\n")



# Print R and package versions

cat("R session info:\n\n")
print(sessionInfo())

cat("\n\n\n")



##*****************************************##
## Read data                             ####
##*****************************************##

# Main MUSE data: three sheets
mfile=paste0(raw_data_dir,"19 October 22 MUSE FEP combined data.xlsx")
muse_baseline=suppressMessages(read_xlsx(mfile,sheet=1))
muse_8week=suppressMessages(read_xlsx(mfile,sheet=2))
muse_12week=suppressMessages(read_xlsx(mfile,sheet=3))
muse_raw=list(muse_baseline,muse_8week,muse_12week)

# Data on number of modes of hallucination
modefile=paste0(raw_data_dir,"number of modalities.xlsx")
modes_raw= suppressMessages(read_xlsx(modefile))

# Randomisation
randfile=paste0(raw_data_dir,"crfRandomisation.csv")
rand_raw= suppressMessages(read.csv(randfile,encoding="UTF-8",stringsAsFactors=FALSE))

# Number of sessions
n_session_file=paste0(raw_data_dir,"Adherence Checklist Dataset 09 22.11.22.xlsx")
n_session_raw=suppressMessages(read_xlsx(n_session_file,sheet=1))

# Demographics
dem_file=paste0(raw_data_dir,"Copy of Combined Demographics v2 for James.xlsx")
dem_raw=suppressMessages(read_xlsx(dem_file,sheet=1))

cat("Loaded data\n\n")

# SF36 mean scores in various populations, used for SF36 processing
sf36file=paste0(working_data_dir,"sf36_distributions.csv")
sf36_dist= suppressMessages(read.csv(sf36file,row.names=1))
if (sf36_use_uk) pref="UK_" else pref="US_"
sf36_dist=sf36_dist[,grep(pref,colnames(sf36_dist))]
colnames(sf36_dist)=gsub(pref,"",colnames(sf36_dist))


##*****************************************##
## Basic processing                      ####
##*****************************************##

muse=list()
for (i in 1:length(muse_raw)) {
  
  ## Get data for corresponding time point
  dat=muse_raw[[i]]
  
  ## Rename columns
  dnames=make.names(muse_raw[[1]][1,])
  col_names=c("site","id","group",
              NA,paste0("psyrats_voices_",dnames[5:15]),"psyrats_voices_total",rep(NA,4),
              NA,paste0("psyrats_delusions_",dnames[22:27]),"psyrats_delusions_total",
              paste0("hpsvq_",dnames[29:37]),"hpsvq_total",NA,
              paste0("dass21_",dnames[40:60]),"dass21_total",NA,NA,NA,
              NA,paste0("choice_",dnames[66:89]),
                NA,paste0("choice_",dnames[90],"_how_do_you_feel"),
                NA,paste0("choice_",dnames[92],"_how_do_you_feel"),NA,
                paste0("choice_",dnames[95:118]),
                NA,paste0("choice_",dnames[119],"_how_satisfied"),
                NA,paste0("choice_",dnames[121],"_how_satisfied"),
              paste0("qpr_",dnames[123:137]),"qpr_total",
              paste0("sf36_",dnames[139:140]),NA,
                paste0("sf36_",dnames[142:151]),NA,              
                paste0("sf36_",dnames[153:156]),NA,
                paste0("sf36_",dnames[158:173]),NA,
                paste0("sf36_",dnames[175:178]),"sf36_total",
              paste0("eq5d_",dnames[180:185]),"eq5d_total",
              paste0("icecap_",dnames[187:191]),"icecap_total")
  if (i %in% 2:3) col_names=col_names[-c(3,17:20,39,62:64)]
  
  colnames(dat)=col_names
    
  ## Remove unnecessary rows
  if (i==1) dat=dat[3:dim(dat)[1],] else dat=dat[2:dim(dat)[1],]
  
  ## Remove unnecessary columns
  dat=dat[,which(!is.na(colnames(dat)))]
  
  ## Convert to data frame
  dat=as.data.frame(dat)
  
  ## Convert all columns to integers
  for (j in 1:dim(dat)[2]) dat[,j]=suppressWarnings(as.numeric(gsub("[^0-9.-]", "", dat[,j])))

  ## Convert 999s etc to NA
  for (j in 1:dim(dat)[2]) dat[which(dat[,j]>998),j]=NA
  
    
  ## Save
  muse[[i]]=dat
}
names(muse)=names(muse_raw)


## Modalities data: affix to main table
modes=as.data.frame(modes_raw[3:84,c(1,2,4,8,11)])
colnames(modes)=c("site","id","modes_baseline","modes_8week","modes_12week")
modes$id=as.numeric(modes$id)
# correct typos
auditory=c("Auditory", "Aufitory", "Auditrory", "Auditry")
visual=c("Visual", "visual")
somatic=c("Somatic", "somatic",  "soamtic")
olfactory=c("Olfactory", "olfactory", "oflactory")
felt_presence=c("Felt presence","Felt sense", "felt presence", 
                "Felt Presence", "felt sense", "fet presence", "felt oresence")
mode_missing=c("999")
for (i in 1:3) {
  mtab=muse[[i]]
  m=match(modes$id,mtab$id)
  mtab$modality_auditory[m]=anygrepl(auditory,modes[[i+2]],miss=mode_missing)
  mtab$modality_visual[m]=anygrepl(visual,modes[[i+2]],miss=mode_missing)
  mtab$modality_somatic[m]=anygrepl(somatic,modes[[i+2]],miss=mode_missing)
  mtab$modality_olfactory[m]=anygrepl(olfactory,modes[[i+2]],miss=mode_missing)
  mtab$modality_felt_presence[m]=anygrepl(felt_presence,modes[[i+2]],miss=mode_missing)
  muse[[i]]=mtab
}

# Randomisation data: affix to main table
rand_id=rand_raw$identifier
muse[[1]]$dateRandomised=rand_raw$dateRandomised[match(rand_id,muse[[1]]$id)]

# Number of sessions: affix to main table
n_sessions_id=n_session_raw$`Participant ID` # ID
n_sessions_muse1=n_session_raw[[65]]        # Number of sessions (recorded in one of two columns)
n_sessions_muse2=n_session_raw[[66]]        # Number of sessions (in one of two columns)
w=which(!is.na(n_sessions_id)); 
n_sessions_id=as.numeric(n_sessions_id[w]); 
n_sessions_muse1=n_sessions_muse1[w]
n_sessions_muse2=n_sessions_muse2[w]

# Combine the two columns in which 'total' is recorded
w1=which(is.na(n_sessions_muse1))
n_sessions_muse=n_sessions_muse1; n_sessions_muse[w1]=n_sessions_muse2[w1]

# To numeric
n_sessions_muse=as.numeric(n_sessions_muse)

# Add to table
nsm=rep(NA,dim(muse[[1]])[1])
nsm[match(n_sessions_id,muse[[1]]$id)]=n_sessions_muse
muse[[1]]$n_sessions_muse=nsm



# Demographic datas: affix to main table
dem=as.data.frame(dem_raw); colnames(dem)=NULL
dem_name=dem_raw[1,]
dem=dem[2:dim(dem)[1],]
colnames(dem)=make.names(dem_name)

# Process IDs
dem_id=dem$Randomisation.ID
dem_id=as.numeric(substring(gsub(" ","",dem_id),3,5))
dem$id=dem_id

dem_sex=rep(NA,dim(dem)[1])
dem_sex[which(dem$Sex=="Male")]="M"
dem_sex[which(dem$Sex=="Female")]="F"
muse[[1]]$sex=dem_sex[match(muse[[1]]$id,dem_id)]

dem_age=as.numeric(dem$Age)
muse[[1]]$age=dem_age[match(muse[[1]]$id,dem_id)]

# Baseline hallucination duration
dem_baseline_duration=rep(NA,dim(dem)[1])
num_only=suppressWarnings(as.numeric(gsub("[ a-zA-Z]","",dem$Duration.of.voices)))
yr=grepl("ear",dem$Duration.of.voices) # 'Years' type answers
ym=grepl("onth",dem$Duration.of.voices) # 'Months' type answers
dem_baseline_duration[which(yr)]=12*num_only[which(yr)]
dem_baseline_duration[which(ym)]=num_only[which(ym)]
muse[[1]]$baseline_hallucination_duration=dem_baseline_duration[match(muse[[1]]$id,dem_id)]

# Date first referred to service
dem_date=rep(as.Date("1900/01/01"),dim(dem)[1])
dem_date1=dem$Open.referred.to.Eip
w1=which(!is.na(suppressWarnings(as.numeric(dem_date1))))
dem_date[w1]=as.Date(as.numeric(dem_date1[w1])-2,origin = '1900-01-01')
w2=setdiff(1:dim(dem)[1],w1)
# Trim free text and standardise format
d2=dem_date1[w2]; d2=gsub("."," ",d2,fixed=TRUE)
wx=unlist(lapply(strsplit(d2," "),function(x) paste(x[1:min(3,length(x))],collapse=" ")))
wx=gsub("nd","",wx); wx=gsub("th","",wx); wx=gsub("st","",wx); 
mth=c("Jan","Feb","Mar","Apr","May","Jun","Jul","Aug","Sep","Oct","Nov","Dec")
for (i in 1:length(mth)) wx=gsub(mth[i],paste0(" ",i," "),wx)
for (i in 1:10) wx=gsub("  "," ",wx); wx=gsub(" ","/",wx)
dem_date[w2]=as.Date(wx,format="%d/%m/%Y")
dem_date[which(dem_date < as.Date("1901/01/01"))]=NA
muse[[1]]$dateService=dem_date[match(muse[[1]]$id,dem_id)]

covid_affected=rep(NA,dim(dem)[1])
covid_affected[which(dem$Treatment.delivered.during.Covid.period.1st.Dec.2021.onwards=="y")]=1
covid_affected[which(dem$Treatment.delivered.during.Covid.period.1st.Dec.2021.onwards=="n")]=0
muse[[1]]$covid_affected=covid_affected[match(muse[[1]]$id,dem_id)]


##*****************************************##
## Collate questionnaires                ####
##*****************************************##

muse_summary=list()
for (i in 1:length(muse)) { # Loop over times
  dat0=muse[[i]]
  dat=dat0[,1:2] # Site and ID
  
  ## Unblind treatment: group 2 is treated
  if (i==1) dat$treatment=as.numeric(muse[[1]]$group==2)

  ## Date of randomisation
  if (i==1) dat$dateRandomised=as.Date(muse[[1]]$dateRandomised)
    
  ## PSYRATS scores
  psyrats = dat0[,grep("psyrats_",colnames(dat0))]
  dat$psyrats_hdis=rsm(psyrats[,c(6,7,8,9,11)])
  dat$psyrats_hfreq=rsm(psyrats[,c(1,2,10)])
  dat$psyrats_hattrcog=rsm(psyrats[,c(3,5)])
  dat$psyrats_hloud=psyrats[,4]
  dat$psyrats_voices_total=rsm(psyrats[,1:11]) # Recalculate this
  dat$psyrats_delusions_total=rsm(psyrats[,13:18]) # Recalculate this
  dat$psyrats_total=rsm(psyrats[,c(1:11,13:18)])
  
  ## HPSVQ scores
  hpsvq = dat0[,grep("hpsvq",colnames(dat0))]
  dat$hpsvq_total = rsm(hpsvq[,c(1:9)]) # Recalculate
  dat$hpsvq_negative = rsm(hpsvq[,c(2,5,6,7)])
  
  ## DASS21 scores
  dass21 = dat0[,grep("dass21",colnames(dat0))]
  dat$dass21_total = rsm(dass21[,c(1:21)]) # Recalculate
  dat$dass21_stress = rsm(dass21[,c(1,6,8,11,12,14,18)])
  dat$dass21_anxiety = rsm(dass21[,c(2,4,7,9,15,19,20)])
  dat$dass21_depression = rsm(dass21[,c(3,5,10,13,16,17,21)])
  
  ## CHOICEs scores  
  choice = dat0[,grep("choice_",colnames(dat0))]
  dat$choice_short_total = rsm(choice[,c(1, 2, 3, 7, 9, 12, 17, 18, 20, 21, 22)])
  
  ## QPR scores
  qpr = dat0[,grep("qpr_",colnames(dat0))]
  dat$qpr_total = rsm(qpr[,1:15])
  
  ## SF36 scores
  sf36 = dat0[,grep("sf36_",colnames(dat0))]
  for (r in c(1, 2, 20, 22, 34, 36)) sf36[,r] = c(100,75,50,25,0)[sf36[,r]]
  for (r in c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12)) sf36[,r] = c(0,50,100)[sf36[,r]]
  for (r in c(13, 14, 15, 16, 17, 18, 19)) sf36[,r] = c(0,100)[sf36[,r]]
  for (r in c(21, 23, 26, 27, 30)) sf36[,r] = c(100,80,60,40,20,0)[sf36[,r]]
  for (r in c(24, 25, 28, 29, 31)) sf36[,r] = c(0,20,40,60,80,100)[sf36[,r]]
  for (r in c(32, 33, 35)) sf36[,r] = c(0,25,50,75,100)[sf36[,r]]
  dat$sf36_PF_physical = rsm(sf36[,c(3, 4, 5, 6, 7, 8, 9, 10, 11, 12)],sf36=TRUE)
  dat$sf36_RP_role_physical = rsm(sf36[,c(13, 14, 15, 16)],sf36=TRUE)
  dat$sf36_BP_pain = rsm(sf36[,c(21, 22)],sf36=TRUE)
  dat$sf36_GH_general = rsm(sf36[,c(1, 33, 34, 35, 36)],sf36=TRUE)
  dat$sf36_VT_vitality = rsm(sf36[,c(23, 27, 29, 31)],sf36=TRUE)
  dat$sf36_SF_social = rsm(sf36[,c(20, 32)],sf36=TRUE)
  dat$sf36_RE_role_emotional = rsm(sf36[,c(17, 18, 19)],sf36=TRUE)
  dat$sf36_MH_mental_health = rsm(sf36[,c(24, 25, 26, 28, 30)],sf36=TRUE)
  subdat = as.matrix(dat[,grep("sf36",colnames(dat))])
  # Rescale to Z-scores according to British/American data
  s_ord=gsub("sf36_","",colnames(subdat)) # Names of columns without 'sf36_'
  sf36_z=t((t(subdat) - sf36_dist[s_ord,]$mean)/sf36_dist[s_ord,]$sd)
  # Generate PCS and MCS scores
  if (sf36_use_z) sc=c(1,0) else sc=c(10,50)
  dat$sf36_pcs = (sf36_z %*% sf36_dist[s_ord,]$PCS_factor)*sc[1] + sc[2]
  dat$sf36_mcs = (sf36_z %*% sf36_dist[s_ord,]$MCS_factor)*sc[1] + sc[2]

  ## EQ5D scores
  eq5d = dat0[,grep("eq5d_",colnames(dat0))]
  dat$eq5d_mobility=rsm(eq5d[,'eq5d_Mobility',drop=F])
  dat$eq5d_self_care=rsm(eq5d[,'eq5d_Self.care',drop=F])
  dat$eq5d_usual_activities=rsm(eq5d[,'eq5d_Usual.activities',drop=F])
  dat$eq5d_pain_discomfort=rsm(eq5d[,'eq5d_Pain.Discomfort',drop=F])
  dat$eq5d_anxiety_depression=rsm(eq5d[,'eq5d_Anxiety.Depression',drop=F])
  dat$eq5d_health_today_vas=rsm(eq5d[,'eq5d_How.good.or.bad.is.your.health.today.',drop=F])
  dat$eq5d_total=rsm(eq5d[,'eq5d_total',drop=F])
  
  ## ICECAP
  icecap = dat0[,grep("icecap_",colnames(dat0))]
  dat$icecap_total = rsm(icecap[,1:5])
  
  ## Hallucination modalities
  modalities = dat0[,grep("modality_",colnames(dat0))]
  dat$modality_total=rsm(modalities)
  dat$modality_total_excl_felt_presence = rsm(modalities[,1:4])
  
  # Save summary table
  muse_summary[[i]]=dat

}

# Add number of sessions
muse_summary[[1]]$n_sessions_muse=muse[[1]]$n_sessions_muse[match(muse[[1]]$id,muse_summary[[1]]$id)]



# Add demographic/treatment information
muse_summary[[1]]$sex=muse[[1]]$sex[match(muse[[1]]$id,muse_summary[[1]]$id)]
muse_summary[[1]]$age=muse[[1]]$age[match(muse[[1]]$id,muse_summary[[1]]$id)]
muse_summary[[1]]$baseline_hallucination_duration=muse[[1]]$baseline_hallucination_duration[match(muse[[1]]$id,muse_summary[[1]]$id)]

time_in_service=as.numeric(as.Date(muse[[1]]$dateRandomised)-muse[[1]]$dateService)
time_in_service=round(time_in_service/365.25)*12 # Convert days to months
time_in_service[which(time_in_service<0)]=NA
muse_summary[[1]]$time_engaged_in_service=time_in_service[match(muse[[1]]$id,muse_summary[[1]]$id)]

muse_summary[[1]]$covid_affected=muse[[1]]$covid_affected[match(muse[[1]]$id,muse_summary[[1]]$id)]



# Correct site discrepancy
for (i in 1:length(muse_summary[[1]]$id)) {
  idx=muse_summary[[1]]$id[i]
  sx=muse_summary[[1]]$site[i]
  w2=which(muse_summary[[2]]$id==idx)
  w3=which(muse_summary[[3]]$id==idx)
  muse_summary[[2]]$site[i]=sx
  muse_summary[[3]]$site[i]=sx
}


## Specify other (non-PPO) outcomes
other_outcomes=c("psyrats_hdis",
                 "psyrats_delusions_total",        ## Added following SAP revision
                 "hpsvq_total","hpsvq_negative",
                 "dass21_total",
                 "dass21_stress",                  ## Added following SAP revision
                 "dass21_anxiety",                 ## Added following SAP revision
                 "dass21_depression",              ## Added following SAP revision
                 "eq5d_mobility",                  ## Added following SAP revision
                 "eq5d_self_care",                 ## Added following SAP revision
                 "eq5d_usual_activities",          ## Added following SAP revision
                 "eq5d_pain_discomfort",           ## Added following SAP revision
                 "eq5d_anxiety_depression",        ## Added following SAP revision
                 "eq5d_health_today_vas",          ## Added following SAP revision
                 "qpr_total",
                 "choice_short_total",
                 "sf36_pcs",
                 "sf36_mcs")

# Long names of outcomes
longnames=c(
  psyrats_voices_total =  "PSYRATs voices score, total",
  psyrats_hdis         =  "PSYRATs voice impact subscale",
  hpsvq_total          =  "HPSVQ total",
  hpsvq_negative       =  "HPSVQ negative aspect subscale",
  dass21_total         =  "DASS21 total",
  qpr_total            =  "QPR total",
  choice_short_total   =  "CHOICEs short form total",
  sf36_pcs             =  "SF-36 physical component score",
  sf36_mcs             =  "SF-36 mental component score"
)
sdiff=setdiff(other_outcomes,names(longnames))
nlongnames=names(longnames)
longnames=c(longnames,sdiff)
names(longnames)=c(nlongnames,sdiff)

# Outcomes for which a higher score is better
swap=c("choice_short_total","qpr_total","icecap_total")



cat("Completed initial data processing\n\n")



##*****************************************##
## Checks and summaries                  ####
##*****************************************##

cat("General data summary\n\n")

cat(paste0("Total number of participants in data: ",
           length(unique(muse_summary[[1]]$id)),"\n"))
cat(paste0("Total number of participants in treatment-as-usual: ",
           length(unique(muse_summary[[1]]$id[which(muse_summary[[1]]$treatment==0)])),"\n"))
cat(paste0("Total number of participants in MUSE treatment group ",
           length(unique(muse_summary[[1]]$id[which(muse_summary[[1]]$treatment==1)])),"\n"))
cat("\n")
cat(paste0("Total number of participants at site 1: ",
           length(unique(muse_summary[[1]]$id[which(muse_summary[[1]]$site==1)])),"\n"))
cat(paste0("Total number of participants at site 2: ",
           length(unique(muse_summary[[1]]$id[which(muse_summary[[1]]$site==2)])),"\n"))
cat("\n")
cat(paste0("Total number of participants with PSYRATS-voices total data at baseline: ",
           length(which(is.finite(muse_summary[[1]]$psyrats_voices_total))),"\n"))
cat(paste0("Total number of participants with PSYRATS-voices total data at 8 weeks: ",
           length(which(is.finite(muse_summary[[2]]$psyrats_voices_total))),"\n"))
cat(paste0("Total number of participants with PSYRATS-voices total data at 12 weeks: ",
           length(which(is.finite(muse_summary[[3]]$psyrats_voices_total))),"\n"))
cat("\n")




# Table: tabulate mean,median etc across these categories
all_id=muse_summary[[1]]$id
cats=list(
  all=all_id,
  TAU=all_id[which(muse_summary[[1]]$treatment==0)],
  MUSE=all_id[which(muse_summary[[1]]$treatment==1)],
  site1_all=all_id[which(muse_summary[[1]]$site==1)],
  site1_TAU=all_id[which((muse_summary[[1]]$site==1) & (muse_summary[[1]]$treatment==0))],
  site1_MUSE=all_id[which((muse_summary[[1]]$site==1) & (muse_summary[[1]]$treatment==1))],
  site2_all=all_id[which(muse_summary[[1]]$site==2)],
  site2_TAU=all_id[which((muse_summary[[1]]$site==2) & (muse_summary[[1]]$treatment==0))],
  site2_MUSE=all_id[which((muse_summary[[1]]$site==2) & (muse_summary[[1]]$treatment==1))]
)


# Tabulate values at baseline
btab=data.frame(
  id=muse_summary[[1]]$id,
  Age=muse_summary[[1]]$age,
  Sex_male=(muse_summary[[1]]$sex=="M"),
  Baseline_hallucination_duration=muse_summary[[1]]$baseline_hallucination_duration,
  Time_engaged_in_service=muse_summary[[1]]$time_engaged_in_service,
  N_sessions_muse=muse_summary[[1]]$n_sessions_muse
)
btab_format=c("real","real","binary","real","real","count")
tab_muse_baseline=c()
id1=btab$id
for (i in 2:length(colnames(btab))) {
    xcol=c()
    for (g in 1:length(cats)) {
      xid=match(cats[[g]],id1)
      if (length(which(is.na(btab[xid,i])))>5) summ="NA" else 
        summ=muse_format(btab[xid,i],type=btab_format[i],dig=ndig)
      xcol=c(xcol,summ)
    }
    tab_muse_baseline=rbind(tab_muse_baseline,xcol)
}
rownames(tab_muse_baseline)=colnames(btab)[2:dim(btab)[2]]
colnames(tab_muse_baseline)=c(names(cats))

cat("\n\n")
cat(paste0("Table of baseline values values across groups; format `n male/n all` for sex, median (IQR) for N_sessions_muse, `mean (SD) m. missing%` otherwise: ","\n\n"))
print(tab_muse_baseline,quote=FALSE)

cat("\n")
tfile="./Analysis/Results/Tables/baseline_data.csv"
write.csv(tab_muse_baseline,file=tfile,quote=FALSE)
cat(paste0("Table saved in: ",tfile))
cat("\n\n\n")





# Tabulate for all three time points
tab_cols3=setdiff(colnames(muse_summary[[1]]),
                  c("site","id","sex","age","treatment","dateRandomised","n_sessions_muse",
                    "baseline_hallucination_duration","time_engaged_in_service","covid_affected"))
tab_format=rep("real",length(tab_cols3))
tab_format[grep("modality",tab_cols3)]="count"
muse_table=c()
for (i in 1:length(tab_cols3)) {
  for (t in 1:3) {
    mtab=muse_summary[[t]]
    mid=mtab$id
    xcol=c()
    for (g in 1:length(cats)) {
      xid=match(cats[[g]],mid)
      summ=muse_format(mtab[xid,tab_cols3[i]],type=tab_format[i],dig=ndig)
      xcol=c(xcol,summ)
    }
    muse_table=rbind(muse_table,xcol)
  }
}
muse_table=cbind(rep(tab_format,each=3),muse_table)
colnames(muse_table)=c("variable_type",names(cats))
rownames(muse_table)=as.vector(outer(c("baseline_","8week_","12week_"),tab_cols3,paste0))

muse_summary_table=muse_table[
  c(grep("psyrats_total",rownames(muse_table)),
    grep("psyrats_voices_total",rownames(muse_table)),
    grep("psyrats_delusions_total",rownames(muse_table))),
  ]

cat(paste0("Table of summary values across groups; format `mean (SD) m. missing%`: ","\n\n"))
print(muse_summary_table,quote=FALSE)

cat("\n")
tfile="./Analysis/Results/Tables/general_table.csv"
write.csv(muse_table,file=tfile,quote=FALSE)
cat(paste0("Full table saved in: ",tfile))
cat("\n\n\n")



##***********************************************##
## Put together data frame for models          ####
##***********************************************##

# Matching indices for IDs, in case IDs not in order
m2=match(muse_summary[[1]]$id,muse_summary[[2]]$id)
m3=match(muse_summary[[1]]$id,muse_summary[[3]]$id)

# Collect data in one data frame
muse_matrix=data.frame(
  # Individual ID
  id=c(muse_summary[[1]]$id,
       muse_summary[[2]]$id,
       muse_summary[[3]]$id),
  
  # Time (1,2, or 3, corresponding to baseline, 8 weeks, 12 weeks respectively)
  time=as.factor(c(rep(1,dim(muse_summary[[1]])[1]),
                   rep(2,dim(muse_summary[[2]])[1]),
                   rep(3,dim(muse_summary[[3]])[1]))),
  
  # Site
  site=c(muse_summary[[1]]$site,
         muse_summary[[2]]$site,
         muse_summary[[3]]$site),

  # Treatment (only recorded in muse_summary[[1]])
  treatment=c(muse_summary[[1]]$treatment,
          muse_summary[[1]]$treatment[m2],
          muse_summary[[1]]$treatment[m3]),
  
  # Pseudo-primary outcome (PSYRATS total score)
  psyrats_voices_total=c(muse_summary[[1]]$psyrats_voices_total,
                  muse_summary[[2]]$psyrats_voices_total,
                  muse_summary[[3]]$psyrats_voices_total)
  
)

# Covariates
covariates=c("age","sex","modality_total_excl_felt_presence",
             "baseline_hallucination_duration","time_engaged_in_service",
             "psyrats_delusions_total")
for (i in 1:length(covariates)) {
  X=c(muse_summary[[1]][[covariates[i]]],
      muse_summary[[1]][[covariates[i]]][m2],
      muse_summary[[1]][[covariates[i]]][m3])
  muse_matrix[covariates[i]]=X
}

# Special treatment for psyrats_delusions_total at baseline
w_pdt=which(colnames(muse_matrix)=="psyrats_delusions_total")
colnames(muse_matrix)[w_pdt]="psyrats_delusions_total_baseline"



# Other outcomes
for (i in 1:length(other_outcomes)) {
  Y=c(muse_summary[[1]][[other_outcomes[i]]],
      muse_summary[[2]][[other_outcomes[i]]],
      muse_summary[[3]][[other_outcomes[i]]])
  muse_matrix[other_outcomes[i]]=Y
}

# Number of treatments
muse_matrix$n_sessions_muse=c(muse_summary[[1]]$n_sessions_muse,
                        muse_summary[[1]]$n_sessions_muse[m2],
                        muse_summary[[1]]$n_sessions_muse[m3])


# Time columns
muse_matrix$time2 = as.numeric(muse_matrix$time==2) # Time 8 weeks
muse_matrix$time3 = as.numeric(muse_matrix$time==3) # Time 12 weeks



##***********************************************##
## 10.2.1: Mixed linear model analysis PPO     ####
## No covariates included other than site        ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.1: Mixed linear model analysis PPO   **\n")
cat("**  No covariates included other than site   **\n")
cat("***********************************************\n")
cat("\n\n")


# Fit model
muse_ppo_lme=lmer(psyrats_voices_total ~ 
                    time2 + 
                    time3 + 
                    treatment:time2 + 
                    treatment:time3 + 
                    factor(site) + 
                    (1|id),
                  data=muse_matrix)

# Point estimates of fixed and random effects
muse_ppo_fixed = summary(muse_ppo_lme)$coefficients
muse_ppo_random = as.data.frame(VarCorr(muse_ppo_lme))


# Confidence intervals
muse_ppo_ci = suppressMessages(confint(muse_ppo_lme))

# Summary table
muse_est=c(muse_ppo_random[c(which(muse_ppo_random$grp=="id"),
                             which(muse_ppo_random$grp=="Residual")),5],
           muse_ppo_fixed[c("(Intercept)", 
                            "time2", 
                            "time3", 
                            "factor(site)2", 
                            "time2:treatment", 
                            "time3:treatment"),"Estimate"])
muse_name=c("σr","σe","β0","βa2","βa3","β(site)","β2","β3")
muse_description=c("SD of per-indivudal random effect",
                   "SD of residuals",
                   "Overall mean",
                   "Fixed effect for 8-week time point",
                   "Fixed effect for 12-week time point",
                   "Fixed effect for site 2 vs site 1",
                   "Fixed effect for treatment at 8-week time point",
                   "Fixed effect for treatment at 12-week time point")
out_tab=data.frame(signif(muse_est,digits=ndig),signif(muse_ppo_ci,digits=ndig),muse_description)
colnames(out_tab)=c("Estimate","CI95_lower","CI95_upper","Description")
rownames(out_tab)=muse_name

# Keep to use later
muse_psyrats_voices_total=out_tab


# Report
cat(paste0("Estimated beta_2 (mean effect of treatment on PSYRATS voices total score at 8 weeks) ",
           signif(muse_ppo_fixed["time2:treatment","Estimate"],digits=ndig)," points; ",
           " 95% CI: (",signif(muse_ppo_ci["time2:treatment","2.5 %"],digits=ndig),", ",
           signif(muse_ppo_ci["time2:treatment","97.5 %"],digits=ndig),") points"))
cat("\n\n")

cat(paste0("All parameter estimates and 95% CIs: ","\n\n"))
print(out_tab)
cat("\n\n")

ofile="./Analysis/Results/Tables/PSYRATs_voices_total/parameter_table_psyrats_voices_total.csv"
write.csv(out_tab,file=ofile,quote=FALSE)
cat(paste0("Parameter table saved to: ",ofile))
cat("\n\n\n")


##***********************************************##
## 10.2.1b: Mixed linear model analysis PPO    ####
## No covariates other than site                 ##
## Sample size calculation (empirical)           ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.1b: Mixed linear model analysis PPO  **\n")
cat("**  No covariates included other than site   **\n")
cat("**  Sample size calculation (empirical)      **\n")
cat("***********************************************\n")
cat("\n\n")



b2=signif(muse_ppo_fixed["time2:treatment","Estimate"],digits=3)
power_file=paste0("./Analysis/Results/R_objects/power_b2_",b2,"_simulation.RData")
if (redo_power_calculation==TRUE | (!file.exists(power_file))) {
  
  # Collect empirical parameter estimates
  muse_power_par=c(beta_0=muse_ppo_fixed[1,1],
                   beta_a1=0,
                   beta_a2=muse_ppo_fixed["time2",1],
                   beta_a3=muse_ppo_fixed["time3",1],
                   beta_2=muse_ppo_fixed["time2:treatment",1],
                   beta_3=muse_ppo_fixed["time3:treatment",1],
                   beta=muse_ppo_fixed["factor(site)2",1],
                   sigma_r=muse_ppo_random[1,5],
                   sigma_e=muse_ppo_random[2,5])
  
  # Simulate
  cat(paste0("Sample size calculation\n\n"))
  set.seed(39283)
  
  
  # Exponentially increase sample size until power>0.99
  sizes=5 
  power=pcheck(sizes[1])
  while(max(power)<0.99) {
    nsize=round(1.5*max(sizes))
    xpower=pcheck(nsize,npower,muse_power_par,alpha)
    sizes=c(sizes,nsize)
    power=c(power,xpower)
  }
  
  # We now have sample sizes at which power>0.99; now find maximum size for which power<0.2
  lsize=max(sizes[which(power<0.2)])
  
  # Now find minimum sample size at which power>1-beta.
  sseq=round(seq(lsize,max(sizes),length=10))
  pseq=rep(0,length(sseq))
  for (i in 1:length(sseq)) pseq[i]=suppressWarnings(pcheck(sseq[i],npower,muse_power_par,alpha))
  sizes=c(sizes,sseq)
  power=c(power,pseq)
  
  # Sort
  os=order(sizes)
  sizes=sizes[os]
  power=power[os]
  
  save(muse_power_par,sizes,power,file=power_file)
  
} else load(power_file)

# Approximate sample size
ssize=suppressWarnings(approx(power,sizes,xout=1-beta)$y)

# Report
if (is.finite(ssize)) {
  
cat(paste0("To have ",round(100*(1-beta)),"% power to detect an effect of MUSE of ",
           signif(muse_power_par["beta_2"],digits=ndig)," points on PSYRATS-voices total, ",
           "with other parameters as estimated, ",
           "we need a minimum of ",round (ssize), " samples per arm"))
  cat("\n\n\n")
} else {
  cat(paste0("With ",round(100*(1-beta)),"% power, we are unable to detect an effect of MUSE of ",
             signif(muse_power_par["beta_2"],digits=ndig)," points on PSYRATS-voices total ",
             "with ",max(sizes), " samples per arm with other parameters as estimated"))
  cat("\n\n\n")
}




##***********************************************##
## 10.2.2: Mixed linear model analysis PPO+cov ####
## Covariates included                           ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.2: Mixed linear model analysis PPO   **\n")
cat("**  Covariates included                      **\n")
cat("***********************************************\n")
cat("\n\n")



# Fit model
muse_ppo_covariates_lme=lmer(psyrats_voices_total ~ 
                    time2 + 
                    time3 + 
                    treatment:time2 + 
                    treatment:time3 + 
                    (1|id) + 
                    
                    factor(site) + 
                    age + 
                    factor(sex) + 
                    modality_total_excl_felt_presence + 
                    baseline_hallucination_duration +
                    time_engaged_in_service + 
                    psyrats_delusions_total_baseline,
                  data=muse_matrix)

# Point estimates of fixed and random effects
muse_ppo_covariates_fixed = summary(muse_ppo_covariates_lme)$coefficients
muse_ppo_covariates_random = as.data.frame(VarCorr(muse_ppo_covariates_lme))


# Confidence intervals
muse_ppo_covariates_ci = suppressMessages(confint(muse_ppo_covariates_lme))

# Summary table
muse_est_covariates=c(muse_ppo_covariates_random[c(which(muse_ppo_covariates_random$grp=="id"),
                             which(muse_ppo_covariates_random$grp=="Residual")),5],
           muse_ppo_covariates_fixed[c("(Intercept)", 
                                       "time2", 
                                       "time3", 
                                       "factor(site)2", 
                                       "age", 
                                       "factor(sex)M", 
                                       "modality_total_excl_felt_presence", 
                                       "baseline_hallucination_duration", 
                                       "time_engaged_in_service", 
                                       "psyrats_delusions_total_baseline", 
                                       "time2:treatment", 
                                       "time3:treatment"),"Estimate"])
muse_name_covariates=c("σr","σe","β0","βa2","βa3",
            "β(site)","β(age)","β(sex_M)",
            "β(mod_total)","β(base_dur)","β(time_eng)","β(delusions)",
            "β2","β3")
muse_description_covariates=c("SD of per-indivudal random effect",
                   "SD of residuals",
                   "Overall mean",
                   "Fixed effect for 8-week time point",
                   "Fixed effect for 12-week time point",
                   "Fixed effect for site 2 vs site 1",
                   "Fixed effect for age",
                   "Fixed effect for sex M vs F",
                   "Fixed effect for number of hallucinatory modalities, excluding sensed presence",
                   "Fixed effect for baseline duration of symptoms",
                   "Fixed effect for length of time engaged in service",
                   "Fixed effect for PSYRATs total delusion score at baseline",
                   "Fixed effect for treatment at 8-week time point",
                   "Fixed effect for treatment at 12-week time point")
out_tab_covariates=data.frame(signif(muse_est_covariates,digits=ndig),
                              signif(muse_ppo_covariates_ci,digits=ndig),
                              muse_description_covariates)
colnames(out_tab_covariates)=c("Estimate","CI95_lower","CI95_upper","Description")
rownames(out_tab_covariates)=muse_name_covariates

# Keep to use later
muse_psyrats_voices_total_covariates=out_tab_covariates


# Report
cat(paste0("Estimated beta_2 (mean effect of treatment on PSYRATS voices total score at 8 weeks) ",
           "including covariates in linear mixed model: ",
           signif(muse_ppo_covariates_fixed["time2:treatment","Estimate"],digits=ndig)," points; ",
           " 95% CI: (",signif(muse_ppo_covariates_ci["time2:treatment","2.5 %"],digits=ndig),", ",
           signif(muse_ppo_covariates_ci["time2:treatment","97.5 %"],digits=ndig),") points"))
cat("\n\n")

cat(paste0("All parameter estimates and 95% CIs, including covariates: ","\n\n"))
print(out_tab_covariates)
cat("\n\n")

ofile="./Analysis/Results/Tables/PSYRATs_voices_total/parameter_table_psyrats_voices_total_covariates.csv"
write.csv(out_tab_covariates,file=ofile,quote=FALSE)
cat(paste0("Parameter table saved to: ",ofile))
cat("\n\n\n")




##***********************************************##
## 10.2.3: Within-treatment analysis of number ####
##  of treatments                                ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.3: Within-treatment analysis of      **\n")
cat("**  number of treatments                     **\n")
cat("***********************************************\n")
cat("\n\n")



# IDs to include: treatment group, and psyrats_voices_total for all time points
incl_id=intersect(
         intersect(
          muse_summary[[1]]$id[which(muse_summary[[1]]$treatment==1)],
          muse_summary[[1]]$id[which(is.finite(muse_summary[[1]]$psyrats_voices_total))]
         ),intersect(
          muse_summary[[1]]$id[which(is.finite(muse_summary[[1]]$psyrats_voices_total))],
          muse_summary[[1]]$id[which(is.finite(muse_summary[[1]]$psyrats_voices_total))]
         )
        )

# Fit model
muse_ntreatment_lme=lmer(psyrats_voices_total ~ 
                           time2 + 
                           time3 + 
                           factor(site) + 
                           n_sessions_muse:time2 + 
                           n_sessions_muse:time3 + 
                           (1|id),
                         data=muse_matrix,
                         subset=which(muse_matrix$id %in% incl_id))

# Point estimates of fixed and random effects
muse_ntreatment_fixed = summary(muse_ntreatment_lme)$coefficients
muse_ntreatment_random = as.data.frame(VarCorr(muse_ntreatment_lme))


# Confidence intervals
muse_ntreatment_ci = suppressMessages(confint(muse_ntreatment_lme))

# Summary table
muse_est=c(muse_ntreatment_random[c(which(muse_ntreatment_random$grp=="id"),
                             which(muse_ntreatment_random$grp=="Residual")),5],
           muse_ntreatment_fixed[c("(Intercept)", 
                            "time2", 
                            "time3", 
                            "factor(site)2", 
                            "time2:n_sessions_muse", 
                            "time3:n_sessions_muse"),"Estimate"])
muse_name=c("σr","σe","β0","βa2","βa3","β(site)","β2","β3")
muse_description=c("SD of per-indivudal random effect",
                   "SD of residuals",
                   "Overall mean",
                   "Fixed effect for 8-week time point",
                   "Fixed effect for 12-week time point",
                   "Fixed effect for site 2 vs site 1",
                   "Fixed effect for number of treatments at 8-week time point",
                   "Fixed effect for number of treatments treatment at 12-week time point")
out_tab=data.frame(signif(muse_est,digits=ndig),signif(muse_ntreatment_ci,digits=ndig),muse_description)
colnames(out_tab)=c("Estimate","CI95_lower","CI95_upper","Description")
rownames(out_tab)=muse_name

# Keep to use later
muse_ntreatment_tab=out_tab


# Report
cat(paste0("Estimated beta_2 (mean effect of number of MUSE treatments on PSYRATS voices total ",
           "score at 8 weeks; incremental effect of one additional treament) ",
           signif(muse_ntreatment_fixed["time2:n_sessions_muse","Estimate"],digits=ndig)," points; ",
           " 95% CI: (",signif(muse_ntreatment_ci["time2:n_sessions_muse","2.5 %"],digits=ndig),", ",
           signif(muse_ntreatment_ci["time2:n_sessions_muse","97.5 %"],digits=ndig),") points"))
cat("\n\n")

cat(paste0("All parameter estimates and 95% CIs: ","\n\n"))
print(out_tab)
cat("\n\n")

ofile="./Analysis/Results/Tables/PSYRATs_voices_total/parameter_table_n_treatments_psyrats_voices_total.csv"
write.csv(out_tab,file=ofile,quote=FALSE)
cat(paste0("Parameter table saved to: ",ofile))
cat("\n\n\n")


##***********************************************##
## 10.2.4: Assessment of missingness,          ####
##  adherence, and other variables               ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.4: Assessment of missingness,        **\n")
cat("**  adherence, and other variables           **\n")
cat("***********************************************\n")
cat("\n\n")



# Do missingness rates in psyrats_voices_total differ at 8 weeks between
#  treated and untreated samples?
id_treated=muse_summary[[1]]$id[which(muse_summary[[1]]$treatment==1)]
id_nontreated=muse_summary[[1]]$id[which(muse_summary[[1]]$treatment==0)]
n_missing_treated=length(which(is.na(muse_summary[[2]]$psyrats_voices_total) & 
                                 (muse_summary[[2]]$id %in% id_treated)))
n_missing_nontreated=length(which(is.na(muse_summary[[2]]$psyrats_voices_total) & 
                                 (muse_summary[[2]]$id %in% id_nontreated)))
cont_matrix_treatment=rbind(c(n_missing_treated,length(id_treated)-n_missing_treated),
                  c(n_missing_nontreated,length(id_nontreated)-n_missing_nontreated))

# Fishers test against null that true odds ratio is 1:
p_missing_psyrats_voices_total_treatment=fisher.test(cont_matrix_treatment)

# Report
cat("\n\n")
cat(paste0("Assessment of missingness in psyrats_voices_total between treated and untreated samples","\n\n"))
cat(paste0("N missing amongst treated samples","\n"))
cat(paste0(n_missing_treated,"/",length(id_treated)))
cat("\n\n")
cat(paste0("N missing at site 2","\n"))
cat(paste0(n_missing_nontreated,"/",length(id_nontreated)))
cat("\n\n")
cat(paste0("Fisher test for difference in proportion against null hypothesis that odds ratio is 1","\n\n"))
print(p_missing_psyrats_voices_total_treatment)
cat("\n\n\n\n")



# Do missingness rates in psyrats_voices_total differ at 8 weeks between sites?
id_site1=muse_summary[[1]]$id[which(muse_summary[[1]]$site==1)]
id_site2=muse_summary[[1]]$id[which(muse_summary[[1]]$site==2)]
n_missing_site1=length(which(is.na(muse_summary[[2]]$psyrats_voices_total) & 
                                 (muse_summary[[2]]$id %in% id_site1)))
n_missing_site2=length(which(is.na(muse_summary[[2]]$psyrats_voices_total) & 
                                    (muse_summary[[2]]$id %in% id_site2)))
cont_matrix_site=rbind(c(n_missing_site1,length(id_site1)-n_missing_site1),
                  c(n_missing_site2,length(id_site2)-n_missing_site2))

# Fishers test against null that true odds ratio is 1:
p_missing_psyrats_voices_total_site=fisher.test(cont_matrix_site)

# Report
cat("\n\n")
cat(paste0("Assessment of missingness in psyrats_voices_total between sites","\n\n"))
cat(paste0("N missing at site 1","\n"))
cat(paste0(n_missing_site1,"/",length(id_site1)))
cat("\n\n")
cat(paste0("N missing at site 2","\n"))
cat(paste0(n_missing_site2,"/",length(id_site2)))
cat("\n\n")
cat(paste0("Fisher test for difference in proportion against null hypothesis that odds ratio is 1","\n\n"))
print(p_missing_psyrats_voices_total_site)
cat("\n\n\n\n")



# Do baseline psyrats_voices_total values differ between missing and non-missing samples?
id_missing=muse_summary[[2]]$id[which(is.na(muse_summary[[2]]$psyrats_voices_total))]
id_nonmissing=muse_summary[[2]]$id[which(!is.na(muse_summary[[2]]$psyrats_voices_total))]
psyrats_voices_total_baseline=muse_summary[[1]]$psyrats_voices_total

pvt_missing=psyrats_voices_total_baseline[which(muse_summary[[1]]$id %in% id_missing)]
pvt_nonmissing=psyrats_voices_total_baseline[which(muse_summary[[1]]$id %in% id_nonmissing)]

# Wilcoxon test against equality of medians
test_psyrats_voices_total_missing=wilcox.test(pvt_missing,pvt_nonmissing)

cat("\n\n")
cat(paste0("Assessment of difference in baseline psyrats_voices_total between missing and non-missing samples","\n\n"))
cat(paste0("Mean (SD) baseline psyrats_voices_total amongst non-missing samples","\n"))
cat(paste0(signif(mean(pvt_nonmissing),digits=ndig)," (",signif(sd(pvt_nonmissing),digits=ndig),")"))
cat("\n\n")
cat(paste0("Mean baseline psyrats_voices_total amongst missing samples","\n"))
cat(paste0(signif(mean(pvt_missing),digits=ndig)))
cat("\n\n")
cat(paste0("Wilcoxon rank-sum test against null hypothesis of equal medians\n\n"))
print(test_psyrats_voices_total_missing)
cat("\n\n\n\n")



# Do baseline HPSVQ_total values differ between missing and non-missing samples?
id_missing=muse_summary[[2]]$id[which(is.na(muse_summary[[2]]$psyrats_voices_total))]
id_nonmissing=muse_summary[[2]]$id[which(!is.na(muse_summary[[2]]$psyrats_voices_total))]
hpsvq_baseline=muse_summary[[1]]$hpsvq_total

ht_missing=hpsvq_baseline[which(muse_summary[[1]]$id %in% id_missing)]
ht_nonmissing=hpsvq_baseline[which(muse_summary[[1]]$id %in% id_nonmissing)]

# Wilcoxon test against equality of medians
test_hpsvq_total_missing=wilcox.test(ht_missing,ht_nonmissing)

cat("\n\n")
cat(paste0("Assessment of difference in baseline HPSVQ_total between missing and non-missing samples","\n\n"))
cat(paste0("Mean (SD) baseline HPSVQ_total amongst non-missing samples","\n"))
cat(paste0(signif(mean(ht_nonmissing),digits=ndig)," (",signif(sd(ht_nonmissing),digits=ndig),")"))
cat("\n\n")
cat(paste0("Mean baseline HPSVQ_total amongst missing samples","\n"))
cat(paste0(signif(mean(ht_missing),digits=ndig)))
cat("\n\n")
cat(paste0("Wilcoxon rank-sum test against null hypothesis of equal medians\n\n"))
print(test_hpsvq_total_missing)
cat("\n\n\n\n")


##***********************************************##
## 10.2.4b: Repeat analysis 10.2.1 on          ####
##  per-protocol population                      ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.4b: Repeat analysis 10.2.1 on        **\n")
cat("**  per-protocol population                  **\n")
cat("***********************************************\n")
cat("\n\n")



# IDs to include:  and psyrats_voices_total for all time points
incl_id_perprotocol=intersect(
  muse_summary[[1]]$id[which(is.finite(muse_summary[[1]]$psyrats_voices_total))],
  intersect(
    muse_summary[[1]]$id[which(is.finite(muse_summary[[1]]$psyrats_voices_total))],
    muse_summary[[1]]$id[which(is.finite(muse_summary[[1]]$psyrats_voices_total))]
  )
)


# Fit model
muse_ppo_perprotocol_lme=lmer(psyrats_voices_total ~ 
                    time2 + 
                    time3 + 
                    treatment:time2 + 
                    treatment:time3 + 
                    factor(site) + 
                    (1|id),
                  data=muse_matrix,
                  subset=which(muse_matrix$id %in% incl_id_perprotocol))

# Point estimates of fixed and random effects
muse_ppo_fixed_perprotocol = summary(muse_ppo_perprotocol_lme)$coefficients
muse_ppo_random_perprotocol = as.data.frame(VarCorr(muse_ppo_perprotocol_lme))


# Confidence intervals
muse_ppo_perprotocol_ci = suppressMessages(confint(muse_ppo_perprotocol_lme))

# Summary table
muse_est_perprotocol=c(muse_ppo_random_perprotocol[c(which(muse_ppo_random_perprotocol$grp=="id"),
                             which(muse_ppo_random_perprotocol$grp=="Residual")),5],
           muse_ppo_fixed_perprotocol[c("(Intercept)", 
                            "time2", 
                            "time3", 
                            "factor(site)2", 
                            "time2:treatment", 
                            "time3:treatment"),"Estimate"])
muse_name=c("σr","σe","β0","βa2","βa3","β(site)","β2","β3")
muse_description=c("SD of per-indivudal random effect",
                   "SD of residuals",
                   "Overall mean",
                   "Fixed effect for 8-week time point",
                   "Fixed effect for 12-week time point",
                   "Fixed effect for site 2 vs site 1",
                   "Fixed effect for treatment at 8-week time point",
                   "Fixed effect for treatment at 12-week time point")
out_tab=data.frame(signif(muse_est_perprotocol,digits=ndig),signif(muse_ppo_perprotocol_ci,digits=ndig),muse_description)
colnames(out_tab)=c("Estimate","CI95_lower","CI95_upper","Description")
rownames(out_tab)=muse_name

# Keep to use later
muse_psyrats_voices_total_perprotocol=out_tab


# Report
cat(paste0("Estimated beta_2 (mean effect of treatment on PSYRATS voices total score at 8 weeks) ",
           "in per-protocol population: ",
           signif(muse_ppo_fixed_perprotocol["time2:treatment","Estimate"],digits=ndig)," points; ",
           " 95% CI: (",signif(muse_ppo_perprotocol_ci["time2:treatment","2.5 %"],digits=ndig),", ",
           signif(muse_ppo_perprotocol_ci["time2:treatment","97.5 %"],digits=ndig),") points"))
cat("\n\n")

cat(paste0("All parameter estimates and 95% CIs (per-protocol): ","\n\n"))
print(out_tab)
cat("\n\n")

ofile="./Analysis/Results/Tables/PSYRATs_voices_total/parameter_table_psyrats_voices_total_perprotocol.csv"
write.csv(out_tab,file=ofile,quote=FALSE)
cat(paste0("Parameter table saved to: ",ofile))
cat("\n\n\n")


##***********************************************##
## 10.2.4c: Effect of COVID                    ####
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.4c: Effect of COVID                  **\n")
cat("***********************************************\n")
cat("\n\n")


# Indices of samples randomised pre- and post- COVID in muse_summary[[1]]
ind_postcovid=which(muse_summary[[1]]$covid_affected==1)
ind_precovid=which(muse_summary[[1]]$covid_affected==0)

# Number of sessions
nmuse=muse_summary[[1]]$n_sessions_muse
  
# Test 
test_n_muse_sessions_covid=suppressWarnings(wilcox.test(
  nmuse[ind_precovid],
  nmuse[ind_postcovid]
))

cat("\n\n")
cat(paste0("Difference in number of MUSE sessions for pre- and post-covid samples","\n\n"))
cat(paste0("Mean (SD) number of sessions for pre-covid samples","\n"))
cat(paste0(signif(mean(nmuse[ind_precovid],na.rm=TRUE),digits=ndig),
           " (",signif(sd(nmuse[ind_precovid],na.rm=TRUE),digits=ndig),")"))
cat("\n\n")
cat(paste0("Mean (SD) number of sessions for post-covid samples","\n"))
cat(paste0(signif(mean(nmuse[ind_postcovid],na.rm=TRUE),digits=ndig),
           " (",signif(sd(nmuse[ind_postcovid],na.rm=TRUE),digits=ndig),")"))
cat("\n\n")
cat(paste0("Wilcoxon rank-sum test against null hypothesis of equal medians\n\n"))
print(test_n_muse_sessions_covid)
cat("\n\n\n\n")


covid_outcomes=c("psyrats_voices_total",
                 "psyrats_hdis",
                 "hpsvq_total",
                 "qpr_total")

for (i in 1:length(covid_outcomes)) {
  c_out=covid_outcomes[[i]]
  
  Y=muse_summary[[1]][[c_out]]
  
  test_Y=suppressWarnings(wilcox.test(
    Y[ind_precovid],
    Y[ind_postcovid]
  ))
  
  cat("\n\n")
  cat(paste0("Difference in ",c_out," for pre- and post-covid samples","\n\n"))
  cat(paste0("Mean (SD) ",c_out," for pre-covid samples","\n"))
  cat(paste0(signif(mean(Y[ind_precovid]),digits=ndig),
             " (",signif(sd(Y[ind_precovid]),digits=ndig),")"))
  cat("\n\n")
  cat(paste0("Mean (SD) ",c_out," for post-covid samples","\n"))
  cat(paste0(signif(mean(Y[ind_postcovid]),digits=ndig),
             " (",signif(sd(Y[ind_postcovid]),digits=ndig),")"))
  cat("\n\n")
  cat(paste0("Wilcoxon rank-sum test against null hypothesis of equal medians\n\n"))
  print(test_Y)
  cat("\n\n\n\n")
  
}


##***********************************************##
## 10.2.5: Mixed linear model analysis for     ####
## other outcomes, no covariates other than      ##
## site                                          ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.5: Mixed linear model analysis for   **\n")
cat("**  other outcomes, no covariates other than **\n")
cat("**  site                                     **\n")
cat("***********************************************\n")
cat("\n\n")



## Most of data frame will remain the same as for analysis of PPO

for (o in 1:length(other_outcomes)) {
  
  # Specify outcome
  outcome=other_outcomes[o]
  
  # Fit model
  xform=as.formula(paste0(outcome," ~ 
                          time2 + 
                          time3 + 
                          treatment:time2 + 
                          treatment:time3 + 
                          factor(site) + 
                          (1|id)"))
  muse_outcome_lme=lmer(xform,data=muse_matrix)
  
  # Point estimates of fixed and random effects
  muse_outcome_fixed = summary(muse_outcome_lme)$coefficients
  muse_outcome_random = as.data.frame(VarCorr(muse_outcome_lme))
  
  
  # Confidence intervals
  muse_outcome_ci = suppressMessages(confint(muse_outcome_lme))
  
  # Summary table
  muse_est=c(muse_outcome_random[c(which(muse_outcome_random$grp=="id"),
                               which(muse_outcome_random$grp=="Residual")),5],
             muse_outcome_fixed[c("(Intercept)", 
                              "time2", 
                              "time3", 
                              "factor(site)2", 
                              "time2:treatment", 
                              "time3:treatment"),"Estimate"])
  out_tab=data.frame(signif(muse_est,digits=ndig),signif(muse_outcome_ci,digits=ndig),muse_description)
  colnames(out_tab)=c("Estimate","CI95_lower","CI95_upper","Description")
  rownames(out_tab)=muse_name
  
  
  # Report
  cat(paste0("**** Analysis for outcome: ",outcome," ****"))
  cat("\n\n")
  
  cat(paste0("Estimated beta_2 (mean effect of treatment on outcome ",outcome," at 8 weeks) ",
             signif(muse_outcome_fixed["time2:treatment","Estimate"],digits=ndig)," points; ",
             " 95% CI: (",signif(muse_outcome_ci["time2:treatment","2.5 %"],digits=ndig),", ",
             signif(muse_outcome_ci["time2:treatment","97.5 %"],digits=ndig),") points"))
  cat("\n\n")
  
  cat(paste0("All parameter estimates and 95% CIs: ","\n\n"))
  print(out_tab)
  cat("\n")
  
  ofile=paste0("./Analysis/Results/Tables/Other_outcomes/parameter_table_",outcome,".csv")
  write.csv(out_tab,file=ofile,quote=FALSE)
  cat(paste0("Parameter table saved to: ",ofile))
  
  cat("\n\n")
  cat(paste0("*********************"))
  cat("\n\n")
  
  # Keep table to use later
  assign(paste0("muse_",outcome),out_tab)
  
}

##*****************************************##
## 10.2.5b: Draw figures                 ####
##*****************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.5b: Draw figures                     **\n")
cat("***********************************************\n")
cat("\n\n")




# Draw table/plot of combined effects
for (wk in c(8,12)) {
  if (wk==8) w_ind="β2" else w_ind="β3"
  
  ptab=muse_psyrats_voices_total[w_ind,1:3]
  for (i in 1:length(other_outcomes)) {
    ptab=rbind(ptab,get(paste0("muse_",other_outcomes[i]))[w_ind,1:3])
  }
  rownames(ptab)=c("psyrats_voices_total",other_outcomes)
  
  # Manage variables for which higher=better
  w_swap=which(rownames(ptab) %in% swap)
  ptab[w_swap,]=-ptab[w_swap,]
  
  xtab=data.frame(Outcome=longnames[rownames(ptab)],
                  Value = paste0(signif(ptab$Estimate,digits=ndig),
                                 " (",signif(ptab$CI95_lower,digits=ndig),", ",
                                 signif(ptab$CI95_upper,digits=ndig),")"))
  xtab[w_swap,1]=paste0(xtab[w_swap,1],"*")
  
  fig_file=paste0(save_dir,"forest_plot_week_",wk,".pdf")
  pdf(fig_file,width=6,height=3)
  par(mar=c(3,1,1,1)); tcex=0.5; rsize=0.2
  xmin=min(ptab[,2:3]); xmax=max(ptab[,2:3])
  yc=dim(xtab)[1]:1
  plot(0,type="n",bty="n",ann=FALSE,xaxt="n",yaxt="n",
       xlim=c(2*xmin-xmax,xmax),ylim=c(0,dim(ptab)[1]))
  text(xmin - 1.1*(xmax-xmin)/2,yc,xtab[,1],adj=c(1,0.5),cex=tcex)
  text(xmin - 0.9*(xmax-xmin)/2,yc,xtab[,2],adj=c(0,0.5),cex=tcex)
  axis(1,at=pretty(range(ptab[,2:3])))
  abline(v=0,lwd=2)
  segments(ptab[,2],yc,ptab[,3],yc,lty=2)
  points(ptab[,1],yc,col="red",pch=16,cex=1.5)
  dev.off()
  
  cat(paste0("\n\nFigure for changes at ",wk," weeks saved to: ",fig_file,"\n\n"))
  
}


##***********************************************##
## 10.2.5c: Tables of mean values              ####
##***********************************************##


# Loop across outcomes
ppos=c("psyrats_voices_total",other_outcomes)
for (i in 1:length(ppos)) {
  dt=data.frame(
    id=muse_summary[[1]]$id,
    treat=muse_summary[[1]]$treatment,
    base=muse_summary[[1]][,ppos[i]],
    week8=muse_summary[[2]][match(muse_summary[[1]]$id,muse_summary[[2]]$id),ppos[i]],
    week12=muse_summary[[3]][match(muse_summary[[1]]$id,muse_summary[[3]]$id),ppos[i]])
  dt1=dt[which(dt$treat==0),]; dt2=dt[which(dt$treat==1),]
  tab0=c(
    length(which(is.finite(dt1$base))),
      msd(dt1$base),
      "-",
      length(which(is.finite(dt2$base))),
      msd(dt2$base),
      "-",
      length(which(is.finite(c(dt1$base,dt2$base)))),
      msdd(dt2$base,dt1$base),
      "-")
  d1=dt1$week8-dt1$base; d2=dt2$week8-dt2$base
  d1=d1[which(is.finite(d1))]; d2=d2[which(is.finite(d2))]
  tab1=c(
    length(which(is.finite(dt1$week8))),
      msd(dt1$week8),
      muse_format(d1,dig=ndig),
      length(which(is.finite(dt2$week8))),
      msd(dt2$week8),
      muse_format(d2,dig=ndig),
      length(which(is.finite(c(dt1$week8,dt2$week8)))),
      msdd(dt2$week8,dt1$week8),
      difci(d2,d1,length(d2),length(d1))
    )
  d1=dt1$week12-dt1$base; d2=dt2$week12-dt2$base
  d1=d1[which(is.finite(d1))]; d2=d2[which(is.finite(d2))]
  tab2=c(
    length(which(is.finite(dt1$week12))),
    msd(dt1$week12),
    muse_format(d1,dig=ndig),
    length(which(is.finite(dt2$week12))),
    msd(dt2$week12),
    muse_format(d2,dig=ndig),
    length(which(is.finite(c(dt1$week12,dt2$week12)))),
    msdd(dt2$week12,dt1$week12),
    difci(d2,d1,length(d2),length(d1))
  )
  tab=rbind(tab0,tab1,tab2)
  rownames(tab)=c("Baseline","Week 8","Week 12")
  colnames(tab)=c("TAU: number of usable outcome scores",
                  "TAU: outcome mean (SD)",
                  "TAU: change from baseline (95% CI)",
                  "MUSE: number of usable outcome scores",
                  "MUSE: outcome mean (SD)",
                  "MUSE: change from baseline (95% CI)",
                  "MUSE-TAU: number of usable outcome scores",
                  "MUSE-TAU: outcome mean (SD)",
                  "MUSE-TAU: change from baseline (95% CI)"
                  )
  
  cat(paste0("Table for outcome: ",ppos[i],": \n\n"))
  print(tab)
  cat("\n\n")
  
  if (i==1) ofile=paste0("./Analysis/Results/Tables/PSYRATs_voices_total/difference_table",ppos[i],".csv")
  if (i>1) ofile=paste0("./Analysis/Results/Tables/Other_outcomes/differences_table_",ppos[i],".csv")
  write.csv(tab,file=ofile,quote=FALSE)
  cat(paste0("Parameter table saved to: ",ofile))
  cat("\n\n")
  
}



##***********************************************##
## 10.2.6: Mixed linear model analysis for     ####
## other outcomes, covariates included           ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.6: Mixed linear model analysis for   **\n")
cat("**  other outcomes, covariates included      **\n")
cat("***********************************************\n")
cat("\n\n")



## Most of data frame will remain the same as for analysis of PPO

for (o in 1:length(other_outcomes)) {
  
  # Specify outcome
  outcome=other_outcomes[o]
  
  if (!(outcome=="psyrats_delusions_total")) {
  
  # Fit model
  xform=as.formula(paste0(outcome," ~ 
                          time2 + 
                          time3 + 
                          treatment:time2 + 
                          treatment:time3 + 
                          (1|id) +
                    factor(site) + 
                    age + 
                    factor(sex) + 
                    modality_total_excl_felt_presence + 
                    baseline_hallucination_duration +
                    time_engaged_in_service + 
                    psyrats_delusions_total_baseline"))
  muse_outcome_lme=lmer(xform,data=muse_matrix)
  
  # Point estimates of fixed and random effects
  muse_outcome_fixed = summary(muse_outcome_lme)$coefficients
  muse_outcome_random = as.data.frame(VarCorr(muse_outcome_lme))
  
    # Confidence intervals
  muse_outcome_ci = suppressMessages(confint(muse_outcome_lme))
  
  # Summary table
  muse_est_covariates=c(muse_outcome_random[c(which(muse_outcome_random$grp=="id"),
                                              which(muse_outcome_random$grp=="Residual")),5],
                        muse_outcome_fixed[c("(Intercept)", 
                                             "time2", 
                                             "time3", 
                                             "factor(site)2", 
                                             "age", 
                                             "factor(sex)M", 
                                             "modality_total_excl_felt_presence", 
                                             "baseline_hallucination_duration", 
                                             "time_engaged_in_service", 
                                             "psyrats_delusions_total_baseline", 
                                             "time2:treatment", 
                                             "time3:treatment"),"Estimate"])
  muse_name_covariates=c("σr","σe","β0","βa2","βa3",
                         "β(site)","β(age)","β(sex_M)",
                         "β(mod_total)","β(base_dur)","β(time_eng)","β(delusions)",
                         "β2","β3")
  muse_description_covariates=c("SD of per-indivudal random effect",
                                "SD of residuals",
                                "Overall mean",
                                "Fixed effect for 8-week time point",
                                "Fixed effect for 12-week time point",
                                "Fixed effect for site 2 vs site 1",
                                "Fixed effect for age",
                                "Fixed effect for sex M vs F",
                                "Fixed effect for number of hallucinatory modalities, excluding sensed presence",
                                "Fixed effect for baseline duration of symptoms",
                                "Fixed effect for length of time engaged in service",
                                "Fixed effect for PSYRATs total delusion score at baseline",
                                "Fixed effect for treatment at 8-week time point",
                                "Fixed effect for treatment at 12-week time point")
  out_tab_covariates=data.frame(signif(muse_est_covariates,digits=ndig),
                                signif(muse_outcome_ci,digits=ndig),
                                muse_description_covariates)
  colnames(out_tab_covariates)=c("Estimate","CI95_lower","CI95_upper","Description")
  rownames(out_tab_covariates)=muse_name_covariates
  
  } else {
    
    
    # Fit model
    xform=as.formula(paste0(outcome," ~ 
                          time2 + 
                          time3 + 
                          treatment:time2 + 
                          treatment:time3 + 
                          (1|id) +
                    factor(site) + 
                    age + 
                    factor(sex) + 
                    modality_total_excl_felt_presence + 
                    baseline_hallucination_duration +
                    time_engaged_in_service"))
    muse_outcome_lme=lmer(xform,data=muse_matrix)
    
    # Point estimates of fixed and random effects
    muse_outcome_fixed = summary(muse_outcome_lme)$coefficients
    muse_outcome_random = as.data.frame(VarCorr(muse_outcome_lme))
    
    # Confidence intervals
    muse_outcome_ci = suppressMessages(confint(muse_outcome_lme))
    
    # Summary table
    muse_est_covariates=c(muse_outcome_random[c(which(muse_outcome_random$grp=="id"),
                                                which(muse_outcome_random$grp=="Residual")),5],
                          muse_outcome_fixed[c("(Intercept)", 
                                               "time2", 
                                               "time3", 
                                               "factor(site)2", 
                                               "age", 
                                               "factor(sex)M", 
                                               "modality_total_excl_felt_presence", 
                                               "baseline_hallucination_duration", 
                                               "time_engaged_in_service",  
                                               "time2:treatment", 
                                               "time3:treatment"),"Estimate"])
    muse_name_covariates=c("σr","σe","β0","βa2","βa3",
                           "β(site)","β(age)","β(sex_M)",
                           "β(mod_total)","β(base_dur)","β(time_eng)",
                           "β2","β3")
    muse_description_covariates=c("SD of per-indivudal random effect",
                                  "SD of residuals",
                                  "Overall mean",
                                  "Fixed effect for 8-week time point",
                                  "Fixed effect for 12-week time point",
                                  "Fixed effect for site 2 vs site 1",
                                  "Fixed effect for age",
                                  "Fixed effect for sex M vs F",
                                  "Fixed effect for number of hallucinatory modalities, excluding sensed presence",
                                  "Fixed effect for baseline duration of symptoms",
                                  "Fixed effect for length of time engaged in service",
                                  "Fixed effect for treatment at 8-week time point",
                                  "Fixed effect for treatment at 12-week time point")
    out_tab_covariates=data.frame(signif(muse_est_covariates,digits=ndig),
                                  signif(muse_outcome_ci,digits=ndig),
                                  muse_description_covariates)
    colnames(out_tab_covariates)=c("Estimate","CI95_lower","CI95_upper","Description")
    rownames(out_tab_covariates)=muse_name_covariates
    
  }
  
  
  
  # Report
  cat(paste0("**** Analysis for outcome: ",outcome,", covariates included ****"))
  cat("\n\n")
  
  cat(paste0("Estimated beta_2 (mean effect of treatment on outcome ",outcome," at 8 weeks) ",
             signif(muse_outcome_fixed["time2:treatment","Estimate"],digits=ndig)," points; ",
             " 95% CI: (",signif(muse_outcome_ci["time2:treatment","2.5 %"],digits=ndig),", ",
             signif(muse_outcome_ci["time2:treatment","97.5 %"],digits=ndig),") points"))
  cat("\n\n")
  
  cat(paste0("All parameter estimates and 95% CIs: ","\n\n"))
  print(out_tab_covariates)
  cat("\n")
  
  ofile=paste0("./Analysis/Results/Tables/Other_outcomes/parameter_table_",outcome,"_covariates.csv")
  write.csv(out_tab_covariates,file=ofile,quote=FALSE)
  cat(paste0("Parameter table saved to: ",ofile))
  
  cat("\n\n")
  cat(paste0("*********************"))
  cat("\n\n")
  
  # Keep table to use later
  assign(paste0("muse_",outcome,"_covariates"),out_tab_covariates)
  
}


##***********************************************##
## 10.2.7: Assessment of potential primary       ##
##  outcomes                                     ##
##***********************************************##

cat("\n\n")
cat("***********************************************\n")
cat("** 10.2.7: Assessment of potential primary   **\n")
cat("**  outcomes                                 **\n")
cat("***********************************************\n")
cat("\n\n")


adherence_tab=c()
ppos=c("psyrats_voices_total",other_outcomes)

for (i in 1:length(ppos)) {
  outcome=ppos[i]
  vec=c(
    length(which(is.na(muse_summary[[1]][[outcome]]))),
    length(which(is.na(muse_summary[[2]][[outcome]]))),
    length(which(is.na(muse_summary[[3]][[outcome]])))
  )
  adherence_tab=rbind(adherence_tab,vec)
}
rownames(adherence_tab)=ppos
colnames(adherence_tab)=c("Baseline","8week","12week")

cat(paste0("Adherence to potential primary outcomes: number of missing values","\n\n"))
print(adherence_tab)
cat("\n")

ofile=paste0("./Analysis/Results/Tables/Comparison/adherence.csv")
write.csv(out_tab_covariates,file=ofile,quote=FALSE)
cat(paste0("Table saved to: ",ofile))
cat("\n\n\n")


# Rank correlation between PPOs, NAs removed
times=c("baseline","8weeks", "12weeks")
for (t in 1:length(times)) {
  subtab=muse_summary[[t]][ppos]
  cor_tab=cor(subtab,use="pairwise.complete.obs",method="spearman")
  
  cat(paste0("Spearman (rank) correlation between potential primary outcomes at time ",times[t],": \n\n"))
  print(cor_tab)
  cat("\n")
  
  ofile=paste0("./Analysis/Results/Tables/Comparison/correlation_",times[t],".csv")
  write.csv(cor_tab,file=ofile,quote=FALSE)
  cat(paste0("Table saved to: ",ofile))
  cat("\n\n\n")
  
}


cat("\n\n")
cat("***********************************************\n")
cat("** Main analysis complete                    **\n")
cat("***********************************************\n")
cat("\n\n")



##*****************************************##
## Check confidence interval cover       ####
##*****************************************##



cat("\n\n\n")
cat(paste0("Checking empirical cover of confidence intervals"))
cat("\n\n")

size=50 # Simulate for this many samples per group
ncheck=2000 # Simulate this many times


conf_int_file=paste0("./Analysis/Results/R_objects/conf_int_cover_simulation.RData")
if (recheck_conf_int | (!file.exists(conf_int_file))) {
  
  set.seed(38729)
  
  
  # Parameters for simulation; correspond to design in SAP
  beta_0=2;
  beta_a1=0
  beta_a2=2.5
  beta_a3=-0.5
  beta_2=3
  beta_3=1
  beta=-1
  sigma_r=2
  sigma_e=1
  
  
  # Parameters in order returned by confint()
  par_check=c(sigma_r,sigma_e,beta_0,beta_a2,beta_a3,beta_2,beta_3,beta)
  names(par_check)=c("sigma_r","sigma_e","beta_0","beta_a2","beta_a3","beta_2","beta_3","beta")
  
  
  cat(paste0("Computing empirical confidence intervals for study size ",size," per group"))
  cat("\n\n")
  
  rec=matrix(0,8,ncheck)
  
  for (i in 1:ncheck) {
    dat=ngen(size,par=c(beta_0,beta_a1,beta_a2,beta_a3,beta_2,beta_3,beta,sigma_r,sigma_e))
    m1=lmer(Y~time + treat_time2 + treat_time3 + site + (1|id),data=dat)
    c1=suppressMessages(confint(m1))
    rec[,i]=(c1[,1]<par_check & c1[,2]>par_check)
    if ((i%%100)==0) cat(paste0("Completed ",i," of ",ncheck,"\n"))
  }
  
save(rec,par_check,file=conf_int_file)
  
} else load(conf_int_file)
  
cat("\n\n")
  
coverage=rowMeans(rec) # Empirical cover of each confidence interval
p_cover=pbinom(coverage*ncheck,ncheck,1-alpha) 

out=data.frame(cover=coverage,quantile=p_cover); rownames(out)=names(par_check)

cat(paste0("Confidence interval cover (nominally 95%) and quantile in Bin(n,0.95):","\n\n"))

prmatrix(out)





##*****************************************##
## Clean up                              ####
##*****************************************##

if (!is.null(text_dir)) sink()
