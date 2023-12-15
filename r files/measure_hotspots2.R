###
#Description of the measure statistics
#a=units inside both true and detected clusters
#b=units inside true and outside detected
#c=units outside true and inside detected
#d=units outside true and outside detected

rm(list=ls())

user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
setwd(paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov"))

# Create file names for all replicate folders
file_names = list.files(path = paste0(getwd()), pattern = "rep*", full.names = FALSE)
file_names=file_names[file_names!="measure_hotspots.R" & file_names!="rep0" & file_names!="repold"
                      & file_names!="replicates.txt" & file_names!="measure_hotspots2.R"]
file_names=file_names[file_names!="rep33" & file_names!="rep63" & file_names!="rep34"
                      & file_names!="rep62" & file_names!="rep93" & file_names!="rep32"]
length(file_names)
file_names

#write.table(file_names, paste0(getwd(), "/file_names.txt"), row.names=FALSE, col.names=FALSE)

type_hotspot = c("single", "multiple","decay", "river", "short", "long", "season", "dry")

############################
for (j in type_hotspot){
  
for (i in 1:length(file_names)){
  
  # User 
  user = strsplit(getwd(), "/", fixed = FALSE, perl = FALSE, useBytes = FALSE)[[1]][5]
  
  # Working directory
  setwd(paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov"))
  
  experiment <- paste0(file_names[i]);experiment
  sed <- 42
  set.seed(sed)
  
  setwd(paste0(getwd(),"/",experiment))
  getwd()
  
  
hotnew<-data.frame(read.table(paste0("hotspots_", j, ".txt"), header = T))
hotnew$kullT03 <- as.numeric(hotnew$kullT03)


##kull
hotnew$kullT03_a[hotnew$hotspot==1 & hotnew$kullT03==1]<-1
hotnew$kullT03_b[hotnew$hotspot==1 & hotnew$kullT03==0]<-1
hotnew$kullT03_c[hotnew$hotspot==0 & hotnew$kullT03==1]<-1
hotnew$kullT03_d[hotnew$hotspot==0 & hotnew$kullT03==0]<-1

hotnew$kullT_a[hotnew$hotspot==1 & hotnew$kullT==1]<-1
hotnew$kullT_b[hotnew$hotspot==1 & hotnew$kullT==0]<-1
hotnew$kullT_c[hotnew$hotspot==0 & hotnew$kullT==1]<-1
hotnew$kullT_d[hotnew$hotspot==0 & hotnew$kullT==0]<-1

hotnew$kullT3y_a[hotnew$hotspot==1 & hotnew$kullT3y==1]<-1
hotnew$kullT3y_b[hotnew$hotspot==1 & hotnew$kullT3y==0]<-1
hotnew$kullT3y_c[hotnew$hotspot==0 & hotnew$kullT3y==1]<-1
hotnew$kullT3y_d[hotnew$hotspot==0 & hotnew$kullT3y==0]<-1

hotnew[is.na(hotnew)]<-0

#misclassification (proportion of mistakenly identified population)
#inside true and outside detect
#outside true and inside detect
#kullT03_Miss<-100 * (sum(hotnew$kullT03_b+hotnew$kullT03_c)/sum(hotnew$kullT03_a+hotnew$kullT03_b+hotnew$kullT03_c+hotnew$kullT03_d))

#sensitivity (proportion of true clusters correctly identified)
kullT03_Sens<-100 * sum(hotnew$kullT03_a)/sum(hotnew$kullT03_a+hotnew$kullT03_b)
#specificity (proportion outside true clusters correctly identified)
kullT03_Spec<-100 * sum(hotnew$kullT03_d)/sum(hotnew$kullT03_c+hotnew$kullT03_d)
#positive predictive values (prop of pop in detected clusters which are true)
kullT03_Ppv<-100 * sum(hotnew$kullT03_a)/sum(hotnew$kullT03_a+hotnew$kullT03_c)

kullT_Sens<-100 * sum(hotnew$kullT_a)/sum(hotnew$kullT_a+hotnew$kullT_b)
kullT_Spec<-100 * sum(hotnew$kullT_d)/sum(hotnew$kullT_c+hotnew$kullT_d)
kullT_Ppv<-100 * sum(hotnew$kullT_a)/sum(hotnew$kullT_a+hotnew$kullT_c)

kullT3y_Sens<-100 * sum(hotnew$kullT3y_a)/sum(hotnew$kullT3y_a+hotnew$kullT3y_b)
kullT3y_Spec<-100 * sum(hotnew$kullT3y_d)/sum(hotnew$kullT3y_c+hotnew$kullT3y_d)
kullT3y_Ppv<-100 * sum(hotnew$kullT3y_a)/sum(hotnew$kullT3y_a+hotnew$kullT3y_c)

results1<-rbind(kullT03_Sens,kullT03_Spec,kullT_Sens,kullT_Spec,kullT3y_Sens,kullT3y_Spec)
colnames(results1)<-experiment

#write.csv(results,"results_si.csv", row.names = T)

#Read file from previous replicate
if (i==1){
  m <- i-1
  exp_old <- paste0("rep",m)
  exp_old
  results2 <- read.csv(paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov/",exp_old,"/results_si.csv"), header = T)
}

if (i>1){
m <- i-1
exp_old <- paste0(file_names[m])
exp_old
results2 <- read.csv(paste0("/scicore/home/penny/",user,"/malaria_hotspots/simulations_nov/",exp_old,"/r_", j, ".csv"), header = T)
}

results3<-data.frame(results2, results1[,1])
write.csv(results3,paste0(getwd(),"/r_", j, ".csv"), row.names = T)

  }
}
  
getwd()

if (i==length(file_names)){
  for (j in type_hotspot){
    # #to calcutate the confidence intervals
    results3 <- read.csv(paste0(getwd(),"/r_", j, ".csv"), header = T)
    
    results3 <- results3[,-c(1:m)]
    results3 <- results3[,-c(2,3)]
    
    results3$median = round(apply(results3[,-c(1)], 1, median), 2)
    results3$LCI = round(apply(results3[,-c(1)], 1, quantile, probs = c(0.25),  na.rm = TRUE), 2)
    results3$UCI = round(apply(results3[,-c(1)], 1, quantile, probs = c(0.75),  na.rm = TRUE), 2)
    
    results3$Median_IQR<-paste(results3$median,"(",results3$LCI,",",results3$UCI,")")
    names(results3)[1] <- "ScanStat"
    
    results4<-results3[,c("ScanStat","median","LCI","UCI", "Median_IQR")]
    
    write.csv(results4,paste0("results_", j, ".csv"), row.names = F)
  }
 }



# save a file with all types of hotspots

# results4$type = "single"
# results5<-rbind(results5,results4)
# write.csv(results5,"results_alltypes.csv", row.names = F)
# 
# 
# # to plot sensitivity and specificity results 
# 
# res1 <- results5
# res1$ScanStat <- as.factor(res1$ScanStat)
# res1$type <- factor(res1$type, levels = c("single", "multiple", "decay", "river","short","long", "season", "dry"))
# 
# res2<-res1[res1$ScanStat=="kullT03_Sens" | res1$ScanStat=="kullT_Sens" | res1$ScanStat=="kullT3y_Sens", ]
# res2$ScanStat <- factor(res2$ScanStat, levels = c("kullT03_Sens","kullT_Sens","kullT3y_Sens"))
# 
# res3<-res1[res1$ScanStat=="kullT03_Spec" | res1$ScanStat=="kullT_Spec" | res1$ScanStat=="kullT3y_Spec", ]
# res3$ScanStat <- factor(res3$ScanStat, levels = c("kullT03_Spec","kullT_Spec","kullT3y_Spec"))
# 
# 
# ggplot() +
#   geom_pointrange(data=res2, mapping=aes(x=ScanStat, y=median, ymin=LCI, ymax=UCI, color=type), width=0.2, size=1) + 
#   facet_wrap(.~type)+
#   labs(x = "PfPR (2-10)", y="incidence reduction (%)",
#        title = "Sensitivity by type of hotspot", subtitle = "method: kulldorf") +
#   theme(axis.text=element_text(size=8),axis.title=element_text(size=12),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         text = element_text(size = 12,face="bold"),legend.text = element_text(size=8,face="bold"))+
#   coord_cartesian(ylim = c(0.0,60.0))
# 
# ggplot() +
#   geom_pointrange(data=res3, mapping=aes(x=ScanStat, y=median, ymin=LCI, ymax=UCI, color=type), width=0.2, size=1) + 
#   facet_wrap(.~type)+
#   labs(x = "PfPR (2-10)", y="incidence reduction (%)",
#        title = "Specificity by type of hotspot", subtitle = "method: kulldorf") +
#   theme(axis.text=element_text(size=8),axis.title=element_text(size=12),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         text = element_text(size = 12,face="bold"),legend.text = element_text(size=8,face="bold"))+
#   coord_cartesian(ylim = c(80.0,100.0))
# 
# # plot the simulated and detected hotspots on grids
# 
# hsingle<-data.frame(read.table("hotspots_single.txt", header = T))
# hsingle<-data.frame(read.table("hotspots_multiple.txt", header = T))
# hsingle<-data.frame(read.table("hotspots_decay.txt", header = T))
# hsingle<-data.frame(read.table("hotspots_river.txt", header = T))
# hsingle<-data.frame(read.table("hotspots_long.txt", header = T))
# hsingle<-data.frame(read.table("hotspots_short.txt", header = T))
# hsingle<-data.frame(read.table("hotspots_season.txt", header = T))
# hsingle<-data.frame(read.table("hotspots_dry.txt", header = T))
# 
# hsingle_long <- gather(hsingle, measures, value, hotspot:kullT3y, factor_key=TRUE)
# head(hsingle_long)
# #hsingle_long$value<-as.factor(hsingle_long$value)
# hsingle_long$measures<-factor(hsingle_long$measures, levels = c("hotspot","kullT03", "kullT", "kullT3y"))
# levels(hsingle_long$measures) = c("hotspot","3 months", "1 year", "3 years")
# 
# ggplot() +
#   geom_point(data=hsingle_long, mapping=aes(x=locx, y=locy, color = value), width=0.2, size=2) + 
#   facet_wrap(.~measures)+
#   labs(x = "Longitude", y="Latitude",
#        title = "Decay gradient", subtitle = "method: kulldorf") +
#   theme_pubclean()+
#   theme(axis.text=element_text(size=8),axis.title=element_text(size=12),axis.text.x = element_text(angle = 45, vjust = 1, hjust=1),
#         text = element_text(size = 12,face="bold"),legend.position = "none")+
#   #scale_fill_manual(values = coloure) +
#   coord_cartesian(ylim = c(0.0,10.0))
# 


