#processing lab files that contain fluoroprobe data and/or toxin results
#Alene Onion
#3/12/2019

##################################################################################################################################
#read raw data files
##################################################################################################################################
#capture current working directory so I can return to it
current_wd<-getwd()
#set working directory to raw data folder
setwd("sections/data/raw.data")

files<-list.files()
nfiles<-length(files)

#read in templates for comparison
template1<-read.csv("template1.csv",stringsAsFactors = FALSE)
template2<-read.csv("template2.csv",stringsAsFactors = FALSE)

#return working directory to original location
setwd(current_wd)
rm(current_wd)

#capture the column names to check formatting of submitted data
colnames1<-colnames(template1)
colnames2<-colnames(template2)

#create a data file with all the appropriate column headers
lab<-template1
lab<-merge(lab,template2)

#remove template files from memory
rm(list=c('template1','template2'))

#set a function to detect if a vector is empty
vector.isnot.empty <- function(x) return(length(x) !=0)

#loop through all the files
#first check if the column names match the template
#then merge them together in one large file called 'lab'
for(i in 1:nfiles){
  temp<-read.csv(paste("sections/data/raw.data/",files[i],sep=""),stringsAsFactors=FALSE)
  #checking to see whether the column names match the template
  #this prints the differences if there are any for the user to check
    colnamestemp<-colnames(temp)
    diff1<-setdiff(colnames1,colnamestemp)
    diff2<-setdiff(colnames2,colnamestemp)
    if(vector.isnot.empty(diff1) && vector.isnot.empty(diff2)){
      print(paste("The file names ",files[i]," does not match the template. These columns are missing ",diff1,sep=""))
      print(paste("The file names ",files[i]," does not match the template. These columns are missing ",diff2,sep=""))
      stop("this script stopped because you need to go fix the lab file column names")
    }
    lab<-merge(lab,temp,all=TRUE)
    rm(list=c('colnamestemp','diff1','diff2','temp'))
}
rm(list=c('colnames1','colnames2','files','nfiles','i','vector.isnot.empty'))

##################################################################################################################################
#manipulate file before assess hab status
##################################################################################################################################

#remove bottom samples from CSLAP samples
lab<-lab[lab$info_type!="BS",]

#reduce all microscopy to lower case
lab$Microscopy<-tolower(lab$Microscopy)

#create HABSTATUS column and HABSTATUS_date column
lab$HABSTATUS<-NA
lab$STATUS_DATE<-Sys.Date()

#create hab status column
lab$HABSTATUS<-ifelse(lab$Bluegreen<25,"No Bloom",lab$HABSTATUS)

##################################################################################################################################
#assess hab status
##################################################################################################################################

lab$HABSTATUS<-ifelse(lab$info_type=="OW",
                                ifelse(lab$Bluegreen>25,
                                      ifelse(grepl("aphanizomenon",lab$Microscopy)|grepl("anabaena",lab$Microscopy)|grepl("dolichospermum",lab$Microscopy)|grepl("woronchinia",lab$Microscopy),
                                             "Open Water", 
                                             "Needs Evaluation"),
                                      lab$HABSTATUS),
                                lab$HABSTATUS)
lab$HABSTATUS<-ifelse(lab$info_type=="SB",
                       ifelse(lab$Bluegreen>25,
                              ifelse(grepl("aphanizomenon",lab$Microscopy)|grepl("anabaena",lab$Microscopy)|grepl("dolichospermum",lab$Microscopy)|grepl("woronchinia",lab$Microscopy),
                                     "Confirmed", 
                                     "Needs Evaluation"),
                              lab$HABSTATUS),
                       lab$HABSTATUS)

#to respond to the microcystin results, we need to convert this to a numeric field
#first create a new column so we can remove ND values without removing them in the final dataset
lab$microcystincopy<-lab$microcystin
#convert ND microcystin results to NA
lab$microcystincopy<-ifelse(lab$microcystincopy=="ND"|lab$microcystincopy=="Nd"|lab$microcystincopy=="nD"|lab$microcystincopy=="nd",NA,lab$microcystincopy)
#remove commas from numbers before converting to numeric
lab$microcystincopy<-gsub(",","",lab$microcystincopy)
#convert microcystin values to numeric
lab$microcystincopy<-as.numeric(lab$microcystincopy,na.rm = TRUE)
#use these values to change hab status
lab$HABSTATUS<-ifelse(lab$info_type=="OW",
                      ifelse(lab$microcystin>=10 & (!is.na(lab$microcystincopy)),
                             "Confirmed with High Toxins",
                             lab$HABSTATUS),
                      ifelse(lab$microcystin>=20 & (!is.na(lab$microcystincopy)),
                              "Confirmed with High Toxins",
                              lab$HABSTATUS))
#remove the microcystincopy column
lab$microcystincopy<-NULL

##################################################################################################################################
#pull in human evaluations
##################################################################################################################################

#capture current working directory so I can return to it
current_wd<-getwd()
#set working directory to raw data folder
setwd("sections/data/needs_evaluation")

files<-list.files()
nfiles<-length(files)

#read in templates to create starting data file
template<-read.csv("template.csv",stringsAsFactors = FALSE)
#create a data file with all the appropriate column headers
needs<-template
#remove template files from memory
rm(list=c('template'))

#return working directory to original location
setwd(current_wd)
rm(current_wd)

#loop through all the files
#first check if the column names match the template
#then merge them together in one large file called 'lab'
for(i in 1:nfiles){
  temp<-read.csv(paste("sections/data/needs_evaluation/",files[i],sep=""),stringsAsFactors=FALSE)
  needs<-merge(needs,temp,all=TRUE)
  rm(temp)
}
rm(list=c('files','nfiles','i'))

#remove hab status date or it won't overlap
needs$STATUS_DATE<-NULL
#remove rows where the evaulation is NA - therefore a human hasn't looked at it
needs$evaluation<-ifelse(needs$evaluation=="",NA,needs$evaluation)
needs<-needs[!is.na(needs$evaluation),]
#make sure the hab status has the proper capitalization
needs$evaluation<-stringr::str_to_title(needs$evaluation)

#merge the two files
lab<-merge(lab,needs,all.x=TRUE)
lab$HABSTATUS<-ifelse(is.na(lab$evaluation),lab$HABSTATUS,lab$evaluation)
lab$evaluation<-NULL
rm(list=c('needs'))

##################################################################################################################################
#subset the list of rows that are needs evaluation - meaning they need human eyes to evaluate them
##################################################################################################################################
#write function to detect if a vector is empty
vector.isnot.empty <- function(x) return(length(x) !=0)

#generate table of needs evaluation
needs<-lab[lab$HABSTATUS=="Needs Evaluation",]
needs<-needs[!is.na(needs$HABSTATUS),]
if(vector.isnot.empty(needs$HABSTATUS)){
  needs$evaluation<-NA
  needs<-needs[!is.na(needs$HABSTATUS),]
  needsfile<-paste("sections/data/needs_evaluation/needs.evaluation_",format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),".csv",sep="")
  write.csv(needs,file=needsfile,row.names = FALSE)
  print(paste("these are the number of records that need evaluating:",length(needs$sample_id),sep=" "))
  rm(list=c('needs','needsfile'))
  stop("this script stopped because you need evaluate some records in the needs_evaluation folder")
}
rm(list=c('vector.isnot.empty','needs'))

##################################################################################################################################
#generate the narrative
##################################################################################################################################

#create the narrative column
lab$narrative<-NA
#create a list of info types in long form
lab$info_type_long<-NA
lab$info_type_long<-ifelse(lab$info_type=="OW","open water","shore bloom")
#generate the narratives
lab$narrative<-ifelse(lab$HABSTATUS=="No Bloom",
                      paste("The results from the ",lab$date,", ",lab$info_type_long," sample submitted showed blue-green chlorophyll a levels below the DEC Confirmed Bloom threshold of 25 µg/L.",sep=""),
                      lab$narrative)
lab$narrative<-ifelse(lab$HABSTATUS=="Confirmed"|lab$HABSTATUS=="Open Water",
                      paste("The results from the ",lab$date,", ",lab$info_type_long," confirm the presence of a cyanobacteria HAB based on blue-green chlorophyll a levels above the DEC Confirmed Bloom threshold of 25 µg/L and a microscopic analysis of cyanobacteria. The size of the bloom according to the sampler is ",lab$extent,".",sep=""),
                      lab$narrative)
lab$narrative<-ifelse(lab$HABSTATUS=="Confirmed with High Toxins",
                      paste("The toxin results from the ",lab$date,", ",lab$info_type_long," confirmed elevated toxin levels. The microcystin concentration was above the DEC Confirmed with High Toxins Bloom threshold of 20 µg/L. The size of the bloom according to the sampler is ",lab$extent,".",sep=""),
                      lab$narrative)
lab$info_type_long<-NULL

##################################################################################################################################
#print the file for upload
##################################################################################################################################
#divide into toxin results and fluroprobe results
toxin<-lab[!is.na(lab$microcystin),]
other<-lab[!is.na(lab$Microscopy),]

filename<-paste("sections/data/for_upload/for_upload_",format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"_toxins.csv",sep="")
write.csv(toxin,filename,row.names=FALSE)
filename<-paste("sections/data/for_upload/for_upload_",format(Sys.time(),"%Y-%m-%d_%H.%M.%S"),"_other.csv",sep="")
write.csv(other,filename,row.names=FALSE)
print(paste("the lab data were successfully assessed and are ready for upload. See ",filename," and _toxins",sep=""))
rm(list=c('toxin','other','lab','filename'))