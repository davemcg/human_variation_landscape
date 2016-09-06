library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
source('analysis/theme_pub.R')
library(caret)
library(ROCR)
parse_vcf_info <- function(info_column, field) {
  # grabs VCF info field
  # first splits by ;
  # then grabs section by grep
  sapply(info_column, function(x) {
    split_up = strsplit(x,';')[[1]]
    split_up[grep(field,split_up)]
    })
}

parse_vcf_info_extended <- function(info_column, field, col_num) {
  # grabs VCF info field
  # first splits by ;
  # then grabs section by grep
  # THEN grabs the n column split by '|'
  sapply(info_column, function(x) {
    split_up = strsplit(x,';')[[1]]
    section <- split_up[grep(field,split_up)]
    strsplit(section,'\\|')[[1]][col_num]
  })
}


grab_value <- function(info,delim) {
  # VCF commonly uses INFO=213412 format
  # in this example '=' would be the delim
  # this function will grab the value and
  # convert to numeric
  sapply(info, function(x) {
    as.numeric(strsplit(x,delim)[[1]][2])
  })
}
# Homo_sapiens_incl_consequences__codingOnly.VLNum.VEPnoPick.GRCh38.20160906.vcf processed by a custom script:
# ~/git/human_variation_landscape/scripts/vcf_to_tab.py Homo_sapiens_incl_consequences__codingOnly.VLNum.VEPnoPick.GRCh38.20160906.vcf &
# grabs selected columns and outputs as tab separated file
clin_vcf <- fread('~/Desktop/Homo_sapiens_incl_consequences__codingOnly.VLNum.VEPnoPick.GRCh38.20160906.my.tab')


# change NA to 0 in GMAF and exac
GMAF0 <- clin_vcf$GMAF
GMAF0[is.na(GMAF0)] <- 0
clin_vcf$GMAF0 <- GMAF0
exac0 <- clin_vcf$ExAC
exac0[is.na(exac0)] <- 0
clin_vcf$ExAC0 <- exac0

# grab the clin pathogenic set and then a matched size non pathogenic set*
# * defined as a maf > 0.05
clin_Path <- clin_vcf %>% filter(Clin_Sig=='CLIN_pathogenic',GMAF0<0.01) %>% mutate(Status='Pathogenic')
clin_Benign <- clin_vcf %>% filter(Clin_Sig=='CLIN_benign') %>% mutate(Status='ClinVar_Benign')
set.seed(95342)
clin_nonPath_LowMAF <- clin_vcf %>% filter(!grepl('pathogenic',Clin_Sig),GMAF0<0.01,ExAC0<0.01) %>% sample_n(nrow(clin_Path)) %>% mutate(Status='Benign_lowMAF')
clin_nonPath_HighMAF <- clin_vcf %>% filter(!grepl('pathogenic',Clin_Sig),(GMAF0>0.05|ExAC0>0.05)) %>% sample_n(nrow(clin_Path)) %>% mutate(Status='Benign_highMAF')
clin_set <- rbind(clin_Path,clin_Benign, clin_nonPath_LowMAF,clin_nonPath_HighMAF)



# grab 60% of the data (rows) for the training set, 20% for query set, and the last 20% for the final test set (only use once!!!!!)
row_index <- seq(1,nrow(clin_set))
set.seed(95342)
train60_index <- sample(row_index, 0.6*nrow(clin_set), replace=FALSE)
remainder = row_index[!(row_index %in% train60_index)]
query20_index <- sample(row_index, 0.5*length(remainder), replace=FALSE)
test20_index <- remainder[!(remainder %in% query20_index)]
# now get the train60set for our playground
training_set <- clin_set[train60_index, ]
# the rest, for later
query20_set <- clin_set[query20_index, ]
test20_set <- clin_set[test20_index, ]


##########
# exploratory analysis, quick stats

# gather/melt
training_set.m <- training_set %>% gather(.,VLN,Value,-Chr,-Position,-ID,-Ref,-Alt,-Variant_Class,-GMAF,-ExAC,-CADD_phred, -CADD_raw, -CADD_raw_rankscore, -PolyPhen, -SIFT, -Clin_Sig, -GMAF0, -ExAC0, -Status)
# density plot
training_set %>% ggplot(aes(x=VLN100all, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,100))
# Uh....everything looks about the same
# summary stats by VLNtype and Pathogenic/Benign
training_set.m %>% group_by(VLN,Status) %>% summarise(mean(Value), median(Value),quantile(Value,0.25),quantile(Value,0.75),Total=n())
# So....pathogenic seems to generally have a higher VLN than the 'benign'. Which is opposite my initial finding. And our expectation. 
# What is going on here? The answer is very simple, but it took me a while to figure out.
# Let's pull in ensembl's clinvar_only variant set, which was used to generate the plot that we've seen before
vcf_file = '~/Desktop/Homo_sapiens_clinically_associated.VLNum.vcf.gz'
rows_to_skip <- as.integer(system( paste(paste('gzcat ',vcf_file),' | head -n 1000 | grep ^## | wc -l'), intern=T))
clinOnly_vcf <- fread(paste('gzcat ',vcf_file),skip=rows_to_skip)
orig_clin <- clinOnly_vcf
# now process to get VLN and CLINSIG info
clinOnly_vcf <- clinOnly_vcf  %>% filter(grepl('VLN100all=',V8))
clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('VLN100all',V8)) %>% mutate(VLN100all=grab_value(parse_vcf_info(V8,'VLN100all'),'='))
clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('VLN100high',V8)) %>% mutate(VLN100high=grab_value(parse_vcf_info(V8,'VLN100high'),'='))
clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('VLN100moderate',V8)) %>% mutate(VLN100moderate=grab_value(parse_vcf_info(V8,'VLN100moderate'),'='))
clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('VLN100low',V8)) %>% mutate(VLN100low=grab_value(parse_vcf_info(V8,'VLN100low'),'='))
clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('CLIN',V8)) %>% mutate(ClinSig=parse_vcf_info(V8,'CLIN'))
# now let's label the clin sig as either benign, pathogenic, or conflicting (dual labelled)
clinOnly_vcf <- clinOnly_vcf %>% mutate(Status=ifelse(grepl('(?=.*CLIN_pathogenic.*)(?=.*CLIN_benign)', V8, perl=TRUE),'Conflicting', ifelse(grepl('pathogenic',V8),'Pathogenic',ifelse(grepl('benign',V8),'Benign','Other'))))
# recreate the density plot we've seen before:
clinOnly_vcf %>% filter(Status!='Conflicting', Status!='Other') %>% 
  ggplot(aes(x=VLN100all, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,100))
# summary stats, as above:
clinOnly_vcf %>% gather(.,VLN,Value,-V1,-V2,-V3,-V4,-V5,-V6,-V7,-V8,-ClinSig,-Status) %>% group_by(VLN, Status) %>% filter(Status!='Conflicting', Status!='Other') %>% summarise(mean(Value), median(Value),quantile(Value,0.25),quantile(Value,0.75),Total=n())
# what immediately pops out is how FEW benign variants there are. Did I screw up the filtering?
# checking number of benign (of ANY kind) in original vcf
orig_clin %>% filter(grepl('benign',V8)) %>% nrow()
# 1238....so for some reason there are VERY VERY few benign's

# I then downloaded the latest clinvar vcf from NCBI's ftp site
# ftp://ftp.ncbi.nlm.nih.gov/pub/clinvar/ 
# clinvar_20160831.vcf.gz
# used vcfanno to add my VLN stats to the vcf then will process the same as above
vcf_file = '~/Desktop/clinvar_20160831.VLN.vcf.gz'
rows_to_skip <- as.integer(system( paste(paste('gzcat ',vcf_file),' | head -n 1000 | grep ^## | wc -l'), intern=T))
clinOnly_vcf_latest <- fread(paste('gzcat ',vcf_file),skip=rows_to_skip)
clinOnly_vcf_latest <- clinOnly_vcf_latest  %>% filter(grepl('VLN100all=',V8))
clinOnly_vcf_latest <- clinOnly_vcf_latest %>% filter(grepl('VLN100all',V8)) %>% mutate(VLN100all=grab_value(parse_vcf_info(V8,'VLN100all'),'='))
clinOnly_vcf_latest <- clinOnly_vcf_latest %>% filter(grepl('VLN100high',V8)) %>% mutate(VLN100high=grab_value(parse_vcf_info(V8,'VLN100high'),'='))
clinOnly_vcf_latest <- clinOnly_vcf_latest %>% filter(grepl('VLN100moderate',V8)) %>% mutate(VLN100moderate=grab_value(parse_vcf_info(V8,'VLN100moderate'),'='))
clinOnly_vcf_latest <- clinOnly_vcf_latest %>% filter(grepl('VLN100low',V8)) %>% mutate(VLN100low=grab_value(parse_vcf_info(V8,'VLN100low'),'='))
clinOnly_vcf_latest <- clinOnly_vcf_latest %>% filter(grepl('CLNSIG',V8)) %>% mutate(ClinSig=grab_value(parse_vcf_info(V8,'CLNSIG'),'='))
# clinvar's release uses '2' as benign and '5' as pathogenic
# let's make a new density plot
clinOnly_vcf_latest %>% filter(ClinSig %in% c(2,5)) %>% mutate(Status=ifelse(ClinSig==2,'Benign','Pathogenic')) %>% 
  ggplot(aes(x=VLN100all, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,100))
# summary stats
clinOnly_vcf_latest %>% filter(ClinSig %in% c(2,5)) %>% mutate(Status=ifelse(ClinSig==2,'Benign','Pathogenic')) %>% 
  gather(.,VLN,Value,-V1,-V2,-V3,-V4,-V5,-V6,-V7,-V8,-ClinSig,-Status) %>% 
  group_by(VLN, Status) %>% summarise(mean(Value), median(Value),quantile(Value,0.25),quantile(Value,0.75),Total=n())
# OK....now we have 30,000 benign variants
# So, unfortunately I used a poor data set for my first pass analysis

# well, let's proceed since I've done 95% of the hard work and see whether the VLN numbers have any ability to distinguish
# benign from pathogenic variants
# keep your expectations low. very low.
# first let's see how well CADD does
# since CADD did a ROC test on higher MAF variants vs clinvar pathogenic, I'll do the same
training_set_sub <- training_set %>% filter(Status=='Pathogenic'|Status=='Benign_highMAF') %>% filter(!is.na(CADD_raw_rankscore))
# train on 60% of the data
bayesglm_Fit_CADDonly <- train(factor(Status) ~ CADD_raw_rankscore, data=training_set_sub,  preProc = c("center", "scale"), method='bayesglm',)
# predict on query20 set
query20_set_sub <- query_20_set %>% filter(Status=='Pathogenic'|Status=='Benign_highMAF') %>% filter(!is.na(CADD_raw_rankscore))
bayesglm_predict_CADDonly <- predict(bayesglm_Fit_CADDonly, newdata=query20_set_sub , type='prob')


pred <- prediction(bayesglm_predict_CADDonly$Benign_highMAF, factor(query20_set_sub$Status))
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col='green')
auc.tmp <- performance(pred,"auc"); as.numeric(auc.tmp@y.values)



# density plot
ggplot(data=training_set.m , aes(x=Value, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,50)) + facet_wrap(~VLN)









####clinvar only vcf######
#vcf_file = '~/Desktop/Homo_sapiens_clinically_associated.VLNum.vcf.gz'
#vcf_file = '~/Desktop/Homo_sapiens_incl_consequences__codingOnly.VEP.GRCh38.VLN.CLINpath_benign.vcf.gz'
vcf_file = '~/Desktop/clinvar_20160831.VLN.vcf.gz'
rows_to_skip <- as.integer(system( paste(paste('gzcat ',vcf_file),' | head -n 1000 | grep ^## | wc -l'), intern=T))
clinOnly_vcf <- fread(paste('gzcat ',vcf_file),skip=rows_to_skip)

clinOnly_vcf <- clinOnly_vcf  %>% filter(grepl('VLN100all=',V8))

clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('VLN100all',V8)) %>% mutate(VLN100all=grab_value(parse_vcf_info(V8,'VLN100all'),'='))
clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('VLN100high',V8)) %>% mutate(VLN100high=grab_value(parse_vcf_info(V8,'VLN100high'),'='))
clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('VLN100moderate',V8)) %>% mutate(VLN100moderate=grab_value(parse_vcf_info(V8,'VLN100moderate'),'='))
clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('VLN100low',V8)) %>% mutate(VLN100low=grab_value(parse_vcf_info(V8,'VLN100low'),'='))

clinOnly_vcf <- clinOnly_vcf %>% filter(grepl('CLNSIG',V8)) %>% mutate(ClinSig=grab_value(parse_vcf_info(V8,'CLNSIG'),'='))

clinOnly_vcf %>% mutate(Status=ifelse(grepl('(?=.*CLIN_pathogenic.*)(?=.*CLIN_benign)', V8, perl=TRUE),'Conflicting', ifelse(grepl('CLIN_pathogenic',V8),'Pathogenic',ifelse(grepl('CLIN_benign',V8),'Benign','Absent')))) %>% 
  ggplot(aes(x=VLN100all, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,100))


######## ROC AUC for CADD and path/benign##########
set <- clin_vcf %>% filter(grepl('CLIN_pathogenic',Clin_Sig)|GMAF>0.05|ExAC>0.05)  %>% filter(CADD_raw_rankscore>0) %>% mutate(GMAF0=ifelse(is.na(GMAF),0,GMAF), ExAC0=ifelse(is.na(ExAC),0,ExAC), Class=ifelse(GMAF0>0.05|ExAC0>0.05,0,1))


row_index <- seq(1,nrow(set))
set.seed(95342)
train60_index <- sample(row_index, 0.6*nrow(set), replace=FALSE)
remainder = row_index[!(row_index %in% train60_index)]
query20_index <- sample(row_index, 0.5*length(remainder), replace=FALSE)
test20_index <- remainder[!(remainder %in% query20_index)]
# now get the train60set for our playground
training_set <- set[train60_index, ]
training_set_sub <- training_set %>% dplyr::select(Class, CADD_raw_rankscore,VLN100high,VLN100moderate,VLN100low)
# scale data
temp <- training_set_sub  %>% dplyr::select(-Class)  %>% scale() %>% data.frame() 
temp$Class <- training_set_sub$Class
training_set_sub <- temp

# the rest, for later
query20_set <- set[query20_index, ]
query20_set_sub <- query20_set %>% dplyr::select(Class, CADD_raw_rankscore,VLN100high,VLN100moderate,VLN100low)
temp <- query20_set_sub  %>% dplyr::select(-Class)  %>% scale() %>% data.frame() 
temp$Class <- query20_set_sub$Class
query20_set_sub <- temp


test20_set <- set[test20_index, ]


library(ROCR) # http://rocr.bioinf.mpi-sb.mpg.de
pred <- prediction(set$CADD_raw_rankscore, set$Class)
perf <- performance(pred, measure = "tpr", x.measure = "fpr") 
plot(perf, col='green',add=TRUE)
auc.tmp <- performance(pred,"auc"); as.numeric(auc.tmp@y.values)

library(e1071) #svm
svm.model <- svm(Class ~ ., data=training_set_sub) # build model on training_set (60% of data)
svm.pred <- predict(svm.model, query20_set_sub[,-5]) # run predictions on query subset, remove the classification
perf.svm <- performance(svm.pred,  'tpr' ,  'fpr' )
auc.tmp <- performance(svm.pred,"auc"); as.numeric(auc.tmp@y.values)
plot(perf.svm, col=rainbow(10))










#######old code##########



clin_vcf_p <- clin_vcf  %>% filter(grepl('VLN100high=',V8))

clin_vcf_p$VLN100high <- grab_value(parse_vcf_info(clin_vcf_p$V8,'VLN100high='))
clin_vcf_p %>% mutate(Status=ifelse(grepl('benign', V8),'Benign','Pathogenic')) %>% 
  mutate(Variants='Missense Only') %>% 
ggplot(aes(x=VLN100high, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,100))

temp <- clin_vcf_p %>% mutate(Status=ifelse(grepl('benign', V8),'Benign','Pathogenic')) %>% 
  mutate(Variants='Missense Only')
#####

vcf_file = '~/Desktop/Homo_sapiens_clinically_associated.VLNum.vcf.gz'

rows_to_skip <- as.integer(system( paste(paste('gzcat ',vcf_file),' | head -n 1000 | grep ^## | wc -l'), intern=T))
clin_vcfAll <- fread(paste('gzcat ',vcf_file),skip=rows_to_skip)

clin_vcfAll_p <- clin_vcfAll  %>% filter(grepl('VLN100=',V8))

clin_vcfAll_p$VLN100 <- grab_value(parse_vcf_info(clin_vcfAll_p$V8,'VLN100='))
clin_vcfAll_p %>% mutate(Status=ifelse(grepl('benign', V8),'Benign','Pathogenic')) %>% 
  mutate(Variants='All') %>% rbind(.,temp) %>% 
  ggplot(aes(x=VLN100, colour=Status))+geom_density() + 
    theme_Publication() + 
    coord_cartesian(xlim=c(0,100)) + 
    facet_wrap(~Variants,nrow=1)


clin_vcfAll_p %>% mutate(Status=ifelse(grepl('benign', V8),'Benign','Pathogenic')) %>% 
  mutate(Variants='All') %>% rbind(.,temp) %>% filter(Variants=='Missense Only') %>% 
  ggplot(aes(x=VLN100, colour=Status))+geom_density() + 
  theme_Publication() + 
  coord_cartesian(xlim=c(0,100)) + 
  facet_wrap(~Variants,nrow=1)

#abca4  1	93992834	94121132
#pnpla6 19	7534003	7561764
#ush2a  1	215622893	216423396
