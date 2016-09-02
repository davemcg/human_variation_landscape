library(data.table)
library(dplyr)
library(ggplot2)
library(tidyr)
source('analysis/theme_pub.R')
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

vcf_file = '~/Desktop/Homo_sapiens_incl_consequences__codingOnly.VEP.GRCh38.VLN.vcf.gz'

rows_to_skip <- as.integer(system( paste(paste('gzcat ',vcf_file),' | head -n 1000 | grep ^## | wc -l'), intern=T))
clin_vcf <- fread(paste('gzcat ',vcf_file),skip=rows_to_skip)

#label variants as benign (not also labelled pathogenic), pathogenic (not also labelled benign), conflicting (double labelled), and not labelled at all (most of them)
clin_vcf <- clin_vcf %>% mutate(ClinVar_Status=ifelse(grepl('(?=.*CLIN_pathogenic.*)(?=.*CLIN_benign)', V8, perl=TRUE),'Conflicting', ifelse(grepl('CLIN_pathogenic',V8),'Pathogenic',ifelse(grepl('CLIN_benign',V8),'Benign','Absent'))))

# pull gmaf (28), only for variants with this info (some are missing???)
clin_vcf <- clin_vcf %>% filter(grepl('CSQ=',V8)) %>% filter(grepl('VLN',V8)) %>% mutate(GMAF=grab_value(parse_vcf_info_extended(V8,'CSQ',28),':'))
# some issue with exac
#clin_vcf <- clin_vcf %>% filter(grepl('CSQ=',V8)) %>% mutate(ExAC_Adj=parse_vcf_info_extended(V8,'CSQ',37))
# get variant status, toss variants without VLN scores
clin_vcf <- clin_vcf %>% filter(grepl('CSQ=',V8)) %>% filter(grepl('VLN',V8)) %>% mutate(GMAF=parse_vcf_info_extended(V8,'CSQ',2),':')
# pull VLN scores
clin_vcf <- clin_vcf %>% mutate(VLN100all=grab_value(parse_vcf_info(V8,'VLN100all'),'='))
clin_vcf <- clin_vcf %>% mutate(VLN100high=grab_value(parse_vcf_info(V8,'VLN100high'),'='))
clin_vcf <- clin_vcf %>% mutate(VLN100moderate=grab_value(parse_vcf_info(V8,'VLN100moderate'),'='))
clin_vcf <- clin_vcf %>% mutate(VLN100low=grab_value(parse_vcf_info(V8,'VLN100low'),'='))
# save, as the above takes OVERNIGHT to calculate
save(clin_vcf,file='clin_vcf.Rdata')

# change NA to 0 in GMAF
GMAF0 <- clin_vcf$GMAF
GMAF0[is.na(GMAF0)] <- 0
clin_vcf$GMAF0 <- GMAF0

# grab the clin pathogenic set and then a matched size non pathogenic set*
# * defined as a maf > 0.05
clin_Path <- clin_vcf %>% filter(ClinVar_Status=='Pathogenic',GMAF0<0.01) %>% mutate(Status='Pathogenic')
clin_Benign <- clin_vcf %>% filter(ClinVar_Status=='Benign') %>% mutate(Status='ClinVar_Benign')
set.seed(95342)
clin_nonPath_LowMAF <- clin_vcf %>% filter(ClinVar_Status!='Pathogenic',ClinVar_Status!='Conflicting',GMAF0<0.01) %>% sample_n(nrow(clin_Path)) %>% mutate(Status='Benign_lowMAF')
clin_nonPath_HighMAF <- clin_vcf %>% filter(ClinVar_Status!='Pathogenic',ClinVar_Status!='Conflicting',GMAF>0.05) %>% sample_n(nrow(clin_Path)) %>% mutate(Status='Benign_highMAF')
clin_set <- rbind(clin_Path,clin_Benign, clin_nonPath_LowMAF,clin_nonPath_HighMAF)



# grab 60% of the data (rows) for the training set, 20% for query set, and the last 20% for the test set (only use once!!!!!)
# but first toss variants that don't have my VLNumbers and are conflicint
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
# exploratory analysis


# quick correct for synonymous (baseline?) rate
# training_clin_vcf$VLNmoderate_over_low <- (training_clin_vcf$VLN100moderate+1) / (training_clin_vcf$VLN100low+1)

# gather/melt
training_set.m <- training_set %>% gather(.,VLN,Value,-V1,-V2,-V3,-V4,-V5,-V6,-V7,-V8,-ClinVar_Status,-Status,-GMAF, -GMAF0, -Variant_Type)
# summary stats by VLNtype and Pathogenic/Benign
training_set.m %>% group_by(VLN,Status) %>% summarise(mean(Value), median(Value),quantile(Value,0.25),quantile(Value,0.75))
# density plot
ggplot(data=training_set.m , aes(x=Value, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,50)) + facet_wrap(~VLN)









####clinvar only vcf######
#vcf_file = '~/Desktop/Homo_sapiens_clinically_associated.VLNum.vcf.gz'
vcf_file = '~/Desktop/Homo_sapiens_incl_consequences__codingOnly.VEP.GRCh38.VLN.CLINpath_benign.vcf.gz'
rows_to_skip <- as.integer(system( paste(paste('gzcat ',vcf_file),' | head -n 1000 | grep ^## | wc -l'), intern=T))
clinOnly_vcf <- fread(paste('gzcat ',vcf_file),skip=rows_to_skip)

clinOnly_vcf <- clinOnly_vcf  %>% filter(grepl('VLN100all=',V8))

clinOnly_vcf$VLN100all <- (grab_value(parse_vcf_info(clinOnly_vcf$V8,'VLN100all='),'='))
clinOnly_vcf %>% mutate(Status=ifelse(grepl('(?=.*CLIN_pathogenic.*)(?=.*CLIN_benign)', V8, perl=TRUE),'Conflicting', ifelse(grepl('CLIN_pathogenic',V8),'Pathogenic',ifelse(grepl('CLIN_benign',V8),'Benign','Absent')))) %>% 
  ggplot(aes(x=VLN100all, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,100))



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