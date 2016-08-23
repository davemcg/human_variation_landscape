
parse_vcf_info <- function(info_column, field) {
  # grabs VCF info field
  # first splits by ;
  # then grabs section by grep
  sapply(info_column, function(x) {
    split_up = strsplit(x,';')[[1]]
    split_up[grep(field,split_up)]
    })
}

grab_value <- function(info) {
  # VCF commonly uses INFO=213412 format
  # this function will grab the value and
  # convert to numeric
  sapply(info, function(x) {
    as.numeric(strsplit(x,'=')[[1]][2])
  })
}

vcf_file = '~/Desktop/Homo_sapiens_clinically_associated.VLNumbers.vcf.gz'

rows_to_skip <- as.integer(system( paste(paste('gzcat ',vcf_file),' | head -n 1000 | grep ^## | wc -l'), intern=T))
clin_vcf <- fread(paste('gzcat ',vcf_file),skip=rows_to_skip)

clin_vcf_p <- clin_vcf  %>% filter(grepl('VLN100=',V8))

clin_vcf_p$VLN100 <- grab_value(parse_vcf_info(clin_vcf_p$V8,'VLN100='))
clin_vcf_p %>% mutate(Status=ifelse(grepl('benign', V8),'Benign','Pathogenic')) %>% 
  mutate(Variants='Missense Only') %>% 
  ggplot(aes(x=VLN100, colour=Status))+geom_density() + theme_Publication() + coord_cartesian(xlim=c(0,100))

temp <- clin_vcf_p %>% mutate(Status=ifelse(grepl('benign', V8),'Benign','Pathogenic')) %>% 
  mutate(Variants='Missense Only')
#####

vcf_file = '~/Desktop/Homo_sapiens_clinically_associated.VLNumbers_all.vcf.gz'

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