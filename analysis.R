# Analysis of V3V4 16S rRNA data from longitudinal cohort of healthy adult males. 
# OTUs and Phylotypes were constructed in mothur code availabile on github repo mothur.script 
# Library load 
library(tidyverse)
library(vegan)
library(ggrepel)
library(reshape2)
# Set working direcotry 
#setwd("/media/uporabnik/Data1/projekti/SMNS2024")

# Determine seed, theme and colors used in the analysis 
set.seed(1996)
# Colors for microbiota (colm) and sporobiota (cols)
colm = c('#47B66A')
cols = c('#D9A534')
colsm = c('#47B66A', '#D9A534')
# Colors for extreme events
col_extreme = c("#8A05EC",
                "#8d1d80",
                "#96aa00",
                "#0054bd",
                "#ffab2d",
                "#ff7d01",
                "#02c2ef",
                "#c10065",
                "#007d36",
                "#505d00",
                "#ffa7dc")


# Import OTU table and taxonomy table and remove rare OTUs from taxtab
  # Use a rarefied table! 
shared_pre = read_tsv('data/mothur/final-48.shared') %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  mutate(total_sample=sum(value)) %>%
  filter(total_sample > 100000) %>%
  ungroup() %>%
  select(-total_sample)

otutab = dcast(shared_pre, Group ~ name, value.var = 'value') %>%
  column_to_rownames('Group')

# Rarefy the data once -as we sequenced so deep that for the first analysis this is not crutial! 
otu_rare = rrarefy(otutab, sample=100000) %>% as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols=2:3002) %>% 
  group_by(name) %>%
  mutate(total_otu=sum(value)) %>%
  filter(total_otu > 4) %>%
  ungroup() %>%
  select(-total_otu) %>%
  pivot_wider(names_from = 'name', values_from = 'value')

rm(otutab)
rm(shared_pre)

otus=colnames(otu_rare %>% column_to_rownames('Group'))
# Import taxonomy table 
taxtab = read_tsv('data/mothur/final.opti_mcc.0.03.cons.taxonomy') %>%
  filter(OTU %in% otus) %>%
  select("OTU", "Taxonomy") %>%
  rename(name = OTU, taxonomy = Taxonomy) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";")

# Import metadata 
metadataS = read.csv('data/metadata.csv', sep=',') %>% 
  column_to_rownames('sporobiota') %>%
  select(-microbiota) %>%
  rownames_to_column('samples') %>%
  add_column(biota = 'Sporobiota')

metadataM = read.delim('data/metadata.csv', sep=',') %>%
  column_to_rownames('microbiota') %>%
  select(-sporobiota) %>%
  rownames_to_column('samples') %>%
  add_column(biota = 'Microbiota')

sample_names = shared$Group
metadata = rbind(metadataM, metadataS) %>%
  filter(samples %in% sample_names) %>%
  as_tibble() %>%
  mutate(person=sub("^(\\D).*$", "\\1", .$original_sample)) %>%
  mutate(date=ymd(date)) 

metadata[metadata == ''] = NA
metadataM = filter(metadata, biota == 'Microbiota')
metadataS = filter(metadata, biota == 'Sporobiota')

# Import calculcations of alpha divrsity (subsampled to 100 000 reads per sample)
alpha = read.delim('data/mothur/final.opti_mcc.0.03.pick.groups.ave-std.summary') %>%
  filter(method=='ave')

# Add metadata
alpha_meta = alpha %>% left_join(metadata, by=join_by('group' == 'samples'))

ggplot(alpha_meta, aes(x=person, y=sobs, color=person)) +
  geom_boxplot() +
  geom_point() +
  labs(x='Person',y='Oberved OTUs') +
  ylim(0,1250) +
  facet_wrap(~biota, nrow=2) +
  theme_bw()
ggsave('plots/mothur/alpha_boxplot.png', dpi=600)

# Changes in observed alpha diversity in a person through time
alpha_meta_tmp = mutate(alpha_meta, person2=person)

ggplot(alpha_meta_tmp, aes(x=day, y=sobs)) +
  geom_line(data=alpha_meta_tmp %>% 
              dplyr::select(-person) %>% 
              filter(biota == 'Microbiota'), 
            aes(group=person2),  color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=alpha_meta_tmp %>% 
              dplyr::select(-person) %>% 
              filter(biota == 'Sporobiota'), 
            aes(group=person2),  color= cols, linewidth=0.5, alpha=0.5) + 
  geom_line(data=alpha_meta_tmp %>% 
              filter(biota == 'Microbiota'),
            aes(color=person), color= colm, linewidth=1.2 )+
  geom_line(data=alpha_meta_tmp %>% 
              filter(biota == 'Sporobiota'),
            aes(color=person), color=cols, linewidth=1.2 )+
  geom_text_repel(aes(label= group), 
                  size=3, max.overlaps = 20) +
  ylim(0, 1250) +
  geom_point(aes(color=event14_type), size=2) +
  scale_color_manual(values= col_extreme, na.value='#B8B9B6') +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Observed OTUs', color= 'Type of extreme event')
ggsave('plots/mothur/alpha_throught_time_Extreme.png', dpi = 600)

# Correlations with alpha diversity 

#Filter each table for Microbiota and sporobiota 
alphaM = filter(alpha_meta, biota == 'Microbiota')
alphaS = filter(alpha_meta, biota == 'Sporobiota')
# Make 1 table for ploting 
corr=left_join(alphaM, alphaS, by=c('original_sample', 'sample_type',
'person', 'sex', 'age', 'height', 'weight', 'place', 'diet', 'alergy', 'food_supplements', 'general_level_activity', 'smoking', 'antibiotics_in_3weeks', 'cecum', 'gallstones', 'operation',
'hospitalization','infection', 'disease', 'which', 'digestion', 'household_members', 'living', 'environment', 'date', 'day', 'time_point', 'diet14', 'antibiotics14', 'antibiotics14_type',
'antibiotics14_duration', 'probiotics14','probiotics14_type', 'prebiotics14', 'prebiotics14_type', 'medication14', 'medication14_type', 'supplements14', 'supplements14_type', 'stress',
'extreme_event14', 'event14_type', 'event14_what','moderate14', 'active14', 'dentist14', 'dentist14_date', 'operation14', 'operation14_date', 'vaccination14', 'vaccination14_type',
'bristol'))

ggplot(corr, aes(x=chao.x, y=chao.y)) +
  geom_point() +
  geom_smooth()

# Correlation betwwen chao and other variables 
corr_data = corr %>%
  select(chao.x, chao.y, person, day, moderate14, active14, stress, bristol)
library(GGally)
ggpairs(corr_data, ggplot2::aes(colour=person))

ggsave('plots/mothur/correlations_ggpairs.png', dpi=600)

# Unique OTUs accumulation curve (how many new unique OTUs were acqured each day)
# Join OTUtable with metadata 
otuPA = otu_rare %>% column_to_rownames('Group')
otuPA[otuPA > 0] = 1
otuPA = rownames_to_column(otuPA, 'Group')
otu_meta <- left_join(otuPA, metadata, by=join_by('Group' == 'samples'))

# Graph that takes into account that some OTUs are present than not and than back again - so all that are completely new
uniqueM = otu_meta %>% 
  filter(biota == 'Microbiota') %>%
    # Select only the columns I need and trasform into longer format
  select(person, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = 3:3002) %>%
  # Group the dataframe by person and otu (OTUs)
  group_by(person, name) %>% 
  # Arrange by day
  arrange(day, .by_group = TRUE) %>% 
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
   # If otu_sum is 1 or more than 1, that means that OTU was present on this day and days before
  # If otu_sum is more than 1, it means it was present in the provious days, so turn that into 0 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  group_by(person, day) %>%
  # count unique OTUs in 1 day, for each person
  summarise(., unique= sum(otu_unique)) %>%
  mutate(biota = 'Microbiota')

uniqueS =otu_meta %>% 
  filter(biota == 'Sporobiota') %>%
    # Select only the columns I need and trasform into longer format
  select(person, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = 3:3002) %>%
  # Group the dataframe by person and otu (OTUs)
  group_by(person, name) %>% 
  # Arrange by day
  arrange(day, .by_group = TRUE) %>% 
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
   # If otu_sum is 1 or more than 1, that means that OTU was present on this day and days before
  # If otu_sum is more than 1, it means it was present in the provious days, so turn that into 0 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  group_by(person, day) %>%
  # count unique OTUs in 1 day, for each person
  summarise(., unique= sum(otu_unique)) %>%
  mutate(biota = 'Sporobiota')

# Check to see if the sum of all OTUs of one person does not exceed the total number of OTUs present in my dataset
uniqueM %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Microbiota') %>% 
  rbind(uniqueS %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Sporobiota')) 

# rbind uniqueM and S 
unique_otus = rbind(uniqueM, uniqueS)

ggplot(unique_otus, aes(x=day, y=unique)) +
  geom_point(aes(color=person), linewidth=1) +
  scale_y_log10() +
  xlim(1,200) +
  #ylim(0, 400) +
  geom_smooth() +
  facet_wrap(~biota) +
  labs(x='Day', y='log(Accumulation of new unique OTUs per time point)', color='Individual') +
  theme_bw()

ggsave('plots/mothur/accumulation_unique_log_point_excl_first.png', dpi=600)

# If I don't separate microbiota and sporobiota 
unique_all = otu_meta %>% 
  select(person, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = 3:3002) %>%
  group_by(person, name) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  group_by(person, day) %>% 
  summarise(., unique= sum(otu_unique))

#unique_all %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Both') %>% 
  #rbind(uniqueM %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Microbiota')) %>% 
uniqueM %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Microbiota') %>%
  left_join(uniqueS %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Sporobiota'), by='person') %>%
  ggplot() +
  geom_segment(aes(x=person, xend=person,y=sum.x, yend=sum.y ), color='grey') +
  geom_point(aes(x=person, y=sum.x), color=colm, size=3) +
  geom_point(aes(x=person, y=sum.y), color=cols, size=3) +
  theme_bw() +
  labs(y='Individual', x='Sum of unique OTUs through duration of the study') +
  ylim(0,2200)
  
ggsave('plots/mothur/unique_otus_perPerson.png', dpi=600)

# What is a persons average new OTU accumulation
unique_all %>%
  filter(day != '0') %>%
  group_by(person) %>%
  summarise(average_new_otu = mean(unique), sd=sd(unique))

# What is the taxonomic determination of the OTUs that are new in later time-points
unique_phylum = otu_meta %>% 
  select(person, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = 3:3002) %>%
  group_by(person, name) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  filter(day != 0) %>%
  filter(otu_unique != 0) %>%
  left_join(taxtab, by=join_by('name' == 'otu')) %>%
  group_by(Phylum) %>%
  summarise(unique_phylum=sum(otu_unique)) %>%
  filter(unique_phylum > 4)

unique_phylum %>%
  ggplot(aes(x=unique_phylum, y=Phylum, fill=Phylum)) +
  geom_bar(stat = 'identity') +
  labs(x='log(Unique acquired OTUs from each phylum)') +
  scale_x_log10() +
  #facet_wrap(~person) +
  theme_bw()
ggsave('plots/mothur/taxonomy_uniques_log.png', dpi=600)

# What is the total abundance of OTUs that are newly acquired in that sample?
otu_meta %>% 
  #filter(biota == 'Microbiota') %>%
  select(person, day, Group, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = 4:3003) %>%
  group_by(person, name) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  left_join(otu_rare %>% pivot_longer(names_to = 'name', values_to = 'value', cols=2:3001)) %>%
  # remove day 0, o get newly acquired OTUs
  filter(day != 0) %>%
  # removes the ones that are not seen later
  filter(otu_unique != 0) %>%
  # Sum the total abundance of this OTUs
  group_by(name) %>%
  summarise(value=sum(value)) %>%
  ungroup() %>%
  ggplot(aes(x=value)) +
  geom_histogram(binwidth = 5) +
  #scale_x_log10()+
  coord_cartesian(xlim = c(0,200)) +
  labs(x='log(Abundance)', y='Number of OTUs') +
  theme_bw()
ggsave('plots/mothur/uniqueOTU_abundance.png', dpi=600)

############# I will come bactk to this 
# Transient or stationary OTUs in microbiota and sporobiota. 
common_otus = otu_rare %>% pivot_longer(names_to = 'name', values_to = 'value', cols=2:3001) %>%
  group_by(name) %>%
  summarize(value=sum(value)) %>%
  filter(value > 10000) %>%
  pull(name)

otu_rel = decostand(otu_rare %>% column_to_rownames('Group'), method='total', MARGIN = 1)
rowSums(otu_rel)

otu_rel %>%
  rownames_to_column('Group') %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols=2:3001) %>%
  left_join(metadata, by=join_by('Group' == 'samples')) %>%
  group_by(name) %>%
  mutate(relabund=log(value/sum(value))) %>%
  filter(name %in% common_otus) %>%
  ungroup() %>%
  # Should only choose OTUs that are present in every time point of all individuals? 
  ggplot(aes(name, time_point, fill=relabund)) +
  geom_tile() +
  facet_grid(biota~person)
# This does not work as a graphical representation
##############

# Residency of OTUS (always present in 1 person vs always present in all individuals?)
always_personM = otu_rare %>% left_join(metadata, by=join_by('Group' == 'samples')) %>% 
  # filter only microbiota samples and regular samples without extreme
  filter(biota== 'Microbiota' & sample_type == 'regular') %>%
  select(Group, person, day, date, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'count', cols = starts_with('Otu')) %>%
  # Make column PA presence/absence
  mutate(PA =ifelse(count > 0, 1, 0)) %>%
  # Group and arange by day and count all the presence of each OTU
  group_by(person, name) %>% 
  arrange(date, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA)) %>%
  ungroup() %>% 
  # Filter otu_sum with 12 counts of PA, which means they were present in all 12 time points! 
  filter(otu_sum >= 12) %>%
  left_join(taxtab, by='name') 

ggplot(always_personM, aes(y=Phylum, fill= Phylum)) +
  geom_bar(stat='count') +
  facet_grid(~person)


# Import beta diverity Bray Curtis with 1000 resaplings from mothur 
readDist = function(phylip_file) {

 # Read the first line of the phylip file to find out how many sequences/samples it contains
    temp_connection = file(phylip_file, 'r')
    len = readLines(temp_connection, n=1)
    len = as.numeric(len)
    close(temp_connection)
    
 
    phylip_data = read.table(phylip_file, fill=T, row.names=1, skip=1, col.names=1:len)
    phylip_matrix = as.dist(phylip_data)
    return(phylip_matrix)
}
dist = readDist('data/mothur/final.opti_mcc.0.03.pick.braycurtis.0.03.lt.ave.dist') 

# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
distM = dist %>% as.matrix() %>% 
  as_tibble(rownames= 'sample') %>%
  pivot_longer(-sample) %>%
  # Remove the distances of the same sample 
  filter(sample != name) %>%
  # Add metadata
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  # Filter so that I have only inter-person comparisons!
  filter(person.x == person.y & biota.x == 'Microbiota' & biota.y == 'Microbiota') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(diff, person.x) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup() %>%
  mutate(biota = 'Microbiota')

# Same as above, but for sporobiota
distS =  dist %>% as.matrix() %>% 
  as_tibble(rownames= 'sample') %>%
  pivot_longer(-sample) %>%
  filter(sample != name) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  filter(person.x == person.y & biota.x == 'Sporobiota' & biota.y == 'Sporobiota') %>%
  mutate(diff=abs(date.x-date.y)) %>%
  group_by(diff, person.x) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup() %>%
  mutate(biota = 'Sporobiota')

rbind(distM, distS) %>% ggplot(aes(x=diff, y=median, color=biota)) +
  geom_point() +
  scale_color_manual(values=c(colm, cols)) +
  geom_smooth(aes(group=biota)) +
  labs(x='Days between sampling points', y='Average Bray-Curtis distance', color='Type of biota') +
  theme_bw()
ggsave('plots/mothur/beta_throughTime.png', dpi=600)
