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
  filter(total_sample > 50000) %>%
  ungroup() %>% 
  select(-total_sample)

otutab = dcast(shared_pre, Group ~ name, value.var = 'value') %>%
  column_to_rownames('Group')

# Rarefy the data once -as we sequenced so deep that for the first analysis this is not crutial! 
otu_rare = rrarefy(otutab, sample=50000) %>% as.data.frame() %>%
  rownames_to_column('Group') %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols=starts_with('Otu')) %>% 
  group_by(name) %>%
  mutate(total_otu=sum(value)) %>%
  filter(total_otu > 4) %>%
  ungroup() %>%
  select(-total_otu) %>%
  pivot_wider(names_from = 'name', values_from = 'value')
saveRDS(otu_rare, 'data/otu_rare.RDS')
# Extract OTUs that are present in otu_rare
otus=colnames(otu_rare %>% column_to_rownames('Group'))

rm(otutab)
rm(shared_pre)

# Make relative abundance tables 
otu_rel = decostand(otu_rare %>% column_to_rownames('Group'), method='total', MARGIN = 1)
rowSums(otu_rel)

otu_rel_long = rownames_to_column(otu_rel, 'Group') %>% pivot_longer(starts_with('Otu')) %>% left_join(metadata, by=join_by('Group' == 'samples'))
# Import taxonomy table 
taxtab = read_tsv('data/mothur/final.opti_mcc.0.03.cons.taxonomy') %>%
  filter(OTU %in% otus) %>%
  select("OTU", "Taxonomy") %>%
  rename(name = OTU, taxonomy = Taxonomy) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\\\|\\\"|\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
           sep=";")
saveRDS(taxtab, 'data/taxtab.RDS')
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

metadata = rbind(metadataM, metadataS) %>%
  as_tibble() %>%
  mutate(person=sub("^(\\D).*$", "\\1", .$original_sample)) %>%
  mutate(date=ymd(date)) 

metadata[metadata == ''] = NA
metadataM = filter(metadata, biota == 'Microbiota')
metadataS = filter(metadata, biota == 'Sporobiota')
saveRDS(metadata, 'data/metadata.RDS')
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

# What percentage of sporobiota is also in microbiota ? 
# How much of sporobiota is found in microbiota in each sample ? 
# in the code for SMNS


# Aleksander's code for determining how relative abundant were OTUs that were present in 1, 2, 3 etc time-points and how many of them are there for each point. 
merged = left_join(metadata, otu_rel %>% rownames_to_column('samples'), by='samples') %>%
  filter(sample_type == 'regular') %>%
  column_to_rownames('samples') %>%
  select(person, biota, starts_with('Otu'))

final <- data.frame()
for (persona in unique(merged$person)) {
  for (bioti in unique(merged$biota)) {
  merged_sub <- merged[merged$biota == bioti & merged$person == persona,]
  merged_sub <- as.data.frame(t(merged_sub[, 3:2996]))
  
  merged_sub2 <- merged_sub
  merged_sub2[merged_sub2>0] <- 1
  merged_sub2$prevalence <- rowSums(merged_sub2)
  merged_sub$person <- persona
  merged_sub$biota <- bioti
  merged_sub$name <- rownames(merged_sub)
  merged_sub2$name <- rownames(merged_sub2)

  melt1 <- melt(merged_sub2, id.vars = c('prevalence', 'name'))
  melt2 <- melt(merged_sub, id.vars = c('person', 'biota', 'name'))
  melt2$prevalence <- melt1$prevalence
    
  final <- rbind(final, melt2)
  }
}

p1 =ggplot(na.omit(final[final$value != 0,]), aes(x = as.factor(prevalence), y = value, color=biota)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_color_manual(values = c(colm, cols)) +
  labs(x='Number of times an OTU was observed', y= 'Relative abundandance of these OTUs', color='Type of biota') +
  theme_bw() +
  coord_flip()
ggsave('plots/mothur/prevalence_relabund.png', dpi=600)

final$count <- 1/12  
# Group_by biota, person and prevalence and summarize count (because each OTU could be present in 12 time points 1/12)
final_agg <- aggregate(count ~ biota + person + prevalence, data = final, FUN = sum)

final_agg_mean = filter(final_agg, prevalence != 0) %>% group_by(biota, prevalence) %>%
  summarise(mean = median(count), sd=sd(count))

p2 = ggplot(final_agg[final_agg$prevalence != 0,], aes(x = prevalence, y = count, color=biota)) +
  geom_point() +
  geom_line(final_agg_mean, mapping=aes(y=mean, color=biota)) +
  scale_color_manual(values = c(colm, cols)) +
  scale_x_continuous(breaks = seq(0,12, by=1))+
  labs(x='Number of times an OTU was observed', y= 'Count of OTUs observed', color='Type of biota') +
  theme_bw() +
  coord_flip()
ggsave('plots/mothur/prevalence_count.png', dpi=600)

# combine plots 
require(ggpubr)
ggarrange(p1, p2, 
          common.legend = TRUE)
ggsave('plots/mothur/prevalence_both.png', dpi=600)
###### End of Aleksander's code 


# Unique OTUs accumulation curve (how many new unique OTUs were acqured each day)
# Join OTUtable with metadata 
otuPA = otu_rare %>% column_to_rownames('Group')
otuPA[otuPA > 0] = 1
otuPA = rownames_to_column(otuPA, 'Group')
otuPA_meta <- left_join(otuPA, metadata, by=join_by('Group' == 'samples'))

# Graph that takes into account that some OTUs are present than not and than back again - so all that are completely new
unique = otuPA_meta %>% 
    # Select only the columns I need and trasform into longer format
  select(person, time_point, biota, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  # Group the dataframe by person and otu (OTUs)
  group_by(person, name, biota) %>% 
  # Arrange by day
  arrange(time_point, .by_group = TRUE) %>% 
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
   # If otu_sum is 1 or more than 1, that means that OTU was present on this day and days before
  # If otu_sum is more than 1, it means it was present in the provious days, so turn that into 0 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  #filter(date != min(date)) %>%
  ungroup() %>%
  group_by(person, time_point, biota) %>%
  # count unique OTUs in 1 day, for each person
  summarise(., unique= sum(otu_unique)) %>%
  ungroup() %>%
  group_by(time_point, biota) %>%
  mutate(percentages_unique = unique/sum(unique), 
          median =median(unique), 
          sd=sd(unique))

# Check to see if the sum of all OTUs of one person does not exceed the total number of OTUs present in my dataset
unique%>% group_by(person) %>% filter(biota == 'Microbiota') %>% summarise(sum=sum(unique))
unique%>% group_by(person) %>% filter(biota == 'Sporobiota') %>% summarise(sum=sum(unique))

ggplot(unique, aes(x=time_point)) +
  geom_point(aes(y=unique, color=person), size=2) +
  geom_line(mapping=aes(y=median), color='gray', size=1.5) +
  geom_ribbon(mapping= aes(ymin=median-sd, ymax=median +sd), fill='grey', alpha=.3) +
  scale_y_log10() +
  facet_wrap(~biota) +
  labs(x='Day', y='Number of new OTUs at each time point', color='Individual') +
  theme_bw()
ggsave('plots/mothur/newOTUs_log.png', dpi=600)

# If I don't separate microbiota and sporobiota 
unique_all = otuPA_meta %>% 
  select(person, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  group_by(person, name) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  group_by(person, day) %>% 
  summarise(., unique= sum(otu_unique))
# What is a persons average new OTU accumulation
unique_all %>%
  filter(day != '0') %>%
  group_by(person) %>%
  summarise(average_new_otu = mean(unique), sd=sd(unique))

# #unique_all %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Both') %>% 
#   #rbind(uniqueM %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Microbiota')) %>% 
# uniqueM %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Microbiota') %>%
#   left_join(uniqueS %>% group_by(person) %>% summarise(sum=sum(unique)) %>% mutate(biota='Sporobiota'), by='person') %>%
#   ggplot() +
#   geom_segment(aes(x=person, xend=person,y=sum.x, yend=sum.y ), color='grey') +
#   geom_point(aes(x=person, y=sum.x), color=colm, size=3) +
#   geom_point(aes(x=person, y=sum.y), color=cols, size=3) +
#   theme_bw() +
#   labs(y='Individual', x='Sum of unique OTUs through duration of the study') +
#   ylim(0,2200)
#   
# ggsave('plots/mothur/unique_otus_perPerson.png', dpi=600)

# What is the taxonomic determination of the OTUs that are new in later time-points
unique_phylum = otuPA_meta %>% 
  select(person, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  group_by(person, name) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  filter(day != 0) %>%
  filter(otu_unique != 0) %>%
  left_join(taxtab, by='name') %>%
  group_by(Phylum) %>%
  summarise(unique_phylum=sum(otu_unique)) %>%
  filter(unique_phylum > 4)

unique_phylum %>%
  ggplot(aes(x=unique_phylum, y=Phylum, fill=Phylum)) +
  geom_bar(stat = 'identity') +
  labs(x='log(Unique acquired OTUs from each phylum)') +
  #scale_x_log10() +
  #facet_wrap(~person) +
  theme_bw()
ggsave('plots/mothur/taxonomy_uniques.png', dpi=600)

# What is the total abundance of OTUs that are newly acquired in that sample?
otuPA_meta %>% 
  #filter(biota == 'Microbiota') %>%
  select(person, day, Group, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  group_by(person, name) %>% 
  arrange(day, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         otu_unique = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  left_join(otu_rare %>% pivot_longer(names_to = 'name', values_to = 'value', cols=starts_with('Otu'))) %>%
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

################################## THIS PART does not work completly as it should, after SMNS
# How many/which OTUs are present in the first sample and than remain in the rest of the samples, how long to they stay in the microbiota ? 
time_otu = otuR_meta %>% 
  filter(biota == 'Microbiota') %>%
    # Select only the columns I need and trasform into longer format
  select(person, date, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols = starts_with('Otu')) %>%
  # Group the dataframe by person and otu (OTUs)
  group_by(person, name) %>% 
  # Arrange by day
  arrange(date, .by_group = TRUE) %>% 
  # Create new column otu_sum is 1 if the OTU is present (PA > 0) on the current day and was not present on any of the previous days
   # If otu_sum is 1 or more than 1, that means that OTU was present on this day and days before
  # If otu_sum is more than 1, it means it was present in the provious days, so turn that into 0 
  mutate(otu_sum = cumsum(value),
         persister = if_else(date == min(date) & value > 0 , 1, 
                             if_else(otu_sum > lag(otu_sum) & lag(otu_sum, default = 0) > 0, 1, 0)), 
         situation = if_else(date == min(date) & value > 0, 'First observation',
                             if_else(persister == '1' & lag(persister) == '1', 'Persister', 
                                     if_else(persister == '1' & lag(persister) == '0', 'New observation',
                                             if_else(persister == '0' & lag(persister) == '1', 'Gone', 
                                                     if_else(persister == '0' & lag(persister) == '0', 'Missing', 'NA')))))) %>% 
  ungroup()

time_otu  %>% mutate(occurance= if_else(situation != 'NA', 1, 0)) %>% 
  #group_by(person, name, situation) %>%
  #summarize(no_situations= sum(occurance)) %>% 
  filter(situation != 'NA')%>%
  #ungroup() %>%
  ggplot(aes(x=as.factor(date), y=name, fill=situation, na.rm=TRUE)) +
  geom_tile() +
  facet_wrap(~person, nrow=3, ncol=3, scales = 'free') +
  theme(axis.text = element_blank(), 
        axis.ticks = element_blank(), 
        axis.title.x = )
ggsave('plots/mothur/time_ofOTUs.png', dpi=600)
#######################################################################################

#### TRANSIENT AND STATIONARY OTUs
# Transient or stationary OTUs in microbiota and sporobiota.
common_otus = otu_rare %>% pivot_longer(names_to = 'name', values_to = 'value', cols=starts_with('Otu')) %>%
  group_by(name) %>%
  summarize(value=sum(value)) %>%
  filter(value > 10000) %>%
  pull(name)

otu_rel %>%
  rownames_to_column('Group') %>%
  pivot_longer(names_to = 'name', values_to = 'value', cols=starts_with('Otu')) %>%
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

##################
# BETA DIVERSITY #
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
distbc = readDist('data/mothur/final.opti_mcc.0.03.pick.braycurtis.0.03.lt.ave.dist') 

# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
distM = distbc %>% as.matrix() %>% 
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
distS =  distbc %>% as.matrix() %>% 
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
ggsave('plots/mothur/beta_BrayCurtis_throughTime.png', dpi=600)

distbc %>% as.matrix() %>% 
  as_tibble(rownames= 'sample') %>%
  pivot_longer(-sample) %>%
  filter(sample != name) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  filter(person.x == person.y) %>%
  group_by(person.x) %>%
  arrange(as.Date(date, '%Y-%m-%d'), .by_group = TRUE) %>%
  mutate(diff=value -lag(value, default = 0)) %>%
  ggplot(aes(x=date, y=diff, color=person)) +
  geom_line() +
  facet_grid(~biota.x)


# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
# With Jaccard distances - which take into account only presence/absence! Not abundance such as Bray-Curtis 
distj = readDist('data/mothur/final.opti_mcc.0.03.pick.jclass.0.03.lt.ave.dist') 

# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
distM = distj %>% as.matrix() %>% 
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
distS =  distj %>% as.matrix() %>% 
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
  labs(x='Days between sampling points', y='Median Jaccard distance', color='Type of biota') +
  theme_bw()
ggsave('plots/mothur/beta_Jaccard_throughTime.png', dpi=600)

# Classical ordination plot - Bray-Curtis 
# Calculate NMDS ordinates
nmds_bc = metaMDS(distbc)
nmds_positions= as.data.frame(scores(nmds_bc, display='sites')) %>%
  rownames_to_column('Group')
# Join with metadata
distbc_meta = nmds_positions %>% left_join(metadata, by=join_by('Group' == 'samples'))
# Ordination plot wih person/biota 
distbc_meta %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person, shape=biota)) +
  geom_point(size=3) +
  #geom_text_repel(aes(label=sample), size= 4, colour='black', max.overlaps = 20) +
  scale_size_continuous(range = c(3,6)) +
  labs(x='', y='', shape='Type of biota', color='Individual') +
  theme_bw()
ggsave('plots/mothur/nmds_braycurtis.png', dpi=600)

# Statistical testing of distances between people/biota of one person/biotas of all samples... 
# Calculate multivariate dispersions of microbiota and sporobiota 
mod = betadisper(distbc, distbc_meta$person, type= 'median')
anova(mod)
# p-value is significant, which means persons are not dispersed homogenous

# tukeys test tells us which groups differ in relation to their variances - NONE
TukeyHSD(mod)
plot(mod)
boxplot(mod)

# Tukeys test, if we separate biotas! 
mod2 = betadisper(distbc, distbc_meta$biota, type='median')
anova(mod2)
TukeyHSD(mod2)
boxplot(mod2)

# Do individuals different from each other?
adonis2(distbc~person, data=distbc_meta, method='bray', permutations = 9999) # yes; p-value= 1e-04; person explains 23,4% variability

# Distances between individual time-points of individual samples
adonis2(distbc~person*time_point, data=distbc_meta, method='bray', permutations = 9999 ) # time_points explain 0,16% % variability. And are not significant

# Distances between individuals, stress, bristol, time-points
adonis2(distbc~person*stress*bristol, data=distbc_meta, method='bray', permutations = 9999, na.action = na.exclude) # Some samples do not have data for physical activity, so na.action = na.exclude 
 # nothing is significant or explains more than 0,1 varability 

# Classical ordination plot - Jaccard
nmds_j = metaMDS(distj)
nmds_positions= as.data.frame(scores(nmds_j, display='sites')) %>%
  rownames_to_column('Group')

distj_meta = nmds_positions %>% left_join(metadata, by=join_by('Group' == 'samples'))
distj_meta %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person, shape=biota)) +
  geom_point(size=3) +
  #geom_text_repel(aes(label=sample), size= 4, colour='black', max.overlaps = 20) +
  scale_size_continuous(range = c(3,6)) +
  labs(x='', y='', shape='Type of biota', color='Individual') +
  theme_bw()
ggsave('plots/mothur/nmds_jaccard.png', dpi=600)
# Statistical testing of distances between people/biota of one person/biotas of all samples... 
# Calculate multivariate dispersions of microbiota and sporobiota 
mod = betadisper(distj, distj_meta$person, type= 'median')
anova(mod)
# p-value is significant, which means persons are not dispersed homogenous

# Tukeys test tells us which groups differ in relation to their variances - NONE
TukeyHSD(mod)
plot(mod)
boxplot(mod)
# Tukeys test, if we separate biotas! YES!
mod2 = betadisper(distj, distj_meta$biota, type='median')
anova(mod2)
TukeyHSD(mod2)
boxplot(mod2)

# Do individuals different from each other?
adonis2(distj~person, data=distj_meta, method='jaccard', permutations = 9999) # yes; p-value= 1e-04; person explains 12,5% variability

# Distances between individual time-points of individual samples
adonis2(distj~person*time_point, data=distj_meta, method='jaccard', permutations = 9999 ) # time_points explain 0,04% % variability. And are not significant!

# Distances between individuals, stress, bristol, time-points
adonis2(distj~person*biota*stress*bristol, data=distj_meta, method='jaccard', permutations = 9999, na.action = na.exclude) # Some samples do not have data for physical activity, so na.action = na.exclude 
# nothing is significant or explains more than 0,2% variability in data 

# Violin plot of distances between samples of 1 individual vs. between individuals Bray-Curtis
# Turn distance matrix into data.frame 
distbc_df = as.data.frame(as.matrix(distbc)) %>%
  rownames_to_column('sample') %>%
  pivot_longer(-sample) %>%
  filter(sample != name) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Sporobiota' & biota.y == 'Sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both')

distbc_df %>% ggplot(aes(x=same_person, y=value, color=which_biota)) +
  #geom_boxplot(outlier.colour="black", outlier.shape=8, outlier.size=1) +
  geom_violin() +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 0.5)) +
  labs(y="Bray-Curtis distance", x="", color='Type of biota') +
  theme_bw()
ggsave('plots/mothur/violin_braycurtis.png', dpi=600)

# Jaccard distance
distj_df = as.data.frame(as.matrix(distj)) %>%
  rownames_to_column('sample') %>%
  pivot_longer(-sample) %>%
  filter(sample != name) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(metadata %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Sporobiota' & biota.y == 'Sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both')

distj_df %>% ggplot(aes(x=same_person, y=value, color=which_biota)) +
  #geom_boxplot(outlier.colour="black", outlier.shape=8, outlier.size=1) +
  geom_violin() +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  theme(axis.text.x=element_text(angle = 45, hjust = 1, vjust = 0.5)) +
  labs(y="Jaccard distance", x="", color='Type of biota') +
  theme_bw()
ggsave('plots/mothur/violin_jaccard.png', dpi=600)
