# Library
require(tidyverse)
require(ggtext)
require(vegan)
require(ggpubr)
# Determine seed, theme and colors used in the analysis 
set.seed(1996)
# Data
# Rarefied data
otu_rare = readRDS('data/otu_rare.RDS')
otu_rel <- apply(decostand(otu_rare %>% column_to_rownames('Group') %>% t(), method="total", MARGIN=2),1, mean)
# Import taxonomy table 
taxtab = readRDS('data/taxtab.RDS')
# Import metadata 
metadata = readRDS('data/metadata.RDS')

# Colors 
colm = c('#47B66A')
cols = c('#D9A534')

# Figure 1 
# In which 'biota' was an OTU present and what is its taxonomy+relative abundance?
otutabM = otu_rare %>% filter(str_detect(Group, '^M')) %>%
  column_to_rownames('Group') %>% t()
otu_PA <- 1*((otutabM>0)==1)
otu_occM <- rowSums(otu_PA)/ncol(otu_PA) 

otutabS = otu_rare %>% filter(str_detect(Group, '^S')) %>%
  column_to_rownames('Group') %>% t()
otu_PA <- 1*((otutabS>0)==1)
otu_occS <- rowSums(otu_PA)/ncol(otu_PA) 

otu_where = data.frame(otu_occM, otu_occS) %>%
  mutate(where = ifelse(otu_occS == 0 & otu_occM > 0, 'Microbiota', 
                        ifelse(otu_occM == 0 & otu_occS > 0, 'Sporobiota', 
                               ifelse(otu_occM > 0 & otu_occS > 0, 'Both', 'None'))))
   
tax_occ = rownames_to_column(otu_where, 'name') %>%
  left_join(taxtab, by='name') %>%
  mutate(Phylum = str_replace(Phylum, 
                             '(.*)_unclassified', 'Unclassified *\\1*'), 
         Phylum = str_replace(Phylum, 
                             '^(\\S*)$', '*\\1*'))

ggplot(tax_occ, aes(x=where, fill=Phylum)) +
  geom_bar(stat='count') +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        legend.title = element_blank(), 
        legend.text = element_markdown(), 
        legend.key.size = unit(18, 'pt')) +
  labs(y='Number of OTUs')
ggsave('plots/where_taxonomy.png', dpi=600)
  
# This does not work as I 
tax_occ_rel = data.frame(tax_occ, otu_rel)
ggplot(tax_occ_rel, aes(x=where, y=otu_rel, fill=where)) +
  geom_boxplot() +
  scale_y_log10() +
  scale_fill_manual(values = c('#1dbeda', colm , cols)) +
  theme_bw() +
  theme(axis.title.x = element_blank(), 
        legend.position = 'none') +
  labs(y='Relative abundance of OTUs')
ggsave('plots/where_relabund.png', dpi=600) 

# Figure 2 
otuPA = otu_rare %>% column_to_rownames('Group')
otuPA[otuPA > 0] = 1
otuPA = rownames_to_column(otuPA, 'Group')
otuPA_meta <- left_join(otuPA, metadata, by=join_by('Group' == 'samples'))

# Graph that takes into account that some OTUs are present than not and than back again - so all that are completely new
new_otus = otuPA_meta %>% 
  filter(time_point < 13) %>%
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
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  #filter(date != min(date)) %>%
  ungroup() %>%
  group_by(person, time_point, biota) %>%
  # percentage of new OTUs in 1 day, for each person
  summarise(., new= sum(new_otu)/sum(PA)*100) %>%
  ungroup() %>%
  group_by(time_point, biota) %>%
  mutate(mean =median(new), 
         sd=sd(new))

ggplot(new_otus, aes(x=time_point)) +
  geom_point(aes(y=new, color=biota),size=2) +
  geom_line(mapping=aes(y=mean, color=biota), linewidth=1) +
  #geom_ribbon(mapping= aes(ymin=mean-sd, ymax=mean +sd),  fill ='grey',  alpha=.2) +
  scale_color_manual(values = c(colm, cols)) +
  scale_x_continuous(position ='top', breaks = seq(1, 14)) +
  labs(x='Sampling point', y='New OTUs (%)', color='Type of sample') +
  theme_bw()
ggsave('plots/newOTUs_percent.png', dpi=600)

# What is the taxonomic determination and relative abundance of the OTUs that are new in later time-points?
new_tax = otuPA_meta %>% 
  filter(time_point < 13) %>%
  select(person, time_point, biota, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  group_by(person, name, biota) %>% 
  arrange(time_point, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  filter(time_point != 1) %>%
  filter(new_otu == 1) %>%
  left_join(taxtab, by='name') %>%
  left_join(data.frame(otu_rel) %>% rownames_to_column('name'), by='name') %>%
  group_by(Phylum, biota) %>%
  summarise(sum_relabund =sum(otu_rel)*100) %>%
  mutate(Phylum = if_else(sum_relabund < 0.005, 'Other (< 0.5%)', Phylum), 
         Phylum = str_replace(Phylum, 
                             '(.*)_unclassified', 'Unclassified \\1'))

ggplot(new_tax, aes(x=Phylum, y=sum_relabund, fill=Phylum)) +
  geom_col() +
  scale_y_log10() +
  labs(x='', y='Mean relative abundance (%)', fill='') +
  coord_flip() +
  theme_bw() +
    theme(legend.position = 'none', 
          axis.text.x = element_markdown(), 
          legend.key.size = unit(18, 'pt'))
ggsave('plots/new_taxonomy.png', dpi=600)

# What is the total abundance of OTUs that are newly acquired in that sample?
new_rel = otuPA_meta %>% 
  filter(time_point < 13) %>%
  select(person, time_point, biota, starts_with('Otu')) %>%
  pivot_longer(names_to = 'name', values_to = 'PA', cols = starts_with('Otu')) %>%
  group_by(person, name, biota) %>% 
  arrange(time_point, .by_group = TRUE) %>% 
  mutate(otu_sum = cumsum(PA), 
         new_otu = ifelse(otu_sum == 1 & lag(otu_sum, default = 0) == 0, 1, 0)) %>%
  filter(time_point != 1) %>%
  filter(new_otu == 1) %>%
  left_join(data.frame(otu_rel) %>% rownames_to_column('name'), by='name') %>%
  # Sum the total abundance of this OTUs
  group_by(name, otu_rel) %>%
  summarise(sum_otus=sum(new_otu)) %>%
  ungroup()

ggplot(new_rel, aes(x=otu_rel)) +
  geom_histogram() +
  scale_x_log10()+
  #coord_cartesian(xlim = c(0,200)) +
  labs(x='log10(Relative abundance)', y='Number of OTUs') +
  theme_bw()
ggsave('plots/new_relabund_log10.png', dpi=600)

# Figure 3 
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

# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
# With Jaccard distances - which take into account only presence/absence! Not abundance such as Bray-Curtis 
distj = readDist('data/mothur/final.opti_mcc.0.03.pick.jclass.0.03.lt.ave.dist') 

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

distj_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c(colm, cols)) +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="Jaccard distance", x="", fill='Type of sample') +
  theme_bw()
ggsave('plots/violin_jaccard.png', dpi=600)

# Figure 4
# Occupancy plot (code by Shade and Stopnisek)
otutab = otu_rare %>% filter(str_detect(Group, '^M')) %>%
  column_to_rownames('Group') %>% t() 

otu_PA <- 1*((otutab>0)==1)                                              # presence-absence data 
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                                # occupancy calculation
otu_rel <- apply(decostand(otutab, method="total", MARGIN=2),1, mean)     # mean relative abundance
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%           # combining occupancy and abundance data frame
  rownames_to_column('name') %>%
  left_join(taxtab, by='name')

# Occupancy abundance plot:
plot1 = ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, fill=colm) +
  labs(x="log10(mean relative abundance)", y="Occupancy") +
  theme_bw()

# For sporobiota
otutab = otu_rare %>% filter(str_detect(Group, '^S')) %>%
  column_to_rownames('Group') %>% t() 

otu_PA <- 1*((otutab>0)==1)                                             
otu_occ <- rowSums(otu_PA)/ncol(otu_PA)                               
otu_rel <- apply(decostand(otutab, method="total", MARGIN=2),1, mean)     
occ_abun <- data.frame(otu_occ=otu_occ, otu_rel=otu_rel) %>%        
  rownames_to_column('name') %>%
  left_join(taxtab, by='name')

# Occupancy abundance plot:
plot2 = ggplot(data=occ_abun, aes(x=log10(otu_rel), y=otu_occ)) +
  geom_point(pch=21, fill=cols) +
  labs(x="log10(mean relative abundance)", y="Occupancy") +
  theme_bw()
  #facet_wrap(~Phylum, ncol=3, nrow = 4)

ggarrange(plot1, plot2, 
          ncol=1)
ggsave('plots/occupancy_mean_relabund.png', dpi=600)
