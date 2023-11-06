# Analysis of V3V4 16S rRNA data from longitudinal cohort of healthy adult males. 
# OTUs and Phylotypes were constructed in mothur code availabile on github repo mothur.script 

# Library load 
library(tidyverse)
library(vegan)
library(ggrepel)
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
otutab= read.delim('data/mothur/final.opti_mcc.0.03.pick.shared')
otus= names(otutab)[-c(1:3)]

taxtab = read_tsv('data/mothur/final.opti_mcc.0.03.cons.taxonomy') %>%
  filter(OTU %in% otus) %>%
  select("OTU", "Taxonomy") %>%
  rename_all(tolower) %>%
  mutate(taxonomy = str_replace_all(taxonomy, "\\(\\d+\\)", ""),
         taxonomy = str_replace(taxonomy, ";$", "")) %>%
  separate(taxonomy,
           into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
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

sample_names = otutab$Group
metadata = rbind(metadataM, metadataS) %>%
  filter(samples %in% sample_names) %>%
  as_tibble() %>%
  mutate(person=sub("^(\\D).*$", "\\1", .$original_sample)) %>%
  mutate(date=ymd(date)) 

metadata[metadata == ''] = NA
metadataM = filter(metadata, biota == 'microbiota')
metadataS = filter(metadata, biota == 'sporobiota')

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
  #geom_point(aes(color=event14_type), size=2) +
  #scale_color_manual(values= col_extreme, na.value='#B8B9B6') +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Observed OTUs', color= 'Type of extreme event')
ggsave('plots/mothur/alpha_throught_time.png', dpi = 600)

# Number of OTUs that have increased or decreased by day in each person
# Calculate how sobs (observed number of OTUs increases or decreases in the next day)
# 
tmp1 = alpha_meta %>%
  # Filter only samples from microbiota
  filter(biota == 'Microbiota')%>%
  # Group_by person and arrange samples by person
  group_by(person) %>%
  arrange(date, .by_group = TRUE) %>%
  # Calculate how sobs observed number of OTUs increases or decreases in the next day and save into sobs_diff
  mutate(sobs_diff= sobs -lag(sobs, default = first(sobs))) %>%
  # Add the value of the first sobs in each group to the column sobs_diff 
  mutate(sobs_diff = ifelse(row_number() == 1, first(sobs), sobs_diff))

# Do the same for sporobiota 
tmp2 = alpha_meta %>%
  filter(biota == 'Sporobiota')%>%
  group_by(person) %>%
  arrange(date, .by_group = TRUE) %>%
  mutate(sobs_diff= sobs - lag(sobs, default = first(sobs))) %>%
  mutate(sobs_diff = ifelse(row_number() == 1, first(sobs), sobs_diff)) 

# Combine both data.frames and turn into tibble
tmp3= rbind(tmp1, tmp2) %>%
  as_tibble()

ggplot(tmp3, aes(x=day, y=sobs_diff)) +
  geom_line(data=tmp3 %>% filter(biota == 'Microbiota'), color= colm, linewidth=1.2) +
  geom_line(data=tmp3 %>% filter(biota == 'Sporobiota'), color= cols, linewidth=1.2) +
  geom_point(aes(y=sobs_diff, color=event14_type), size=2) +
  scale_color_manual(values=col_extreme, na.value='#B8B9B6') +
  facet_wrap(~person) +
  labs(x='Day', y= 'Observed OTUs', color= 'Type of event')

ggsave('plots/mothur/increase_decrease_daily.png', dpi = 600)

# Unique OTUs accumulation curve (how many new unique OTUs were aqured each day)
# Join OTUtable with metadata 
otuPA = otu_rare
otuPA[otuPA > 0] = 1
otuPA = rownames_to_column(otuPA, 'samples')
otu_meta <- left_join(metadata, otuPA)

# Calculate the number of unique new OTUs acquired each day from all previous days
uniqueM = otu_meta %>%
  filter(biota == 'Microbiota') %>%
  select(person, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'otu', values_to = 'PA', cols = 3:765) %>%
  group_by(person) %>%
  dplyr::arrange(otu, day, .by_group = TRUE) %>%
  #filter(PA != 0) %>%  # IS this necessary?
  # Create new column new_otus is 1 if the OTU is present (PA > 0) on the current date and was not present on the previous date (lag(PA))
  mutate(new_otus = ifelse(PA > 0 & lag(PA, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  group_by(person, day) %>%
  summarise(., unique= sum(new_otus)) %>%
  mutate(., unique_old_new = cumsum(unique)) %>%
  mutate(biota = 'Microbiota')

uniqueS = otu_meta %>%
  filter(biota == 'Sporobiota') %>%
  select(person, day, starts_with('Otu')) %>%
  pivot_longer(names_to = 'otu', values_to = 'PA', cols = 3:765) %>%
  group_by(person) %>%
  dplyr::arrange(otu, day, .by_group = TRUE) %>%
  #filter(PA != 0) %>%  # IS this necessary?
  # Create new column new_otus is 1 if the OTU is present (PA > 0) on the current date and was not present on the previous date (lag(PA))
  mutate(new_otus = ifelse(PA > 0 & lag(PA, default = 0) == 0, 1, 0)) %>%
  ungroup() %>%
  group_by(person, day) %>%
  summarise(., unique= sum(new_otus)) %>%
  mutate(., unique_old_new = cumsum(unique)) %>%
  mutate(biota = 'Sporobiota')

unique_otus = rbind(uniqueM, uniqueS)

ggplot(unique_otus, aes(x=day, y=unique_old_new)) +
  geom_line(aes(color=person)) +
  ylim(0, 1000) +
  geom_smooth() +
  facet_wrap(~biota) +
  labs(x='Day', y='Accumulation of new unique OTUs') +
  theme_bw()

ggsave('plots/mothur/unique_otus_accumulation.png', dpi=600)

# The graphs and code above do not take into account the fact that some OTUs may disepear and than become apperent again becouse of the method of sequnecing 
# Graph that takes into account that some otus are present than not and than back again - so all that are completly new

# A list of OTUs in the first sample of person A 
otuA1 = otutab %>%
  as_tibble() %>%
  pivot_longer(names_to = 'otu', values_to = 'value', cols=4:766) %>%
  left_join(metadata, by=join_by('Group' == 'samples')) %>%
  select(Group, person, day, otu, value) %>%
  filter(day==0 & value > 0 & person == 'A') %>%
  select(otu) %>%
  as.list()

otuPA = otu_rare
otuPA[otuPA > 0] = 1
otuPA = rownames_to_column(otuPA, 'samples')

otutab %>% left_join(metadata, by=join_by('Group' == 'samples')) %>%
  as_tibble() %>%
  filter(biota== 'Microbiota') %>%
  select(person, day,starts_with('Otu')) %>%
  group_by(person) %>%
  arrange(day, .by_group = TRUE)

