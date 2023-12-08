# Library load 
library(tidyverse)
library(vegan)
library(ggrepel)
library(reshape2)
# Determine seed, theme and colors used in the analysis 
set.seed(1996)
# Colors for microbiota (colm) and sporobiota (cols)
colm = c('#47B66A')
cols = c('#D9A534')
colsm = c('#47B66A', '#D9A534')

# Data
# Rarefied data
otu_rare = readRDS('data/otu_rare.RDS')
otu_rel <- apply(decostand(otu_rare %>% column_to_rownames('Group') %>% t(), method="total", MARGIN=2),1, mean)
# Import taxonomy table 
taxtab = readRDS('data/taxtab.RDS')
# Import metadata 
metadata = readRDS('data/metadata.RDS')

# How many/which OTUs are present in the first sample and than remain in the rest of the samples, how long to they stay in the microbiota ? 
time_otu = otu_rare %>% pivot_longer(names_to = 'name', values_to = 'value', cols = starts_with('Otu')) %>% left_join(metadata, by=join_by('Group' == 'samples')) %>%
  filter(biota == 'Microbiota') %>%
  # Select only the columns I need and trasform into longer format
  select(person, date, name, value, starts_with('Otu')) %>%
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
        axis.title.x = '')