#### Does relative abundance in microbiota or sporobiota sample of the same original sample correlate ? 

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
# Import taxonomy table 
taxtab = readRDS('data/taxtab.RDS')
# Import metadata 
metadata = readRDS('data/metadata.RDS')

otu_rel <- decostand(column_to_rownames(otu_rare, 'Group') %>% t(), method="total", MARGIN=2) %>% t()
  
rel_meta = as.data.frame(otu_rel) %>%
  rownames_to_column('samples') %>% 
  left_join(metadata %>% select(samples, original_sample), by='samples') %>%
  filter(grepl('^M', samples)) %>%
  pivot_longer(names_to = 'otu', values_to = 'relabund_microbiota', cols = starts_with('Otu')) %>%
  left_join(as.data.frame(otu_rel) %>%
              rownames_to_column('samples') %>% 
              left_join(metadata %>% select(samples, original_sample), by='samples') %>%
              filter(grepl('^S', samples)) %>%
              pivot_longer(names_to = 'otu', values_to = 'relabund_sporobiota', cols = starts_with('Otu')),
            by=join_by('original_sample' == 'original_sample', 'otu' == 'otu')) %>%
  select(-samples.x, -samples.y)

rel_meta %>% left_join(metadata %>% filter(grepl('^M', samples)), by='original_sample') %>%
  left_join(taxtab, by=join_by('otu' == 'name')) %>%
  ggplot(aes(x=relabund_microbiota, y=relabund_sporobiota, color=Phylum)) +
  geom_point() +
  # xlim(0,1) +
  # ylim(0,1) 
  scale_x_log10() +
  scale_y_log10()

corr_otu = rel_meta %>% group_by(otu) %>%
  summarize(R = cor.test(relabund_microbiota, relabund_sporobiota, method = 'pearson')$estimate, 
         pvalue = cor.test(relabund_microbiota, relabund_sporobiota, method = 'pearson')$p.value) 

positive_corr = filter(corr_otu, R > 0 & pvalue < 0.05)
negative_corr =filter(corr_otu, R < 0 & pvalue < 0.05)

# Not quite sure where this is going 
rel_meta_where = rel_meta %>% mutate(where = ifelse(relabund_microbiota == 0 & relabund_sporobiota == 0, 'Not present in this sample', 
                                   ifelse(relabund_microbiota == 0 & relabund_sporobiota > 0, 'Only sporobiota', 
                                          ifelse(relabund_sporobiota == 0 & relabund_microbiota > 0, 'Only microbiota', 
                                                 ifelse(relabund_microbiota > 0 & relabund_sporobiota > 0 & relabund_microbiota > relabund_sporobiota, 'More in microbiota', 
                                                        ifelse(relabund_microbiota > 0 & relabund_sporobiota > 0 & relabund_microbiota < relabund_sporobiota, 'More in sporobiota', 
                                                               ifelse(relabund_microbiota > 0 & relabund_sporobiota > 0 & relabund_microbiota == relabund_sporobiota, 'Same amount', 'NA'))))))) %>%
  #filter(where != 'Not present in this sample') %>%
  group_by(otu) %>%
  mutate(R = cor.test(relabund_microbiota, relabund_sporobiota, method = 'pearson')$estimate, 
         pvalue = cor.test(relabund_microbiota, relabund_sporobiota, method = 'pearson')$p.value) 
  ggplot(aes(x=where)) +
  geom_histogram(stat = 'count')

