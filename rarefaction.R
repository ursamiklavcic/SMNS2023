# rarefaction/normalization analysis of longitudinal data 
library(tidyverse)
library(vegan)
library(reshape2)

set.seed(1996)

shared = read_tsv('data/mothur/final-48.shared') %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  # Calculate the total number of seqs in each sample and remove samples with less than 150 000
  mutate(total_sample=sum(value)) %>%
  filter(total_sample > 150000) %>%
  group_by(name) %>%
  # Calculate the total number of seqs for each otu and remove otus with 0 reads
  mutate(total_otu=sum(value)) %>%
  filter(total_otu != 0) %>%
  ungroup() %>%
  select(-total_sample, -total_otu)

# 1. Comparison of exact vs empirical approach 
# If we choose a sampling depth how many OTUs will be observed at that point?
# A function to subsample data once 
subsample = function(data, sample_size){
  data %>%
    group_by(Group) %>%
    uncount(value) %>%
    sample_n(size=sample_size) %>%
    summarise(noOTU = n_distinct(name))
}

# Do multiple subsamplings
subsamplings = map_dfr(1:100, ~subsample(shared, 150000), .id = 'iters')

# Calculate from our function the empirical number of OTUs we would get at a choosen sequencing depth
empirical = subsamplings %>% 
  group_by(Group) %>%
  summarise(noOTUs =mean(noOTU))

# Use VEGANS exact number of rarefaction! 
exact = shared %>%
  group_by(Group) %>%
  summarise(noOTUs =rarefy(value, 150000))

bind_rows(empirical=empirical, exact=exact, .id='approach')%>%
  ggplot(aes(x=approach, y=noOTUs, group=Group)) +
  geom_line() +
  theme_bw()

ggsave('plots/mothur/empirical_exact_line.png', dpi=600)

inner_join(empirical, exact, by='Group') %>%
  ggplot(aes(x=noOTUs.x, y=noOTUs.y)) +
  geom_point() +
  geom_abline(intercept = 0, slope = 1, color='grey') +
  theme_bw()

ggsave('plots/mothur/empirical_exact_slope.png', dpi=600)

# What is the procent error between empirical and exact
inner_join(empirical, exact, by='Group') %>%
  mutate(error=100*(noOTUs.y-noOTUs.x)/noOTUs.y) %>%
  ggplot(aes(x=error)) +
  geom_density()+
  theme_bw()
ggsave('plots/mothur/exact_empirical_density.png', dpi=600)

# calculcate the exact mean and sd of error 
inner_join(empirical, exact, by='Group') %>%
  mutate(error=100*(noOTUs.y-noOTUs.x)/noOTUs.y) %>%
  summarise(mean= mean(error), sd=sd(error))

# 2. Comparison of rarefaction vs normalization vs relative abundance vs no rarefaction EFFECT ON DISTANCES (Bray-Curtis)
# This matrix has samples that are not statisticly different from eachother
rand = shared %>%
  uncount(value) %>%
  mutate(rand_name = sample(name)) %>%
  select(-name) %>%
  count(Group, rand_name)
# Turn ran into matrix so it can be used in distance calculations
rand_df = rand %>%
  pivot_wider(names_from="rand_name", values_from="n", values_fill = 0) %>%
  column_to_rownames('Group') %>%
  as.data.frame()

rand_matrix <- rand_df %>%
  as.matrix()

# From rand matrix! 
# Calculate one itteration of distance matrix and avgdist - default is 100 
norare_dist_matrix <- vegdist(rand_matrix, method="bray")
rare_dist_matrix <- avgdist(rand_matrix, dmethod="bray", sample=150000, iterations=9)

norare_dist_tibble <- norare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

rare_dist_tibble <- rare_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

# Calculate avgdist from relative abundance data
relabund_matrix = rand %>%
  group_by(Group) %>%
  mutate(rel_abund = n/sum(n)) %>%
  ungroup() %>%
  select(-n) %>%
  pivot_wider(names_from = 'rand_name', values_from = 'rel_abund', values_fill = 0) %>%
  as.data.frame() %>%
  column_to_rownames('Group') %>%
  as.matrix()

relabund_dist_matrix = vegdist(relabund_matrix, method='bray')

relabund_dist_tibble <- relabund_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample)

# Normalization 
# What is the smallest sample (the minimum number of sequences in a sample?)
rand_group_count = rand %>% 
  group_by(Group) %>%
  summarise(n=sum(n))
min_group = min(rand_group_count$n)

# Make a normalizied matrix
# Do normalization by hand

# rand %>%
#   group_by(Group) %>%
#   mutate(rel_abund = n/sum(n)) %>%
#   ungroup() %>%
#   select(-n) %>%
#   mutate(scaled=round(rel_abund*min_group, 0)) %>%
#   group_by(Group) %>%
#   summarise(n_scaled=sum(scaled))

# I see that not all samples were normalized to the exact same number! 
# package SRS has built in functions for normalization of ecology data: 
library(SRS)  # https://pubmed.ncbi.nlm.nih.gov/32832266/

# SRS needs a otutab that has OTUs as rows and samples by columns
# Normalized to the min_group number! 
normalized = rand_df %>%
  t() %>%
  as.data.frame() %>%
  SRS(Cmin=min_group)
# Check if all samples have min_group! 
normalized %>% as_tibble(rownames = 'otu') %>%
  pivot_longer(-otu) %>%
  group_by(name) %>%
  summarise(n=sum(value))

normalized_dist_matrix = vegdist(t(normalized), method='bray') 

normalized_dist_tibble <- normalized_dist_matrix %>%
  as.matrix() %>%
  as_tibble(rownames="sample") %>%
  pivot_longer(-sample) %>%
  filter(name < sample) 

comparison <- inner_join(norare_dist_tibble, rare_dist_tibble, by=c("sample", "name")) %>%
  inner_join(., relabund_dist_tibble, by=c('sample', 'name')) %>%
  inner_join(., normalized_dist_tibble, by=c('sample', 'name')) %>%
  select(sample, name, norarefied=value.x, rarefied=value.y, relabund=value.x.x, normalized=value.y.y) %>%
  inner_join(., rand_group_count, by=c("sample" = "Group")) %>%
  inner_join(., rand_group_count, by=c("name" = "Group")) %>%
  mutate(n_diff = abs(n.x-n.y)) %>%
  select(-n.x, -n.y)

comparison %>%
  pivot_longer(cols=c("norarefied", "rarefied", 'relabund', 'normalized'), names_to="type", values_to="dist") %>%
  ggplot(aes(x=n_diff,  y=dist)) +
  geom_point(size=0.25, alpha=0.25) +
  facet_wrap(~type, nrow=4, scales = 'free_y')

ggsave('plots/mothur/norare_rare-relabund-normal.png', dpi=600)

# 3. How rarefaction/normalization/no rarefaction / relative abundance EFFECTS ALPHA DIVERSITY MATRICES (Observed/Shannon/Chao)
shared %>%
  # Make a randomized otutab, without subsampling !
  uncount(value) %>%                                      
  mutate(name = sample(name)) %>%
  count(Group, name, name='value') %>%
  group_by(Group) %>%
  # Calculate alpha diversity metrics ! 
  summarise(observed = specnumber(value),               
            shannon = diversity(value, index = 'shannon'), 
            simpson = diversity(value, index='simpson'), 
            invsimpson = diversity(value, index = 'invsimpson'), 
            n=sum(value)) %>%
  pivot_longer(cols = c(observed, shannon, simpson, invsimpson), 
               names_to = 'metric') %>%
  ggplot(aes(x=n, y=value)) +
  geom_point() +
  geom_smooth()+
  facet_wrap(~metric, nrow=4, scales='free_y')

ggsave('plots/mothur/alpha_samplingEffect.png', dpi = 600)

# What is the number to which I should rarefy my data? 
shared1 = read_tsv('data/mothur/final-48.shared') %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group) %>%
  # Calculate the total number of seqs for each otu and remove otus with 0 reads
  group_by(name) %>%
  mutate(total_otu=sum(value)) %>%
  filter(total_otu != 0) %>%
  ungroup() %>%
  select(-total_otu)

# rand = shared %>%
#   uncount(value) %>%
#   mutate(name=sample(name)) %>%
#   count(Group, name, name = 'value')

# Calculate the sampling coverage for each sample 
sampling_covergae = shared1 %>%
  group_by(Group) %>%
  summarise(n_seqs=sum(value))
# Plot sampling covergae with a histogram plot 
sampling_covergae %>% ggplot(aes(x=n_seqs)) +
  geom_histogram(binwidth = 10000) +
  coord_cartesian(xlim=c(0,200000))

ggsave('plots/mothur/sampling_coverage_histo_partial.png', dpi=600)

sampling_covergae %>% ggplot(aes(x=1, y=n_seqs)) +
  geom_jitter() +
  scale_y_log10()
ggsave('plots/mothur/sampling_coverage_jitter.png', dpi=600)

sampling_covergae %>%
  arrange(n_seqs) %>%
  ggplot(aes(x=1:nrow(.), y=n_seqs))+
  geom_line() + 
  coord_cartesian(xlim=c(0,50), ylim=c(0, 200000))
ggsave('plots/mothur/sampling_coverage_line_partial.png', dpi=600)

# make a list of the samples that I will lose: 
sampling_covergae %>%
  arrange(n_seqs)  %>%
  print(n=20)

# What can I do with that sampling depth that I have left?
# Good's covergae = fraction of sequences that appear in an OTU that has been seen more than 1 
coverage_stats = shared %>%
  group_by(Group) %>%
  summarise(n_seqs=sum(value), 
            n_singletons = sum(value==1), 
            goods= 100*(1- n_singletons/n_seqs)) 

coverage_stats%>%
  ggplot(aes(x=n_seqs, y=goods)) +
  geom_point()
ggsave('plots/mothur/goods_coverage_150000.png', dpi=600)

coverage_stats%>%
  arrange(goods)

#######################
# 4. How to rarefy in R 
# To do rarefaction in R, you have to do multiple rrarefy and calculate the metric you wish and than continue
# That's why I'm doing rarefaction in mothur ! 
shared = read_tsv('data/mothur/final-48.shared') %>%
  select(Group, starts_with('Otu')) %>%
  pivot_longer(-Group) %>%
  group_by(Group) %>%
  # Calculate the total number of seqs in each sample and remove samples with less than 150 000
  mutate(total_sample=sum(value)) %>%
  filter(total_sample > 100000) %>%
  group_by(name) %>%
  # Calculate the total number of seqs for each otu and remove otus with 0 reads
  mutate(total_otu=sum(value)) %>%
  filter(total_otu != 0) %>%
  ungroup() %>%
  select(-total_sample, -total_otu)



# Rarefaction 
# Transform otutab into shape for vegan::rarecurve 
oturare = shared %>% dcast(Group ~name , value.var = 'value') %>%
  column_to_rownames('Group')
# Calculate rarefaction curve 
rarecurve(oturare, step = 100, xlab= 'Sample Size', ylab='OTUs')
# Rarefaction to 150 000 reads per sample with 999 iterattions. Which samples have less than 150000 redas?
insuff = rownames(oturare)[rowSums(oturare) < 100000]
# Exclude this samples from analysis. 
oturar = oturare %>%
  as.data.frame() %>%
  rownames_to_column('samples') %>%
  filter(., !(samples %in% insuff)) %>%
  column_to_rownames('samples') %>%
  as.data.frame()
# Rrarefy otutable, to get a multiple subsampled table
rrarefy(oturar, sample= 100000)

rarefactions = map_dfr(1:3, ~rrarefy(oturar, sample = 100000), .id='iters')

subsamplings = map_dfr(1:100, ~subsample(shared, 150000), .id = 'iters')

# Calculate from our function the empirical number of OTUs we would get at a choosen sequencing depth
empirical = subsamplings %>% 
  group_by(Group) %>%
  summarise(noOTUs =mean(noOTU))