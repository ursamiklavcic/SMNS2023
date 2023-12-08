# Phylogenetic analysis of the data 
# Library 
require(tidyverse)
library(vegan)
library(ape)
require(phyloseq)
require(ggrepel)
require(ggpubr)
library(glue)
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

# Load tree file from mafft 
tree = read.tree('data/mothur/fasttree_mafft_align2.tree')

# Final count table (before making OTUs)
count_table = read.table('data/mothur/final.full.count_table', header = TRUE, sep = "\t") %>%
  filter(total > 1)
# Table of only sequences and read counts
PAtable = select(count_table, seq = Representative_Sequence, 3:236) %>% column_to_rownames('seq')
# Get the name of samples that have less than 100 000 reads per sample!
reads= 100000
remove_samples = as.data.frame(t(PAtable)) %>% mutate(sum=rowSums(select(., starts_with('V')))) %>% filter(sum > reads) %>% rownames_to_column('name') %>% pull(name)

# Remove samples that have less than 100 000 reads per sample! 
PAfilter = as.data.frame(t(PAtable)) %>% rownames_to_column('samples') %>%filter(samples %in% remove_samples) %>% column_to_rownames('samples') %>% as.data.frame()
# rarefy the table to 100 000 reads per sample
rare_table =  rrarefy(PAfilter, sample=reads) %>% as.data.frame()

# Taxonomy 
taxonomy = read.csv('data/mothur/final.taxonomy', header=FALSE, sep='\t') %>%
  rename(seq = V1, taxa = V2) %>%
   mutate(taxa = str_replace_all(taxa, "\\\\|\\\"|\\(\\d+\\)", ""),
          taxa = str_replace(taxa, ";$", "")) %>%
   separate(taxa, into=c("Domain", "Phylum", "Class", "Order", "Family", "Genus"),
            sep=";") %>%
  filter(seq %in% colnames(rare_table)) %>%
  column_to_rownames('seq')
  

# Import metadata 
metadata = readRDS('data/metadata.RDS') %>% 
  filter(samples %in% remove_samples) %>%
  column_to_rownames('samples')

# Make sure that sequences that we filtered (only have 1 read or less) are removed from the tree as well 
clean_tree = drop.tip(phy=tree, 
                      tip=setdiff(tree$tip.label, rownames(t(rare_table))))

# Calculate Faith's index from picante
PD = pd(rare_table, clean_tree, include.root = FALSE)
PD_meta = PD %>% rownames_to_column(('samples')) %>% left_join(metadata %>% rownames_to_column('samples'), by='samples')

ggplot(PD_meta, aes(x=person, y=PD, color=person)) +
  geom_boxplot() +
  geom_point() +
  labs(x='Person',y='Phylogenetic diversity') +
  facet_wrap(~biota, nrow=2) +
  theme_bw()
ggsave('plots/PD_individual.png', dpi=600)
# PD through time 
PD_meta_tmp = mutate(PD_meta, person2=person)

ggplot(PD_meta_tmp, aes(x=day, y=PD)) +
  geom_line(data=PD_meta_tmp %>% 
              dplyr::select(-person) %>% 
              filter(biota == 'Microbiota'), 
            aes(group=person2),  color= colm, linewidth=0.5, alpha=0.5) +
  geom_line(data=PD_meta_tmp %>% 
              dplyr::select(-person) %>% 
              filter(biota == 'Sporobiota'), 
            aes(group=person2),  color= cols, linewidth=0.5, alpha=0.5) + 
  geom_line(data=PD_meta_tmp %>% 
              filter(biota == 'Microbiota'),
            aes(color=person), color= colm, linewidth=1.2 )+
  geom_line(data=PD_meta_tmp %>% 
              filter(biota == 'Sporobiota'),
            aes(color=person), color=cols, linewidth=1.2 )+
  geom_text_repel(aes(label= samples), 
                  size=3, max.overlaps = 20) +
  geom_point(aes(color=event14_type), size=2) +
  scale_color_manual(values= col_extreme, na.value='#B8B9B6') +
  facet_wrap(~person, scales = 'free') +
  labs(x='Day', y= 'Phylogenetic diversity', color= 'Type of extreme event') +
  theme_bw()
ggsave('plots/PD_through_time.png', dpi=600)

# How does PD correlate with SR (species richness)?
ggplot(PD_meta, aes(x=PD, y=SR, color=biota)) +
  geom_point()
cor.test(PD_meta$PD, PD_meta$SR, method = 'pearson') # Strongly positive (0,96) significant (p-value=2.2e-16) correlation between phylogenetic diversity and Species richness

# BETA DIVERSITY # 
# Calculate UniFrac (unweight/weight)
# Construct phyloseq object for easier manipulation 
ps = phyloseq(otu_table(as.matrix(rare_table), taxa_are_rows = FALSE), sample_data(metadata), tax_table(as.matrix(taxonomy)), phy_tree(clean_tree))

unifrac_w = UniFrac(ps, weighted = TRUE, normalized = FALSE, parallel = TRUE)
unifrac_u = UniFrac(ps, weighted = FALSE, parallel = TRUE)

# PCoA  weighted 
pcoa = cmdscale(unifrac_w, k=2, eig=TRUE, add=TRUE)
positions = pcoa$points
colnames(positions) <- c('pcoa1', 'pcoa2')

percent_explained = 100 * pcoa$eig/ sum(pcoa$eig)
pretty_pe = format(round(percent_explained[1:2], digits = 1), nsmall = 1, trim=TRUE)
labs = c(glue('PCo 1 ({pretty_pe[1]}%)'), 
         glue('PCo 1 ({pretty_pe[2]}%)'))

positions %>%
  as_tibble(rownames = 'samples') %>%
  left_join(rownames_to_column(metadata, 'samples'), by='samples') %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color=person)) +
  geom_point(size=2) +
  #scale_color_manual(values = colsm) + 
  labs(x=labs[1], y=labs[2])

tibble(pe=cumsum(percent_explained), 
       axis=1:length(percent_explained)) %>%
  ggplot(aes(x=axis, y=pe)) +
  geom_line() +
  coord_cartesian(xlim=c(1,10))
# PCoA only explains in 2 axis 20% of variation of my data 

# PCoA unweighted 
pcoa = cmdscale(unifrac_u, k=2, eig=TRUE, add=TRUE)
positions = pcoa$points
colnames(positions) <- c('pcoa1', 'pcoa2')

percent_explained = 100 * pcoa$eig/ sum(pcoa$eig)
pretty_pe = format(round(percent_explained[1:2], digits = 1), nsmall = 1, trim=TRUE)
labs = c(glue('PCo 1 ({pretty_pe[1]}%)'), 
         glue('PCo 1 ({pretty_pe[2]}%)'))

positions %>%
  as_tibble(rownames = 'samples') %>%
  left_join(rownames_to_column(metadata, 'samples'), by='samples') %>%
  ggplot(aes(x=pcoa1, y=pcoa2, color=person)) +
  geom_point(size=2) +
  #scale_color_manual(values = colsm) + 
  labs(x=labs[1], y=labs[2])

tibble(pe=cumsum(percent_explained), 
       axis=1:length(percent_explained)) %>%
  ggplot(aes(x=axis, y=pe)) +
  geom_line() +
  coord_cartesian(xlim=c(1,10))

# Basic NMDS
nmds_wuf = metaMDS(unifrac_w)
nmds_positions= as.data.frame(scores(nmds_wuf, display='sites')) %>%
  rownames_to_column('samples')
# Join with metadata
distwuf_meta = nmds_positions %>% left_join(rownames_to_column(metadata, 'samples'), by='samples')
# Ordination plot wih person/biota 
distwuf_meta %>%
  ggplot(aes(x=NMDS1, y=NMDS2, color=person, shape=biota)) +
  geom_point(size=3) +
  #geom_text_repel(aes(label=sample), size= 4, colour='black', max.overlaps = 20) +
  scale_size_continuous(range = c(3,6)) +
  labs(x='', y='', shape='Type of biota', color='Individual') +
  theme_bw()

# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
unifracW_df = as.data.frame(as.matrix(unifrac_w)) %>%
  rownames_to_column('sample') %>%
  pivot_longer(-sample) %>%
  filter(sample != name) %>%
  left_join(rownames_to_column(metadata, 'samples') %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(rownames_to_column(metadata, 'samples')  %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Sporobiota' & biota.y == 'Sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracW_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  geom_boxplot() +
  #geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c(colm, cols)) +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="weighted UniFrac distance", x="", fill='Type of sample') +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))
ggsave('plots/UniFrac_W_boxplot.png', dpi=600)

# Unweighted UniFrac distances 
unifracU_df = as.data.frame(as.matrix(unifrac_u)) %>%
  rownames_to_column('sample') %>%
  pivot_longer(-sample) %>%
  filter(sample != name) %>%
  left_join(rownames_to_column(metadata, 'samples') %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(rownames_to_column(metadata, 'samples')  %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Sporobiota' & biota.y == 'Sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both')

unifracU_df %>% ggplot(aes(x=same_person, y=value, fill=which_biota)) +
  #geom_boxplot() +
  geom_violin(draw_quantiles = c(0.5)) +
  scale_fill_manual(values = c(colm, cols)) +
  #geom_point(size=0.5, alpha=0.3) +
  stat_compare_means(aes(group=paste0(same_person, which_biota))) +
  labs(y="unweighted UniFrac distance", x="", fill='Type of sample') +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))
ggsave('plots/UniFrac_U.png', dpi=600)

# Weighted 
# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
unifrac = as.matrix(unifrac_w) %>% 
  as_tibble(rownames= 'sample') %>%
  pivot_longer(-sample) %>%
  # Remove the distances of the same sample 
  filter(sample != name) %>%
  # Add metadata
  left_join(rownames_to_column(metadata, 'samples') %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(rownames_to_column(metadata, 'samples') %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Sporobiota' & biota.y == 'Sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both' & same_person != 'Different individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(which_biota, diff, person.x) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup()

ggplot(unifrac, aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  scale_color_manual(values=c(colm, cols)) +
  geom_smooth(aes(group=which_biota), method = 'lm') +
  labs(x='Days between sampling points', y='Median Weighted UniFrac distance', color='Type of sample') +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))
ggsave('plots/UniFrac_W_time.png', dpi=600)
# Pearsons correlation between median of distance between samples and time
unifracM = filter(unifrac, which_biota == 'Microbiota')
cor.test(as.numeric(unifracM$diff), unifracM$median, method='pearson') # Negative non significant correlation
unifracS = filter(unifrac, which_biota == 'Sporobiota')
cor.test(as.numeric(unifracS$diff), unifracS$median, method='pearson') # Positive non significant correlation

# Unweighted 
# Calculate if there is a difference in the distance between samples of individuals if they were sampled closer together and more appart between microbiota and sporobiota. 
unifrac = as.matrix(unifrac_u) %>% 
  as_tibble(rownames= 'sample') %>%
  pivot_longer(-sample) %>%
  # Remove the distances of the same sample 
  filter(sample != name) %>%
  # Add metadata
  left_join(rownames_to_column(metadata, 'samples') %>% select(samples, person, date, biota), by=join_by('sample' == 'samples')) %>%
  left_join(rownames_to_column(metadata, 'samples') %>% select(samples, person, date, biota), by=join_by('name' == 'samples')) %>%
  # Filter so that I have only inter-person comparisons!
  mutate(same_person= ifelse(person.x==person.y, 'Same individual', 'Different individual'), 
         which_biota= ifelse(biota.x == 'Microbiota' & biota.y == 'Microbiota', 'Microbiota',
                             ifelse(biota.x == 'Sporobiota' & biota.y == 'Sporobiota', 'Sporobiota', 'Both'))) %>%
  filter(which_biota != 'Both' & same_person != 'Different individual') %>%
  # Calculate the difference between sampling times
  mutate(diff=abs(date.x-date.y)) %>%
  # group by difference between days and person
  group_by(which_biota, diff, person.x) %>%
  summarise(median=median(value), sd= sd(value)) %>%
  ungroup()

ggplot(unifrac, aes(x=diff, y=median, color=which_biota)) +
  geom_point() +
  scale_color_manual(values=c(colm, cols)) +
  geom_smooth(aes(group=which_biota), method = 'lm') +
  labs(x='Days between sampling points', y='Median Unweighted UniFrac distance', color='Type of sample') +
  theme_bw(base_size = 18) +
  theme( text = element_text(family = "Calibri"))
ggsave('plots/UniFrac_U_time.png', dpi=600)
# Pearsons correlation between median of distance between samples and time
unifracM = filter(unifrac, which_biota == 'Microbiota')
cor.test(as.numeric(unifracM$diff), unifracM$median, method='pearson') # Positive (0,19) correlation that is very significant (0,00005039); 95 percent confidence interval: 0.1038342 0.2897334
unifracS = filter(unifrac, which_biota == 'Sporobiota')
cor.test(as.numeric(unifracS$diff), unifracS$median, method='pearson') # Positive (0,14) correlation that is just significant (0,0089); 95 percent confidence interval: 0.03583862 0.24524096

# 