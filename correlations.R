library(tidyverse)
library(vegan)
library(ggrepel)
library(reshape2)
library(ggpubr)
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
ggsave('plots/mothur/comparisons/correlations_ggpairs.png', dpi=600)

# Comparison between xy variable and Chao 
comparisons1 = list(c('2', '3'), 
                   c('2', '4'), 
                   c('2', '5'), 
                   c('2', '6'), 
                   c('3', '4'), 
                   c('3', '5'), 
                   c('3', '6'),
                   c('4', '5'), 
                   c('4', '6'), 
                   c('5', '6'))

# Level of activity
alpha_meta %>% mutate(activity = active14+moderate14) %>%
  ggplot(aes(x=activity, y=chao, group = activity, color=person)) +
  geom_boxplot() +
  geom_point() +
  facet_grid(~biota) +
  labs(x='Chao1', y='Level of activity', color= 'Individual')
ggsave('plots/mothur/comparisons/chao_activity.png', dpi=600)

# Stress
comparisons2 = list(c(1,2), 
                    c(1,3), 
                    c(1,4), 
                    c(1,5), 
                    c(2,3), 
                    c(2,4), 
                    c(2,5), 
                    c(3,4), 
                    c(3,5), 
                    c(4,5))
alpha_meta %>%
  filter(stress != '0' & stress != '8') %>%
  ggplot(aes(x=stress, y=chao, group = stress, color=person)) +
  geom_boxplot() + 
  geom_point() +
  facet_wrap(~biota) +
  stat_compare_means(method='t.test', p.adjust.methods = 'BH', comparisons = comparisons2) +
  labs(x='Chao1', y='Level of stress', color= 'Individual')
ggsave('plots/mothur/comparisons/chao_stress.png', dpi=600)

# Bristol stool scale 
comparisons3 = list(c(2,3), 
                    c(2,4), 
                    c(2,5), 
                    c(3,4), 
                    c(3,5), 
                    c(4,5))
alpha_meta %>%
  ggplot(aes(x=bristol, y=chao, group = bristol, color=person)) +
  geom_boxplot() +
  geom_point() +
  facet_wrap(~biota) +
  #stat_compare_means(comparisons = comparisons3, method = 't.test', p.adjust.methods = 'BH') +
  labs(x='Chao1', y='Bristol stool scale', color= 'Individual')
ggsave('plots/mothur/comparisons/chao_bristol.png', dpi=600)

# 
alpha_meta %>%
  ggplot(aes(x=day, y=chao, color=person)) +
  geom_point() +
  geom_smooth(se=FALSE) +
  facet_wrap(~biota, ncol =1)

alpha_meta %>%
  ggplot(aes(x=diet14, y=chao)) +
  geom_point() +
  facet_wrap(~biota, ncol =1) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

alpha_meta %>%
  filter(antibiotics14 != '') %>%
  ggplot(aes(x=antibiotics14, y=chao)) +
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(~biota, ncol =1)

alpha_meta %>%
  filter(prebiotics14 != '') %>%
  ggplot(aes(x=prebiotics14, y=chao)) +
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(~biota, ncol =1) +
  theme(axis.text.x = element_text(angle=45, hjust=1))

alpha_meta %>%
  filter(probiotics14 != '') %>%
  ggplot(aes(x=probiotics14, y=chao)) +
  geom_boxplot() +
  stat_compare_means() +
  facet_wrap(~biota, ncol =1) +
  theme(axis.text.x = element_text(angle=45, hjust=1))


