# Graph x=Sample Size y= number of OTUs, but for all microbiota and sporobiota samples togehter with SD
# Dont know how this makes sence from rarecurve 
otutab= read.delim('data/mothur/final.opti_mcc.0.03.pick.shared') %>%
  select(-label, -numOtus)

otuM = otutab %>%
  filter(grepl('^M', Group)) %>%
  as_tibble() %>%
  pivot_longer(names_to = 'otu', values_to = 'value', cols=2:764) %>%
  group_by(otu) %>%
  summarise(otu_value=sum(value)) %>%
  mutate(biota= 'Microbiota')

otuS = otutab %>%
  filter(grepl('^S', Group)) %>%
  as_tibble() %>%
  pivot_longer(names_to = 'otu', values_to = 'value', cols=2:764) %>%
  group_by(otu) %>%
  summarise(otu_value=sum(value)) %>%
  mutate(biota= 'Sporobiota')

library(reshape2)
otus =rbind(otuM, otuS) %>%
  dcast(biota~otu, value.var = 'otu_value') %>%
  column_to_rownames('biota')

rarecurve(otus, step=100)

# Graph x= Otus, y=mean relabund from most to least abundand
otutab= read.delim('data/mothur/final.opti_mcc.0.03.pick.shared') %>%
  select(-label, -numOtus) %>%
  column_to_rownames('Group')
oturel = decostand(otutab, MARGIN=1, method='total')
rowSums(oturel)

otu_log10 = oturel %>%
  pivot_longer(names_to = 'otu', values_to = 'value', cols=1:763) %>%
  group_by(otu) %>%
  summarize(value=mean(value)) %>%
  ungroup() 

sum(otu_log10$value)

ggplot(otu_log10, aes(x=otu, y=log10(value))) +
  geom_point()
ggsave('plots/mothur/abundanceOTUs.png', dpi=600)
