library(dplyr)
library(readr)
library(ggplot2)
library(tidyr)

Human_Inserts <- read.csv('human_insert.csv')
Chimp_Inserts <- read.csv('chimp_insert.csv')
Bonobo_Inserts <- read.csv('bonobo_insert.csv')

Human_Inserts %>% 
  mutate(Species='Human') %>% 
  separate(col=insertion_name, into=c('repName','specifics'), sep='_range=') -> Human_Inserts

Chimp_Inserts %>% 
  mutate(Species='Chimpanzee') %>% 
  separate(col=insertion_name, into=c('repName','specifics'), sep='_range=') -> Chimp_Inserts

Bonobo_Inserts %>% 
  mutate(Species='Bonobo') %>% 
  separate(col=insertion_name, into=c('repName','specifics'), sep='_range=') -> Bonobo_Inserts


anti_join(Human_Inserts,Chimp_Inserts, by=c('repName','gene_name','chrom','instrand','genstrand','classification')) %>% 
  anti_join(Bonobo_Inserts, by=c('repName','gene_name','chrom','instrand','genstrand','classification')) %>% 
  group_by(gene_name,classification) %>% 
  summarize(class_count = n()) -> Homo2Chim2Bono

anti_join(Human_Inserts,Bonobo_Inserts, by=c('repName','gene_name','chrom','instrand','genstrand','classification')) %>% 
  anti_join(Chimp_Inserts, by=c('repName','gene_name','chrom','instrand','genstrand','classification')) %>% 
  group_by(gene_name,classification) %>% 
  summarize(class_count = n()) -> Homo2Bono2Chim

anti_join(Chimp_Inserts,Bonobo_Inserts, by=c('repName','gene_name','chrom','instrand','genstrand','classification')) %>% 
  anti_join(Human_Inserts, by=c('repName','gene_name','chrom','instrand','genstrand','classification')) %>% 
  group_by(gene_name,classification) %>% 
  summarize(class_count = n()) -> Chim2Bono2Homo
  