ov_annotations <- read_csv("/home/dwarriaz/data/lab_work/OverwriTE_annotations.csv")

ov_annotations %>% 
  group_by(name2) %>% 
  select(-name,-classification) %>% 
  distinct() %>% 
  write_csv(.,file = '/home/dwarriaz/data/lab_work/data_for_queenie.csv')

batch14 <- read_csv("/home/dwarriaz/data/lab_work/output/Batch14_output.csv")

group_by('name2') %>% 
  select(-name,-classification) %>% 
  distinct() %>% 
  write_csv(.,file = '/home/dwarriaz/data/lab_work/updated_batch14_output.csv')


library(tidyverse)


filenames <- list.files("/home/dwarriaz/data/lab_work/1vsAll/", pattern  = "*.results",full.names = TRUE)

ldf <- lapply(filenames, function(x){
  denom <- str_split(x, pattern = '_vs')[[1]][1]
  denom <- str_split(denom, pattern = '//')[[1]][2]
  read_csv(x,skip = 7) %>% 
  mutate(tissue=denom) %>% 
  filter(name %in% ov_annotations$name2) %>% 
  filter(padj <= 0.01,abs(log2FoldChange)>3)
  }) %>% bind_rows()

ldf %>% 
  group_by(name) %>%
  summarize(count=n()) %>% 
  arrange(-count) 
  
ldf %>% 
  merge(ov_annotations,by.x = 'name',by.y = 'name2') -> vall_ov

vall_ov %>%
  group_by(name, name.y) %>% 
  summarise(count = n()) %>% 
  group_by(name) %>% 
  summarise(isocount = n()) -> isocount

vall_ov %>% 
  filter(classification == "Transcription Start Site") %>% 
  group_by(tissue, name, repFamily) %>% 
  summarise(count = n(),log2fc = mean(log2FoldChange)) %>% 
  merge(isocount, by = 'name') %>% 
  mutate(count=count/isocount) -> beertime

ggplot(beertime, aes(x=count,y=log2fc, color=repFamily)) + 
  geom_point() +
  facet_wrap(~tissue) + 
  geom_hline(yintercept = 0)
  

