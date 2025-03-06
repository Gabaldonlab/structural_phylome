library(tidyverse)
source("workflow/scripts/functions.R")
theme_set(theme_classic())

df_ml <- read_delim("results/test/trees/UP000005640_mltrees.txt") %>% 
  mutate(model=factor(gsub("3DI", "3Di", gsub("GTR20", "GTR", gsub("\\+.*", "", Model))), levels=models), 
         rate=sub("\\+", "", gsub("GTR20|3DI|LG|AF|LLM", "", Model))) %>%
  group_by(gene, targets) %>% 
  filter(model!="LG")


df_ml %>% 
  mutate(rank_LL=dense_rank(desc(LogL)),
         rank_BIC=dense_rank(BIC),
         rank_AIC=dense_rank(AIC),
         rank_AICc=dense_rank(AICc)) %>% 
  select(targets, model, contains("rank")) %>% 
  pivot_longer(contains("rank")) %>% 
  ggplot(aes(value, fill=model)) + 
  geom_bar(position = "dodge") + 
  scale_fill_manual(values = palettes_model) + 
  scale_y_continuous(expand = expansion(0,0)) + 
  facet_grid(name~.) + 
  labs(x="Rank", y="N")
