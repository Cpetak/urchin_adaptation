# R script to visualise decay of LD over distances

# system("plink --bfile afterQC --cow --r2 --ld-window 1000 --ld-window-kb 1000  --ld-window-r2 0 --out LdExample")

library(tidyverse)

# read in LD results file
LdValues <- read_table("resultLDall2.ld")

# calculate LD in 20 kb bins to display the trendline
averageLD <- LdValues %>%
  mutate(markerDistance = abs(BP_B - BP_A)/1000) %>%
  dplyr::filter(markerDistance < 5000) %>%
  mutate(intervals = cut_width(markerDistance, 20, boundary = 0)) %>%
  group_by(intervals) %>%
  summarise_at(vars(R2),funs(mean(., na.rm=TRUE))) %>%
  rename(averageR2=R2)

# calculate inter marker distances
fullLD <- LdValues %>%
  mutate(markerDistance = abs(BP_B - BP_A)/1000) %>%
  dplyr::filter(markerDistance < 5000) %>%
  mutate(intervals = cut_width(markerDistance, 20, boundary = 0))

#merge the two data sets (full LD info and average per bin)
mergedLD <- full_join(fullLD,averageLD, by = "intervals")

# visualize LD decay
  ggplot(mergedLD) +
  geom_point(aes(x=markerDistance, y=R2)) +
    geom_line(aes(x=markerDistance, y=averageR2), color="red", size=2)
