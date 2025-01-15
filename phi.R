library(tidyverse)

phi <- read_csv("phi.csv")  # phi(y,b) calls for x = 10^12

# b = C = 8 has many calls, decreasing until b = 431, 
# then b = a = 10453 has many calls
table(phi$b)


acbrt <- 11 * 1e4
bsize <- 2^20
blocks <- seq(1, 1e7, bsize)

# plot calls by block for selected b values
phi %>% filter(b %in% c(8, 100, 10453)) %>%
  ggplot(aes(y, fill = factor(b))) +
  geom_histogram(closed="right", breaks = blocks, position="dodge") + 
  geom_vline(xintercept = acbrt) + 
  scale_x_continuous(breaks = blocks, labels = \(x) as.hexmode(x)) +
  theme(axis.text.x=element_text(angle=45, hjust=1))
