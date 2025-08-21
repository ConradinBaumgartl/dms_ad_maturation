source("00_init.r")

GC_curves <- viruses_fasta %>%
  map_df(get_GC_sliding, .id="species") %>%
  mutate(species = factor(species, levels = c("HAdC5", "HAdC2", "HAdB7", "HAdD26", "HAdE4", "HAdF41", "Bat Ad2", "Avian Ad celo")))

ggplot(data=GC_curves, aes(x=abs_loc, y=GC, color=species)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("absolute position") + 
  ylab("GC") +
  coord_cartesian(xlim=c(0,1), ylim=c(0.3, 0.8), expand=F) +
  scale_color_manual(values=ad_colors)
ggsave("plots/GC_curves_absolute.png", dpi=300, width=8, height = 3)
ggsave("plots/GC_curves_absolute.svg", dpi=300, width=8, height = 3)


ggplot(data=GC_curves, aes(x=loc, y=GC, color=species)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0.3, 0.8), expand=c(0,0)) +
  scale_color_manual(values=ad_colors) +
  facet_wrap(~species)
ggsave("plots/GC_curves_seperate.png", dpi=300, width=8, height = 3)
ggsave("plots/GC_curves_seperate.svg", dpi=300, width=8, height = 3)


