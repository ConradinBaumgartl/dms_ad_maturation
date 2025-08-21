source("01_GC_curves.r")
source("02_hexon_fiber_zoom.r")


# read mutation frequency from python script
hexon_mut_freq <- read_delim("../protein_identity_GC/mut_freq_.NEW.hexons.txt", delim = "\n", col_names = c("mutation_frequency"))
hexon_mut_freq$abs_loc <- seq(0, 1, length.out=nrow(hexon_mut_freq))
# smooth the mutation frequency with a loess and check the squigglyness
hexon_mut_freq$loess <- predict(loess(hexon_mut_freq$mutation_frequency ~ hexon_mut_freq$abs_loc, span = 0.1))
ggplot(data=hexon_mut_freq, aes(x=abs_loc, y=loess)) +
  geom_line()

# plot the hexons
hexon_plot <- L3_L5_entireSeq_GC %>% filter(product=="hexon")
ggplot() +
  geom_line(data=hexon_plot, aes(x=abs_loc, y=GC, color=species), size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0, 1), expand=F) +
  scale_color_manual(values=ad_colors) +
  geom_line(data=hexon_mut_freq, aes(x=abs_loc, y=loess), size=1) + 
  ggtitle("hexon - start 0 end 1") +
  scale_y_continuous(sec.axis = sec_axis(~.*100, name = "mutation frequency [%]"))
ggsave("plots/hexon_together_centered_0start1end_proteinLVLmutfreq.png", dpi=300, height = 5, width = 8)
ggsave("plots/hexon_together_centered_0start1end_proteinLVLmutfreq.svg", dpi=300, height = 5, width = 8)
ggplot() +
  geom_line(data=hexon_plot, aes(x=abs_loc, y=GC, color=species), size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0, 1), expand=F) +
  scale_color_manual(values=ad_colors) +
  geom_line(data=hexon_mut_freq, aes(x=abs_loc, y=loess), size=1.5) + 
  ggtitle("hexon - start 0 end 1") +
  scale_y_continuous(sec.axis = sec_axis(~.*100, name = "mutation frequency [%]")) +
  facet_wrap(~species)
ggsave("plots/hexon_seperate_centered_0start1end_proteinLVLmutfreq.png", dpi=300, height = 8, width = 8)
ggsave("plots/hexon_seperate_centered_0start1end_proteinLVLmutfreq.svg", dpi=300, height = 8, width = 8)

#### Fibers ####
# read mutation frequency from python script
fiber_mut_freq <- read_delim("../protein_identity_GC/mut_freq_.NEW.fibers.txt", delim = "\n", col_names = c("mutation_frequency"))
fiber_mut_freq$abs_loc <- seq(0, 1, length.out=nrow(fiber_mut_freq))
# smooth the mutation frequency with a loess and check the squigglyness
fiber_mut_freq$loess <- predict(loess(fiber_mut_freq$mutation_frequency ~ fiber_mut_freq$abs_loc, span = 0.1))
ggplot(data=fiber_mut_freq, aes(x=abs_loc, y=loess)) +
  geom_line()

fiber_plot <- L3_L5_entireSeq_GC%>%filter(product!="hexon")%>%filter(product!="short fiber")
ggplot() +
  geom_line(data=fiber_plot, aes(x=abs_loc, y=GC, color=species), size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0, 1), expand=F) +
  scale_color_manual(values=ad_colors) +
  geom_line(data=fiber_mut_freq, aes(x=abs_loc, y=loess), size=1) + 
  ggtitle("fiber - start 0 end 1") +
  scale_y_continuous(sec.axis = sec_axis(~.*100, name = "mutation frequency [%]"))
ggsave("plots/fiber_together_centered_0start1end_proteinLVLmutfreq.png", dpi=300, height = 5, width = 8)
ggsave("plots/fiber_together_centered_0start1end_proteinLVLmutfreq.svg", dpi=300, height = 5, width = 8)
ggplot() +
  geom_line(data=fiber_plot, aes(x=abs_loc, y=GC, color=species), size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0, 1), expand=F) +
  scale_color_manual(values=ad_colors) +
  geom_line(data=fiber_mut_freq, aes(x=abs_loc, y=loess), size=1.5) + 
  ggtitle("fiber - start 0 end 1") +
  scale_y_continuous(sec.axis = sec_axis(~.*100, name = "mutation frequency [%]")) +
  facet_wrap(~species)
ggsave("plots/fiber_seperate_centered_0start1end_proteinLVLmutfreq.png", dpi=300, height = 8, width = 8)
ggsave("plots/fiber_seperate_centered_0start1end_proteinLVLmutfreq.svg", dpi=300, height = 8, width = 8)

