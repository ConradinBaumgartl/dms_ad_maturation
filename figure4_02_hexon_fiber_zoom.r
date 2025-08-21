source("00_init.r")

viruses_CDS <- viruses_gff3 %>%
  map(filter, type=="CDS")


# Incredibly tedious way of combining all hexon and fiber genes into one dataframe
# but I did not find a better way, since the Hexon and fiber genes are all named slightly differently in the gff3 files
c5_hexon <- dplyr::filter(viruses_CDS$HAdC5, grepl("L3 pII", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdC5", product="hexon")
c5_fiber <- dplyr::filter(viruses_CDS$HAdC5, grepl("L5", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdC5", product="fiber")
c2_hexon <- dplyr::filter(viruses_CDS$HAdC2, grepl("hexon", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdC2", product="hexon")
c2_fiber <- dplyr::filter(viruses_CDS$HAdC2, grepl("fiber", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdC2", product="fiber")
b7_hexon <- dplyr::filter(viruses_CDS$HAdB7, grepl("hexon", product)) %>%
  dplyr::filter(grepl("L3", gene)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdB7", product="hexon")
b7_fiber <- dplyr::filter(viruses_CDS$HAdB7, grepl("fiber", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdB7", product="fiber")
d26_hexon <- dplyr::filter(viruses_CDS$HAdD26, grepl("hexon", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdD26", product="hexon")
d26_fiber <- dplyr::filter(viruses_CDS$HAdD26, grepl("fiber", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdD26", product="fiber")
e4_hexon <- dplyr::filter(viruses_CDS$HAdE4, grepl("hexon", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdE4", product="hexon")
e4_fiber <- dplyr::filter(viruses_CDS$HAdE4, grepl("fiber", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdE4", product="fiber")
f41_hexon <- dplyr::filter(viruses_CDS$HAdF41, grepl("hexon", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdF41", product="hexon")
f41_fiber <- dplyr::filter(viruses_CDS$HAdF41, grepl("fiber", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="HAdF41", product="fiber")
bat_hexon <- dplyr::filter(viruses_CDS$`Bat Ad2`, grepl("hexon", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="Bat Ad2", product="hexon")
bat_fiber <- dplyr::filter(viruses_CDS$`Bat Ad2`, grepl("fiber", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="Bat Ad2", product="fiber")
avian_hexon <- dplyr::filter(viruses_CDS$`Avian Ad celo`, grepl("hexon", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="Avian Ad celo", product="hexon")
avian_fiber <- dplyr::filter(viruses_CDS$`Avian Ad celo`, grepl("fiber", product)) %>%
  as_tibble() %>%
  select(c(seqnames, start, end, strand, product)) %>%
  mutate(species="Avian Ad celo", product="fiber")

L3_L5 <- bind_rows(
  c5_hexon,
  c5_fiber,
  c2_hexon,
  c2_fiber,
  b7_hexon,
  b7_fiber,
  d26_hexon,
  d26_fiber,
  e4_hexon,
  e4_fiber,
  f41_hexon,
  f41_fiber,
  bat_hexon,
  bat_fiber,
  avian_hexon,
  avian_fiber
  )

# get the lengths of every gene
L3_L5 %>% 
  filter(product == "hexon") %>%
  mutate(width = end-start) %>%
  select(c(product, species, width)) %>%
  write_delim("hexon_lengths.tab", delim="\t")
L3_L5 %>% 
  filter(product == "fiber") %>%
  mutate(width = end-start) %>%
  select(c(product, species, width)) %>%
  write_delim("fiber_lengths.tab", delim="\t")

# center every gene on the start point. 500 bp before and 2500 after
after_start <- 3000
before_start <- 1000
total_extension <- after_start+before_start
L3_L5_centered <- L3_L5 %>%
  mutate(end = start+after_start, start=start-before_start) %>%
  rowwise() %>%
  mutate(sequence = toString(viruses_fasta[[species]][[1]][start:end])) # get the sequence of every section

wind_s <- 250
step_s <- 10
number_of_windows <- (total_extension - wind_s) / step_s
test <- L3_L5_centered$sequence %>%
  map(get_GC_sliding, windowsize=wind_s, stepsize=step_s, nchar=T) %>%
  map_df(as_tibble) %>%
  mutate(species = c(rep("HAdC5", 2*number_of_windows),
                     rep("HAdC2", 2*number_of_windows),
                     rep("HAdB7", 2*number_of_windows),
                     rep("HAdD26", 2*number_of_windows),
                     rep("HAdE4", 2*number_of_windows),
                     rep("HAdF41", 3*number_of_windows),
                     rep("Bat Ad2", 2*number_of_windows),
                     rep("Avian Ad celo", 3*number_of_windows))) %>%
  mutate(product = c(rep(c(rep("hexon", number_of_windows),
                           rep("fiber", number_of_windows)), 5),
                     rep("hexon", number_of_windows),
                     rep("short fiber", number_of_windows),
                     rep("long fiber", number_of_windows),
                     rep("hexon", number_of_windows),
                     rep("fiber", number_of_windows),
                     rep("hexon", number_of_windows),
                     rep("long fiber", number_of_windows),
                     rep("short fiber", number_of_windows))) %>%
  mutate(species = factor(species, levels = c("HAdC5", "HAdC2", "HAdB7", "HAdD26", "HAdE4", "HAdF41", "Bat Ad2", "Avian Ad celo")))


ggplot(data=test%>%filter(product=="hexon"), aes(x=loc, y=GC, color=species)) +
  geom_line(size=1) +
  geom_vline(xintercept =before_start) + 
  ggtitle("hexon start codon indicated") +
  facet_wrap(~species) +
  theme_bw() +
  scale_color_manual(values=ad_colors) + 
  coord_cartesian(ylim=c(0.3, 0.8), expand=F)
ggsave("plots/hexon_centered.png", dpi=300)
ggsave("plots/hexon_centered.svg", dpi=300)


ggplot(data=test%>%filter(product!="hexon")%>%filter(product!="short fiber"), aes(x=loc, y=GC, color=species)) +
  geom_line(size=1) +
  geom_vline(xintercept =before_start) +
  ggtitle("fiber start codon indicated") +
  facet_wrap(~species) +
  theme_bw() +
  scale_color_manual(values=ad_colors)
ggsave("plots/fiber_centered.png", dpi=300)
ggsave("plots/fiber_centered.svg", dpi=300)

# add hexon fiber positions to the GC-curves
L3_L5_plot <- L3_L5 %>%
  select(product, species, start) %>%
  rename(start = "loc")

ggplot(data=GC_curves, aes(x=loc, y=GC, color=species)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0.3, 0.8), expand=c(0,0)) +
  scale_color_manual(values=ad_colors) +
  facet_wrap(~species) +
  geom_vline(data=L3_L5_plot%>%filter(product=="hexon"), aes(xintercept=loc)) +
  ggtitle("Hexon start position indicated")
ggsave("plots/GC_curves_seperate_hexon.png", dpi=300, width=8, height = 8)
ggsave("plots/GC_curves_seperate_hexon.svg", dpi=300, width=8, height = 8)

ggplot(data=GC_curves, aes(x=loc, y=GC, color=species)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0.3, 0.8), expand=c(0,0)) +
  scale_color_manual(values=ad_colors) +
  facet_wrap(~species) +
  geom_vline(data=L3_L5_plot%>%filter(product=="fiber"), aes(xintercept=loc)) +
  ggtitle("Fiber start position indicated")
ggsave("plots/GC_curves_seperate_fiber.png", dpi=300, width=8, height = 8)
ggsave("plots/GC_curves_seperate_fiber.svg", dpi=300, width=8, height = 8)

#### Do the GC-content over the entire, relative Gene sequences
TTS_ext <- 0
TSS_ext <- 0
L3_L5_entireSeq <- L3_L5 %>%
  mutate(end = end + TTS_ext, start = start - TSS_ext) %>%
  rowwise() %>%
  mutate(sequence = toString(viruses_fasta[[species]][[1]][start:end])) # get the sequence of every section
# calc GC
wind_s <- 100
step_s <- 10
L3_L5_entireSeq_GC <- L3_L5_entireSeq$sequence %>%
  map(get_GC_sliding, windowsize=wind_s, stepsize=step_s, nchar=T)

# get the correct number of windows for every gene
number_of_windows <- L3_L5_entireSeq_GC %>%
  map(nrow) %>%
  unlist()
# define vectors of species and product according to L3_L5
had_species_v <- c(
  "HAdC5","HAdC5", "HAdC2","HAdC2", "HAdB7","HAdB7", "HAdD26","HAdD26", "HAdE4","HAdE4", "HAdF41","HAdF41","HAdF41", "Bat Ad2","Bat Ad2", "Avian Ad celo", "Avian Ad celo", "Avian Ad celo"
)
product_v <- c(
  "hexon", "fiber","hexon", "fiber","hexon", "fiber","hexon", "fiber","hexon", "fiber","hexon", "short fiber", "long fiber", "hexon", "fiber", "hexon", "long fiber", "short fiber"
)

# add the correct species and product name to the GC content df
L3_L5_entireSeq_GC <- L3_L5_entireSeq_GC %>%
  map_df(as_tibble) %>%
  mutate(species = repeat_vector(had_species_v, number_of_windows)) %>%
  mutate(product = repeat_vector(product_v, number_of_windows)) %>%
  mutate(species = factor(species, levels = c("HAdC5", "HAdC2", "HAdB7", "HAdD26", "HAdE4", "HAdF41", "Bat Ad2", "Avian Ad celo")))

ggplot(data=L3_L5_entireSeq_GC%>%filter(product=="hexon"), aes(x=abs_loc, y=GC, color=species)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0.3, 0.8), expand=F) +
  scale_color_manual(values=ad_colors) +
  facet_wrap(~species) +
  ggtitle("hexon - start 0 end 1")
ggsave("plots/hexon_seperate_centered_0start1end.png", dpi=300, height = 8, width = 8)
ggsave("plots/hexon_seperate_centered_0start1end.svg", dpi=300, height = 8, width = 8)
ggplot(data=L3_L5_entireSeq_GC%>%filter(product=="hexon"), aes(x=abs_loc, y=GC, color=species)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0.3, 0.8), expand=F) +
  scale_color_manual(values=ad_colors) +
  ggtitle("hexon - start 0 end 1")
ggsave("plots/hexon_centered_0start1end.png", dpi=300)
ggsave("plots/hexon_centered_0start1end.svg", dpi=300)

ggplot(data=L3_L5_entireSeq_GC%>%filter(product!="hexon")%>%filter(product!="short fiber"), aes(x=abs_loc, y=GC, color=species)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0.2, 0.7), expand=F) +
  scale_color_manual(values=ad_colors) +
  facet_wrap(~species) +
  ggtitle("fiber - start 0 end 1")
ggsave("plots/fiber_seperate_centered_0start1end.png", dpi=300, height = 8, width = 8)
ggsave("plots/fiber_seperate_centered_0start1end.svg", dpi=300, height = 8, width = 8)
ggplot(data=L3_L5_entireSeq_GC%>%filter(product!="hexon")%>%filter(product!="short fiber"), aes(x=abs_loc, y=GC, color=species)) +
  geom_line(size=1) +
  theme_bw() +
  xlab("position") + 
  ylab("GC") +
  coord_cartesian(ylim=c(0.2, 0.7), expand=F) +
  scale_color_manual(values=ad_colors) +
  ggtitle("fiber - start 0 end 1")
ggsave("plots/fiber_centered_0start1end.png", dpi=300)
ggsave("plots/fiber_centered_0start1end.svg", dpi=300)

# write gene sequences
L3_L5 %>%
  rowwise() %>%
  mutate(sequence = toString(viruses_fasta[[species]][[1]][start:end])) %>% # get the sequence of every section
  write_delim("hexon_fiber_positions_sequence.txt")

