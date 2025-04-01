# FIG3

# SNVs related: Richness, AF densities and hotspots
# MJ OLMO-UCEDA
# 25/03/14
################################################################################
tumv.proteins <- read.csv("TuMV_Anc.gff",
                          sep = "\t",
                          skip = 2) %>%
  filter(region == "mature_protein_region_of_CDS") %>%
  dplyr::select(c("X1", "X9832", "ID.TuMV_Anc.1..9832.Dbxref.taxon.12230.gbkey.Src.genome.genomic.mol_type.genomic.RNA")) %>%
  rename("start" = "X1",
         "end" = "X9832",
         "description" = "ID.TuMV_Anc.1..9832.Dbxref.taxon.12230.gbkey.Src.genome.genomic.mol_type.genomic.RNA") %>%
  rowwise() %>%
  mutate(protein = str_split(description, "product=")[[1]][2]) %>%
  dplyr::select(c("protein", "start", "end", "description"))
  

allSNVs <- readxl::read_xlsx("data/All variants LoFreq.xlsx",
                             sheet = "All data")
allSNVs$treatment %>% unique()
allSNVs$original_genotype %>% unique()

mut.fixed <-  read_sav("data/Number of mutations.sav")
mut.fixed %>%
  filter(timing %in% c(2,4)) %>%
  ggplot(.,
         aes(x = paseevo,
             y = nums,
             group = interaction(trans, timing),
             colour = as.factor(timing),
             #linetype = as.factor(phenotype)
         )) +
  geom_line(aes(linetype = as.factor(trans)),
            linewidth = 1,
            alpha = .8) +
  #facet_grid(timing~.) +
  scale_color_manual(values = c("2" = "#5A87BB",
                                "4" = "#BFBC25")) +
  labs(x = "Passage",
       y = "Viral reads per million base") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1,10,1)) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

################################################################################
# ONLY SELECTED COMBIANTIONS 
snv <- allSNVs %>%
  filter(DP > 20,
         AF > 0.01,
         treatment %in% c("gradual", "intermediate"),
         original_genotype %in% c("Gy-0", "jin1")
  )


richness <- snv %>%
  group_by(treatment,
           original_genotype,
           passage) %>%
  summarise(snv.richness = n()) %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  ggplot(.,
         aes(x = passage,
             y = snv.richness,
             group = interaction(treatment, original_genotype),
         colour = as.factor(treatment),
         #linetype = as.factor(phenotype)
  )) +
  # NUMBER OF TOTAL SNV
  geom_line(#aes(linetype = original_genotype),
            linewidth = 1,
            alpha = .8) +
  # NUMBER OF SNVs NEW TO EACH PASSAGE
  geom_line(data = snv.transmitted %>%
              group_by(treatment,
                       original_genotype,
                       passage,
                       transmitted) %>%
              summarise(snv.transmitted = n()) %>%
              filter(transmitted == "no") %>%
              mutate(treatment = ifelse(treatment == "gradual",
                                                                   "Gradual",
                                                                   "Sudden"),
                                                host = ifelse(original_genotype == "Gy-0",
                                                              "Gy-0 to Oy-0",
                                                              "jin1 to eds8-1")),
            aes(x = passage,
                y = snv.transmitted,
                group = as.factor(treatment),
                colour = as.factor(treatment)),
            linetype = "dashed",
            inherit.aes = T) +
  facet_grid(host~.) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  labs(x = "Passage",
       y = "Richness") +
  theme_bw() +
  scale_x_continuous(breaks = c(1,4,7,10)) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

richness.transmitted <- snv %>%
  group_by(treatment,
           original_genotype,
           passage) %>%
  summarise(snv.richness = n()) %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  ggplot(.,
         aes(x = passage,
             y = snv.richness,
             group = interaction(treatment, original_genotype),
             colour = as.factor(treatment),
             #linetype = as.factor(phenotype)
         )) +
  # NUMBER OF TOTAL SNV
  geom_line(#aes(linetype = original_genotype),
    linewidth = 1,
    alpha = .8) +
  # NUMBER OF SNVs NEW TO EACH PASSAGE
  geom_line(data = snv.transmitted %>%
              group_by(treatment,
                       original_genotype,
                       transmitted) %>%
              summarise(snv.transmitted = n()) %>%
              filter(transmitted == "yes") %>%
              mutate(treatment = ifelse(treatment == "gradual",
                                        "Gradual",
                                        "Sudden"),
                     host = ifelse(original_genotype == "Gy-0",
                                   "Gy-0 to Oy-0",
                                   "jin1 to eds8-1")),
            aes(x = passage,
                y = snv.transmitted,
                group = as.factor(treatment),
                colour = as.factor(treatment)),
            linetype = "solid",
            linewidth = .5,
            inherit.aes = T) +
  geom_line(data = snv.transmitted %>%
              group_by(treatment,
                       original_genotype,
                       passage,
                       transmitted) %>%
              summarise(snv.transmitted = n()) %>%
              filter(transmitted == "no") %>%
              mutate(treatment = ifelse(treatment == "gradual",
                                        "Gradual",
                                        "Sudden"),
                     host = ifelse(original_genotype == "Gy-0",
                                   "Gy-0 to Oy-0",
                                   "jin1 to eds8-1")),
            aes(x = passage,
                y = snv.transmitted,
                group = as.factor(treatment),
                colour = as.factor(treatment)),
            linetype = "dashed",
            linewidth = .5,
            inherit.aes = T) +
  facet_grid(host~.) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  labs(x = "Passage",
       y = "Richness") +
  theme_bw() +
  scale_x_continuous(breaks = c(1,4,7,10)) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )
  


AF.density <- snv %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  #filter(AF > 0.05) %>%
  ggplot(.,
         aes(x = AF,
             y = as.factor(passage),
             fill = treatment
             #group = interaction(host, treatment)
             )) +
  ggridges::geom_density_ridges(alpha = .6) +
  facet_grid(host ~treatment) +
  #coord_flip() +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  scale_shape_manual(values = c("synonymous" = 20,
                                "nonsynonymous" = 18)) +
  labs(y = "Passage",
       x = "Alelle Frequency") +
  theme_bw() +
  theme(legend.position = "none", 
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "transparent"),
    strip.text = element_text(color = "black"),
    panel.grid = element_line(colour = "black",
                              linewidth = 0.05),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )

richness + AF.density +
  plot_layout(widths = c(1,3))

richness + snv.per.cistron +
  plot_layout(widths = c(1,3))

richness + AF.transmitted +
  plot_layout(widths = c(1,3))

richness.transmitted + AF.density +
  plot_layout(widths = c(1,3))

richness.transmitted + AF.transmitted +
  plot_layout(widths = c(1,3))

################################################################################
(snv %>%
    mutate(treatment = ifelse(treatment == "gradual",
                              "Gradual",
                              "Sudden"),
           host = ifelse(original_genotype == "Gy-0",
                         "Gy-0 to Oy-0",
                         "jin1 to eds8-1")) %>%
    filter(AF > 0.05) %>%
    ggplot(.,
           aes(x = POS,
               y = (passage),
               fill = products_nod,
               shape = aa_effect)) +
   # geom_line(aes(group = id_snv_nod),
   #           linewidth = .1) +
   geom_point(aes(size = AF),
               alpha = .7) +
    scale_y_continuous(breaks = c(1,4,7,10),
                       expand = c(0.1, 0.1)) +
    facet_grid(host ~ treatment) +
    #coord_flip() +
    scale_fill_manual(values = color.proteins) +
    scale_shape_manual(values = c("synonymous" = 21,
                                  "nonsynonymous" = 23)) +
   labs(y = "Passage",
        x = "Genomic position") +
    theme_bw() +
    theme(legend.position = "bottom", 
      axis.text = element_text(color = "black"),
      strip.background = element_rect(fill = "transparent"),
      strip.text = element_text(color = "black"),
      panel.grid = element_line(colour = "black",
                                linewidth = 0.05),
      panel.grid.minor.x = element_blank(),
      panel.grid.major.x = element_blank(),
      panel.grid.minor.y = element_blank()
      
    ))
################################################################################
color.proteins <- c("5'UTR" = "white",
                    "P1" = "#abdfcf",
                    "HC-Pro" = "#5e6b6d",
                    "P3" = "#364649",
                    "6K1" = "#a66899",
                    "CI" = "#8ad4eb",
                    "6K2" = "#ff9665",
                    "NIa-VPg" = "#cccc30", 
                    "NIa-Pro" = "#f3c812",
                    "NIb" = "#fd6236",
                    "coat" = "#8c439a",
                    "3'UTR" = "white")

################################################################################
# MUTATIONS
################################################################################
snvs.per.combination <- list()


for (t in unique(snv$treatment)){
  for (g in unique(snv$original_genotype)){
    for (p in unique(snv$passage)){
      print(t)
      print(g)
      print(p)
      temp <- snv %>%
        filter(treatment == t,
               original_genotype == g,
               passage == p)
      print(length(temp$id_snv_nod))
      snvs.per.combination[[paste0(t, ".", g, "_", p)]] <- temp$id_snv_nod
    }
  }
}

nSNV <-  as.data.frame(array(lapply(snvs.per.combination, length),
      dim = c(4, 4))) %>%
  mutate(passage = c(1, 4, 7, 10))

colnames(nSNV) <- c("gradual.Gy-0", "gradual.jin1", "sudden.Gy-0", "sudden.jin1", "passage")


nSNV %>%
  pivot_longer(cols = -passage, 
               names_to = c("treatment", "genotype"), 
               names_sep = "\\.") %>%  
  unnest(value) %>%
  ggplot(.,
         aes(x = passage,
             y = value,
             group = interaction(treatment, genotype))) +
  geom_point()

for (i in seq(1, length(snvs.per.combination), 4)){
  print(names(snvs.per.combination[(i): (i+3)]))
  print(Reduce(intersect, snvs.per.combination[(i): (i+3)]))
}

Reduce(intersect, snvs.per.combination)

UpSetR::upset(fromList(snvs.per.combination),
              nsets = 16,
              order.by = "freq")

Heatmap(ComplexHeatmap::list_to_matrix(snvs.per.combination),
        column_order = NULL,
        col = c("white", "black"))


################################################################################
# NOVEL MUTATIONS
################################################################################
in_previous_passages <- function(geno, treat, pass, id_snv){
  # P4,7,10
  if (pass > 1){
    df.sub <- snv %>%
      filter(original_genotype == geno,
             treatment == treat,
             passage < pass)
    previous = ifelse(id_snv %in% df.sub$id_snv_nod,
           "yes",
           "no")
  } else {previous = "no"}
  return(previous)
}


snv.transmitted <- snv %>%
  rowwise() %>%
  mutate(transmitted = in_previous_passages(original_genotype, treatment, passage, id_snv_nod))


# snv.transmitted %>%
#   group_by(original_genotype, 
#            treatment,
#            passage) %>%
#   summarise(numberSNVs = n()) %>%
#   mutate(treatment = ifelse(treatment == "gradual",
#                             "Gradual",
#                             "Sudden"),
#          host = ifelse(original_genotype == "Gy-0",
#                        "Gy-0 to Oy-0",
#                        "jin1 to eds8-1")) %>%
#   writexl::write_xlsx(.,
#                       "SNVs_up5pc_totalRichness.xlsx")
# 
# snv.transmitted %>%
#   group_by(original_genotype, 
#            treatment,
#            passage,
#            transmitted) %>%
#   summarise(numberSNVs = n()) %>%
#   mutate(treatment = ifelse(treatment == "gradual",
#                             "Gradual",
#                             "Sudden"),
#          host = ifelse(original_genotype == "Gy-0",
#                        "Gy-0 to Oy-0",
#                        "jin1 to eds8-1")) %>%
#   writexl::write_xlsx(.,
#                       "SNVs_up5pc_richness.xlsx")


snv.transmitted %>%
  group_by(original_genotype, 
           treatment,
           passage,
           transmitted) %>%
  summarise(n = n()) %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  ggplot(.,
         aes(x = as.factor(passage),
             y = n,
             group = interaction(treatment, transmitted),
             color = as.factor(treatment),
         )) +
  geom_line(aes(linetype = transmitted)) +
  facet_grid(host~treatment) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  scale_linetype_manual(values = c("yes" = "solid",
                                   "no" = "dashed")) +
  labs(x = "Passage",
       y = "Richness") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )


snv.transmitted %>%
  group_by(original_genotype, 
           treatment,
           passage,
           transmitted) %>%
  summarise(n = n()) %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  ggplot(.,
         aes(x = as.factor(transmitted),
             y = n,
             group = interaction(treatment, passage),
             #group = interaction(original_genotype, transmitted),
             color = as.factor(treatment),
             
         )) +
  ggbeeswarm::geom_quasirandom(aes(fill = as.factor(passage)),
                               width = .1) +
  geom_line(linetype = "dotted") +
  stat_summary(aes(group = interaction(treatment, transmitted),
                   shape = as.factor(transmitted))) +
  facet_grid(host~treatment) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  labs(x = "Passage",
       y = "Richness") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

snv.transmitted %>%
  group_by(original_genotype, 
           treatment,
           passage,
           transmitted) %>%
  summarise(n = n()) %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  aov(n ~ as.factor(transmitted) * as.factor(treatment)* as.factor(original_genotype),
      data = .) %>%
  summary()


## NUMBER PER CISTRON
snv.per.cistron <- snv.transmitted %>%
  group_by(original_genotype, 
           treatment,
           passage,
           products_nod,
           aa_effect) %>%
  summarise(n = n()) %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1"),
         protein = factor(products_nod, 
                           levels = names(color.proteins))) %>%
  ggplot(.,
         aes(x = protein,
             #y = n,
             group = passage,
             fill = protein)) +
  geom_hline(yintercept = 0) +
  geom_line(aes(x = protein,
                y = ifelse(aa_effect == "nonsynonymous",
                           n,
                           -n),
                group = interaction(passage, aa_effect),
                linetype = as.factor(passage)),
            color = "black",
            linewidth = .2) +
  geom_segment(aes(x = protein,
                   xend = protein,
                   y = 0,
                   yend = ifelse(aa_effect == "nonsynonymous",
                                 n,
                                 -n),
                   color = protein)) +
  geom_point(aes(fill = protein,
                 y = ifelse(aa_effect == "nonsynonymous",
                            n,
                            -n),
                 shape = aa_effect)) +
  facet_grid(host~treatment) +
  scale_fill_manual(values = color.proteins) +
  scale_color_manual(values = color.proteins) +
  scale_shape_manual(values = c("synonymous" = 21,
                                "nonsynonymous" = 23)) +
  labs(x = "Cistron",
       y = "Nunber of mutations") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

snv$products_nod %in% names(color.proteins)

###################3
snv.effect <- snv.transmitted %>%
  group_by(original_genotype, 
           treatment,
           passage,
           aa_effect) %>%
  summarise(n = n()) %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  ggplot(.,
         aes(x = aa_effect,
             y = n,
             group = interaction(treatment, passage),
             #group = interaction(original_genotype, transmitted),
             color = as.factor(treatment),
             
         )) +
  geom_line(linetype = "dotted") +
  ggbeeswarm::geom_quasirandom(aes(shape = aa_effect,
                                   fill = as.factor(treatment)),
                               alpha = .8,
                               color = "black",
                               width = .1) +
  stat_summary(aes(group = interaction(aa_effect,
                                       treatment),
                   shape = aa_effect,
                   fill = treatment),
               color = "black") +
  stat_summary(aes(group = interaction(treatment)),
               geom = "line") +
  facet_grid(host~.) +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  scale_shape_manual(values = c("synonymous" = 21,
                                "nonsynonymous" = 23)) +
  labs(x = "Mutation effect",
       y = "Number of mutations") +
  scale_x_discrete(expand = c(0.1,0.1)) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

  
snv.effect + snv.per.cistron +
  plot_layout(widths = c(1,4))


{snv.effect + snv.per.cistron +
  plot_layout(widths = c(1, 4))} / snv.effect.cistron +
  plot_layout(heights = c(2,1))

snv.effect.cistron <- snv.transmitted %>%
  filter(in_CDS_nod == 1) %>%
  mutate(protein = factor(products_nod,
                          levels = names(color.proteins))) %>%
  group_by(original_genotype, 
           treatment,
           passage,
           protein,
           aa_effect) %>%
  summarise(n = n()) %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  ggplot(.,
         aes(x = aa_effect,
             y = n,
             group = interaction(treatment, passage),
             #group = interaction(original_genotype, transmitted),
             color = as.factor(treatment),
             
         )) +
  geom_line(linetype = "dotted") +
  ggbeeswarm::geom_quasirandom(aes(shape = aa_effect,
                                   fill = as.factor(treatment)),
                               size = .8,
                               color = "black",
                               alpha = .8,
                               width = .1) +
  stat_summary(aes(group = interaction(aa_effect,
                                       treatment),
                   shape = aa_effect,
                   fill = treatment),
               color = "black") +
  stat_summary(aes(group = interaction(treatment)),
               geom = "line") +
  facet_grid(host~protein) +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  scale_shape_manual(values = c("synonymous" = 21,
                                "nonsynonymous" = 23)) +
  labs(x = "Mutation effect",
       y = "Number of mutations") +
  scale_x_discrete(expand = c(0.1,0.1)) +
  scale_y_continuous(breaks = seq(1, 12, by = 2)) +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 90,
                                   vjust = 0.5),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

#############

AF.transmitted <- snv.transmitted %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  #filter(AF > 0.05) %>%
  ggplot(.,
         aes(x = AF,
             y = as.factor(passage),
             fill = treatment
         )) +
  ggridges::geom_density_ridges(aes(linetype = transmitted,
                                    point_shape = as.factor(transmitted)
    ),
                                 alpha = .4,
                                jittered_points = F,
                                # position = ggridges::position_points_jitter(width = 0.01, 
                                #                                             height = 0.1,
                                #                                             yoffset = ifelse(snv.transmitted$treatment == "Gradual", 
                                #                                                              0.1, 
                                #                                                              -0.1)),
                                #point_shape = 21, 
                                point_size = 2, 
                                point_alpha = .5, 
                                alpha = 0.5) +
  facet_grid(host ~ treatment) +
  #coord_flip() +x
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  scale_shape_manual(values = c("yes" = 3,
                                "no" = 25)) +
  
  scale_linetype_manual(values = c("yes" = "solid",
                                   "no" = "dashed")) +
  labs(y = "Passage",
       x = "Alelle Frequency") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()
  )



{snv.effect + snv.effect.cistron +
    plot_layout(widths = c(1, 6))}

# snv.transmitted %>%
#   filter(in_CDS_nod == 1) %>%
#   mutate(protein = factor(products_nod,
#                           levels = names(color.proteins))) %>%
#   group_by(original_genotype, 
#            treatment,
#            passage,
#            protein,
#            aa_effect) %>%
#   summarise(n = n()) %>%
#   mutate(treatment = ifelse(treatment == "gradual",
#                             "Gradual",
#                             "Sudden"),
#          host = ifelse(original_genotype == "Gy-0",
#                        "Gy-0 to Oy-0",
#                        "jin1 to eds8-1")) %>%
#   writexl::write_xlsx(.,
#                       "Mutations_per_cistron.xlsx")
# 
# 
# snv.transmitted %>%
#   mutate(protein = factor(products_nod,
#                           levels = names(color.proteins))) %>%
#   group_by(original_genotype, 
#            treatment,
#            passage,
#            aa_effect) %>%
#   summarise(n = n()) %>%
#   mutate(treatment = ifelse(treatment == "gradual",
#                             "Gradual",
#                             "Sudden"),
#          host = ifelse(original_genotype == "Gy-0",
#                        "Gy-0 to Oy-0",
#                        "jin1 to eds8-1")) %>%
#   writexl::write_xlsx(.,
#                       "Mutations_aaefect.xlsx")

# 
# snv.transmitted %>%
#   mutate(treatment = ifelse(treatment == "gradual",
#                             "Gradual",
#                             "Sudden"),
#          host = ifelse(original_genotype == "Gy-0",
#                        "Gy-0 to Oy-0",
#                        "jin1 to eds8-1")) %>%
#   dplyr::select(c("id_snv_nod", "AF", "products_nod", "nonsynonymous_nod", "aa_effect", "original_genotype", "passage", 
#                   "treatment", "transmitted")) %>%
#   writexl::write_xlsx("allSNVs.with.transmitted.info.xlsx")
