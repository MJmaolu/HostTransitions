# FIG 6

# BENEFICIAL MUTATIONS
# MJ OLMO-UCEDA
# 25/03/16
################################################################################
# LOAD S STIMATES
jin1.gradual <- readxl::read_xlsx("data/Fitness estimates SNP.xlsx",
                                  sheet = "jin1->eds81 gradual")
jin1.sudden <- readxl::read_xlsx("data/Fitness estimates SNP.xlsx",
                                  sheet = "jin1->eds81 intermediate")
Gy0.gradual <- readxl::read_xlsx("data/Fitness estimates SNP.xlsx",
                                  sheet = "Gy->Oy gradual")
Gy0.sudden <- readxl::read_xlsx("data/Fitness estimates SNP.xlsx",
                                 sheet = "Gy->Oy intermediate")

all.s.annot <- rbind(jin1.gradual %>%
  mutate(treatment = "Gradual",
         host = "jin1 to eds8-1") %>%
    left_join(.,
              snv.transmitted %>%
                mutate(Locus = id_snv_nod) %>%
                dplyr::select(c("Locus", 
                                "products_nod", 
                                "in_CDS_nod", 
                                "nonsynonymous_nod", 
                                "aa_effect")) %>%
                distinct(),
              by = "Locus") %>%
    mutate(protein = factor(products_nod, levels = names(color.proteins))),
Gy0.gradual %>%
  mutate(treatment = "Gradual",
         host = "Gy-0 to Oy-0") %>%
  left_join(.,
            snv.transmitted %>%
              mutate(Locus = id_snv_nod) %>%
              dplyr::select(c("Locus", 
                              "products_nod", 
                              "in_CDS_nod", 
                              "nonsynonymous_nod", 
                              "aa_effect")) %>%
              distinct(),
            by = "Locus") %>%
  mutate(protein = factor(products_nod, levels = names(color.proteins))),
jin1.sudden %>%
  mutate(treatment = "Sudden",
         host = "jin1 to eds8-1") %>%
  left_join(.,
            snv.transmitted %>%
              mutate(Locus = id_snv_nod) %>%
              dplyr::select(c("Locus", 
                              "products_nod", 
                              "in_CDS_nod", 
                              "nonsynonymous_nod", 
                              "aa_effect")) %>%
              distinct(),
            by = "Locus") %>%
  mutate(protein = factor(products_nod, levels = names(color.proteins))),
Gy0.sudden %>%
  mutate(treatment = "Sudden",
         host = "Gy-0 to Oy-0") %>%
           left_join(.,
                     snv.transmitted %>%
                       mutate(Locus = id_snv_nod) %>%
                       dplyr::select(c("Locus", 
                                       "products_nod", 
                                       "in_CDS_nod", 
                                       "nonsynonymous_nod", 
                                       "aa_effect")) %>%
                       distinct(),
                     by = "Locus") %>%
           mutate(protein = factor(products_nod, levels = names(color.proteins))))



# S density distributions
s.distribution <- all.s.annot %>%
  filter(!is.na(nonsynonymous_nod)) %>%
  ggplot(.,
         aes(x = median_s,
             y = as.factor(aa_effect),
             fill = treatment,
             shape = as.factor(aa_effect)
         )) +
  ggridges::geom_density_ridges(aes(#linetype = aa_effect,
                                    point_shape = as.factor(aa_effect)                                    ),
        alpha = .4,
        jittered_points = T,
        point_size = 1, 
        point_alpha = .5, 
        alpha = 0.5) +
  facet_grid(host ~ treatment,
             #scales = "free_x"
             ) +
  geom_text(data = . %>%
              filter(median_s > 0.05),
            aes(label = ifelse(aa_effect == "synonymous",
                               Locus,
                               nonsynonymous_nod),
                color = as.factor(protein))) +
  
  geom_text(data = . %>%
              filter(median_s < - 0.05),
            aes(label = ifelse(aa_effect == "synonymous",
                               Locus,
                               nonsynonymous_nod),
                color = as.factor(protein))) +
  #coord_flip() +
  scale_discrete_manual(aesthetics = "point_shape", values = c("synonymous" = 21,
                                                               "nonsynonymous" = 23)) +
  scale_color_manual(values = color.proteins) +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  scale_shape_manual(values = c("synonymous" = 21,
                                "nonsynonymous" = 23)) +
  scale_linetype_manual(values = c("synonymous" = "solid",
                                   "nonsynonymous" = "dashed")) +
  labs(y = "Mutation effect",
       x = "Selection coefficient (s)") +
  
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

all.s.annot %>%
  filter(!is.na(nonsynonymous_nod)) %>%
  ggplot(.,
         aes(y = median_s,
             x = as.factor(aa_effect),
             fill = treatment,
             shape = as.factor(aa_effect)
         )) +
  # ggridges::geom_density_ridges(
  #   aes(linetype = aa_effect),
  #   alpha = .4,
  #   jittered_points = F,
  #   #point_size = 2,
  # ) +
  
  geom_violin(aes(fill = as.factor(treatment)),
               width = .4,
              alpha = .6) +
  # geom_boxplot(aes(fill = as.factor(treatment)),
  #              width = .4) +
  geom_jitter(aes(),
                  width = .1) +
  geom_text(data = . %>%
              filter(median_s < -0.1),
            aes(label = ifelse(aa_effect == "synonymous",
                               Locus,
                               nonsynonymous_nod),
                color = as.factor(protein))) +
  # ggbeeswarm::geom_quasirandom(aes(),
  #                              #height = .4,
  #                              width = .2,
  #                              size = 1) +  
  facet_grid(host ~ treatment,
             scales = "free") +
  scale_color_manual(values = color.proteins) +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  scale_shape_manual(values = c("synonymous" = 21,
                                "nonsynonymous" = 23)) +
  scale_linetype_manual(values = c("synonymous" = "solid",
                                   "nonsynonymous" = "dashed")) +
  labs(y = "Selection coefficient (s)",
       x = "") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black", linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank())



(AF.trajectories +
  geom_text(data = snv.transmitted %>%
              mutate(treatment = ifelse(treatment == "gradual",
                                        "Gradual",
                                        "Sudden"),
                     host = ifelse(original_genotype == "Gy-0",
                                   "Gy-0 to Oy-0",
                                   "jin1 to eds8-1")) %>%
              filter(transmitted == "yes",
                     AF > 0.2),
            aes(label = ifelse(aa_effect == "synonymous",
                               id_snv_nod,
                               nonsynonymous_nod),
                color = products_nod))) / s.distribution +
  plot_layout(heights = c(2, 1))

################################################################################
# AF TRAJECTORIES
################################################################################
AF.trajectories <- snv %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  #filter(AF > 0.05) %>%
  ggplot(.,
         aes(y = AF,
             x = as.factor(passage),
             fill = products_nod
             #group = interaction(host, treatment)
         )) +
  ggbeeswarm::geom_quasirandom(aes(shape = aa_effect,
                                   #fill = products_nod
                                   #color = aa_effect
  ),
  size = 2,
  ) +
  geom_line(aes(group = id_snv_nod,
                color = products_nod
  )) +
  scale_color_manual(values = color.proteins) +
  scale_fill_manual(values = color.proteins) +
  facet_grid(host ~treatment) +
  #coord_flip() +
  scale_shape_manual(values = c("synonymous" = 21,
                                "nonsynonymous" = 23)) +
  labs(x = "Passage",
       y = "Alelle Frequency") +
  theme_bw() +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black",
                                  size = 14),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()
  )


AF.trajectories +
  geom_text(data = snv.transmitted %>%
              mutate(treatment = ifelse(treatment == "gradual",
                                        "Gradual",
                                        "Sudden"),
                     host = ifelse(original_genotype == "Gy-0",
                                   "Gy-0 to Oy-0",
                                   "jin1 to eds8-1")) %>%
              filter(transmitted == "yes",
                     AF > 0.2),
            aes(label = ifelse(aa_effect == "synonymous",
                               id_snv_nod,
                               nonsynonymous_nod),
                color = products_nod))
# AF increasing
snv.transmitted %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  filter(AF > 0.05,
         #transmitted == "yes"
  ) %>%
  group_by(treatment,
           host,
           passage) %>%
  summarise(sumAF = sum(AF)) %>% #writexl::write_xlsx("sumAF.xlsx")
  ggplot(.,
         aes(x = as.factor(passage),
             y = sumAF,
             fill = as.factor(treatment),
             #linetype = as.factor(aa_effect)
         )) +
  geom_point(aes(shape = as.factor(host),
                 fill = as.factor(treatment)),
             size = 2,
             alpha = .6) +
  geom_line(aes(group = interaction(treatment, host),
                color = as.factor(treatment),
                linetype = as.factor(host)
  )
  ) +
  #facet_grid(host ~ .) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  scale_shape_manual(values = c("Gy-0 to Oy-0" = 24,
                                "jin1 to eds8-1" = 22)) +
  # scale_linetype_manual(values = c("synonymous" = "solid",
  #                                  "nonsynonymous" = "dashed")) +
  labs(x = "Passage",
       y = "Sum of Alelle Frequencies") +
  theme_bw() +
  theme(#legend.position = "none", 
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "transparent"),
    strip.text = element_text(color = "black"),
    panel.grid = element_line(colour = "black",
                              linewidth = 0.05),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )


snv.transmitted %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  filter(AF > 0.05,
         #transmitted == "yes"
  ) %>%
  group_by(treatment,
           host,
           passage,
           transmitted) %>%
  summarise(sumAF = sum(AF)) %>%
  ggplot(.,
         aes(x = as.factor(passage),
             y = sumAF,
             fill = as.factor(treatment),
             linetype = as.factor(transmitted)
         )) +
  geom_point(aes(shape = as.factor(host)),
             size = 2,
             alpha = .6) +
  geom_line(aes(group = interaction(treatment, transmitted),
                color = as.factor(treatment),
                linetype = as.factor(transmitted)
  )
  ) +
  facet_grid(host ~ .) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                "Gradual" = "#BFBC25")) +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  # scale_shape_manual(values = c("synonymous" = 21,
  #                               "nonsynonymous" = 23)) +
  # scale_linetype_manual(values = c("synonymous" = "solid",
  #                                  "nonsynonymous" = "dashed")) +
  labs(x = "Passage",
       sumAF = "Alelle Frequency") +
  theme_bw() +
  theme(#legend.position = "none", 
    axis.text = element_text(color = "black"),
    strip.background = element_rect(fill = "transparent"),
    strip.text = element_text(color = "black"),
    panel.grid = element_line(colour = "black",
                              linewidth = 0.05),
    panel.grid.minor.x = element_blank(),
    panel.grid.major.x = element_blank(),
    panel.grid.minor.y = element_blank()
  )



snv.transmitted %>%
  mutate(treatment = ifelse(treatment == "gradual",
                            "Gradual",
                            "Sudden"),
         host = ifelse(original_genotype == "Gy-0",
                       "Gy-0 to Oy-0",
                       "jin1 to eds8-1")) %>%
  filter(AF > 0.05,
         #transmitted == "yes"
         ) %>%
  group_by(treatment,
           host,
           passage,
           aa_effect) %>%
  summarise(sumAF = sum(AF)) %>%
  ggplot(.,
         aes(x = as.factor(passage),
             y = sumAF,
             fill = as.factor(treatment),
             #linetype = as.factor(aa_effect)
             )) +
  geom_point(aes(shape = as.factor(aa_effect)),
             size = 2,
             alpha = .6) +
  geom_line(aes(group = interaction(aa_effect, treatment),
                color = as.factor(treatment),
                linetype = as.factor(aa_effect)
                )
            ) +
  facet_grid(host ~ .) +
  scale_color_manual(values = c("Sudden" = "#5A87BB",
                                                     "Gradual" = "#BFBC25")) +
  scale_fill_manual(values = c("Sudden" = "#5A87BB",
                               "Gradual" = "#BFBC25")) +
  scale_shape_manual(values = c("synonymous" = 21,
                                "nonsynonymous" = 23)) +
  scale_linetype_manual(values = c("synonymous" = "solid",
                                   "nonsynonymous" = "dashed")) +
  labs(x = "Passage",
       sumAF = "Alelle Frequency") +
  theme_bw() +
  theme(#legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.minor.y = element_blank()
  )
  
