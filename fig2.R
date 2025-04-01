# fig 2: Rates of disease-related traits and VL
# MJ OLMO-UCEDA
################################################################################

treatment.colors <- c("Gradual" = "#BFBC25",
                      "gradual" = "#BFBC25",
                      "4" = "#BFBC25",
                      "Sudden" = "#5A87BB",
                      "sudden" = "#5A87BB",
                      "2" = "#5A87BB")
## LOAD DATA
infectivity <- read_sav("data/Infectivity.sav")
symptoms <- read_sav("data/Symptoms 12 dpi.sav")
VL <- read_sav("data/Evolution viral load.sav")

vl <- read_xlsx("data/gradualTransitions_withQuantification_241215.xlsx")

## INFECTIVITY
infectivity$Trans
# value               label
# 0            Constant
# 1        Sudden early
# 2 Sudden intermediate
# 3         Sudden late
# 4             Gradual
infectivity$Order
# 1   Gy-0 -> Oy-0
# 2   Oy-0 -> Gy-0
# 3 jin1 -> eds8-1
# 4 eds8-1 -> jin1
# 5           Gy-0
# 6           Oy-0
# 7           jin1
# 8         eds8-1

p1 <- infectivity %>%
  filter(Trans %in% c(2, 4),
         Order %in% c(3, 1)) %>%
  ggplot(.,
         aes(x = as.numeric(pase),
             y = as.numeric(NumInfec),
             group = interaction(Trans, Order),
             colour = as.factor(Trans),
             #linetype = as.factor(phenotype)
         )) +
  geom_line(linewidth = 1,
            alpha = .8) +
  # stat_summary(aes(group = Trans),
  #              geom = "line") +
  facet_grid(Order~.) +
  scale_color_manual(values = treatment.colors) +
  labs(x = "Passage",
       y = "Number of infected plants 12 dpi") +
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

p2 <- infectivity %>%
  filter(Trans %in% c(2, 4),
         Order %in% c(3, 1)) %>%
  ggplot(.,
         aes(x = as.numeric(pase),
             y = as.numeric(AUDPS),
             group = interaction(Trans, Order),
             colour = as.factor(Trans),
             #linetype = as.factor(phenotype)
         )) +
  geom_line(linewidth = 1,
            alpha = .8) +
  # stat_summary(aes(group = Trans),
  #              geom = "line") +
  facet_grid(Order~.) +
  scale_color_manual(values = treatment.colors) +
  labs(x = "Passage",
       y = "AUDPS") +
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


## SYMPTOMS
symptoms$Trans
# value               label
# 1        Sudden early
# 2 Sudden intermediate
# 3         Sudden late
# 4             Gradual
symptoms$Order
# value          label
# 1   Gy-0 -> Oy-0
# 2   Oy-0 -> Gy-0
# 3 jin1 -> eds8-1
# 4 eds8-1 -> jin1


p3 <- symptoms %>%
  filter(Trans %in% c(2, 4),
         Order %in% c(3, 1)) %>%
  ggplot(.,
         aes(x = as.numeric(Passage),
             y = as.numeric(Symptoms),
             group = interaction(Trans, Order),
             colour = as.factor(Trans),
             #linetype = as.factor(phenotype)
         )) +
  geom_jitter(size = .2,
              height = .1,
              width = .1,
              alpha = .8) +
  stat_summary(aes(group = Trans),
               geom = "line",
               linewidth = 1,
               alpha = .8) +
  stat_summary(aes(group = Trans),
               geom = "errorbar",
               width = .1) +
  facet_grid(Order~.) +
  scale_color_manual(values = treatment.colors) +
  scale_fill_manual(values = treatment.colors) +
  labs(x = "Passage",
       y = "Mean symptoms 12 dpi") +
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

## VL
p4 <- VL %>%
  filter(treatment %in% c("2", "3"),
         genotype %in% c("2", "3")
         ) %>%
  ggplot(.,
       aes(x = as.numeric(passage),
           y = as.numeric(VPMB),
           group = interaction(genotype, phenotype),
           colour = as.factor(treatment),
           #linetype = as.factor(phenotype)
           )) +
  geom_line(aes(group = treatment),
            linewidth = 1,
            alpha = .8) +
  facet_grid(genotype~.) +
  scale_color_manual(values = c("3" = "#5A87BB",
                                "2" = "#BFBC25")) +
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

VL$treatment
# Labels:
#   value        label
# 1        early
# 2      gradual
# 3 intermediate
# 4         late

VL$phenotype
# value     label
# 1 Resistant
# 2 Sensitive

VL$genotype
# Labels:
#   value  label
# 1 eds8-1
# 2   Gy-0
# 3  jin-1
# 4   Oy-0
################################################################################
(p1 + p2 ) / (p3 + p4)


################################################################################
allTraits <- merge(infectivity  %>% 
        filter(Trans %in% c(2, 4),
               Order %in% c(3, 1)) %>%
        mutate(comb = paste0(ifelse(Trans == 2,
                                    "sudden",
                                    "gradual"),
                             ".",
                             ifelse(Order == 1,
                                    "Gy0", 
                                    "jin1"),
                             "_",
                             pase)),
      symptoms %>%
        filter(Trans %in% c(2, 4),
               Order %in% c(3, 1)) %>%
            group_by(Trans, Order, Passage) %>%
            summarise(meanSymptoms = mean(Symptoms)) %>%
        mutate(comb = paste0(ifelse(Trans == 2,
                                    "sudden",
                                    "gradual"),
                             ".",
                             ifelse(Order == 1,
                                    "Gy0", 
                                    "jin1"),
                             "_",
                             Passage)),
      by = "comb") %>%
  left_join(.,
            VL %>%
              filter(treatment %in% c("2", "3"),
                     genotype %in% c("2", "3")
              ) %>%
              mutate(comb = paste0(ifelse(treatment == 2,
                                          "gradual",
                                          "sudden"),
                                   ".",
                                   ifelse(genotype == 2,
                                          "Gy0", 
                                          "jin1"),
                                   "_",
                                   passage))) %>%
  arrange(Passage) 



t1 <- allTraits %>%
  ggplot(.,
         aes(x = as.numeric(AUDPS),
             y = as.numeric(NumInfec),
             group = as.factor(Trans.x),
             color = as.factor(Trans.x))) +
  geom_point(aes(size = pase),
             alpha = .8) +
  geom_smooth(method = "lm",
              se = F,
              linetype = "dotted") +
  geom_path(arrow = arrow(type = "closed",
                          angle = 15)) +
  scale_color_manual(values = treatment.colors) +
  theme_bw() +
  labs(x = "AUDPS",
       y = "Number infected plants") +
  facet_grid(Order.x ~.) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

t2 <- allTraits %>%
  filter(!is.na(VPMB)) %>%
  ggplot(.,
         aes(x = as.numeric(VPMB),
             y = as.numeric(meanSymptoms),
             group = as.factor(Trans.x),
             color = as.factor(Trans.x))) +
  geom_point(#aes(size = pase),
    size = 2.5,
    alpha = .8) +
  
  # geom_smooth(method = "lm",
  #             se = F,
  #             linetype = "dotted") +
  geom_path(arrow = arrow(type = "open",
                          angle = 15,
                          length = unit(0.3, "cm"),
                          ends = "last")) +
  scale_color_manual(values = treatment.colors) +
  theme_bw() +
  labs(x = "VPMB",
       y = "Mean symptoms") +
  facet_grid(Order.x ~.) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

t1 + t2


t3 <- allTraits %>%
  ggplot(.,
         aes(x = as.numeric(meanSymptoms),
             y = as.numeric(AUDPS),
             group = as.factor(Trans.x),
             color = as.factor(Trans.x))) +
  geom_point(aes(size = pase),
             alpha = .8) +
  geom_smooth(method = "lm",
              se = F,
              linetype = "dotted") +
  geom_path(arrow = arrow(type = "closed",
                          angle = 15)) +
  scale_color_manual(values = treatment.colors) +
  theme_bw() +
  labs(x = "Mean Symptoms",
       y = "AUDPS") +
  facet_grid(Order.x ~.) +
  #scale_x_continuous(breaks = seq(1,10,1)) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

t4 <- allTraits %>%
  ggplot(.,
         aes(x = as.numeric(meanSymptoms),
             y = as.numeric(NumInfec),
             group = as.factor(Trans.x),
             color = as.factor(Trans.x))) +
  geom_point(aes(size = pase),
             alpha = .8) +
  geom_smooth(method = "lm",
              se = F,
              linetype = "dotted") +
  geom_path(arrow = arrow(type = "closed",
                          angle = 15)) +
  scale_color_manual(values = treatment.colors) +
  theme_bw() +
  labs(x = "Mean Symptoms",
       y = "Number infected plants") +
  facet_grid(Order.x ~.) +
  #scale_x_continuous(breaks = seq(1,10,1)) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

t5 <- allTraits %>%
  filter(!is.na(VPMB)) %>%
  ggplot(.,
         aes(x = as.numeric(VPMB),
             y = as.numeric(NumInfec),
             group = as.factor(Trans.x),
             color = as.factor(Trans.x))) +
  geom_point(#aes(size = pase),
    size = 2.5,
    alpha = .8) +

  # geom_smooth(method = "lm",
  #             se = F,
  #             linetype = "dotted") +
  geom_path(arrow = arrow(type = "open",
                          angle = 15,
                          length = unit(0.3, "cm"),
                          ends = "last")) +
  scale_color_manual(values = treatment.colors) +
  theme_bw() +
  labs(x = "VPMB",
       y = "Number infected plants") +
  facet_grid(Order.x ~.) +
  #scale_x_continuous(breaks = seq(1,10,1)) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

t6 <- allTraits %>%
  filter(!is.na(VPMB)) %>%
  ggplot(.,
         aes(x = as.numeric(VPMB),
             y = as.numeric(AUDPS),
             group = as.factor(Trans.x),
             color = as.factor(Trans.x))) +
  geom_point(#aes(size = pase),
    size = 2.5,
    alpha = .8) +
  
  # geom_smooth(method = "lm",
  #             se = F,
  #             linetype = "dotted") +
  geom_path(arrow = arrow(type = "open",
                          angle = 15,
                          length = unit(0.3, "cm"),
                          ends = "last")) +
  scale_color_manual(values = treatment.colors) +
  theme_bw() +
  labs(x = "VPMB",
       y = "AUDPS") +
  facet_grid(Order.x ~.) +
  #scale_x_continuous(breaks = seq(1,10,1)) +
  theme(legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

(t1 + t2) / (t3 + t4) / (t5 + t6)

t2 + t5 + t6

################################################################################
VL %>%
  filter(treatment %in% c("1", "2", "3"),
          genotype %in% c("2", "3")
  ) %>%
  ggplot(.,
         aes(x = as.numeric(passage),
             y = as.numeric(VPMB),
             group = interaction(genotype, phenotype),
             colour = as.factor(treatment),
             #linetype = as.factor(phenotype)
         )) +
  geom_line(aes(group = treatment,
                linetype = as.factor(genotype)),
            linewidth = 1,
            alpha = .8) +
  facet_grid(genotype~.) +
  scale_color_manual(values = c("3" = "#5A87BB",
                                "2" = "#BFBC25")) +
  labs(x = "Passage",
       y = "Viral reads per million base") +
  theme_bw() +
  scale_x_continuous(breaks = seq(1,10,1)) +
  theme(#legend.position = "none", 
        axis.text = element_text(color = "black"),
        strip.background = element_rect(fill = "transparent"),
        strip.text = element_text(color = "black"),
        panel.grid = element_line(colour = "black",
                                  linewidth = 0.05),
        panel.grid.minor.x = element_blank(),
        panel.grid.major.x = element_blank(),
  )

symptoms %>%
  filter(Trans %in% c(1, 2, 4),
         Order %in% c(3, 1)) %>%
  ggplot(.,
         aes(x = as.numeric(Passage),
             y = as.numeric(Symptoms),
             group = interaction(Trans, Order),
             colour = as.factor(Trans),
             #linetype = as.factor(phenotype)
         )) +
  geom_jitter(size = .2,
              height = .1,
              width = .1,
              alpha = .8) +
  stat_summary(aes(group = Trans),
               geom = "line",
               linewidth = 1,
               alpha = .8) +
  stat_summary(aes(group = Trans),
               geom = "errorbar",
               width = .1) +
  facet_grid(Order~.) +
  scale_color_manual(values = treatment.colors) +
  scale_fill_manual(values = treatment.colors) +
  labs(x = "Passage",
       y = "Mean symptoms 12 dpi") +
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

infectivity %>%
  filter(Trans %in% c(1, 2, 4),
         Order %in% c(3, 1)) %>%
  ggplot(.,
         aes(x = as.numeric(pase),
             y = as.numeric(AUDPS),
             group = interaction(Trans, Order),
             colour = as.factor(Trans),
             #linetype = as.factor(phenotype)
         )) +
  geom_line(linewidth = 1,
            alpha = .8) +
  # stat_summary(aes(group = Trans),
  #              geom = "line") +
  facet_grid(Order~.) +
  scale_color_manual(values = treatment.colors) +
  labs(x = "Passage",
       y = "Number of infected plants 12 dpi") +
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
# Fig2
