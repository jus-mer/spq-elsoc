# Predicted values and interaction plot for the within-person (WE) three-way
# interaction of subjective social status: own x family of origin x children,
# with children's status at its 25th and 75th percentiles.
# Run from the project root: Rscript output/interactions/int-triple-we-p25p75.R
#
# Y = market justice preferences (MJP) index
# X = own SSS (WE);  lines = origin SSS (WE);  facets = children's SSS (WE)
# Own and origin: Low / Mean / High = -1 SD / 0 / +1 SD of the within-person
# deviation. Children: P25 and P75 of the within-person deviation.

library(dplyr)
library(lme4)
library(lmerTest)
library(ggplot2)

load("input/data/proc/df_final.RData")

# Same analysis sample as the paper ------------------------------------------
d <- df_study1 %>%
  dplyr::select(-dplyr::any_of(c("prev_salud", "prev_salud_fill", "cambio_previ_salud", "inc_origen", "inc_destino", "traj_ingreso", "dir_ingreso"))) %>%
  na.omit()

# Between-person terms centred on the grand mean (person level)
p  <- d %>% distinct(idencuesta, ess_og_mean, ess_mean, ess_hijos_mean)
gm <- colMeans(p[, -1])

d <- d %>%
  mutate(og_be  = ess_og_mean    - gm["ess_og_mean"],
         own_be = ess_mean       - gm["ess_mean"],
         hij_be = ess_hijos_mean - gm["ess_hijos_mean"],
         og_we  = ess_og_cwc,
         own_we = ess_cwc,
         hij_we = ess_hijos_cwc)

# Model: three-way interactions (WE and BE) with the paper's controls --------
m <- lmer(mjp_index ~ ola_num + ola_2 + ideo + sex + age +
            dummy_decile10_mean + dummy_educ_mean +
            og_we * own_we * hij_we +
            og_be * own_be * hij_be +
            (1 + ola_num | idencuesta),
          data = d, REML = FALSE,
          control = lmerControl(optimizer = "bobyqa",
                                optCtrl = list(maxfun = 2e5)))

# Predicted values -----------------------------------------------------------
# Controls at their sample means, BE status terms at the grand mean (0).
# 95% CI from the fixed-effects covariance.
sd_we <- sapply(d[, c("own_we", "og_we")], sd)
lv    <- c(Low = -1, Mean = 0, High = 1)
q_hij <- quantile(d$hij_we, c(P25 = 0.25, P75 = 0.75), names = FALSE)
names(q_hij) <- c("P25", "P75")
print(round(q_hij, 3))

grid <- expand.grid(own = names(lv), og = names(lv), hij = names(q_hij),
                    stringsAsFactors = FALSE) %>%
  mutate(own_we = lv[own] * sd_we["own_we"],
         og_we  = lv[og]  * sd_we["og_we"],
         hij_we = q_hij[hij],
         og_be = 0, own_be = 0, hij_be = 0)

b  <- fixef(m)
V  <- as.matrix(vcov(m))
mm <- colMeans(model.matrix(m))
status <- c("own_we", "og_we", "hij_we", "own_be", "og_be", "hij_be")

X <- t(sapply(seq_len(nrow(grid)), function(i) {
  sapply(names(b), function(term) {
    parts <- strsplit(term, ":")[[1]]
    if (all(parts %in% status)) prod(unlist(grid[i, parts])) else mm[[term]]
  })
}))

grid$pred <- as.numeric(X %*% b)
grid$se   <- sqrt(rowSums((X %*% V) * X))
grid$lwr  <- grid$pred - qnorm(0.975) * grid$se
grid$upr  <- grid$pred + qnorm(0.975) * grid$se

lab     <- c(Low = "Low (-1 SD)", Mean = "Mean", High = "High (+1 SD)")
lab_hij <- setNames(sprintf("%s (%+.2f)", names(q_hij), q_hij), names(q_hij))

pred <- grid %>%
  transmute(own_we_value = round(own_we, 3),
            origin_we_value = round(og_we, 3),
            children_we_value = round(hij_we, 3),
            children_we = factor(lab_hij[hij], levels = lab_hij),
            origin_we   = factor(lab[og],  levels = lab),
            own_we      = factor(lab[own], levels = lab),
            predicted = round(pred, 3), se = round(se, 3),
            ci_low = round(lwr, 3), ci_high = round(upr, 3)) %>%
  arrange(children_we, origin_we, own_we)

write.csv(pred, "output/interactions/pred-int-triple-we-p25p75.csv", row.names = FALSE)
print(as.data.frame(pred))

# Plot -----------------------------------------------------------------------
# Origin status: blue = low, neutral grey = mean, brick = high, plus point
# shape so the lines are not told apart by colour alone.
ramp <- c("#2a78d6", "#7a7974", "#b5432e")
ink  <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid_col <- "#e1e0d9"

plot_df <- pred %>%
  mutate(children_we = factor(paste0("Children's SSS (WE): ", children_we),
                              levels = paste0("Children's SSS (WE): ", lab_hij)))
dodge <- position_dodge(width = 0.18)

g <- ggplot(plot_df, aes(x = own_we, y = predicted, group = origin_we,
                         colour = origin_we, shape = origin_we)) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = dodge,
                 linewidth = 0.4, alpha = 0.55, show.legend = FALSE) +
  geom_line(position = dodge, linewidth = 0.8) +
  geom_point(position = dodge, size = 2.6, stroke = 0.6, fill = "white") +
  facet_wrap(~ children_we, nrow = 1) +
  scale_colour_manual(values = setNames(ramp, lab), name = "Origin SSS (WE)") +
  scale_shape_manual(values = setNames(c(21, 15, 17), lab), name = "Origin SSS (WE)") +
  scale_y_continuous(breaks = seq(1.9, 2.4, 0.05)) +
  labs(x = "Own SSS (WE)", y = "Predicted MJP index (1-5)",
       caption = paste0("Predicted values with 95% CI. Controls at sample means; between-person status at the grand mean.\n",
                        "Own and origin: Low / High = -1 / +1 SD of the within-person deviation. Children: 25th and 75th\n",
                        "percentiles of the within-person deviation. The y axis does not start at 1.")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
        axis.text = element_text(colour = ink2),
        axis.title = element_text(colour = ink),
        strip.text = element_text(colour = ink, face = "bold", hjust = 0),
        legend.position = "top", legend.justification = "left",
        legend.title = element_text(colour = ink), legend.text = element_text(colour = ink2),
        plot.caption = element_text(colour = muted, hjust = 0, size = 8),
        plot.caption.position = "plot",
        panel.spacing.x = unit(1.4, "lines"),
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("output/interactions/fig-int-triple-we-p25p75.png", g, width = 8, height = 4.6, dpi = 300)
ggsave("output/interactions/fig-int-triple-we-p25p75.pdf", g, width = 8, height = 4.6)
