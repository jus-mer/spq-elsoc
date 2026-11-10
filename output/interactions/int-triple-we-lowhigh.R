# Predicted values and interaction plot for the within-person (WE) three-way
# interaction of subjective social status: own x family of origin x children,
# showing only low and high levels of each (no mean).
# Run from the project root: Rscript output/interactions/int-triple-we-lowhigh.R
#
# Y = market justice preferences (MJP) index
# X = own SSS (WE);  lines = origin SSS (WE);  facets = children's SSS (WE)
# Low / High = 25th / 75th percentile of each within-person deviation, i.e.
# points below / above the respondent's own average on the 0-10 scale.

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
# Controls at their sample means, BE status terms at the grand mean (0), WE
# status terms at their 25th and 75th percentiles.
# 95% CI from the fixed-effects covariance.
qs <- sapply(d[, c("own_we", "og_we", "hij_we")], quantile, probs = c(0.25, 0.75))
rownames(qs) <- c("P25", "P75")
print(round(qs, 3))

grid <- expand.grid(own = rownames(qs), og = rownames(qs), hij = rownames(qs),
                    stringsAsFactors = FALSE) %>%
  mutate(own_we = qs[own, "own_we"],
         og_we  = qs[og,  "og_we"],
         hij_we = qs[hij, "hij_we"],
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

# Labels carry the percentile and its value in scale points
mk_lab <- function(v) setNames(sprintf("%s (%s, %+.2f)", c("Low", "High"), rownames(qs), qs[, v]),
                               rownames(qs))
lab_own <- mk_lab("own_we"); lab_og <- mk_lab("og_we"); lab_hij <- mk_lab("hij_we")

pred <- grid %>%
  transmute(own_we_value = round(own_we, 3),
            origin_we_value = round(og_we, 3),
            children_we_value = round(hij_we, 3),
            children_we = factor(lab_hij[hij], levels = lab_hij),
            origin_we   = factor(lab_og[og],   levels = lab_og),
            own_we      = factor(lab_own[own], levels = lab_own),
            predicted = round(pred, 3), se = round(se, 3),
            ci_low = round(lwr, 3), ci_high = round(upr, 3)) %>%
  arrange(children_we, origin_we, own_we)

write.csv(pred, "output/interactions/pred-int-triple-we-lowhigh.csv", row.names = FALSE)
print(as.data.frame(pred))

# Plot -----------------------------------------------------------------------
# Origin status: blue = low, brick = high, plus point shape so the lines are
# not told apart by colour alone.
ramp <- c("#2a78d6", "#b5432e")
ink  <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid_col <- "#e1e0d9"

plot_df <- pred %>%
  mutate(children_we = factor(paste0("Children's SSS (WE): ", children_we),
                              levels = paste0("Children's SSS (WE): ", lab_hij)))
dodge <- position_dodge(width = 0.12)

g <- ggplot(plot_df, aes(x = own_we, y = predicted, group = origin_we,
                         colour = origin_we, shape = origin_we)) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = dodge,
                 linewidth = 0.4, alpha = 0.55, show.legend = FALSE) +
  geom_line(position = dodge, linewidth = 0.8) +
  geom_point(position = dodge, size = 2.6, stroke = 0.6, fill = "white") +
  facet_wrap(~ children_we, nrow = 1) +
  scale_colour_manual(values = setNames(ramp, lab_og), name = "Origin SSS (WE)") +
  scale_shape_manual(values = setNames(c(21, 17), lab_og), name = "Origin SSS (WE)") +
  scale_y_continuous(breaks = seq(1.9, 2.4, 0.05)) +
  labs(x = "Own SSS (WE)", y = "Predicted MJP index (1-5)",
       caption = paste0("Predicted values with 95% CI. Controls at sample means; between-person status at the grand mean.\n",
                        "Low / High = 25th / 75th percentile of the within-person deviation, in points of the 0-10 scale\n",
                        "below / above the respondent's own average. The y axis does not start at 1.")) +
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

ggsave("output/interactions/fig-int-triple-we-lowhigh.png", g, width = 7.5, height = 4.6, dpi = 300)
ggsave("output/interactions/fig-int-triple-we-lowhigh.pdf", g, width = 7.5, height = 4.6)
