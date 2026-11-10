# Predicted values and interaction plot for the within-person (WE) two-way
# interaction of subjective social status: own x family of origin.
# Based on M2: all two-way interactions among origin, own and children's
# status (WE and BE), no three-way terms.
# Run from the project root: Rscript output/interactions/int-origin-own-we.R
#
# Y = market justice preferences (MJP) index
# X = own SSS (WE);  lines = origin SSS (WE);  children's SSS (WE) at its mean
# Low / Mean / High = -1 SD / 0 / +1 SD of each within-person deviation.

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

# M2: two-way interactions (WE and BE) with the paper's controls -------------
m <- lmer(mjp_index ~ ola_num + ola_2 + ideo + sex + age +
            dummy_decile10_mean + dummy_educ_mean +
            (og_we + own_we + hij_we)^2 +
            (og_be + own_be + hij_be)^2 +
            (1 + ola_num | idencuesta),
          data = d, REML = FALSE,
          control = lmerControl(optimizer = "bobyqa",
                                optCtrl = list(maxfun = 2e5)))

print(round(coef(summary(m))[c("og_we", "own_we", "og_we:own_we"), ], 4))

# Predicted values -----------------------------------------------------------
# Controls at their sample means, children's WE status and all BE status terms
# at their means (0), own and origin WE status at -1 SD / 0 / +1 SD.
# 95% CI from the fixed-effects covariance.
sd_we <- sapply(d[, c("own_we", "og_we")], sd)
lv    <- c(Low = -1, Mean = 0, High = 1)

grid <- expand.grid(own = names(lv), og = names(lv), stringsAsFactors = FALSE) %>%
  mutate(own_we = lv[own] * sd_we["own_we"],
         og_we  = lv[og]  * sd_we["og_we"],
         hij_we = 0, og_be = 0, own_be = 0, hij_be = 0)

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

# Simple slopes of own SSS (WE) at each level of origin SSS (WE)
slopes <- data.frame(origin_we = names(lv), og_value = lv * sd_we["og_we"]) %>%
  mutate(slope = b["own_we"] + b["og_we:own_we"] * og_value,
         se = sqrt(V["own_we", "own_we"] + og_value^2 * V["og_we:own_we", "og_we:own_we"] +
                     2 * og_value * V["own_we", "og_we:own_we"]),
         p = 2 * pnorm(-abs(slope / se)))
print(slopes, digits = 3)

lab <- c(Low = "Low (-1 SD)", Mean = "Mean", High = "High (+1 SD)")
pred <- grid %>%
  transmute(own_we_value = round(own_we, 3),
            origin_we_value = round(og_we, 3),
            origin_we = factor(lab[og],  levels = lab),
            own_we    = factor(lab[own], levels = lab),
            predicted = round(pred, 3), se = round(se, 3),
            ci_low = round(lwr, 3), ci_high = round(upr, 3)) %>%
  arrange(origin_we, own_we)

write.csv(pred, "output/interactions/pred-int-origin-own-we.csv", row.names = FALSE)
print(as.data.frame(pred))

# Plot -----------------------------------------------------------------------
# Origin status: blue = low, neutral grey = mean, brick = high, plus point
# shape so the lines are not told apart by colour alone.
ramp <- c("#2a78d6", "#7a7974", "#b5432e")
ink  <- "#0b0b0b"; ink2 <- "#52514e"; muted <- "#898781"; grid_col <- "#e1e0d9"
dodge <- position_dodge(width = 0.18)

g <- ggplot(pred, aes(x = own_we, y = predicted, group = origin_we,
                      colour = origin_we, shape = origin_we)) +
  geom_linerange(aes(ymin = ci_low, ymax = ci_high), position = dodge,
                 linewidth = 0.4, alpha = 0.55, show.legend = FALSE) +
  geom_line(position = dodge, linewidth = 0.8) +
  geom_point(position = dodge, size = 2.6, stroke = 0.6, fill = "white") +
  scale_colour_manual(values = setNames(ramp, lab), name = "Origin SSS (WE)") +
  scale_shape_manual(values = setNames(c(21, 15, 17), lab), name = "Origin SSS (WE)") +
  scale_y_continuous(breaks = seq(1.9, 2.4, 0.05)) +
  labs(x = "Own SSS (WE)", y = "Predicted MJP index (1-5)",
       caption = paste0("Predicted values with 95% CI. Controls at sample means; children's status (WE)\n",
                        "and between-person status at their means. Low / High = -1 / +1 SD of the\n",
                        "within-person deviation. The y axis does not start at 1.")) +
  theme_minimal(base_size = 11) +
  theme(panel.grid.minor = element_blank(),
        panel.grid.major.x = element_blank(),
        panel.grid.major.y = element_line(colour = grid_col, linewidth = 0.3),
        axis.text = element_text(colour = ink2),
        axis.title = element_text(colour = ink),
        legend.position = "top", legend.justification = "left",
        legend.title = element_text(colour = ink), legend.text = element_text(colour = ink2),
        plot.caption = element_text(colour = muted, hjust = 0, size = 8),
        plot.caption.position = "plot",
        plot.background = element_rect(fill = "white", colour = NA))

ggsave("output/interactions/fig-int-origin-own-we.png", g, width = 6, height = 4.6, dpi = 300)
ggsave("output/interactions/fig-int-origin-own-we.pdf", g, width = 6, height = 4.6)
