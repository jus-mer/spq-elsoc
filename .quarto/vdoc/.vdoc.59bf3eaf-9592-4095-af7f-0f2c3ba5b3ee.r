#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: packages
#| include: true
#| warning: false

if (! require("pacman")) install.packages("pacman")

pacman::p_load(tidyverse, 
               sjmisc, 
               sjPlot, 
               lme4, 
               here, 
               performance,
               influence.ME, 
               srvyr,
               ordinal,
               texreg, 
               ggdist,
               misty,
               kableExtra,
               ggalluvial, 
               shadowtext,
               MetBrewer,
               patchwork,
               sjlabelled,
               corrplot)


options(scipen=999)
#
#
#
#| label: data
#| echo: true

load(file = here("input/data/proc/df_final.RData"))

df_model <- df_study1 %>%
  select(-c(prev_salud, dir_ingreso, inc_origen, inc_destino)) %>%
  mutate(traj_ingreso = forcats::fct_na_value_to_level(traj_ingreso, level = "Missings")) %>%
  na.omit()
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#| label: fig-mjp-likert
#| fig-cap: "Distribution of responses to the three market justice items (analytical sample, all person-wave observations)"
#| fig-width: 9
#| fig-height: 3.5
#| fig-align: center

items_lab <- c(just_educ = "Education", just_salud = "Health care", just_pension = "Pensions")

df_items <- df_model %>%
  ungroup() %>%
  select(all_of(names(items_lab)))

df_likert <- df_items %>%
  pivot_longer(everything(), names_to = "item", values_to = "resp") %>%
  mutate(item = factor(item, levels = rev(names(items_lab)), labels = rev(items_lab)),
         resp = factor(resp, levels = levels(df_model$just_educ))) %>%
  count(item, resp) %>%
  group_by(item) %>%
  mutate(p = n / sum(n)) %>%
  ungroup()

ggplot(df_likert, aes(x = p, y = item, fill = resp)) +
  geom_col(position = position_stack(reverse = TRUE), width = 0.7) +
  geom_text(aes(label = ifelse(p >= 0.04, scales::percent(p, accuracy = 1), ""),
                group = resp, colour = resp %in% levels(resp)[c(1, 5)]),
            position = position_stack(vjust = 0.5, reverse = TRUE),
            size = 3.2, show.legend = FALSE) +
  scale_colour_manual(values = c("TRUE" = "white", "FALSE" = "grey15")) +
  scale_fill_manual(values = met.brewer("Hiroshige", 5), name = NULL) +
  scale_x_continuous(labels = scales::percent_format(), expand = c(0, 0)) +
  labs(x = NULL, y = NULL,
       caption = paste0("Person-wave observations in the analytical sample: N = ", format(nrow(df_items), big.mark = ","), ".")) +
  guides(fill = guide_legend(nrow = 1)) +
  theme_minimal(base_size = 12) +
  theme(legend.position = "bottom", panel.grid = element_blank())
#
#
#
#| label: fig-mjp-corr
#| fig-cap: "Spearman correlations between the three market justice items (analytical sample)"
#| fig-width: 5
#| fig-height: 5
#| fig-align: center

mjp_num <- df_model %>%
  ungroup() %>%
  transmute(Education   = as.numeric(just_educ),
            `Health care` = as.numeric(just_salud),
            Pensions    = as.numeric(just_pension))

cor_mjp <- cor(mjp_num, method = "spearman", use = "pairwise.complete.obs")

corrplot(cor_mjp,
         method = "color",
         type = "lower",
         diag = FALSE,
         addCoef.col = "black",
         number.cex = 1,
         number.digits = 2,
         tl.col = "black",
         tl.srt = 45,
         col = colorRampPalette(c("#B2182B", "white", "#2166AC"))(200),
         col.lim = c(-1, 1),
         is.corr = TRUE,
         mar = c(0, 0, 1, 0))
#
#
#
#| label: mjp-alpha
#| echo: true

mjp_items <- df_model %>%
  ungroup() %>%
  transmute(just_educ    = as.numeric(just_educ),
            just_salud   = as.numeric(just_salud),
            just_pension = as.numeric(just_pension))

alpha_mjp <- performance::cronbachs_alpha(mjp_items)
round(alpha_mjp, 2)
#
#
#
#
#
#| label: fig-mjp-dist
#| fig-cap: "Distribution of the market justice preferences (MJP) index by wave"
#| fig-width: 9
#| fig-height: 5.5
#| fig-align: center

df_mjp <- df_model %>%
  ungroup() %>%
  filter(!is.na(mjp_index)) %>%
  mutate(mjp_bin = cut(mjp_index, breaks = seq(1, 5, by = 0.5),
                       include.lowest = TRUE, right = FALSE))

df_mjp_mean <- df_mjp %>%
  group_by(ola) %>%
  summarise(m = mean(mjp_index), .groups = "drop") %>%
  mutate(lab = paste0("Mean = ", sprintf("%.2f", m)))

ggplot(df_mjp, aes(x = mjp_index)) +
  geom_histogram(aes(y = after_stat(count / tapply(count, PANEL, sum)[PANEL])),
                 breaks = seq(0.75, 5.25, by = 0.5), closed = "right",
                 fill = met.brewer("Hiroshige", 1, direction = -1)[1],
                 colour = "white") +
  geom_vline(data = df_mjp_mean, aes(xintercept = m),
             linetype = "dashed", colour = "grey20") +
  geom_text(data = df_mjp_mean, aes(x = 5, y = Inf, label = lab),
            hjust = 1, vjust = 1.5, size = 3.2, inherit.aes = FALSE) +
  scale_x_continuous(breaks = 1:5) +
  scale_y_continuous(labels = scales::percent_format(accuracy = 1),
                     expand = expansion(mult = c(0, 0.1))) +
  facet_wrap(~ ola, nrow = 2) +
  labs(x = "MJP index (1 = strongly disagree, 5 = strongly agree)", y = NULL,
       caption = "Bins of width 0.5 centred on 1, 1.5, ..., 5. Dashed line: wave mean.") +
  theme_minimal(base_size = 12) +
  theme(strip.text = element_text(face = "bold"),
        panel.grid.minor = element_blank())
#
#
#
#
#
#
#
#
#
#
#
#
#| label: tbl-mobility-vars
#| tbl-cap: "Subjective social mobility variables constructed from the subjective social status (SSS) items"
#| tbl-cap-location: top
#| echo: false

mov_vars_tab <- tibble::tribble(
  ~Variable, ~Short_label, ~Operacionalization, ~Concept,
  "ess_og", "Origin SSS",
  "Subjective social status of the family of origin (0-10)",
  "Perceived position of the family of origin on the social ladder",
  "ess", "Own SSS",
  "Subjective social status of the respondent (0-10)",
  "Perceived current position of the respondent on the social ladder",
  "ess_hijos", "Children's SSS",
  "Expected subjective social status of the respondent's children (0-10)",
  "Position where respondents expect their children to be located",
  "mov_exp_og", "Children - origin",
  "Children's SSS minus family-of-origin SSS",
  "Children expected to be better off than the family of origin"
)

mov_vars_tab %>%
  dplyr::select(Concept, Operacionalization, Short_label) %>%
  kbl(format = if (knitr::is_latex_output()) "latex" else "html",
      booktabs = TRUE, caption = NULL, align = "lll") %>%
  kable_styling(full_width = TRUE, font_size = 12,
                latex_options = "scale_down") %>%
  footnote(general = "Positive values indicates upward subjective status or mobility.",
           general_title = "Note:", footnote_as_chunk = TRUE)
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
#
