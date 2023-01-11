# Conference Paper
# Pandemic Effects on Nigerian Informal Sector: Evidence from Robust Bayesian Meta-Analysis
# Objs
# 1. What is the impact of the pandemic on the informal sector?
# 2. Are there any improvement post-pandemic?
# 10 Jan 2023

# Loading Packages
library(tidyverse)
library(data.table)

# loading dataset
setwd('./src/')
isector <- read.csv('./data/conf_econdept_unn.csv') %>%
  setnames('s.error', 's_error')
str(isector)
# Filtering and generating standard errors and deviations
isector <- isector %>% select(c(authors, sd, coeficient, s_error, escat, outcome_category, sample,
								sector_groups, company_loc, method))
str(isector)
sum(is.na(isector))

# standard deviation. Given binary predictors, wéused (sqrt(sample)*s_error)*1/2
isector$sd[5:12] <- (sqrt(isector$sample[5:12]) * isector$s_error[5:12]) * (1 / 2)
isector$sd[15:26] <- (sqrt(isector$sample[15:26]) * isector$s_error[15:26]) * (1 / 2)
isector$sd[29:37] <- (sqrt(isector$sample[29:37]) * isector$s_error[29:37]) * (1 / 2)
isector$sd[45:109] <- (sqrt(isector$sample[45:109]) * isector$s_error[45:109]) * (1 / 2)
isector$sd[111] <- (sqrt(isector$sample[111]) * isector$s_error[111]) * (1 / 2)
isector$sd[116:121] <- (sqrt(isector$sample[116:121]) * isector$s_error[116:121]) * (1 / 2)
isector$sd[125:138] <- (sqrt(isector$sample[125:138]) * isector$s_error[125:138]) * (1 / 2)
isector$sd[143:147] <- (sqrt(isector$sample[143:147]) * isector$s_error[143:147]) * (1 / 2)

# Fill Standard errors
isector$s_error[1:4] <- isector$sd[1:4] / sqrt(isector$sample[1:4])
isector$s_error[13:14] <- isector$sd[13:14] / sqrt(isector$sample[13:14])
isector$s_error[27:28] <- isector$sd[27:28] / sqrt(isector$sample[27:28])
isector$s_error[38:86] <- isector$sd[38:86] / sqrt(isector$sample[38:86])
isector$s_error[110] <- isector$sd[110] / sqrt(isector$sample[110])
isector$s_error[112:118] <- isector$sd[112:118] / sqrt(isector$sample[112:118])

sdv <- sqrt(1 / 2 * (1 - 1 / 2)) # assuming 1/2 chance of the having significant beta coefficient. In the absence
# of standard errors
isector$sd[45:86] <- sdv
isector$sd[116:118] <- sdv
isector$sd[125:138] <- sdv
isector$s_error[45:86] <- isector$sd[45:86] / sqrt(isector$sample[45:86])
isector$s_error[116:118] <- isector$sd[116:118] / sqrt(isector$sample[116:118])
isector$s_error[125:138] <- isector$sd[125:138] / sqrt(isector$sample[125:138])
str(isector)
sum(is.na(isector))

# Cohen's D and Standard Error of Cohen's D (SeD)
isector <- isector %>% mutate(
  t_stat = coeficient / s_error,
  cohenD = coeficient / sd,
  SeD = (sqrt(4 / sample + sd^2 / (2 * sample))) #Se of cohen's D
)

# Covid on Nigerian Overall Informal Sector
unique(isector$escat)
unique(isector$outcome_category)
unique(isector$sector_groups)

# Visualisations -- Appendix
str(isector)
pct_format <- scales::percent_format(accuracy = .1) # for bar percentages
bysectorgroup <- isector %>%
  ggplot(., aes(x = outcome_category, fill = sector_groups)) +
  geom_bar(position = 'dodge', col = 'white') +
  geom_text(aes(label = pct_format(after_stat(count) / sum(after_stat(count)))),
			position = position_dodge(width = 0.9), vjust = -0.25,
			stat = 'count', colour = 'black') +
  facet_grid(~escat) +
  labs(title = 'Performance indicators by sector',
	   x = 'COVID',
	   y = 'count of coefficients') +
  theme(legend.position = 'bottom')
bysectorgroup

# bydesign <- covidisector %>%
#   ggplot(., aes(x = outcome_category, fill = design)) +
#   geom_bar(position = 'dodge', col = 'white') +
#   labs(title = 'COVID on Informal Sectors by Design',
# 	   x = 'COVID',
# 	   y = 'count of coefficients') +
#   theme_bw()
# bydesign

bymethod <- isector %>%
  ggplot(., aes(x = outcome_category, fill = method)) +
  geom_bar(position = 'dodge', col = 'white') +
  labs(title = 'Performance indicators by location by Method',
	   x = 'COVID',
	   y = 'count of coefficients') +
  facet_grid(~escat) +
  theme_bw()
bymethod

bycomploc <- isector %>%
  ggplot(., aes(x = outcome_category, fill = company_loc)) +
  geom_bar(position = 'dodge', col = 'white') +
  labs(title = 'Performance indicators by location',
	   y = 'count of coefficients') +
  facet_grid(~escat) +
  theme(legend.position = 'bottom')
bycomploc

byoutcome <- isector %>%
  ggplot(., aes(x = outcome_category, y = coeficient)) +
  geom_col() +
  facet_grid(~escat) +
  labs(x = 'Performance') +
  ggtitle('Performance during and post pandemic') +
  theme_bw()
byoutcome

library(patchwork)
pdf(file = './figures/fig73_pandemic_sector_xtics.pdf', width = 18, height = 15)
bysectorgroup + bycomploc + plot_annotation(tag_levels = 'A') + plot_layout(ncol = 1, nrow = 2)
dev.off()

pdf(file = './figures/fig73a_pandemic_sector_xtics.pdf', width = 15, height = 15)
bymethod + byoutcome + plot_annotation(tag_levels = 'A') + plot_layout(ncol = 1, nrow = 2)
dev.off()


# RoBMA-PSMA of Covid on Informal sector performance
library(RoBMA)
library(psych)
library(xtable)
library(ggpubr)
library(patchwork)

sumstat <- isector %>%
  filter(escat == 'covid') %>%
  select(c(coeficient, t_stat, s_error, cohenD, sample))
sumstat <- describe(sumstat)
sumstat
sumstat <- select(sumstat, c(mean, sd, min, max, n))
sumstat

covidsect <- isector %>%
  filter(escat == 'covid') %>%
  select(c(authors, coeficient, s_error, sample, cohenD))
sectperfom <- RoBMA(d = covidsect$cohenD,
					n = covidsect$sample,
					model_type = 'psma', # pp and 2w
					study_names = covidsect$authors,
					parallel = T,
					transformation = 'fishers_z',
					seed = 1)

summary(sectperfom, output_scale = 'r', conditional = T)
interpret(sectperfom, output_scale = 'r')
# Results -- Positive but weak evidence for the negative effect. Inconclusive

# Models Overview
pmp_sectperfom <- summary(sectperfom, output_scale = 'r', type = 'models')
pmp_sectperfom
names(pmp_sectperfom)
pmp_sectperfom <- pmp_sectperfom$summary
pmp_sectperfom <- select(pmp_sectperfom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_sectperfom, type = "latex", label = 'covid-sectorperformance-models-overview', auto = T),
	  # tabular.environment = 'longtable', floating = F,
	  file = "./tables/tab1_covid_sectorperformance_models_overview.tex")

# Plotting Models Overview
models <- plot_models(sectperfom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(sectperfom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig1_covid_sectorperformance_models_overview.pdf', width = 7, height = 22)
models
dev.off()

# Diagnostice Test -- Appendix
diagnos_sectperfom <- summary(sectperfom, output_scale = 'r', type = 'diagnostics')
diagnos_sectperfom
names(diagnos_sectperfom)
diagnos_sectperfom <- diagnos_sectperfom$diagnostics
# view(diagnos_sectperfom)
print(xtable(diagnos_sectperfom, type = "latex", label = 'covid-on-sectorperform-diagnostics', auto = T),
	  # tabular.environment = 'longtable', floating = F,
	  index = F,
	  file = "./tables/tab2_covid_on_sectorperform.tex")

# Mean
chains_mean <- diagnostics(sectperfom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_mean <- diagnostics(sectperfom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_mean <- diagnostics(sectperfom, parameter = "mu", type = "densities", plot_type = 'ggplot', show_models = 36)

# Heterogeneity
chains_tau <- diagnostics(sectperfom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_tau <- diagnostics(sectperfom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_tau <- diagnostics(sectperfom, parameter = "tau", type = "densities", plot_type = 'ggplot', show_models = 36)

pdf(file = './figures/fig2_covid_sectorperform_diagnostics.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots
forest <- forest(sectperfom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(sectperfom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig3_covid_sectorperform_forestplot.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities
esize <- plot(sectperfom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
			  main = 'PMP for the Effect of COVID on ISectors', xlab = 'Mean')
het <- plot(sectperfom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
			main = 'PMP for Heterogeneity of COVID on ISectors', xlab = 'Heterogeneity')
selm <- plot(sectperfom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
			 plot_type = 'ggplot', main = 'PMP for the Selection Models of COVID on ISectors',
			 xlab = 'Selection models publication bias adjustments')
petpeese <- plot(sectperfom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
				 plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of COVID on ISectors',
				 xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig4_covid_sectorperform_posteriors.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
																											  'bold'))
dev.off()

# https://fbartos.github.io/RoBMA/articles/CustomEnsembles.html
#### COVID on ISectors -- custom priors ####
# Model 1 -- Assuming no Priors: effect, heterogeneity and publication bias
sectperform_custom <- RoBMA(d = covidsect$cohenD,
							n = covidsect$sample,
							study_names = covidsect$authors,
							# model_type = 'psma',
							transformation = 'fishers_z',
							parallel = T,
							priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
							priors_effect_null = prior("spike", parameters = list(location = 0)),
							priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
							priors_bias_null = prior_none(),
							seed = 1)

summary(sectperform_custom, output_scale = 'r')
interpret(sectperform_custom, output_scale = 'r')
summary(sectperform_custom, output_scale = 'r', type = "models")

# Before we add the second model to the ensemble, we need to decide on the prior distribution for the mean parameter.
# If a household has electricity access, then, children's study hours will likely increase.
# A positive effect would mean more hours that children's could study.
# This is because any deviation from randomness could be characterized as an effect.
# Therefore, I decided to use the first moment and truncate at the minimum of the distribution. These are the beliefs
# from the surveyed literature on this relationship in the global south. That is mean = 4.84, standard deviation 11.90
# and the minimum = -22. This sets the prior density around the beliefs from the literature.
# To get a better grasp of the prior distribution, we visualize it using the plot() function (the figure can be
# also created using the ggplot2 package by adding plot_type == "ggplot" argument).
sumstat
pdf(file = './figures/fig5_covid_sectperform_beliefs.pdf', width = 7, height = 5)
plot(prior("normal", parameters = list(mean = 56.21, sd = 336.57), truncation = list(lower = -1038)),
	 plot_type = 'ggplot', main = 'Literature Beliefs for COVID on ISectors')
dev.off()

# We add the second model to the ensemble using the update.RoBMA() function. The function can also be used for many
# other purposes - updating settings, prior model weights, and refitting failed models. Here, we supply the fitted
# ensemble object and add an argument specifying the prior distributions of each components for the additional model.
# Since we want to add Model 2 - we set the prior for the μ parameter to be treated as a prior belonging to the
# alternative hypothesis of the effect size component and the remaining priors treated as belonging to the alternative
# hypotheses. If we wanted, we could also specify prior_weights argument, to change the prior probability of the fitted
# model but we do not utilize this option here and keep the default value, which sets the prior weights for the new
# model to 1. (Note that the arguments for specifying prior distributions in update.RoBMA() function are
# prior_X - in singular, in comparison to RoBMA() function that uses priors_X in plural.)

# Model 2
sectperform_custom <- update(sectperform_custom,
							 prior_effect = prior("normal",
												  parameters = list(mean = 56.21, sd = 336.57),
												  truncation = list(lower = -1038)),
							 prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
							 prior_bias_null = prior_none())

# Before we add the remaining models to the ensemble using the update() function, we need to decide on the remaining
# prior distributions. Specifically, on the prior distribution for the heterogeneity parameter τ, and the publication
# bias adjustment parameters ω (for the selection models’ weightfunctions) and PET and PEESE for the PET and PEESE
# adjustment.
#
# For Model 3, we use the usual inverse-gamma(1, .15) prior distribution based on empirical heterogeneity estimates
# (Erp et al., 2017) for the heterogeneity parameter τ. For Models 4.1-4.4 we use the default settings for the
# publication bias adjustments as outlined the Appendix B of (Bartoš, Maier, et al., 2021).
#
# Now, we just need to add the remaining models to the ensemble using the update() function as already illustrated.

# Model 3
sectperform_custom <- update(sectperform_custom,
							 prior_effect = prior("normal",
												  parameters = list(mean = 56.21, sd = 336.57),
												  truncation = list(lower = -1038)),
							 prior_heterogeneity = prior("invgamma", parameters = list(shape = 1, scale = 0.05)),
							 prior_bias_null = prior_none())

# Model 4.1 -- one-sided selection operating on significant p-values (Hf1,pb1)
sectperform_custom <- update(sectperform_custom,
							 prior_effect_null = prior("spike", parameters = list(location = 0)),
							 prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
							 prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1),
																							  steps = 0.05)))

# Model 4.2 -- one-sided selection operating on significant and marginally significant p-values (Hf1,pb2)
sectperform_custom <- update(sectperform_custom,
							 prior_effect_null = prior("spike", parameters = list(location = 0)),
							 prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
							 prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1, 1),
																							  steps = c(0.05, 0.10))))

# Model 4.3 -- PET correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors (Hf1,pb3)
sectperform_custom <- update(sectperform_custom,
							 prior_effect_null = prior("spike", parameters = list(location = 0)),
							 prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
							 prior_bias = prior_PET("Cauchy", parameters = list(0, 1),
													truncation = list(lower = 0)))

# Model 4.4 -- PEESE correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors squared (Hf1,pb4).
sectperform_custom <- update(sectperform_custom,
							 prior_effect_null = prior("spike", parameters = list(location = 0)),
							 prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
							 prior_bias = prior_PEESE("Cauchy", parameters = list(0, 5),
													  truncation = list(lower = 0)))

summary(sectperform_custom, output_scale = 'r', conditional = T)
interpret(sectperform_custom, output_scale = 'r')
# Results -- strong evidence for the positive. Thus, on average, Grid electrification increases children's study hours
# by (0.62*60) = 37.2mins in developing countries.

# Here! Table for the effect (2w, pp, psma, and custom) -- tab4

# Models Overview Custom Beliefs
pmp_sectperform_custom <- summary(sectperform_custom, output_scale = 'r', type = 'models')
pmp_sectperform_custom
names(pmp_sectperform_custom)
pmp_sectperform_custom <- pmp_sectperform_custom$summary
pmp_sectperform_custom <- select(pmp_sectperform_custom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_sectperform_custom, type = "latex", label = 'covid-sectperform-models-overview'),
	  # tabular.environment = 'longtable', floating = F,
	  file = "./tables/tab3_covid_sectperform_models_overview_custom.tex")

# Plotting Models Overview Custom
models <- plot_models(sectperform_custom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(sectperform_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')

pdf(file = './figures/fig6_covid_sectperform_models_overview_custom.pdf', width = 10, height = 7)
ggarrange(models, models_cond, labels = c('A', 'B'), ncol = 2, font.label = list(face = 'bold'))
dev.off()

# Diagnostics Overview Custom Beliefs -- Appendix
diagnos_sectperform_custom <- summary(sectperform_custom, output_scale = 'r', type = 'diagnostics')
diagnos_sectperform_custom
names(diagnos_sectperform_custom)
diagnos_sectperform_custom <- diagnos_sectperform_custom$diagnostics
# view(diagnos_ea_shours_custom)
print(xtable(diagnos_sectperform_custom, type = "latex", label = 'covid-sectperform-diadnostics', auto = T),
	  # tabular.environment = 'longtable', floating = F,
	  file = "./tables/tab4_diagnos_covid_sectperform_custom.tex")

# Mean
chains_mean <- diagnostics(sectperform_custom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_mean <- diagnostics(sectperform_custom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot',
						  show_models = 3)
densities_mean <- diagnostics(sectperform_custom, parameter = "mu", type = "densities", plot_type = 'ggplot',
							  show_models = 3)

# Heterogeneity
chains_tau <- diagnostics(sectperform_custom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_tau <- diagnostics(sectperform_custom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot',
						 show_models = 3)
densities_tau <- diagnostics(sectperform_custom, parameter = "tau", type = "densities", plot_type = 'ggplot',
							 show_models = 3)

pdf(file = './figures/fig7_diagnos_covid_sectperforn_custom.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots Custom
forest <- forest(sectperform_custom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(sectperform_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig8_covid_sectperform_forestplot_custom.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities Custom
esize <- plot(sectperform_custom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
			  main = 'PMP for the Effect of COVID on ISectors', xlab = 'Mean')
het <- plot(sectperform_custom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
			main = 'PMP for Heterogeneity of COVID on ISectors', xlab = 'Heterogeneity')
selm <- plot(sectperform_custom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
			 plot_type = 'ggplot', main = 'PMP Plot for the Selection Models of COVID on ISectors',
			 xlab = 'Selection models publication bias adjustments')
petpeese <- plot(sectperform_custom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
				 plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of COVID on ISectors',
				 xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig9_covid_sectperform_posteriors_custom.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
																											  'bold'))
dev.off()


# RoBMA for Postcovid on informal Sectors
sumstat <- isector %>%
  filter(escat == 'postcovid') %>%
  select(c(coeficient, t_stat, s_error, cohenD, sample))
sumstat <- describe(sumstat)
sumstat
sumstat <- select(sumstat, c(mean, sd, min, max, n))
sumstat

postcovidsect <- isector %>%
  filter(escat == 'postcovid') %>%
  select(c(authors, coeficient, s_error, sample, cohenD))
postcovid_sectperfom <- RoBMA(d = postcovidsect$cohenD,
					n = postcovidsect$sample,
					model_type = 'psma', # pp and 2w
					study_names = postcovidsect$authors,
					parallel = T,
					transformation = 'fishers_z',
					seed = 1)

summary(postcovid_sectperfom, output_scale = 'r', conditional = T)
interpret(postcovid_sectperfom, output_scale = 'r')
# Results -- Positive but weak evidence for the negative effect. Inconclusive

# Models Overview
pmp_postcovid_sectperfom <- summary(postcovid_sectperfom, output_scale = 'r', type = 'models')
pmp_postcovid_sectperfom
names(pmp_postcovid_sectperfom)
pmp_postcovid_sectperfom <- pmp_postcovid_sectperfom$summary
pmp_postcovid_sectperfom <- select(pmp_postcovid_sectperfom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_postcovid_sectperfom, type = "latex", label = 'postcovid-sectorperformance-models-overview', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab5_postcovid_sectorperformance_models_overview.tex")

# Plotting Models Overview
models <- plot_models(postcovid_sectperfom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(postcovid_sectperfom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig10_postcovid_sectorperformance_models_overview.pdf', width = 7, height = 22)
models
dev.off()

# Diagnostice Test -- Appendix
diagnos_postcovid_sectperfom <- summary(postcovid_sectperfom, output_scale = 'r', type = 'diagnostics')
diagnos_postcovid_sectperfom
names(diagnos_postcovid_sectperfom)
diagnos_postcovid_sectperfom <- diagnos_postcovid_sectperfom$diagnostics
# view(diagnos_postcovid_sectperfom)
print(xtable(diagnos_postcovid_sectperfom, type = "latex", label = 'postcovid-on-sectorperform-diagnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     index = F,
     file = "./tables/tab6_postcovid_on_sectorperform.tex")

# Mean
chains_mean <- diagnostics(postcovid_sectperfom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_mean <- diagnostics(postcovid_sectperfom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_mean <- diagnostics(postcovid_sectperfom, parameter = "mu", type = "densities", plot_type = 'ggplot', show_models = 36)

# Heterogeneity
chains_tau <- diagnostics(postcovid_sectperfom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_tau <- diagnostics(postcovid_sectperfom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_tau <- diagnostics(postcovid_sectperfom, parameter = "tau", type = "densities", plot_type = 'ggplot', show_models = 36)

pdf(file = './figures/fig11_postcovid_sectorperform_diagnostics.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots
forest <- forest(postcovid_sectperfom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(postcovid_sectperfom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig12_postcovid_sectorperform_forestplot.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities
esize <- plot(postcovid_sectperfom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of post-COVID on ISectors', xlab = 'Mean')
het <- plot(postcovid_sectperfom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of post-COVID on ISectors', xlab = 'Heterogeneity')
selm <- plot(postcovid_sectperfom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP for the Selection Models of post-COVID on ISectors',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(postcovid_sectperfom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of post-COVID on ISectors',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig13_postcovid_sectorperform_posteriors.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()


# https://fbartos.github.io/RoBMA/articles/CustomEnsembles.html
#### PostCOVID on ISectors -- custom priors ####
# Model 1 -- Assuming no Priors: effect, heterogeneity and publication bias
postcovid_sectperform_custom <- RoBMA(d = postcovidsect$cohenD,
                     n = postcovidsect$sample,
                     study_names = postcovidsect$authors,
                     # model_type = 'psma',
                     transformation = 'fishers_z',
                     parallel = T,
                     priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
                     priors_effect_null = prior("spike", parameters = list(location = 0)),
                     priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                     priors_bias_null = prior_none(),
                     seed = 1)

summary(postcovid_sectperform_custom, output_scale = 'r')
interpret(postcovid_sectperform_custom)
summary(postcovid_sectperform_custom, output_scale = 'r', type = "models")

# Before we add the second model to the ensemble, we need to decide on the prior distribution for the mean parameter.
# If a household has electricity access, then, children's study hours will likely increase.
# A positive effect would mean more hours that children's could study.
# This is because any deviation from randomness could be characterized as an effect.
# Therefore, I decided to use the first moment and truncate at the minimum of the distribution. These are the beliefs
# from the surveyed literature on this relationship in the global south. That is mean = 4.84, standard deviation 11.90
# and the minimum = -22. This sets the prior density around the beliefs from the literature.
# To get a better grasp of the prior distribution, we visualize it using the plot() function (the figure can be
# also created using the ggplot2 package by adding plot_type == "ggplot" argument).
sumstat
pdf(file = './figures/fig14_postcovid_sectperform_beliefs.pdf', width = 7, height = 5)
plot(prior("normal", parameters = list(mean = 2.47, sd = 2.50), truncation = list(lower = 0.16)),
    plot_type = 'ggplot', main = 'Literature Beliefs for post-COVID on ISectors')
dev.off()

# We add the second model to the ensemble using the update.RoBMA() function. The function can also be used for many
# other purposes - updating settings, prior model weights, and refitting failed models. Here, we supply the fitted
# ensemble object and add an argument specifying the prior distributions of each components for the additional model.
# Since we want to add Model 2 - we set the prior for the μ parameter to be treated as a prior belonging to the
# alternative hypothesis of the effect size component and the remaining priors treated as belonging to the alternative
# hypotheses. If we wanted, we could also specify prior_weights argument, to change the prior probability of the fitted
# model but we do not utilize this option here and keep the default value, which sets the prior weights for the new
# model to 1. (Note that the arguments for specifying prior distributions in update.RoBMA() function are
# prior_X - in singular, in comparison to RoBMA() function that uses priors_X in plural.)

# Model 2
postcovid_sectperform_custom <- update(postcovid_sectperform_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 2.47, sd = 2.50),
                                      truncation = list(lower = 0.16)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias_null = prior_none())

# Before we add the remaining models to the ensemble using the update() function, we need to decide on the remaining
# prior distributions. Specifically, on the prior distribution for the heterogeneity parameter τ, and the publication
# bias adjustment parameters ω (for the selection models’ weightfunctions) and PET and PEESE for the PET and PEESE
# adjustment.
#
# For Model 3, we use the usual inverse-gamma(1, .15) prior distribution based on empirical heterogeneity estimates
# (Erp et al., 2017) for the heterogeneity parameter τ. For Models 4.1-4.4 we use the default settings for the
# publication bias adjustments as outlined the Appendix B of (Bartoš, Maier, et al., 2021).
#
# Now, we just need to add the remaining models to the ensemble using the update() function as already illustrated.

# Model 3
postcovid_sectperform_custom <- update(postcovid_sectperform_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 2.47, sd = 2.50),
                                      truncation = list(lower = 0.16)),
                      prior_heterogeneity = prior("invgamma", parameters = list(shape = 1, scale = 0.05)),
                      prior_bias_null = prior_none())

# Model 4.1 -- one-sided selection operating on significant p-values (Hf1,pb1)
postcovid_sectperform_custom <- update(postcovid_sectperform_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1),
                                                                       steps = 0.05)))

# Model 4.2 -- one-sided selection operating on significant and marginally significant p-values (Hf1,pb2)
postcovid_sectperform_custom <- update(postcovid_sectperform_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1, 1),
                                                                       steps = c(0.05, 0.10))))

# Model 4.3 -- PET correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors (Hf1,pb3)
postcovid_sectperform_custom <- update(postcovid_sectperform_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PET("Cauchy", parameters = list(0, 1),
                                       truncation = list(lower = 0)))

# Model 4.4 -- PEESE correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors squared (Hf1,pb4).
postcovid_sectperform_custom <- update(postcovid_sectperform_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PEESE("Cauchy", parameters = list(0, 5),
                                         truncation = list(lower = 0)))

summary(postcovid_sectperform_custom, output_scale = 'r', conditional = T)
interpret(postcovid_sectperform_custom, output_scale = 'r')
# Results -- strong evidence for the positive. Thus, on average, Grid electrification increases children's study hours
# by (0.62*60) = 37.2mins in developing countries.

# Here! Table for the effect (2w, pp, psma, and custom) -- tab4

# Models Overview Custom Beliefs
pmp_postcovid_sectperform_custom <- summary(postcovid_sectperform_custom, output_scale = 'r', type = 'models')
pmp_postcovid_sectperform_custom
names(pmp_postcovid_sectperform_custom)
pmp_postcovid_sectperform_custom <- pmp_postcovid_sectperform_custom$summary
pmp_postcovid_sectperform_custom <- select(pmp_postcovid_sectperform_custom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_postcovid_sectperform_custom, type = "latex", label = 'postcovid-postcovid_sectperform-models-overview'),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab7_postcovid_sectperform_models_overview_custom.tex")

# Plotting Models Overview Custom
models <- plot_models(postcovid_sectperform_custom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(postcovid_sectperform_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')

pdf(file = './figures/fig15_postcovid_sectperform_models_overview_custom.pdf', width = 10, height = 7)
models
dev.off()

# Diagnostics Overview Custom Beliefs -- Appendix
diagnos_postcovid_sectperform_custom <- summary(postcovid_sectperform_custom, output_scale = 'r', type = 'diagnostics')
diagnos_postcovid_sectperform_custom
names(diagnos_postcovid_sectperform_custom)
diagnos_postcovid_sectperform_custom <- diagnos_postcovid_sectperform_custom$diagnostics
# view(diagnos_ea_shours_custom)
print(xtable(diagnos_postcovid_sectperform_custom, type = "latex", label = 'post-covid-postcovid_sectperform-diadnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab8_diagnos_covid_postcovid_sectperform_custom.tex")

# Mean
chains_mean <- diagnostics(postcovid_sectperform_custom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_mean <- diagnostics(postcovid_sectperform_custom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot',
                    show_models = 3)
densities_mean <- diagnostics(postcovid_sectperform_custom, parameter = "mu", type = "densities", plot_type = 'ggplot',
                       show_models = 3)

# Heterogeneity
chains_tau <- diagnostics(postcovid_sectperform_custom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_tau <- diagnostics(postcovid_sectperform_custom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot',
                   show_models = 3)
densities_tau <- diagnostics(postcovid_sectperform_custom, parameter = "tau", type = "densities", plot_type = 'ggplot',
                      show_models = 3)

pdf(file = './figures/fig16_diagnos_postcovid_sectperforn_custom.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots Custom
forest <- forest(postcovid_sectperform_custom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(postcovid_sectperform_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig17_postcovid_sectperform_forestplot_custom.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities Custom
esize <- plot(postcovid_sectperform_custom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of post-COVID on ISectors', xlab = 'Mean')
het <- plot(postcovid_sectperform_custom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of post-COVID on ISectors', xlab = 'Heterogeneity')
selm <- plot(postcovid_sectperform_custom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP Plot for the Selection Models of post-COVID on ISectors',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(postcovid_sectperform_custom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of post-COVID on ISectors',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig18_postcovid_sectperform_posteriors_custom.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()


# RoBMA for covid on company profits
sumstat <- isector %>%
  filter(escat == 'covid' & outcome_category == 'profit') %>%
  select(c(coeficient, t_stat, s_error, cohenD, sample))
sumstat <- describe(sumstat)
sumstat
sumstat <- select(sumstat, c(mean, sd, min, max, n))
sumstat

covidsect_profit <- isector %>%
  filter(escat == 'covid' & outcome_category == 'profit') %>%
  select(c(authors, coeficient, s_error, sample, cohenD))
covid_profit <- RoBMA(d = covidsect_profit$cohenD,
					n = covidsect_profit$sample,
					model_type = 'psma', # pp and 2w
					study_names = covidsect_profit$authors,
					parallel = T,
					transformation = 'fishers_z',
					seed = 1)

summary(covid_profit, output_scale = 'r', conditional = T)
interpret(covid_profit, output_scale = 'r')
# Results -- Positive but weak evidence for the negative effect. Inconclusive

# Models Overview
pmp_covid_profit <- summary(covid_profit, output_scale = 'r', type = 'models')
pmp_covid_profit
names(pmp_covid_profit)
pmp_covid_profit <- pmp_covid_profit$summary
pmp_covid_profit <- select(pmp_covid_profit, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_covid_profit, type = "latex", label = 'covid-profit-models-overview', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab9_covid_profit_models_overview.tex")

# Plotting Models Overview
models <- plot_models(covid_profit, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(covid_profit, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig19_covid_profit_models_overview.pdf', width = 7, height = 22)
models
dev.off()

# Diagnostice Test -- Appendix
diagnos_covid_profit <- summary(covid_profit, output_scale = 'r', type = 'diagnostics')
diagnos_covid_profit
names(diagnos_covid_profit)
diagnos_covid_profit <- diagnos_covid_profit$diagnostics
# view(diagnos_covid_profit)
print(xtable(diagnos_covid_profit, type = "latex", label = 'covid-on-profit-diagnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     index = F,
     file = "./tables/tab10_covid_on_profit.tex")

# Mean
chains_mean <- diagnostics(covid_profit, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_mean <- diagnostics(covid_profit, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_mean <- diagnostics(covid_profit, parameter = "mu", type = "densities", plot_type = 'ggplot', show_models = 36)

# Heterogeneity
chains_tau <- diagnostics(covid_profit, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_tau <- diagnostics(covid_profit, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_tau <- diagnostics(covid_profit, parameter = "tau", type = "densities", plot_type = 'ggplot', show_models = 36)

pdf(file = './figures/fig20_covid_profit_diagnostics.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots
forest <- forest(covid_profit, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(covid_profit, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig21_covid_profit_forestplot.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities
esize <- plot(covid_profit, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of COVID on Profit', xlab = 'Mean')
het <- plot(covid_profit, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of COVID on Profit', xlab = 'Heterogeneity')
selm <- plot(covid_profit, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP for the Selection Models of COVID on Profit',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(covid_profit, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of COVID on Profit',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig22_covid_profit_posteriors.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()

# https://fbartos.github.io/RoBMA/articles/CustomEnsembles.html
#### COVID on Profit -- custom priors ####
# Model 1 -- Assuming no Priors: effect, heterogeneity and publication bias
covid_profit_custom <- RoBMA(d = covidsect_profit$cohenD,
                     n = covidsect_profit$sample,
                     study_names = covidsect_profit$authors,
                     # model_type = 'psma',
                     transformation = 'fishers_z',
                     parallel = T,
                     priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
                     priors_effect_null = prior("spike", parameters = list(location = 0)),
                     priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                     priors_bias_null = prior_none(),
                     seed = 1)

summary(covid_profit_custom, output_scale = 'r')
interpret(covid_profit_custom)
summary(covid_profit_custom, output_scale = 'r', type = "models")

# Before we add the second model to the ensemble, we need to decide on the prior distribution for the mean parameter.
# If a household has electricity access, then, children's study hours will likely increase.
# A positive effect would mean more hours that children's could study.
# This is because any deviation from randomness could be characterized as an effect.
# Therefore, I decided to use the first moment and truncate at the minimum of the distribution. These are the beliefs
# from the surveyed literature on this relationship in the global south. That is mean = 4.84, standard deviation 11.90
# and the minimum = -22. This sets the prior density around the beliefs from the literature.
# To get a better grasp of the prior distribution, we visualize it using the plot() function (the figure can be
# also created using the ggplot2 package by adding plot_type == "ggplot" argument).
sumstat
pdf(file = './figures/fig23_covid_profit_beliefs.pdf', width = 7, height = 5)
plot(prior("normal", parameters = list(mean = 122.82, sd = 506.02), truncation = list(lower = -1038)),
    plot_type = 'ggplot', main = 'Literature Beliefs for COVID on Profit')
dev.off()

# We add the second model to the ensemble using the update.RoBMA() function. The function can also be used for many
# other purposes - updating settings, prior model weights, and refitting failed models. Here, we supply the fitted
# ensemble object and add an argument specifying the prior distributions of each components for the additional model.
# Since we want to add Model 2 - we set the prior for the μ parameter to be treated as a prior belonging to the
# alternative hypothesis of the effect size component and the remaining priors treated as belonging to the alternative
# hypotheses. If we wanted, we could also specify prior_weights argument, to change the prior probability of the fitted
# model but we do not utilize this option here and keep the default value, which sets the prior weights for the new
# model to 1. (Note that the arguments for specifying prior distributions in update.RoBMA() function are
# prior_X - in singular, in comparison to RoBMA() function that uses priors_X in plural.)

# Model 2
covid_profit_custom <- update(covid_profit_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias_null = prior_none())

# Before we add the remaining models to the ensemble using the update() function, we need to decide on the remaining
# prior distributions. Specifically, on the prior distribution for the heterogeneity parameter τ, and the publication
# bias adjustment parameters ω (for the selection models’ weightfunctions) and PET and PEESE for the PET and PEESE
# adjustment.
#
# For Model 3, we use the usual inverse-gamma(1, .15) prior distribution based on empirical heterogeneity estimates
# (Erp et al., 2017) for the heterogeneity parameter τ. For Models 4.1-4.4 we use the default settings for the
# publication bias adjustments as outlined the Appendix B of (Bartoš, Maier, et al., 2021).
#
# Now, we just need to add the remaining models to the ensemble using the update() function as already illustrated.

# Model 3
covid_profit_custom <- update(covid_profit_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity = prior("invgamma", parameters = list(shape = 1, scale = 0.05)),
                      prior_bias_null = prior_none())

# Model 4.1 -- one-sided selection operating on significant p-values (Hf1,pb1)
covid_profit_custom <- update(covid_profit_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1),
                                                                       steps = 0.05)))

# Model 4.2 -- one-sided selection operating on significant and marginally significant p-values (Hf1,pb2)
covid_profit_custom <- update(covid_profit_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1, 1),
                                                                       steps = c(0.05, 0.10))))

# Model 4.3 -- PET correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors (Hf1,pb3)
covid_profit_custom <- update(covid_profit_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PET("Cauchy", parameters = list(0, 1),
                                       truncation = list(lower = 0)))

# Model 4.4 -- PEESE correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors squared (Hf1,pb4).
covid_profit_custom <- update(covid_profit_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PEESE("Cauchy", parameters = list(0, 5),
                                         truncation = list(lower = 0)))

summary(covid_profit_custom, output_scale = 'r', conditional = T)
interpret(covid_profit_custom, output_scale = 'r')
# Results -- strong evidence for the positive. Thus, on average, Grid electrification increases children's study hours
# by (0.62*60) = 37.2mins in developing countries.

# Here! Table for the effect (2w, pp, psma, and custom) -- tab4

# Models Overview Custom Beliefs
pmp_covid_profit_custom <- summary(covid_profit_custom, output_scale = 'r', type = 'models')
pmp_covid_profit_custom
names(pmp_covid_profit_custom)
pmp_covid_profit_custom <- pmp_covid_profit_custom$summary
pmp_covid_profit_custom <- select(pmp_covid_profit_custom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_covid_profit_custom, type = "latex", label = 'covid-profit-models-overview'),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab11_covid_profit_models_overview_custom.tex")

# Plotting Models Overview Custom
models <- plot_models(covid_profit_custom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(covid_profit_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')

pdf(file = './figures/fig24_covid_profit_models_overview_custom.pdf', width = 10, height = 7)
models
dev.off()

# Diagnostics Overview Custom Beliefs -- Appendix
diagnos_covid_profit_custom <- summary(covid_profit_custom, output_scale = 'r', type = 'diagnostics')
diagnos_covid_profit_custom
names(diagnos_covid_profit_custom)
diagnos_covid_profit_custom <- diagnos_covid_profit_custom$diagnostics
# view(diagnos_ea_shours_custom)
print(xtable(diagnos_covid_profit_custom, type = "latex", label = 'covid-profit-diadnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab12_diagnos_covid_profit_custom.tex")

# Mean
chains_mean <- diagnostics(covid_profit_custom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_mean <- diagnostics(covid_profit_custom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot',
                    show_models = 3)
densities_mean <- diagnostics(covid_profit_custom, parameter = "mu", type = "densities", plot_type = 'ggplot',
                       show_models = 3)

# Heterogeneity
chains_tau <- diagnostics(covid_profit_custom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_tau <- diagnostics(covid_profit_custom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot',
                   show_models = 3)
densities_tau <- diagnostics(covid_profit_custom, parameter = "tau", type = "densities", plot_type = 'ggplot',
                      show_models = 3)

pdf(file = './figures/fig25_diagnos_covid_profit_custom.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots Custom
forest <- forest(covid_profit_custom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(covid_profit_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig26_covid_profit_forestplot_custom.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities Custom
esize <- plot(covid_profit_custom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of COVID on Profit', xlab = 'Mean')
het <- plot(covid_profit_custom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of COVID on Profit', xlab = 'Heterogeneity')
selm <- plot(covid_profit_custom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP Plot for the Selection Models of COVID on Profit',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(covid_profit_custom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of COVID on Profit',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig27_covid_profit_posteriors_custom.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()


# RoBMA for covid on company sales
sumstat <- isector %>%
  filter(escat == 'covid' & outcome_category == 'sales') %>%
  select(c(coeficient, t_stat, s_error, cohenD, sample))
sumstat <- describe(sumstat)
sumstat
sumstat <- select(sumstat, c(mean, sd, min, max, n))
sumstat

covidsect_sales <- isector %>%
  filter(escat == 'covid' & outcome_category == 'sales') %>%
  select(c(authors, coeficient, s_error, sample, cohenD))
covid_sales <- RoBMA(d = covidsect_sales$cohenD,
               n = covidsect_sales$sample,
               model_type = 'psma', # pp and 2w
               study_names = covidsect_sales$authors,
               parallel = T,
               transformation = 'fishers_z',
               seed = 1)

summary(covid_sales, output_scale = 'r', conditional = T)
interpret(covid_sales, output_scale = 'r')
# Results -- Positive but weak evidence for the negative effect. Inconclusive

# Models Overview
pmp_covid_sales <- summary(covid_sales, output_scale = 'r', type = 'models')
pmp_covid_sales
names(pmp_covid_sales)
pmp_covid_sales <- pmp_covid_sales$summary
pmp_covid_sales <- select(pmp_covid_sales, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_covid_sales, type = "latex", label = 'covid-sales-models-overview', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab13_covid_sales_models_overview.tex")

# Plotting Models Overview
models <- plot_models(covid_sales, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(covid_sales, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig28_covid_sales_models_overview.pdf', width = 7, height = 22)
models
dev.off()

# Diagnostice Test -- Appendix
diagnos_covid_sales <- summary(covid_sales, output_scale = 'r', type = 'diagnostics')
diagnos_covid_sales
names(diagnos_covid_sales)
diagnos_covid_sales <- diagnos_covid_sales$diagnostics
# view(diagnos_covid_sales)
print(xtable(diagnos_covid_sales, type = "latex", label = 'covid-on-sales-diagnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     index = F,
     file = "./tables/tab14_covid_on_sales.tex")

# Mean
chains_mean <- diagnostics(covid_sales, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_mean <- diagnostics(covid_sales, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_mean <- diagnostics(covid_sales, parameter = "mu", type = "densities", plot_type = 'ggplot', show_models = 36)

# Heterogeneity
chains_tau <- diagnostics(covid_sales, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_tau <- diagnostics(covid_sales, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_tau <- diagnostics(covid_sales, parameter = "tau", type = "densities", plot_type = 'ggplot', show_models = 36)

pdf(file = './figures/fig29_covid_sales_diagnostics.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots
forest <- forest(covid_sales, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(covid_sales, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig30_covid_sales_forestplot.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities
esize <- plot(covid_sales, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of COVID on Sales', xlab = 'Mean')
het <- plot(covid_sales, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of COVID on Sales', xlab = 'Heterogeneity')
selm <- plot(covid_sales, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP for the Selection Models of COVID on Sales',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(covid_sales, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of COVID on Sales',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig31_covid_sales_posteriors.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()

# https://fbartos.github.io/RoBMA/articles/CustomEnsembles.html
#### COVID on Sales -- custom priors ####
# Model 1 -- Assuming no Priors: effect, heterogeneity and publication bias
covid_sales_custom <- RoBMA(d = covidsect_sales$cohenD,
                     n = covidsect_sales$sample,
                     study_names = covidsect_sales$authors,
                     # model_type = 'psma',
                     transformation = 'fishers_z',
                     parallel = T,
                     priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
                     priors_effect_null = prior("spike", parameters = list(location = 0)),
                     priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                     priors_bias_null = prior_none(),
                     seed = 1)

summary(covid_sales_custom, output_scale = 'r')
interpret(covid_sales_custom)
summary(covid_sales_custom, output_scale = 'r', type = "models")

# Before we add the second model to the ensemble, we need to decide on the prior distribution for the mean parameter.
# If a household has electricity access, then, children's study hours will likely increase.
# A positive effect would mean more hours that children's could study.
# This is because any deviation from randomness could be characterized as an effect.
# Therefore, I decided to use the first moment and truncate at the minimum of the distribution. These are the beliefs
# from the surveyed literature on this relationship in the global south. That is mean = 4.84, standard deviation 11.90
# and the minimum = -22. This sets the prior density around the beliefs from the literature.
# To get a better grasp of the prior distribution, we visualize it using the plot() function (the figure can be
# also created using the ggplot2 package by adding plot_type == "ggplot" argument).
sumstat
pdf(file = './figures/fig32_covid_sales_beliefs.pdf', width = 7, height = 5)
plot(prior("normal", parameters = list(mean = 3.89, sd = 16.88), truncation = list(lower = -36.59)),
    plot_type = 'ggplot', main = 'Literature Beliefs for COVID on sales')
dev.off()

# We add the second model to the ensemble using the update.RoBMA() function. The function can also be used for many
# other purposes - updating settings, prior model weights, and refitting failed models. Here, we supply the fitted
# ensemble object and add an argument specifying the prior distributions of each components for the additional model.
# Since we want to add Model 2 - we set the prior for the μ parameter to be treated as a prior belonging to the
# alternative hypothesis of the effect size component and the remaining priors treated as belonging to the alternative
# hypotheses. If we wanted, we could also specify prior_weights argument, to change the prior probability of the fitted
# model but we do not utilize this option here and keep the default value, which sets the prior weights for the new
# model to 1. (Note that the arguments for specifying prior distributions in update.RoBMA() function are
# prior_X - in singular, in comparison to RoBMA() function that uses priors_X in plural.)

# Model 2
covid_sales_custom <- update(covid_sales_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias_null = prior_none())

# Before we add the remaining models to the ensemble using the update() function, we need to decide on the remaining
# prior distributions. Specifically, on the prior distribution for the heterogeneity parameter τ, and the publication
# bias adjustment parameters ω (for the selection models’ weightfunctions) and PET and PEESE for the PET and PEESE
# adjustment.
#
# For Model 3, we use the usual inverse-gamma(1, .15) prior distribution based on empirical heterogeneity estimates
# (Erp et al., 2017) for the heterogeneity parameter τ. For Models 4.1-4.4 we use the default settings for the
# publication bias adjustments as outlined the Appendix B of (Bartoš, Maier, et al., 2021).
#
# Now, we just need to add the remaining models to the ensemble using the update() function as already illustrated.

# Model 3
covid_sales_custom <- update(covid_sales_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity = prior("invgamma", parameters = list(shape = 1, scale = 0.05)),
                      prior_bias_null = prior_none())

# Model 4.1 -- one-sided selection operating on significant p-values (Hf1,pb1)
covid_sales_custom <- update(covid_sales_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1),
                                                                       steps = 0.05)))

# Model 4.2 -- one-sided selection operating on significant and marginally significant p-values (Hf1,pb2)
covid_sales_custom <- update(covid_sales_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1, 1),
                                                                       steps = c(0.05, 0.10))))

# Model 4.3 -- PET correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors (Hf1,pb3)
covid_sales_custom <- update(covid_sales_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PET("Cauchy", parameters = list(0, 1),
                                       truncation = list(lower = 0)))

# Model 4.4 -- PEESE correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors squared (Hf1,pb4).
covid_sales_custom <- update(covid_sales_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PEESE("Cauchy", parameters = list(0, 5),
                                         truncation = list(lower = 0)))

summary(covid_sales_custom, output_scale = 'r', conditional = T)
interpret(covid_sales_custom, output_scale = 'r')
# Results -- strong evidence for the positive. Thus, on average, Grid electrification increases children's study hours
# by (0.62*60) = 37.2mins in developing countries.

# Here! Table for the effect (2w, pp, psma, and custom) -- tab4

# Models Overview Custom Beliefs
pmp_covid_sales_custom <- summary(covid_sales_custom, output_scale = 'r', type = 'models')
pmp_covid_sales_custom
names(pmp_covid_sales_custom)
pmp_covid_sales_custom <- pmp_covid_sales_custom$summary
pmp_covid_sales_custom <- select(pmp_covid_sales_custom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_covid_sales_custom, type = "latex", label = 'covid-sales-models-overview'),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab15_covid_sales_models_overview_custom.tex")

# Plotting Models Overview Custom
models <- plot_models(covid_sales_custom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(covid_sales_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')

pdf(file = './figures/fig33_covid_sales_models_overview_custom.pdf', width = 10, height = 7)
models
dev.off()

# Diagnostics Overview Custom Beliefs -- Appendix
diagnos_covid_sales_custom <- summary(covid_sales_custom, output_scale = 'r', type = 'diagnostics')
diagnos_covid_sales_custom
names(diagnos_covid_sales_custom)
diagnos_covid_sales_custom <- diagnos_covid_sales_custom$diagnostics
# view(diagnos_ea_shours_custom)
print(xtable(diagnos_covid_sales_custom, type = "latex", label = 'covid-sales-diadnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab16_diagnos_covid_sales_custom.tex")

# Mean
chains_mean <- diagnostics(covid_sales_custom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_mean <- diagnostics(covid_sales_custom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot',
                    show_models = 3)
densities_mean <- diagnostics(covid_sales_custom, parameter = "mu", type = "densities", plot_type = 'ggplot',
                       show_models = 3)

# Heterogeneity
chains_tau <- diagnostics(covid_sales_custom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_tau <- diagnostics(covid_sales_custom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot',
                   show_models = 3)
densities_tau <- diagnostics(covid_sales_custom, parameter = "tau", type = "densities", plot_type = 'ggplot',
                      show_models = 3)

pdf(file = './figures/fig34_diagnos_covid_sales_custom.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots Custom
forest <- forest(covid_sales_custom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(covid_sales_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig35_covid_sales_forestplot_custom.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities Custom
esize <- plot(covid_sales_custom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of COVID on Sales', xlab = 'Mean')
het <- plot(covid_sales_custom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of COVID on Sales', xlab = 'Heterogeneity')
selm <- plot(covid_sales_custom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP Plot for the Selection Models of COVID on Sales',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(covid_sales_custom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of COVID on Sales',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig36_covid_sales_posteriors_custom.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()


# RoBMA for covid on company survival
sumstat <- isector %>%
  filter(escat == 'covid' & outcome_category == 'survival') %>%
  select(c(coeficient, t_stat, s_error, cohenD, sample))
sumstat <- describe(sumstat)
sumstat
sumstat <- select(sumstat, c(mean, sd, min, max, n))
sumstat

covidsect_survival <- isector %>%
  filter(escat == 'covid' & outcome_category == 'survival') %>%
  select(c(authors, coeficient, s_error, sample, cohenD))
covid_survival <- RoBMA(d = covidsect_survival$cohenD,
               n = covidsect_survival$sample,
               model_type = 'psma', # pp and 2w
               study_names = covidsect_survival$authors,
               parallel = T,
               transformation = 'fishers_z',
               seed = 1)
summary(covid_survival, output_scale = 'r', conditional = T)
interpret(covid_survival, output_scale = 'r')
# Results -- Positive but weak evidence for the negative effect. Inconclusive

# Models Overview
pmp_covid_survival <- summary(covid_survival, output_scale = 'r', type = 'models')
pmp_covid_survival
names(pmp_covid_survival)
pmp_covid_survival <- pmp_covid_survival$summary
pmp_covid_survival <- select(pmp_covid_survival, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_covid_survival, type = "latex", label = 'covid-survival-models-overview', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab17_covid_survival_models_overview.tex")

# Plotting Models Overview
models <- plot_models(covid_survival, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(covid_survival, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig37_covid_survival_models_overview.pdf', width = 7, height = 22)
models
dev.off()

# Diagnostice Test -- Appendix
diagnos_covid_survival <- summary(covid_survival, output_scale = 'r', type = 'diagnostics')
diagnos_covid_survival
names(diagnos_covid_survival)
diagnos_covid_survival <- diagnos_covid_survival$diagnostics
# view(diagnos_covid_survival)
print(xtable(diagnos_covid_survival, type = "latex", label = 'covid-on-survival-diagnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     index = F,
     file = "./tables/tab18_covid_on_survival.tex")

# Mean
chains_mean <- diagnostics(covid_survival, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_mean <- diagnostics(covid_survival, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_mean <- diagnostics(covid_survival, parameter = "mu", type = "densities", plot_type = 'ggplot', show_models = 36)

# Heterogeneity
chains_tau <- diagnostics(covid_survival, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_tau <- diagnostics(covid_survival, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_tau <- diagnostics(covid_survival, parameter = "tau", type = "densities", plot_type = 'ggplot', show_models = 36)

pdf(file = './figures/fig38_covid_survival_diagnostics.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots
forest <- forest(covid_survival, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(covid_survival, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig39_covid_survival_forestplot.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities
esize <- plot(covid_survival, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of COVID on Survival', xlab = 'Mean')
het <- plot(covid_survival, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of COVID on Survival', xlab = 'Heterogeneity')
selm <- plot(covid_survival, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP for the Selection Models of COVID on Survival',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(covid_survival, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of COVID on Survival',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig40_covid_survival_posteriors.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()

# https://fbartos.github.io/RoBMA/articles/CustomEnsembles.html
#### COVID on survival -- custom priors ####
# Model 1 -- Assuming no Priors: effect, heterogeneity and publication bias
covid_survival_custom <- RoBMA(d = covidsect_survival$cohenD,
                     n = covidsect_survival$sample,
                     study_names = covidsect_survival$authors,
                     # model_type = 'psma',
                     transformation = 'fishers_z',
                     parallel = T,
                     priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
                     priors_effect_null = prior("spike", parameters = list(location = 0)),
                     priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                     priors_bias_null = prior_none(),
                     seed = 1)

summary(covid_survival_custom, output_scale = 'r')
interpret(covid_survival_custom)
summary(covid_survival_custom, output_scale = 'r', type = "models")

# Before we add the second model to the ensemble, we need to decide on the prior distribution for the mean parameter.
# If a household has electricity access, then, children's study hours will likely increase.
# A positive effect would mean more hours that children's could study.
# This is because any deviation from randomness could be characterized as an effect.
# Therefore, I decided to use the first moment and truncate at the minimum of the distribution. These are the beliefs
# from the surveyed literature on this relationship in the global south. That is mean = 4.84, standard deviation 11.90
# and the minimum = -22. This sets the prior density around the beliefs from the literature.
# To get a better grasp of the prior distribution, we visualize it using the plot() function (the figure can be
# also created using the ggplot2 package by adding plot_type == "ggplot" argument).
sumstat
pdf(file = './figures/fig41_covid_survival_beliefs.pdf', width = 7, height = 5)
plot(prior("normal", parameters = list(mean = 3.89, sd = 16.88), truncation = list(lower = -36.59)),
    plot_type = 'ggplot', main = 'Literature Beliefs for COVID on Survival')
dev.off()

# We add the second model to the ensemble using the update.RoBMA() function. The function can also be used for many
# other purposes - updating settings, prior model weights, and refitting failed models. Here, we supply the fitted
# ensemble object and add an argument specifying the prior distributions of each components for the additional model.
# Since we want to add Model 2 - we set the prior for the μ parameter to be treated as a prior belonging to the
# alternative hypothesis of the effect size component and the remaining priors treated as belonging to the alternative
# hypotheses. If we wanted, we could also specify prior_weights argument, to change the prior probability of the fitted
# model but we do not utilize this option here and keep the default value, which sets the prior weights for the new
# model to 1. (Note that the arguments for specifying prior distributions in update.RoBMA() function are
# prior_X - in singular, in comparison to RoBMA() function that uses priors_X in plural.)

# Model 2
covid_survival_custom <- update(covid_survival_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias_null = prior_none())

# Before we add the remaining models to the ensemble using the update() function, we need to decide on the remaining
# prior distributions. Specifically, on the prior distribution for the heterogeneity parameter τ, and the publication
# bias adjustment parameters ω (for the selection models’ weightfunctions) and PET and PEESE for the PET and PEESE
# adjustment.
#
# For Model 3, we use the usual inverse-gamma(1, .15) prior distribution based on empirical heterogeneity estimates
# (Erp et al., 2017) for the heterogeneity parameter τ. For Models 4.1-4.4 we use the default settings for the
# publication bias adjustments as outlined the Appendix B of (Bartoš, Maier, et al., 2021).
#
# Now, we just need to add the remaining models to the ensemble using the update() function as already illustrated.

# Model 3
covid_survival_custom <- update(covid_survival_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity = prior("invgamma", parameters = list(shape = 1, scale = 0.05)),
                      prior_bias_null = prior_none())

# Model 4.1 -- one-sided selection operating on significant p-values (Hf1,pb1)
covid_survival_custom <- update(covid_survival_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1),
                                                                       steps = 0.05)))

# Model 4.2 -- one-sided selection operating on significant and marginally significant p-values (Hf1,pb2)
covid_survival_custom <- update(covid_survival_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1, 1),
                                                                       steps = c(0.05, 0.10))))

# Model 4.3 -- PET correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors (Hf1,pb3)
covid_survival_custom <- update(covid_survival_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PET("Cauchy", parameters = list(0, 1),
                                       truncation = list(lower = 0)))

# Model 4.4 -- PEESE correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors squared (Hf1,pb4).
covid_survival_custom <- update(covid_survival_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PEESE("Cauchy", parameters = list(0, 5),
                                         truncation = list(lower = 0)))

summary(covid_survival_custom, output_scale = 'r', conditional = T)
interpret(covid_survival_custom, output_scale = 'r')
# Results -- strong evidence for the positive. Thus, on average, Grid electrification increases children's study hours
# by (0.62*60) = 37.2mins in developing countries.

# Here! Table for the effect (2w, pp, psma, and custom) -- tab4

# Models Overview Custom Beliefs
pmp_covid_survival_custom <- summary(covid_survival_custom, output_scale = 'r', type = 'models')
pmp_covid_survival_custom
names(pmp_covid_survival_custom)
pmp_covid_survival_custom <- pmp_covid_survival_custom$summary
pmp_covid_survival_custom <- select(pmp_covid_survival_custom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_covid_survival_custom, type = "latex", label = 'covid-survival-models-overview'),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab19_covid_survival_models_overview_custom.tex")

# Plotting Models Overview Custom
models <- plot_models(covid_survival_custom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(covid_survival_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')

pdf(file = './figures/fig42_covid_survival_models_overview_custom.pdf', width = 10, height = 7)
models
dev.off()

# Diagnostics Overview Custom Beliefs -- Appendix
diagnos_covid_survival_custom <- summary(covid_survival_custom, output_scale = 'r', type = 'diagnostics')
diagnos_covid_survival_custom
names(diagnos_covid_survival_custom)
diagnos_covid_survival_custom <- diagnos_covid_survival_custom$diagnostics
# view(diagnos_ea_shours_custom)
print(xtable(diagnos_covid_survival_custom, type = "latex", label = 'covid-survival-diadnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab20_diagnos_covid_survial_custom.tex")

# Mean
chains_mean <- diagnostics(covid_survival_custom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_mean <- diagnostics(covid_survival_custom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot',
                    show_models = 3)
densities_mean <- diagnostics(covid_survival_custom, parameter = "mu", type = "densities", plot_type = 'ggplot',
                       show_models = 3)

# Heterogeneity
chains_tau <- diagnostics(covid_survival_custom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_tau <- diagnostics(covid_survival_custom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot',
                   show_models = 3)
densities_tau <- diagnostics(covid_survival_custom, parameter = "tau", type = "densities", plot_type = 'ggplot',
                      show_models = 3)

pdf(file = './figures/fig43_diagnos_covid_survival_custom.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots Custom
forest <- forest(covid_survival_custom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(covid_survival_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig44_covid_survival_forestplot_custom.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities Custom
esize <- plot(covid_survival_custom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of COVID on Survival', xlab = 'Mean')
het <- plot(covid_survival_custom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of COVID on Survival', xlab = 'Heterogeneity')
selm <- plot(covid_survival_custom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP Plot for the Selection Models of COVID on Survival',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(covid_survival_custom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of COVID on Survival',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig45_covid_survival_posteriors_custom.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()


# RoBMA for postcovid on company profit
sumstat <- isector %>%
  filter(escat == 'postcovid' & outcome_category == 'profit') %>%
  select(c(coeficient, t_stat, s_error, cohenD, sample))
sumstat <- describe(sumstat)
sumstat
sumstat <- select(sumstat, c(mean, sd, min, max, n))
sumstat

postcovidsect_profit <- isector %>%
  filter(escat == 'postcovid' & outcome_category == 'profit') %>%
  select(c(authors, coeficient, s_error, sample, cohenD))
postcovid_profit <- RoBMA(d = postcovidsect_profit$cohenD,
               n = postcovidsect_profit$sample,
               model_type = 'psma', # pp and 2w
               study_names = postcovidsect_profit$authors,
               parallel = T,
               transformation = 'fishers_z',
               seed = 1)
summary(postcovid_profit, output_scale = 'r', conditional = T)
interpret(postcovid_profit, output_scale = 'r')
# Results -- Positive but weak evidence for the negative effect. Inconclusive

# Models Overview
pmp_postcovid_profit <- summary(postcovid_profit, output_scale = 'r', type = 'models')
pmp_postcovid_profit
names(pmp_postcovid_profit)
pmp_postcovid_profit <- pmp_postcovid_profit$summary
pmp_postcovid_profit <- select(pmp_postcovid_profit, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_postcovid_profit, type = "latex", label = 'postcovid-profit-models-overview', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab21_postcovid_profit_models_overview.tex")

# Plotting Models Overview
models <- plot_models(postcovid_profit, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(postcovid_profit, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig46_postcovid_profit_models_overview.pdf', width = 7, height = 22)
models
dev.off()

# Diagnostice Test -- Appendix
diagnos_postcovid_profit <- summary(postcovid_profit, output_scale = 'r', type = 'diagnostics')
diagnos_postcovid_profit
names(diagnos_postcovid_profit)
diagnos_postcovid_profit <- diagnos_postcovid_profit$diagnostics
# view(diagnos_postcovid_profit)
print(xtable(diagnos_postcovid_profit, type = "latex", label = 'postcovid-on-profit-diagnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     index = F,
     file = "./tables/tab22_postcovid_on_profit.tex")

# Mean
chains_mean <- diagnostics(postcovid_profit, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_mean <- diagnostics(postcovid_profit, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_mean <- diagnostics(postcovid_profit, parameter = "mu", type = "densities", plot_type = 'ggplot', show_models = 36)

# Heterogeneity
chains_tau <- diagnostics(postcovid_profit, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_tau <- diagnostics(postcovid_profit, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_tau <- diagnostics(postcovid_profit, parameter = "tau", type = "densities", plot_type = 'ggplot', show_models = 36)

pdf(file = './figures/fig47_postcovid_profit_diagnostics.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots
forest <- forest(postcovid_profit, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(postcovid_profit, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig48_postcovid_profit_forestplot.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities
esize <- plot(postcovid_profit, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of post-COVID on Profit', xlab = 'Mean')
het <- plot(postcovid_profit, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of post-COVID on Profit', xlab = 'Heterogeneity')
selm <- plot(postcovid_profit, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP for the Selection Models of post-COVID on Profit',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(postcovid_profit, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of post-COVID on Profit',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig49_postcovid_profit_posteriors.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()

# https://fbartos.github.io/RoBMA/articles/CustomEnsembles.html
#### PostCOVID on profit -- custom priors ####
# Model 1 -- Assuming no Priors: effect, heterogeneity and publication bias
postcovid_profit_custom <- RoBMA(d = postcovidsect_profit$cohenD,
                     n = postcovidsect_profit$sample,
                     study_names = postcovidsect_profit$authors,
                     # model_type = 'psma',
                     transformation = 'fishers_z',
                     parallel = T,
                     priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
                     priors_effect_null = prior("spike", parameters = list(location = 0)),
                     priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                     priors_bias_null = prior_none(),
                     seed = 1)

summary(postcovid_profit_custom, output_scale = 'r')
interpret(postcovid_profit_custom)
summary(postcovid_profit_custom, output_scale = 'r', type = "models")

# Before we add the second model to the ensemble, we need to decide on the prior distribution for the mean parameter.
# If a household has electricity access, then, children's study hours will likely increase.
# A positive effect would mean more hours that children's could study.
# This is because any deviation from randomness could be characterized as an effect.
# Therefore, I decided to use the first moment and truncate at the minimum of the distribution. These are the beliefs
# from the surveyed literature on this relationship in the global south. That is mean = 4.84, standard deviation 11.90
# and the minimum = -22. This sets the prior density around the beliefs from the literature.
# To get a better grasp of the prior distribution, we visualize it using the plot() function (the figure can be
# also created using the ggplot2 package by adding plot_type == "ggplot" argument).
sumstat
pdf(file = './figures/fig50_postcovid_profit_beliefs.pdf', width = 7, height = 5)
plot(prior("normal", parameters = list(mean = 3.89, sd = 16.88), truncation = list(lower = -36.59)),
    plot_type = 'ggplot', main = 'Literature Beliefs for post-COVID on Profit')
dev.off()

# We add the second model to the ensemble using the update.RoBMA() function. The function can also be used for many
# other purposes - updating settings, prior model weights, and refitting failed models. Here, we supply the fitted
# ensemble object and add an argument specifying the prior distributions of each components for the additional model.
# Since we want to add Model 2 - we set the prior for the μ parameter to be treated as a prior belonging to the
# alternative hypothesis of the effect size component and the remaining priors treated as belonging to the alternative
# hypotheses. If we wanted, we could also specify prior_weights argument, to change the prior probability of the fitted
# model but we do not utilize this option here and keep the default value, which sets the prior weights for the new
# model to 1. (Note that the arguments for specifying prior distributions in update.RoBMA() function are
# prior_X - in singular, in comparison to RoBMA() function that uses priors_X in plural.)

# Model 2
postcovid_profit_custom <- update(postcovid_profit_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias_null = prior_none())

# Before we add the remaining models to the ensemble using the update() function, we need to decide on the remaining
# prior distributions. Specifically, on the prior distribution for the heterogeneity parameter τ, and the publication
# bias adjustment parameters ω (for the selection models’ weightfunctions) and PET and PEESE for the PET and PEESE
# adjustment.
#
# For Model 3, we use the usual inverse-gamma(1, .15) prior distribution based on empirical heterogeneity estimates
# (Erp et al., 2017) for the heterogeneity parameter τ. For Models 4.1-4.4 we use the default settings for the
# publication bias adjustments as outlined the Appendix B of (Bartoš, Maier, et al., 2021).
#
# Now, we just need to add the remaining models to the ensemble using the update() function as already illustrated.

# Model 3
postcovid_profit_custom <- update(postcovid_profit_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity = prior("invgamma", parameters = list(shape = 1, scale = 0.05)),
                      prior_bias_null = prior_none())

# Model 4.1 -- one-sided selection operating on significant p-values (Hf1,pb1)
postcovid_profit_custom <- update(postcovid_profit_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1),
                                                                       steps = 0.05)))

# Model 4.2 -- one-sided selection operating on significant and marginally significant p-values (Hf1,pb2)
postcovid_profit_custom <- update(postcovid_profit_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1, 1),
                                                                       steps = c(0.05, 0.10))))

# Model 4.3 -- PET correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors (Hf1,pb3)
postcovid_profit_custom <- update(postcovid_profit_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PET("Cauchy", parameters = list(0, 1),
                                       truncation = list(lower = 0)))

# Model 4.4 -- PEESE correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors squared (Hf1,pb4).
postcovid_profit_custom <- update(postcovid_profit_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PEESE("Cauchy", parameters = list(0, 5),
                                         truncation = list(lower = 0)))

summary(postcovid_profit_custom, output_scale = 'r', conditional = T)
interpret(postcovid_profit_custom, output_scale = 'r')

# Here! Table for the effect (2w, pp, psma, and custom) -- tab4

# Models Overview Custom Beliefs
pmp_postcovid_profit_custom <- summary(postcovid_profit_custom, output_scale = 'r', type = 'models')
pmp_postcovid_profit_custom
names(pmp_postcovid_profit_custom)
pmp_postcovid_profit_custom <- pmp_postcovid_profit_custom$summary
pmp_postcovid_profit_custom <- select(pmp_postcovid_profit_custom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_postcovid_profit_custom, type = "latex", label = 'postcovid-profit-models-overview'),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab23_postcovid_profit_models_overview_custom.tex")

# Plotting Models Overview Custom
models <- plot_models(postcovid_profit_custom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(postcovid_profit_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')

pdf(file = './figures/fig51_postcovid_profit_models_overview_custom.pdf', width = 10, height = 7)
models
dev.off()

# Diagnostics Overview Custom Beliefs -- Appendix
diagnos_postcovid_profit_custom <- summary(postcovid_profit_custom, output_scale = 'r', type = 'diagnostics')
diagnos_postcovid_profit_custom
names(diagnos_postcovid_profit_custom)
diagnos_postcovid_profit_custom <- diagnos_postcovid_profit_custom$diagnostics
# view(diagnos_ea_shours_custom)
print(xtable(diagnos_postcovid_profit_custom, type = "latex", label = 'postcovid-profit-diadnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab24_diagnos_postcovid_profit_custom.tex")

# Mean
chains_mean <- diagnostics(postcovid_profit_custom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_mean <- diagnostics(postcovid_profit_custom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot',
                    show_models = 3)
densities_mean <- diagnostics(postcovid_profit_custom, parameter = "mu", type = "densities", plot_type = 'ggplot',
                       show_models = 3)

# Heterogeneity
chains_tau <- diagnostics(postcovid_profit_custom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_tau <- diagnostics(postcovid_profit_custom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot',
                   show_models = 3)
densities_tau <- diagnostics(postcovid_profit_custom, parameter = "tau", type = "densities", plot_type = 'ggplot',
                      show_models = 3)

pdf(file = './figures/fig52_diagnos_postcovid_profit_custom.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots Custom
forest <- forest(postcovid_profit_custom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(postcovid_profit_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig53_postcovid_profit_forestplot_custom.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities Custom
esize <- plot(postcovid_profit_custom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of post-COVID on Profit', xlab = 'Mean')
het <- plot(postcovid_profit_custom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of post-COVID on Profit', xlab = 'Heterogeneity')
selm <- plot(postcovid_profit_custom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP Plot for the Selection Models of post-COVID on Profit',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(postcovid_profit_custom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of post-COVID on Profit',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig54_postcovid_profit_posteriors_custom.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face = 'bold'))
dev.off()


# RoBMA for postcovid on company sales
sumstat <- isector %>% filter(escat == 'postcovid' & outcome_category == 'sales') %>%
  select(c(coeficient, t_stat, s_error, cohenD, sample))
sumstat <- describe(sumstat)
sumstat
sumstat <- select(sumstat, c(mean, sd, min, max, n))
sumstat

postcovidsect_sales <- isector %>%
  filter(escat == 'postcovid' & outcome_category == 'sales') %>%
  select(c(authors, coeficient, s_error, sample, cohenD))
postcovid_sales <- RoBMA(d = postcovidsect_sales$cohenD,
               n = postcovidsect_sales$sample,
               model_type = 'psma', # pp and 2w
               study_names = postcovidsect_sales$authors,
               parallel = T,
               transformation = 'fishers_z',
               seed = 1)
summary(postcovid_sales, output_scale = 'r', conditional = T)
interpret(postcovid_sales, output_scale = 'r')
# Results -- Positive but weak evidence for the negative effect. Inconclusive

# Models Overview
pmp_postcovid_sales <- summary(postcovid_sales, output_scale = 'r', type = 'models')
pmp_postcovid_sales
names(pmp_postcovid_sales)
pmp_postcovid_sales <- pmp_postcovid_sales$summary
pmp_postcovid_sales <- select(pmp_postcovid_sales, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_postcovid_sales, type = "latex", label = 'postcovid-sales-models-overview', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab25_postcovid_sales_models_overview.tex")

# Plotting Models Overview
models <- plot_models(postcovid_sales, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(postcovid_sales, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig55_postcovid_sales_models_overview.pdf', width = 7, height = 22)
models
dev.off()

# Diagnostice Test -- Appendix
diagnos_postcovid_sales <- summary(postcovid_sales, output_scale = 'r', type = 'diagnostics')
diagnos_postcovid_sales
names(diagnos_postcovid_sales)
diagnos_postcovid_sales <- diagnos_postcovid_sales$diagnostics
# view(diagnos_postcovid_sales)
print(xtable(diagnos_postcovid_sales, type = "latex", label = 'postcovid-on-sales-diagnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     index = F,
     file = "./tables/tab26_postcovid_on_sales_diagnostics.tex")

# Mean
chains_mean <- diagnostics(postcovid_sales, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_mean <- diagnostics(postcovid_sales, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_mean <- diagnostics(postcovid_sales, parameter = "mu", type = "densities", plot_type = 'ggplot', show_models = 36)

# Heterogeneity
chains_tau <- diagnostics(postcovid_sales, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_tau <- diagnostics(postcovid_sales, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_tau <- diagnostics(postcovid_sales, parameter = "tau", type = "densities", plot_type = 'ggplot', show_models = 36)

pdf(file = './figures/fig56_postcovid_sales_diagnostics.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots
forest <- forest(postcovid_sales, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(postcovid_sales, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig57_postcovid_sales_forestplot.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities
esize <- plot(postcovid_sales, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of post-COVID on Sales', xlab = 'Mean')
het <- plot(postcovid_sales, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of post-COVID on Sales', xlab = 'Heterogeneity')
selm <- plot(postcovid_sales, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP for the Selection Models of post-COVID on Sales',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(postcovid_sales, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of post-COVID on Sales',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig58_postcovid_sales_posteriors.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()

# https://fbartos.github.io/RoBMA/articles/CustomEnsembles.html
#### PostCOVID on sales -- custom priors ####
# Model 1 -- Assuming no Priors: effect, heterogeneity and publication bias
postcovid_sales_custom <- RoBMA(d = postcovidsect_sales$cohenD,
                     n = postcovidsect_sales$sample,
                     study_names = postcovidsect_sales$authors,
                     # model_type = 'psma',
                     transformation = 'fishers_z',
                     parallel = T,
                     priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
                     priors_effect_null = prior("spike", parameters = list(location = 0)),
                     priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                     priors_bias_null = prior_none(),
                     seed = 1)

summary(postcovid_sales_custom, output_scale = 'r')
interpret(postcovid_sales_custom)
summary(postcovid_sales_custom, output_scale = 'r', type = "models")

# Before we add the second model to the ensemble, we need to decide on the prior distribution for the mean parameter.
# If a household has electricity access, then, children's study hours will likely increase.
# A positive effect would mean more hours that children's could study.
# This is because any deviation from randomness could be characterized as an effect.
# Therefore, I decided to use the first moment and truncate at the minimum of the distribution. These are the beliefs
# from the surveyed literature on this relationship in the global south. That is mean = 4.84, standard deviation 11.90
# and the minimum = -22. This sets the prior density around the beliefs from the literature.
# To get a better grasp of the prior distribution, we visualize it using the plot() function (the figure can be
# also created using the ggplot2 package by adding plot_type == "ggplot" argument).
sumstat
pdf(file = './figures/fig59_postcovid_sales_beliefs.pdf', width = 7, height = 5)
plot(prior("normal", parameters = list(mean = 3.89, sd = 16.88), truncation = list(lower = -36.59)),
    plot_type = 'ggplot', main = 'Literature Beliefs for post-COVID on Sales')
dev.off()

# We add the second model to the ensemble using the update.RoBMA() function. The function can also be used for many
# other purposes - updating settings, prior model weights, and refitting failed models. Here, we supply the fitted
# ensemble object and add an argument specifying the prior distributions of each components for the additional model.
# Since we want to add Model 2 - we set the prior for the μ parameter to be treated as a prior belonging to the
# alternative hypothesis of the effect size component and the remaining priors treated as belonging to the alternative
# hypotheses. If we wanted, we could also specify prior_weights argument, to change the prior probability of the fitted
# model but we do not utilize this option here and keep the default value, which sets the prior weights for the new
# model to 1. (Note that the arguments for specifying prior distributions in update.RoBMA() function are
# prior_X - in singular, in comparison to RoBMA() function that uses priors_X in plural.)

# Model 2
postcovid_sales_custom <- update(postcovid_sales_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias_null = prior_none())

# Before we add the remaining models to the ensemble using the update() function, we need to decide on the remaining
# prior distributions. Specifically, on the prior distribution for the heterogeneity parameter τ, and the publication
# bias adjustment parameters ω (for the selection models’ weightfunctions) and PET and PEESE for the PET and PEESE
# adjustment.
#
# For Model 3, we use the usual inverse-gamma(1, .15) prior distribution based on empirical heterogeneity estimates
# (Erp et al., 2017) for the heterogeneity parameter τ. For Models 4.1-4.4 we use the default settings for the
# publication bias adjustments as outlined the Appendix B of (Bartoš, Maier, et al., 2021).
#
# Now, we just need to add the remaining models to the ensemble using the update() function as already illustrated.

# Model 3
postcovid_sales_custom <- update(postcovid_sales_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity = prior("invgamma", parameters = list(shape = 1, scale = 0.05)),
                      prior_bias_null = prior_none())

# Model 4.1 -- one-sided selection operating on significant p-values (Hf1,pb1)
postcovid_sales_custom <- update(postcovid_sales_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1),
                                                                       steps = 0.05)))

# Model 4.2 -- one-sided selection operating on significant and marginally significant p-values (Hf1,pb2)
postcovid_sales_custom <- update(postcovid_sales_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1, 1),
                                                                       steps = c(0.05, 0.10))))

# Model 4.3 -- PET correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors (Hf1,pb3)
postcovid_sales_custom <- update(postcovid_sales_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PET("Cauchy", parameters = list(0, 1),
                                       truncation = list(lower = 0)))

# Model 4.4 -- PEESE correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors squared (Hf1,pb4).
postcovid_sales_custom <- update(postcovid_sales_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PEESE("Cauchy", parameters = list(0, 5),
                                         truncation = list(lower = 0)))

summary(postcovid_sales_custom, output_scale = 'r', conditional = T)
interpret(postcovid_sales_custom, output_scale = 'r')
# Results -- strong evidence for the positive. Thus, on average, Grid electrification increases children's study hours
# by (0.62*60) = 37.2mins in developing countries.

# Here! Table for the effect (2w, pp, psma, and custom) -- tab4

# Models Overview Custom Beliefs
pmp_postcovid_sales_custom <- summary(postcovid_sales_custom, output_scale = 'r', type = 'models')
pmp_postcovid_sales_custom
names(pmp_postcovid_sales_custom)
pmp_postcovid_sales_custom <- pmp_postcovid_sales_custom$summary
pmp_postcovid_sales_custom <- select(pmp_postcovid_sales_custom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_postcovid_sales_custom, type = "latex", label = 'postcovid-sales-models-overview'),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab27_postcovid_sales_models_overview_custom.tex")

# Plotting Models Overview Custom
models <- plot_models(postcovid_sales_custom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(postcovid_sales_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')

pdf(file = './figures/fig60_postcovid_sales_models_overview_custom.pdf', width = 10, height = 7)
models
dev.off()

# Diagnostics Overview Custom Beliefs -- Appendix
diagnos_postcovid_sales_custom <- summary(postcovid_sales_custom, output_scale = 'r', type = 'diagnostics')
diagnos_postcovid_sales_custom
names(diagnos_postcovid_sales_custom)
diagnos_postcovid_sales_custom <- diagnos_postcovid_sales_custom$diagnostics
# view(diagnos_ea_shours_custom)
print(xtable(diagnos_postcovid_sales_custom, type = "latex", label = 'postcovid-sales-diadnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab28_diagnos_postcovid_sales_custom.tex")

# Mean
chains_mean <- diagnostics(postcovid_sales_custom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_mean <- diagnostics(postcovid_sales_custom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot',
                    show_models = 3)
densities_mean <- diagnostics(postcovid_sales_custom, parameter = "mu", type = "densities", plot_type = 'ggplot',
                       show_models = 3)

# Heterogeneity
chains_tau <- diagnostics(postcovid_sales_custom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_tau <- diagnostics(postcovid_sales_custom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot',
                   show_models = 3)
densities_tau <- diagnostics(postcovid_sales_custom, parameter = "tau", type = "densities", plot_type = 'ggplot',
                      show_models = 3)

pdf(file = './figures/fig61_diagnos_postcovid_sales_custom.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots Custom
forest <- forest(postcovid_sales_custom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(postcovid_sales_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig62_postcovid_sales_forestplot_custom.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities Custom
esize <- plot(postcovid_sales_custom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of post-COVID on Sales', xlab = 'Mean')
het <- plot(postcovid_sales_custom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of post-COVID on Sales', xlab = 'Heterogeneity')
selm <- plot(postcovid_sales_custom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP Plot for the Selection Models of post-COVID on Sales',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(postcovid_sales_custom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of post-COVID on Sales',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig63_postcovid_sales_posteriors_custom.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face = 'bold'))
dev.off()


# RoBMA for postcovid on company survival
sumstat <- isector %>% filter(escat == 'postcovid' & outcome_category == 'survival') %>%
  select(c(coeficient, t_stat, s_error, cohenD, sample))
sumstat <- describe(sumstat)
sumstat
sumstat <- select(sumstat, c(mean, sd, min, max, n))
sumstat

postcovidsect_survival <- isector %>%
  filter(escat == 'postcovid' & outcome_category == 'survival') %>%
  select(c(authors, coeficient, s_error, sample, cohenD))
postcovid_survival <- RoBMA(d = postcovidsect_survival$cohenD,
               n = postcovidsect_survival$sample,
               model_type = 'psma', # pp and 2w
               study_names = postcovidsect_survival$authors,
               parallel = T,
               transformation = 'fishers_z',
               seed = 1)
summary(postcovid_survival, output_scale = 'r', conditional = T)
interpret(postcovid_survival, output_scale = 'r')
# Results -- Positive but weak evidence for the negative effect. Inconclusive

# Models Overview
pmp_postcovid_survival <- summary(postcovid_survival, output_scale = 'r', type = 'models')
pmp_postcovid_survival
names(pmp_postcovid_survival)
pmp_postcovid_survival <- pmp_postcovid_survival$summary
pmp_postcovid_survival <- select(pmp_postcovid_survival, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_postcovid_survival, type = "latex", label = 'postcovid-survival-models-overview', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab29_postcovid_surival_models_overview.tex")

# Plotting Models Overview
models <- plot_models(postcovid_survival, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(postcovid_survival, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig64_postcovid_survival_models_overview.pdf', width = 7, height = 22)
models
dev.off()

# Diagnostice Test -- Appendix
diagnos_postcovid_survival <- summary(postcovid_survival, output_scale = 'r', type = 'diagnostics')
diagnos_postcovid_survival
names(diagnos_postcovid_survival)
diagnos_postcovid_survival <- diagnos_postcovid_survival$diagnostics
# view(diagnos_postcovid_survival)
print(xtable(diagnos_postcovid_survival, type = "latex", label = 'postcovid-on-survival-diagnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     index = F,
     file = "./tables/tab30_postcovid_on_survival.tex")

# Mean
chains_mean <- diagnostics(postcovid_survival, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_mean <- diagnostics(postcovid_survival, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_mean <- diagnostics(postcovid_survival, parameter = "mu", type = "densities", plot_type = 'ggplot', show_models = 36)

# Heterogeneity
chains_tau <- diagnostics(postcovid_survival, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 36)
acorr_tau <- diagnostics(postcovid_survival, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot', show_models = 36)
densities_tau <- diagnostics(postcovid_survival, parameter = "tau", type = "densities", plot_type = 'ggplot', show_models = 36)

pdf(file = './figures/fig65_postcovid_survival.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots
forest <- forest(postcovid_survival, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(postcovid_survival, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig66_postcovid_survival_forestplot.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities
esize <- plot(postcovid_survival, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of post-COVID on Survival', xlab = 'Mean')
het <- plot(postcovid_survival, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of post-COVID on Survival', xlab = 'Heterogeneity')
selm <- plot(postcovid_survival, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP for the Selection Models of post-COVID on Survival',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(postcovid_survival, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of post-COVID on Survival',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig67_postcovid_survival_posteriors.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face =
                                                                                   'bold'))
dev.off()

# https://fbartos.github.io/RoBMA/articles/CustomEnsembles.html
#### PostCOVID on survival -- custom priors ####
# Model 1 -- Assuming no Priors: effect, heterogeneity and publication bias
postcovid_survival_custom <- RoBMA(d = postcovidsect_survival$cohenD,
                     n = postcovidsect_survival$sample,
                     study_names = postcovidsect_survival$authors,
                     # model_type = 'psma',
                     transformation = 'fishers_z',
                     parallel = T,
                     priors_effect = NULL, priors_heterogeneity = NULL, priors_bias = NULL,
                     priors_effect_null = prior("spike", parameters = list(location = 0)),
                     priors_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                     priors_bias_null = prior_none(),
                     seed = 1)

summary(postcovid_survival_custom, output_scale = 'r')
interpret(postcovid_survival_custom)
summary(postcovid_survival_custom, output_scale = 'r', type = "models")

# Before we add the second model to the ensemble, we need to decide on the prior distribution for the mean parameter.
# If a household has electricity access, then, children's study hours will likely increase.
# A positive effect would mean more hours that children's could study.
# This is because any deviation from randomness could be characterized as an effect.
# Therefore, I decided to use the first moment and truncate at the minimum of the distribution. These are the beliefs
# from the surveyed literature on this relationship in the global south. That is mean = 4.84, standard deviation 11.90
# and the minimum = -22. This sets the prior density around the beliefs from the literature.
# To get a better grasp of the prior distribution, we visualize it using the plot() function (the figure can be
# also created using the ggplot2 package by adding plot_type == "ggplot" argument).
sumstat
pdf(file = './figures/fig68_postcovid_survival_beliefs.pdf', width = 7, height = 5)
plot(prior("normal", parameters = list(mean = 3.89, sd = 16.88), truncation = list(lower = -36.59)),
    plot_type = 'ggplot', main = 'Literature Beliefs for post-COVID on Survival')
dev.off()

# We add the second model to the ensemble using the update.RoBMA() function. The function can also be used for many
# other purposes - updating settings, prior model weights, and refitting failed models. Here, we supply the fitted
# ensemble object and add an argument specifying the prior distributions of each components for the additional model.
# Since we want to add Model 2 - we set the prior for the μ parameter to be treated as a prior belonging to the
# alternative hypothesis of the effect size component and the remaining priors treated as belonging to the alternative
# hypotheses. If we wanted, we could also specify prior_weights argument, to change the prior probability of the fitted
# model but we do not utilize this option here and keep the default value, which sets the prior weights for the new
# model to 1. (Note that the arguments for specifying prior distributions in update.RoBMA() function are
# prior_X - in singular, in comparison to RoBMA() function that uses priors_X in plural.)

# Model 2
postcovid_survival_custom <- update(postcovid_survival_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias_null = prior_none())

# Before we add the remaining models to the ensemble using the update() function, we need to decide on the remaining
# prior distributions. Specifically, on the prior distribution for the heterogeneity parameter τ, and the publication
# bias adjustment parameters ω (for the selection models’ weightfunctions) and PET and PEESE for the PET and PEESE
# adjustment.
#
# For Model 3, we use the usual inverse-gamma(1, .15) prior distribution based on empirical heterogeneity estimates
# (Erp et al., 2017) for the heterogeneity parameter τ. For Models 4.1-4.4 we use the default settings for the
# publication bias adjustments as outlined the Appendix B of (Bartoš, Maier, et al., 2021).
#
# Now, we just need to add the remaining models to the ensemble using the update() function as already illustrated.

# Model 3
postcovid_survival_custom <- update(postcovid_survival_custom,
                      prior_effect = prior("normal",
                                      parameters = list(mean = 122.82, sd = 506.02),
                                      truncation = list(lower = -1038)),
                      prior_heterogeneity = prior("invgamma", parameters = list(shape = 1, scale = 0.05)),
                      prior_bias_null = prior_none())

# Model 4.1 -- one-sided selection operating on significant p-values (Hf1,pb1)
postcovid_survival_custom <- update(postcovid_survival_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1),
                                                                       steps = 0.05)))

# Model 4.2 -- one-sided selection operating on significant and marginally significant p-values (Hf1,pb2)
postcovid_survival_custom <- update(postcovid_survival_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_weightfunction("one.sided", parameters = list(alpha = c(1, 1, 1),
                                                                       steps = c(0.05, 0.10))))

# Model 4.3 -- PET correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors (Hf1,pb3)
postcovid_survival_custom <- update(postcovid_survival_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PET("Cauchy", parameters = list(0, 1),
                                       truncation = list(lower = 0)))

# Model 4.4 -- PEESE correction for publication bias which adjusts for the relationship between effect sizes
# and standard errors squared (Hf1,pb4).
postcovid_survival_custom <- update(postcovid_survival_custom,
                      prior_effect_null = prior("spike", parameters = list(location = 0)),
                      prior_heterogeneity_null = prior("spike", parameters = list(location = 0)),
                      prior_bias = prior_PEESE("Cauchy", parameters = list(0, 5),
                                         truncation = list(lower = 0)))

summary(postcovid_survival_custom, output_scale = 'r', conditional = T)
interpret(postcovid_survival_custom, output_scale = 'r')
# Results -- strong evidence for the positive. Thus, on average, Grid electrification increases children's study hours
# by (0.62*60) = 37.2mins in developing countries.

# Here! Table for the effect (2w, pp, psma, and custom) -- tab4

# Models Overview Custom Beliefs
pmp_postcovid_survival_custom <- summary(postcovid_survival_custom, output_scale = 'r', type = 'models')
pmp_postcovid_survival_custom
names(pmp_postcovid_survival_custom)
pmp_postcovid_survival_custom <- pmp_postcovid_survival_custom$summary
pmp_postcovid_survival_custom <- select(pmp_postcovid_survival_custom, c(Model:prior_prob, post_prob:inclusion_BF))
print(xtable(pmp_postcovid_survival_custom, type = "latex", label = 'postcovid-survival-models-overview'),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab31_postcovid_survival_models_overview_custom.tex")

# Plotting Models Overview Custom
models <- plot_models(postcovid_survival_custom, output_scale = 'r', plot_type = 'ggplot')
models_cond <- plot_models(postcovid_survival_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')

pdf(file = './figures/fig69_postcovid_survival_models_overview_custom.pdf', width = 10, height = 7)
models
dev.off()

# Diagnostics Overview Custom Beliefs -- Appendix
diagnos_postcovid_survival_custom <- summary(postcovid_survival_custom, output_scale = 'r', type = 'diagnostics')
diagnos_postcovid_survival_custom
names(diagnos_postcovid_survival_custom)
diagnos_postcovid_survival_custom <- diagnos_postcovid_survival_custom$diagnostics
# view(diagnos_ea_shours_custom)
print(xtable(diagnos_postcovid_survival_custom, type = "latex", label = 'postcovid-survival-diadnostics', auto = T),
     # tabular.environment = 'longtable', floating = F,
     file = "./tables/tab32_diagnos_postcovid_survival_custom.tex")

# Mean
chains_mean <- diagnostics(postcovid_survival_custom, parameter = "mu", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_mean <- diagnostics(postcovid_survival_custom, parameter = "mu", type = "autocrrelation", plot_type = 'ggplot',
                    show_models = 3)
densities_mean <- diagnostics(postcovid_survival_custom, parameter = "mu", type = "densities", plot_type = 'ggplot',
                       show_models = 3)

# Heterogeneity
chains_tau <- diagnostics(postcovid_survival_custom, parameter = "tau", type = "chains", plot_type = 'ggplot', show_models = 3)
acorr_tau <- diagnostics(postcovid_survival_custom, parameter = "tau", type = "autocrrelation", plot_type = 'ggplot',
                   show_models = 3)
densities_tau <- diagnostics(postcovid_survival_custom, parameter = "tau", type = "densities", plot_type = 'ggplot',
                      show_models = 3)

pdf(file = './figures/fig70_diagnos_postcovid_survival_custom.pdf', width = 12, height = 6)
chains_mean +
  acorr_mean +
  densities_mean +
  chains_tau +
  acorr_tau +
  densities_tau +
  plot_annotation(tag_levels = 'A')
dev.off()

# Forest Plots Custom
forest <- forest(postcovid_survival_custom, output_scale = 'r', plot_type = 'ggplot')
forest_cond <- forest(postcovid_survival_custom, conditional = T, output_scale = 'r', plot_type = 'ggplot')
pdf(file = './figures/fig71_postcovid_survival_forestplot_custom.pdf', width = 6, height = 25)
forest
dev.off()

# Plotting posterior probabilities Custom
esize <- plot(postcovid_survival_custom, parameter = "mu", output_scale = 'r', prior = T, plot_type = 'ggplot',
           main = 'PMP for the Effect of post-COVID on Survival', xlab = 'Mean')
het <- plot(postcovid_survival_custom, parameter = "tau", output_scale = 'r', prior = T, plot_type = 'ggplot',
         main = 'PMP for Heterogeneity of post-COVID on Survival', xlab = 'Heterogeneity')
selm <- plot(postcovid_survival_custom, parameter = "weightfunction", output_scale = 'r', prior = T, rescale_x = T,
          plot_type = 'ggplot', main = 'PMP Plot for the Selection Models of post-COVID on Survival',
          xlab = 'Selection models publication bias adjustments')
petpeese <- plot(postcovid_survival_custom, parameter = "PET-PEESE", output_scale = 'r', rescale_x = T, prior = T,
             plot_type = 'ggplot', main = 'PMP for the PET-PEESE Models of post-COVID on Survival',
             xlab = 'PET-PEESE publication bias adjustments')

pdf(file = './figures/fig72_postcovid_survival_posteriors_custom.pdf', width = 10, height = 6)
ggarrange(esize, het, selm, petpeese, labels = c('A', 'B', 'C', 'D'), ncol = 2, nrow = 2, font.label = list(face = 'bold'))
dev.off()
