# Clean the workspace and console
closeAllConnections(); rm(list=ls())
cat("\014"); graphics.off()
#R0 plot
require(bdskytools)
#devtools::install_github("iamamutt/ggdistribute", dependencies=TRUE)
require(ggdistribute)
require(coda)
require(dplyr)
require(ggplot2)
require(dplyr)
require(wesanderson)
library(viridis)

#HONG_KONG

setwd("/Users/rhysinward/Documents/Master's_paper/Data/Hong_Kong/BDSKY/violin_plot")

WORKDIR = "/Users/rhysinward/Documents/Master's_paper/Data/Hong_Kong/BDSKY/violin_plot"

#uniform

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Uniform/uniform_sampling_alligned_ORF1AB_removed.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Uniform/datatable_uniform.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Uniform")
deme_labels <- c("Uniform")
params <- paste("reproductiveNumber_BDSKY_Serial.1")                               

logfile_long_uniform <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_uniform$Deme <- factor(
  x = logfile_long_uniform$Deme,
  levels = params,
  labels = deme_labels)

#unsampled

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Unsampled/unsampled_new_Sampling.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Unsampled/datatable_unsampled.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Unsampled")
deme_labels <- c("Unsampled")
params <- paste("reproductiveNumber_BDSKY_Serial.1")                               

logfile_long_unsampled <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_unsampled$Deme <- factor(
  x = logfile_long_unsampled$Deme,
  levels = params,
  labels = deme_labels)

#proportional

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Proportional/Hong_Kong_case_sample_alligned_ORF1AB.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Proportional/datatable_proportional.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Proportional")
deme_labels <- c("Proportional")
params <- paste("reproductiveNumber_BDSKY_Serial.1")                               

logfile_long_proportional <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_proportional$Deme <- factor(
  x = logfile_long_proportional$Deme,
  levels = params,
  labels = deme_labels)

#inverse

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Inverse/inverse_proportion_Hong_Kong.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Inverse/datatable_inverse.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Reciprocal-proportional")
deme_labels <- c("Reciprocal-proportional")
params <- paste("reproductiveNumber_BDSKY_Serial.1")                               

logfile_long_inverse <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_inverse$Deme <- factor(
  x = logfile_long_inverse$Deme,
  levels = params,
  labels = deme_labels)

#merge 

logfile_long <- rbind (logfile_long_inverse,logfile_long_uniform,logfile_long_proportional,logfile_long_unsampled)

#plot

HongKong <- ggplot(logfile_long) +
  aes(x=Re, y=Deme) +
  geom_posterior(
    aes(color=),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="Sampling Scheme", x=expression(R[0])) +
  xlim(0,8)+
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "none")

#Brazil

setwd("/Users/rhysinward/Documents/Master's_paper/Data/Brazil/BDSKY")

WORKDIR = "/Users/rhysinward/Documents/Master's_paper/Data/Brazil/BDSKY"

#uniform

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Uniform/trail_2/uniform_Brazil.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Uniform/trail_2/datatable_uniform.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Uniform")
deme_labels <- c("Uniform")
params <- paste("reproductiveNumber_BDSKY_Serial.1")                               

logfile_long_uniform <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_uniform$Deme <- factor(
  x = logfile_long_uniform$Deme,
  levels = params,
  labels = deme_labels)

#unsampled

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Unsampled/BDSKY_Trail_2/unsampled_Brazil.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Unsampled/BDSKY_Trail_2/datatable_unsampled.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Unsampled")
deme_labels <- c("Unsampled")
params <- paste("reproductiveNumber_BDSKY_Serial.1")                               

logfile_long_unsampled <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_unsampled$Deme <- factor(
  x = logfile_long_unsampled$Deme,
  levels = params,
  labels = deme_labels)

#proportional

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Proportional/trail_2/beast.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Proportional/trail_2/datatable_proportional.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Proportional")
deme_labels <- c("Proportional")
params <- paste("reproductiveNumber_BDSKY_Serial.1")                               

logfile_long_proportional <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_proportional$Deme <- factor(
  x = logfile_long_proportional$Deme,
  levels = params,
  labels = deme_labels)

#inverse

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Inverse/trail_2/inverse_Brazil.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Inverse/trail_2/datatable_inverse.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Reciprocal-proportional")
deme_labels <- c("Reciprocal-proportional")
params <- paste("reproductiveNumber_BDSKY_Serial.1")                               

logfile_long_inverse <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_inverse$Deme <- factor(
  x = logfile_long_inverse$Deme,
  levels = params,
  labels = deme_labels)

#merge 

logfile_long <- rbind (logfile_long_inverse,logfile_long_uniform,logfile_long_proportional,logfile_long_unsampled)

#plot

Brazil <- ggplot(logfile_long) +
  aes(x=Re, y=Deme) +
  geom_posterior(
    aes(color=),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="Sampling Scheme", x=expression(R[0])) +
  scale_x_continuous(breaks=seq(0, 8, 2)) +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "none")


#plot TMRCA 

#HONG_KONG

setwd("/Users/rhysinward/Documents/Master's_paper/Data/Hong_Kong/BDSKY/violin_plot")

WORKDIR = "/Users/rhysinward/Documents/Master's_paper/Data/Hong_Kong/BDSKY/violin_plot"

#uniform

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Uniform/uniform_sampling_alligned_ORF1AB_removed.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Uniform/TMRCA_uniform.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Uniform")
deme_labels <- c("Uniform")
params <- paste("origin_BDSKY_Serial")                               

logfile_long_uniform <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_uniform$Deme <- factor(
  x = logfile_long_uniform$Deme,
  levels = params,
  labels = deme_labels)

#unsampled

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Unsampled/unsampled_new_Sampling.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Unsampled/TMRCA_unsampled.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Unsampled")
deme_labels <- c("Unsampled")
params <- paste("origin_BDSKY_Serial")                               

logfile_long_unsampled <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_unsampled$Deme <- factor(
  x = logfile_long_unsampled$Deme,
  levels = params,
  labels = deme_labels)

#proportional

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Proportional/Hong_Kong_case_sample_alligned_ORF1AB.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Proportional/TMRCA_proportional.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Proportional")
deme_labels <- c("Proportional")
params <- paste("origin_BDSKY_Serial")                               

logfile_long_proportional <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_proportional$Deme <- factor(
  x = logfile_long_proportional$Deme,
  levels = params,
  labels = deme_labels)

#inverse

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Inverse/inverse_proportion_Hong_Kong.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Inverse/TMRCA_inverse.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Reciprocal-proportional")
deme_labels <- c("Reciprocal-proportional")
params <- paste("origin_BDSKY_Serial")                               

logfile_long_inverse <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_inverse$Deme <- factor(
  x = logfile_long_inverse$Deme,
  levels = params,
  labels = deme_labels)

#merge 

logfile_long <- rbind (logfile_long_inverse,logfile_long_uniform,logfile_long_proportional,logfile_long_unsampled)

logfile_long$Re <- 2020.2664 - logfile_long$Re

logfile_long$Re <- date_decimal(logfile_long$Re)

logfile_long$Re <- as.Date(logfile_long$Re)

#plot

HongKong_TMRCA <- ggplot(logfile_long) +
  aes(x=Re, y=Deme) +
  geom_posterior(
    aes(color=),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="Sampling Scheme", x="TMRCA (Time)") +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "none")

#Brazil

setwd("/Users/rhysinward/Documents/Master's_paper/Data/Brazil/BDSKY")

WORKDIR = "/Users/rhysinward/Documents/Master's_paper/Data/Brazil/BDSKY"

#uniform

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Uniform/trail_2/uniform_Brazil.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Uniform/trail_2/TMRCA_uniform.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Uniform")
deme_labels <- c("Uniform")
params <- paste("origin_BDSKY_Serial")                               

logfile_long_uniform <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_uniform$Deme <- factor(
  x = logfile_long_uniform$Deme,
  levels = params,
  labels = deme_labels)

#unsampled

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Unsampled/BDSKY_Trail_2/unsampled_Brazil.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Unsampled/BDSKY_Trail_2/TMRCA_unsampled.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Unsampled")
deme_labels <- c("Unsampled")
params <- paste("origin_BDSKY_Serial")                               

logfile_long_unsampled <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_unsampled$Deme <- factor(
  x = logfile_long_unsampled$Deme,
  levels = params,
  labels = deme_labels)

#proportional

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Proportional/trail_2/beast.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Proportional/trail_2/TMRCA_proportional.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Proportional")
deme_labels <- c("Proportional")
params <- paste("origin_BDSKY_Serial")                               

logfile_long_proportional <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_proportional$Deme <- factor(
  x = logfile_long_proportional$Deme,
  levels = params,
  labels = deme_labels)

#inverse

OUTPUT <- paste(WORKDIR, "figures/r0.pdf", sep = "/")
LOGFILE <- paste(WORKDIR, "/Inverse/trail_2/inverse_Brazil.log", sep = "")
TRACERTABLE <- paste(WORKDIR, "/Inverse/trail_2/TMRCA_inverse.txt", sep = "")
BURNIN <- 0.1  # 10% already excluded by logcombiner
LOGNORM_PRIOR_MEAN <- 0.8
LOGNORM_PRIOR_SD <- 0.5

logfile <- bdskytools::readLogfile(LOGFILE, burnin = BURNIN)
tracer_table <- read.delim(TRACERTABLE)

# Generate samples from under the prior
samples_under_prior <- data.frame(
  prior = rlnorm(
    n = nrow(logfile), meanlog = LOGNORM_PRIOR_MEAN, sdlog = LOGNORM_PRIOR_SD))

# Reformat logfile data
demes <- c("Reciprocal-proportional")
deme_labels <- c("Reciprocal-proportional")
params <- paste("origin_BDSKY_Serial")                               

logfile_long_inverse <- tidyr::pivot_longer(
  data = logfile,
  cols = params,
  names_to = "Deme",
  values_to = "Re")

# Clean up deme names for plotting
logfile_long_inverse$Deme <- factor(
  x = logfile_long_inverse$Deme,
  levels = params,
  labels = deme_labels)

#merge 

logfile_long <- rbind (logfile_long_inverse,logfile_long_uniform,logfile_long_proportional,logfile_long_unsampled)

logfile_long$Re <- 2021.1027 - logfile_long$Re

logfile_long$Re <- date_decimal(logfile_long$Re)

logfile_long$Re <- as.Date(logfile_long$Re)

#plot

Brazil_TRMCA <- ggplot(logfile_long) +
  aes(x=Re, y=Deme) +
  geom_posterior(
    aes(color=),
    midline=NULL,
    mirror=TRUE,
    fill="#FFFFFF",
    draw_sd=FALSE,
    interval_type="hdi",
    vjust=0,
    position=position_spread(
      reverse=TRUE, # order of spreaded groups within panels
      padding=0.6, # shrink heights of distributions
      height=2 # scale by heights within panels
    ),
    adjust=1.5,
  ) +
  labs(y="Sampling Scheme", x="TMRCA (Time)") +
  scale_x_date (date_breaks = "1 month", date_labels = "%Y %b %d") +
  theme_bw() +
  theme(
    text = element_text(size = 15),
    legend.position = "none")

library(ggpubr)
ggarrange(
  HongKong, HongKong_TMRCA, Brazil, Brazil_TRMCA, labels = c("A", "B", 'C','D'),
  common.legend = FALSE
)
  
  

