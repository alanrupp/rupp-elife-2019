library(tidyverse)
library(lme4); library(emmeans); library(lmerTest)

# - Read in data --------------------------------------------------------------
wt <- read_csv("wt.csv")
mko <- read_csv("mko.csv")
mo <- read_csv("mo.csv")

# - Plotting temperature over 48 hrs ------------------------------------------
# function to make plot
two_day_plot <- function(df, color) {
  df %>%
    gather(-Time, key = "Mouse", value = "Temp") %>%
    group_by(Time) %>%
    summarize(mean = mean(Temp),
              sem = sd(Temp)) %>%
    ggplot(aes(x = Time, y = mean)) +
    geom_line(color = color, size = 1) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), 
                alpha = 0.4, fill = color) +
    geom_rect(aes(ymin = 35, ymax = 35.2, xmin = 0, xmax = 2), 
              color = "black", fill = "black") +
    geom_rect(aes(ymin = 35, ymax = 35.2, xmin = 0.25, xmax = 0.75), 
              color = "black", fill = "white") +
    geom_rect(aes(ymin = 35, ymax = 35.2, xmin = 1.25, xmax = 1.75), 
              color = "black", fill = "white") +
    geom_rect(aes(ymin = 35, ymax = 35.2, xmin = 1+20/24, xmax = 1+23/24), 
              color = "black", fill = "white") +
    scale_x_continuous(limits = c(0, 2), 
                       breaks = seq.int(0, 2, length.out = 9),
                       labels = rep(c(18, 0, 6, 12), length.out = 9),
                       expand = c(0, 0)) +
    scale_y_continuous(limits = c(35, 38.6), expand = c(0, 0),
                       breaks = seq(35, 38.5, by = 0.5),
                       labels = c(35, "", 36, "", 37, "", 38, "")) +
    labs(x = "Zeitgeber time (hr)", y = "Body temperature (°C)") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"))
}

# Fig 1B: WT 48 hrs
two_day_plot(wt, "black")

# Fig 1D: MO 48 hrs
two_day_plot(mo, "forestgreen")

# Fig 1E: MKO 48 hrs
two_day_plot(mko, "maroon")


# - Find diurnal amplitude ----------------------------------------------------
calc_amplitude <- function(df, geno) {
  df %>%
    gather(-Time, key = "Mouse", value = "Temp") %>%
    mutate("tod" = case_when(
      Time > 0.25 & Time <= 0.75 ~ "Day",
      Time > 0.75 & Time <= 1.25 ~ "Night"
    )) %>%
    group_by(tod, Mouse) %>%
    summarize(Temp = mean(Temp)) %>%
    group_by(Mouse) %>%
    mutate(Temp = Temp - lag(Temp)) %>%
    filter(!is.na(tod) & !is.na(Temp)) %>%
    select(-tod) %>%
    mutate("Geno" = geno)
}

# calculate all and combine
amplitudes <- bind_rows(
  calc_amplitude(wt, "WT"),
  calc_amplitude(mo, "MO"),
  calc_amplitude(mko, "MKO")
  ) %>%
  mutate(Geno = factor(Geno, levels = c("WT", "MO", "MKO")))

# Fig 1F: diurnal amplitudes
amplitudes %>%
  group_by(Geno) %>%
  summarize(mean = mean(Temp),
            sd = sd(Temp)) %>%
  ggplot(aes(x = Geno, y = mean)) +
  geom_col(fill = "white", color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.4) +
  geom_jitter(data = amplitudes, aes(y = Temp, color = Geno), 
              width = 0.3, show.legend = FALSE) +
  scale_color_manual(values = c("black", "forestgreen", "maroon")) +
  scale_x_discrete(labels = c("Wildtype", "Melanopsin-only", "Melanopsin KO")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.6), 
                     breaks = c(0, 0.5, 1, 1.5)) +
  labs(y = "Diurnal amplitude (°C)", x = element_blank()) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

# ANOVA model comparing amplitude across genotypes
amplitude_model <- aov(Temp ~ Geno, data = amplitudes)
summary(amplitude_model)


# - Plots of relative body temperature during light pulse ---------------------
# find time closest to 20:00
find_pulse_start <- function(df, time0) {
  which(abs(df$Time - time0) == min(abs(df$Time - time0)))
}

make_relative <- function(df, time0) {
  df2 <- apply(df[,2:ncol(df)], 2, function(x) x - x[find_pulse_start(df, time0)])
  df2 <- as.data.frame(df2) %>%
    mutate(Time = df$Time)
}


pulse_plot <- function(df, color) {
  # plot
  make_relative(df, 1+20/24) %>% 
    gather(-Time, key = "Mouse", value = "Temp") %>%
    group_by(Time) %>%
    summarize(mean = mean(Temp),
              sd = sd(Temp)) %>%
    ggplot(aes(x = Time, y = mean)) +
    geom_line(color = color) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_rect(aes(xmin = 1+19/24, xmax = 2, ymin = -1.6, ymax = -1.48),
              color = "black", fill = "black") +
    geom_rect(aes(xmin = 1+20/24, xmax = 1+23/24, ymin = -1.6, ymax = -1.48),
              color = "black", fill = "white") +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.4, fill = color) +
    scale_x_continuous(limits = c(1+19/24, 2), expand = c(0, 0),
                       breaks = seq(1+19/24, 2, by = 1/24),
                       labels = seq(13, 18, by = 1)) +
    scale_y_continuous(limits = c(-1.6, 1),
                       breaks = c(-1, 0, 1),
                       expand = c(0, 0)) +
    labs(y = "Change in body temperature (°C)", x = "Zeitgeber time (hr)") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"))
}

# Fig 1C: WT light pulse
pulse_plot(wt, "black")

# Fig 1G: MO light pulse
pulse_plot(mo, "forestgreen")

# Fig 1I: MKO light pulse
pulse_plot(mko, "maroon")


# - Paired plots of body temperature during control and light pulse -----------
paired_data <- function(df) {
  df %>%
    gather(-Time, key = "Mouse", value = "Temp") %>%
    mutate("Epoch" = case_when(
      Time > 20/24 & Time <= 23/24 ~ "Dark",
      Time > 1+20/24 & Time <= 1+23/24 ~ "Light"
    )) %>%
    filter(!is.na(Epoch)) %>%
    group_by(Mouse, Epoch) %>%
    summarize(Temp = mean(Temp))
}

wt_paired <- paired_data(wt)
mo_paired <- paired_data(mo)
mko_paired <- paired_data(mko)

# plot data
paired_plot <- function(df, color) {
  df %>%
    ggplot(aes(x = Epoch, y = Temp, group = Mouse)) +
    geom_line(color = color) +
    geom_point(color = color) +
    scale_y_continuous(limits = c(35, 38.6), expand = c(0, 0)) +
    labs(y = "Mean body temperature (°C)", x = element_blank()) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

# Fig 1H: MO paired
paired_plot(mo_paired, "forestgreen")

# Fig 1J: MKO paired
paired_plot(mko_paired, "maroon")

# paired t-tests of Dark vs. Light
paired_stats <- function(df) {
  t.test(
    filter(df, Epoch == "Dark")$Temp,
    filter(df, Epoch == "Light")$Temp,
    paired = TRUE
  )$p.value
}

paired_stats(wt_paired)
paired_stats(mo_paired)
paired_stats(mko_paired)


# - Sensitivity --------------------------------------------------------------
lux_files <- list.files(pattern = "*lux")
lux_expts <- str_extract(lux_files, "(.+)(?=lux)")
lux <- map(lux_files, read_csv)
# remove 203 from 100lux2 experiment
lux[[3]]$`203` <- NA
lux <- map(1:length(lux), ~ mutate(lux[[.x]], "Intensity" = lux_expts[.x])) %>%
  bind_rows()

# keep data for 2 days only
lux <- filter(lux, Time >= 0 & Time <= 2)

# make time equal to closest "true" sequence to match data across files
true_time <- function(values) {
  true_seq <- seq(0, 2, by = 1/24/60) %>% round(4)
  best_time <- function(value) {
    differences <- abs(true_seq - value)
    hit <- which(differences == min(differences))
    return(true_seq[hit])
  }
  return(unlist(map(values, best_time)))
}

lux <- mutate(lux, Time = true_time(Time))

sensitivity_plot <- function(lux_df) {
  lux_df <- lux_df %>%
    gather(-Time, -Intensity, key = "Mouse", value = "Temp")
  
  pulses <-
    lux_df %>%
    filter(Time > 1+19/24 & Time <= 1+24/24) %>%
    group_by(Mouse, Time, Intensity) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE)) %>%
    group_by(Time, Intensity) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE))
  
  # use previous night from 1000 lux experiment as control
  controls <-
    lux_df %>%
    filter(Time > 19/24 & Time <= 24/24) %>%
    filter(Intensity == 1000) %>%
    mutate(Intensity = as.character(0)) %>%
    group_by(Time, Intensity) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE)) %>%
    ungroup() %>%
    mutate(Time = Time + 1)
  
  bind_rows(pulses, controls) %>%
    ggplot(aes(x = Time, y = Temp, color = Intensity)) +
    geom_line(show.legend = FALSE) +
    geom_rect(aes(xmin = 1+19/24, xmax = 2, ymin = 35.5, ymax = 35.7),
              color = "black", fill = "black") +
    geom_rect(aes(xmin = 1+20/24, xmax = 1+23/24, ymin = 35.5, ymax = 35.7),
              color = "black", fill = "white") +
    scale_color_manual(values = paste0("gray", seq(70, 0, -(70/5)))) +
    scale_x_continuous(limits = c(1+19/24, 1+24/24),
                       breaks = seq(1+19/24, 1+24/24, 1/24),
                       expand = c(0,0), labels = seq(13, 18, 1)) +
    scale_y_continuous(limits = c(35.5, 38.5), expand = c(0, 0)) +
    labs(x = "Zeitgeber time (hr)", y = "Body temperature (°C)") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"))
}

# Fig 1-S1B: Sensitivity timecourse
sensitivity_plot(lux)

# Find mean temperature during the pulse
find_mean_sensitivity <- function(lux_df) {
  lux_df <- lux_df %>%
    gather(-Time, -Intensity, key = "Mouse", value = "Temp")
  
  pulses <-
    lux_df %>%
    filter(Time > 1+20/24 & Time <= 1+23/24) %>%
    group_by(Mouse, Time, Intensity) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE)) %>%
    group_by(Mouse, Intensity) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE))
  
  # all preceding ZT14-17 averaged
  controls <-
    lux_df %>%
    filter(Time > 20/24 & Time <= 23/24) %>%
    mutate(Intensity = as.character(0)) %>%
    group_by(Mouse, Intensity) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE))
  
  return(bind_rows(pulses, controls))
}
  
  
sensitivity_mean <- find_mean_sensitivity(lux)
sensitivity_mean$Intensity <- as.numeric(sensitivity_mean$Intensity)

sensitivity_mean_plot <- function(sensitivity_mean_df) {
  sensitivity_mean %>%
    group_by(Intensity) %>%
    summarize(mean = mean(Temp),
              sd = sd(Temp)) %>%
    ggplot(aes(x = log10(Intensity + 0.1), y = mean)) +
    geom_point(size = 2.5) +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.1,
                  size = 0.8) +
    geom_smooth(se = FALSE, color = "black", span = 2, size = 1) +
    scale_x_continuous(labels = c(0, 1, 10, 100, 1000)) +
    scale_y_continuous(limits = c(36, 38), expand = c(0, 0)) +
    labs(x = "Light intensity (lux)", y = "Mean body temperature (°C)") +
    theme_bw() +
    theme(panel.grid = element_blank(),
          axis.text = element_text(color = "black"))
}

# Fig 1-S1C
sensitivity_mean_plot(sensitivity_mean)

# sensitivity statistical model
sensitivity_model <- lmer(Temp ~ log10(Intensity + 0.1) + (1|Mouse), 
                          data = sensitivity_mean)
summary(sensitivity_model)


# - Relative body temperature change -----------------------------------------
relative_data <- function(df, time0, duration, geno) {
  relative <- make_relative(df, time0)
  relative <- filter(relative, Time >= round(time0, 3) & 
                       Time <= round(time0,3) + duration / 24)
  relative <- select(relative, -Time)
  relative <- colMeans(relative) %>% as.data.frame() %>%
    set_names("Change") %>%
    rownames_to_column("Mouse") %>%
    mutate("Geno" = geno)
  
  if (time0 < 1) {
    expt <- "Ctrl"
  } else {
    expt <- "Light"
  }
  relative <- mutate(relative, Expt = expt)
  return(relative)
}

# calculate relative change during light pulse
relative <- bind_rows(
  relative_data(wt, 20/24, 3, "WT"),
  relative_data(mo, 20/24, 3, "MO"),
  relative_data(mko, 20/24, 3, "MKO"),
  relative_data(wt, 1+20/24, 3, "WT"),
  relative_data(mo, 1+20/24, 3, "MO"),
  relative_data(mko, 1+20/24, 3, "MKO")
) %>%
  mutate(Geno = factor(Geno, levels = c("WT", "MO", "MKO"),
         labels = c("Wildtype", "Melanopsin-only", "Melanopsin KO")))

# Fig 1-S2: relative change in each geno
relative %>%
  group_by(Geno, Expt) %>%
  summarize(mean = mean(Change),
            sd = sd(Change)) %>%
  ggplot(aes(x = Expt, y = mean)) +
  geom_col(aes(fill = Expt), color = "black", show.legend = FALSE) +
  scale_fill_manual(values = c("white", "gray")) +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd),
                width = 0.4) +
  geom_jitter(data = relative, aes(y = Change, color = Geno),
              show.legend = FALSE, width = 0.2, height = 0) +
  scale_color_manual(values = c("black", "forestgreen", "maroon")) +
  geom_hline(aes(yintercept = 0)) +
  labs(y = "Mean change in body temperature (°C)", x = element_blank()) +
  facet_wrap(~Geno) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 10))

relative_model <- lmer(Change ~ Geno * Expt + (1|Mouse), data = relative)
pairs(emmeans(relative_model, ~ Expt | Geno))
pairs(emmeans(relative_model, ~ Geno | Expt))
