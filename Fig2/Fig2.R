library(tidyverse)
library(lme4); library(emmeans); library(lmerTest)

# - Read in data --------------------------------------------------------------
ctrl <- read_csv("ctrl.csv")
brn3b <- read_csv("brn3b-dta.csv")

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

# Fig 2A: Ctrl 48 hrs
two_day_plot(ctrl, "black")

# Fig 2B: Brn3b-DTA 48 hrs
two_day_plot(brn3b, "navyblue")

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
  calc_amplitude(ctrl, "Ctrl"),
  calc_amplitude(brn3b, "Brn3b-DTA")
) %>%
  mutate(Geno = factor(Geno, levels = c("Ctrl", "Brn3b-DTA")))

# Fig 2C: diurnal amplitudes
amplitudes %>%
  group_by(Geno) %>%
  summarize(mean = mean(Temp),
            sd = sd(Temp)) %>%
  ggplot(aes(x = Geno, y = mean)) +
  geom_col(fill = "white", color = "black") +
  geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.4) +
  geom_jitter(data = amplitudes, aes(y = Temp, color = Geno), 
              width = 0.3, show.legend = FALSE) +
  scale_color_manual(values = c("black", "navyblue")) +
  scale_x_discrete(labels = c("Control", "Brn3b-DTA")) +
  scale_y_continuous(expand = c(0, 0), limits = c(0, 1.6), 
                     breaks = c(0, 0.5, 1, 1.5)) +
  labs(y = "Diurnal amplitude (°C)", x = element_blank()) +
  theme_classic() +
  theme(axis.text = element_text(color = "black"),
        axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))

t.test(
  filter(amplitudes, Geno == "Ctrl")$Temp,
  filter(amplitudes, Geno == "Brn3b-DTA")$Temp
)$p.value


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
    geom_rect(aes(xmin = 1+19/24, xmax = 2, ymin = -1.9, ymax = -1.75),
              color = "black", fill = "black") +
    geom_rect(aes(xmin = 1+20/24, xmax = 1+23/24, ymin = -1.9, ymax = -1.75),
              color = "black", fill = "white") +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.4, fill = color) +
    scale_x_continuous(limits = c(1+19/24, 2), expand = c(0, 0),
                       breaks = seq(1+19/24, 2, by = 1/24),
                       labels = seq(13, 18, by = 1)) +
    scale_y_continuous(limits = c(-1.9, 1),
                       breaks = c(-1, 0, 1),
                       expand = c(0, 0)) +
    labs(y = "Change in body temperature (°C)", x = "Zeitgeber time (hr)") +
    theme_classic() +
    theme(axis.text = element_text(color = "black"))
}

# Fig 2D: Ctrl light pulse
pulse_plot(ctrl, "black")

# Fig 2F: Brn3b-DTA light pulse
pulse_plot(brn3b, "navyblue")


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

ctrl_paired <- paired_data(ctrl)
brn3b_paired <- paired_data(brn3b)

# plot data
paired_plot <- function(df, color) {
  df %>%
    ggplot(aes(x = Epoch, y = Temp, group = Mouse)) +
    geom_line(color = color) +
    geom_point(color = color) +
    scale_y_continuous(limits = c(34.5, 39), expand = c(0, 0)) +
    labs(y = "Mean body temperature (°C)", x = element_blank()) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"),
          axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
}

# Fig 2E: Ctrl paired
paired_plot(ctrl_paired, "black")

# Fig 2G: Brn3-DTA paired
paired_plot(brn3b_paired, "navyblue")

# paired t-tests of Dark vs. Light
paired_stats <- function(df) {
  t.test(
    filter(df, Epoch == "Dark")$Temp,
    filter(df, Epoch == "Light")$Temp,
    paired = TRUE
  )$p.value
}

paired_stats(ctrl_paired)
paired_stats(brn3b_paired)


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
  relative_data(ctrl, 20/24, 3, "Control"),
  relative_data(brn3b, 20/24, 3, "Brn3b-DTA"),
  relative_data(ctrl, 1+20/24, 3, "Control"),
  relative_data(brn3b, 1+20/24, 3, "Brn3b-DTA")
) %>%
  mutate(Geno = factor(Geno, levels = c("Control", "Brn3b-DTA")))

# Fig 2-S1: relative change in each geno
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
  scale_color_manual(values = c("black", "navyblue")) +
  geom_hline(aes(yintercept = 0)) +
  labs(y = "Mean change in body temperature (C)", x = element_blank()) +
  facet_wrap(~Geno) +
  theme_classic() +
  theme(axis.text.x = element_text(color = "black", size = 10))

relative_model <- lmer(Change ~ Geno * Expt + (1|Mouse), data = relative)
pairs(emmeans(relative_model, ~ Expt | Geno))
pairs(emmeans(relative_model, ~ Geno | Expt))

