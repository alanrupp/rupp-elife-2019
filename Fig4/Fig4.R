library(tidyverse)
library(lme4); library(lmerTest); library(emmeans)

# - 24-hr sleep / wake with 1-hr bins -----------------------------------------
# one-hour binned data
one_hr <- readxl::read_xlsx('Sleep Raw Data.xlsx', sheet = 2) %>%
  gather(-Hour, key = "Mouse", value = "Sleep") %>%
  mutate("Genotype" = ifelse(str_detect(Mouse, "^Control"), "Control", "Brn3b-DTA")) %>%
  mutate(Genotype = factor(Genotype, levels = c("Control", "Brn3b-DTA"))) %>%
  mutate(Hour = Hour + 0.5)

one_day_plot <- function(one_hr_bin_df, genotype, color) {
  one_hr_bin_df %>%
    filter(Genotype == genotype) %>%
    group_by(Hour) %>%
    summarize(mean = mean(Sleep), sem = sd(Sleep)/sqrt(n())) %>%
    ggplot(aes(x = Hour, y = mean)) +
    geom_rect(aes(xmin = 0, xmax = 24, ymin = 0, ymax = 4),
              fill = "black", color = "black") +
    geom_rect(aes(xmin = 0, xmax = 12, ymin = 0, ymax = 4),
              fill = "white", color = "black") +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem), 
                alpha = 0.4, show.legend = FALSE, fill = color) +
    geom_line(show.legend = FALSE, color = color) +
    theme_classic() +
    scale_x_continuous(expand = c(0,0), limits = c(0, 24), 
                       breaks = seq(0, 24, 6)) +
    scale_y_continuous(expand = c(0,0), limits = c(0, 100)) + 
    labs(x = 'Zeitgeber time (hr)', y = 'Time spent asleep (%)')
}

# Fig 4A: Control 24 hr
one_day_plot(one_hr, "Control", "black")

# Fig 4B: Brn3b-DTA 24 hr
one_day_plot(one_hr, "Brn3b-DTA", "navyblue")

# statistical model of sleep by genotype by one-hour bin
one_hr_model <- lmer(Sleep ~ Genotype * factor(Hour) + (1|Mouse), data = one_hr)
pairs(emmeans(one_hr_model, ~ Genotype | Hour))


# - Total Wake, Sleep, NREM, and REM -----------------------------------------
wake_sleep <- readxl::read_xlsx('Sleep Raw Data.xlsx', sheet = 1, n_max = 2) %>%
  rename("State" = X__1) %>% 
  gather(-State, key = "Mouse", value = "Percent") %>%
  mutate("Genotype" = ifelse(str_detect(Mouse, "^Control"), "Control", "Brn3b-DTA")) %>%
  mutate(Genotype = factor(Genotype, levels = c("Control", "Brn3b-DTA"))) %>%
  mutate(State = factor(State, levels = c('Wake', 'Sleep')))

state_percent_plot <- function(state_df) {
  state_df %>%
    group_by(Genotype, State) %>%
    summarize(mean = mean(Percent), sem = sd(Percent)/sqrt(n())) %>% 
    ggplot(aes(x = Genotype, y = mean)) +
    geom_col(position = 'dodge', color = 'black', 
             fill = 'white', show.legend = FALSE) +
    geom_errorbar(aes(ymin= mean - sem, ymax = mean + sem), 
                  width = 0.4, color = 'black') +
    geom_jitter(data = state_df, 
                aes(x = Genotype, y = Percent, color = Genotype),
              show.legend = FALSE, width = 0.2, alpha = 0.5, stroke = 0) +
    scale_y_continuous(expand = c(0, 0), 
                       limits = c(0, max(state_df$Percent*1.1))) +
    scale_color_manual(values = c('black', 'navyblue')) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1)) +
    labs(y = 'Total time (%)', x = element_blank()) +
    facet_wrap(~State)
}

# Fig 4C: total time spent wake / sleep
state_percent_plot(wake_sleep)

# stats
t.test(
  filter(wake_sleep, State == 'Wake' & Genotype == 'Control')$Percent,
  filter(wake_sleep, State == 'Wake' & Genotype == 'Brn3b-DTA')$Percent
  )$p.value

# NREM / REM changes
nrem_rem <- readxl::read_xlsx('Sleep Raw Data.xlsx', sheet = 1) %>%
  rename("State" = X__1) %>% 
  filter(State %in% c("NREM", "REM")) %>%
  gather(-State, key = "Mouse", value = "Percent") %>%
  mutate("Genotype" = ifelse(str_detect(Mouse, "^Control"), "Control", "Brn3b-DTA")) %>%
  mutate(Genotype = factor(Genotype, levels = c("Control", "Brn3b-DTA"))) %>%
  mutate(State = factor(State, levels = c('NREM', 'REM')))

# Fig 4-S1A: total time spent nrem / rem
state_percent_plot(nrem_rem) +
  scale_y_continuous(limits = c(0, 100), expand = c(0, 0)) +
  ylab("24-hr sleep time\n(% total sleep)")

# t-test of NREM
t.test(
  filter(nrem_rem, State == 'NREM' & Genotype == 'Control')$Percent,
  filter(nrem_rem, State == 'NREM' & Genotype == 'Brn3b-DTA')$Percent
)$p.value

# t-test of REM (should be same as NREM)
t.test(
  filter(nrem_rem, State == 'REM' & Genotype == 'Control')$Percent,
  filter(nrem_rem, State == 'REM' & Genotype == 'Brn3b-DTA')$Percent
)$p.value


# - 3-hr Light Pulse ----------------------------------------------------------
reshape_pulse <- function(geno, skip, n_max) {
  df <- 
    readxl::read_xlsx('Sleep Raw Data.xlsx', sheet = 4, skip = skip,
                      n_max = n_max) %>%
    rename("Time" = `Time Post Light Pulse`) %>%
    gather(-Time, key = "Mouse", value = "Percent") %>%
    mutate("Treatment" = str_extract(Mouse, "^[A-Za-z]+")) %>%
    mutate(Mouse = ifelse(str_detect(Mouse, "[0-9]$"), 
                          str_extract(Mouse, "[0-9]+$"), 0)) %>%
    mutate("Genotype" = geno) %>%
    mutate(Mouse = paste(Genotype, Mouse, sep = "_")) %>%
    mutate(Time = Time - 15)
  return(df)
}

# read in data and add annotation
pulse_30min <- bind_rows(
  reshape_pulse("Control", 1, 6),
  reshape_pulse("Brn3b-DTA", 9, 6)
) %>%
  mutate(Genotype = factor(Genotype, levels = c("Control", "Brn3b-DTA"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Baseline", "Light")))


# Plot Percent sleep with time
plot_light_pulse <- function(pulse_df, geno, colors) {
  pulse_df %>%
    filter(Genotype == geno) %>%
    group_by(Treatment, Time) %>%
    summarize(mean = mean(Percent),
              sem = sd(Percent)/sqrt(n())) %>%
    ggplot(aes(x = Time, y = mean, group = Treatment)) +
    geom_line(aes(fill = Treatment), show.legend = FALSE) +
    geom_ribbon(aes(ymin = mean - sem, ymax = mean + sem, fill = Treatment),
                alpha = 0.4, show.legend = FALSE) +
    scale_fill_manual(values = colors) +
    scale_color_manual(values = colors) +
    scale_y_continuous(limits = c(0, 100), expand = c(0,0)) +
    labs(x = "Time from lights-on (min)", y = "Time spent asleep (%)") +
    theme_classic()
}

# Fig 4D
plot_light_pulse(pulse_30min, "Control", c("black", "gray50"))

# Fig 4F
plot_light_pulse(pulse_30min, "Brn3b-DTA", c("navyblue", "lightblue"))

# Light pulse model
pulse_model <- lmer(Percent ~ Time * sqrt(Time) * Genotype * Treatment + (1|Mouse), 
                    data = pulse_30min)
pairs(emmeans(
  ref_grid(pulse_model, cov.reduce = FALSE), ~ Treatment | Genotype + Time)
  )



# - Paired summary -----------------------------------------------------------
reshape_paired <- function(geno, skip, n_max) {
  df <- 
    readxl::read_xlsx('Sleep Raw Data.xlsx', sheet = 3, skip = skip,
                      n_max = n_max)
  colnames(df)[1] <- "State"
  df <- df %>%
    gather(-State, key = "Mouse", value = "Percent") %>%
    mutate("Treatment" = str_extract(Mouse, "^[A-Za-z]+")) %>%
    mutate(Mouse = ifelse(str_detect(Mouse, "[0-9]$"), 
                          str_extract(Mouse, "[0-9]+$"), 0)) %>%
    mutate("Genotype" = geno) %>%
    mutate(Mouse = paste(Genotype, Mouse, sep = "_"))
  return(df)
}

paired <- bind_rows(
  reshape_paired("Control", 0, 4),
  reshape_paired("Brn3b-DTA", 5, 4)
) %>%
  mutate(Genotype = factor(Genotype, levels = c("Control", "Brn3b-DTA"))) %>%
  mutate(Treatment = factor(Treatment, levels = c("Baseline", "Light"))) %>%
  filter(!is.na(Percent))

paired_plot <- function(paired_df, geno, states, color) {
  paired_df <- paired_df %>%
    filter(State %in% states)
  
  # set plot attributes
  if (states == "Sleep") {
    y_title = "Time spent asleep (%)"
    y_max = max(paired_df$Percent) * 1.1
  } else if (states == c("NREM", "REM")) {
    y_title = "ZT14-17 sleep time\n(% total sleep)"
    y_max = 100
  }
  
  # plot
  plt <-
    paired_df %>%
    filter(Genotype == geno) %>%
    ggplot(aes(x = Treatment, y = Percent, group = Mouse)) +
    geom_line(color = color) +
    geom_point(color = color) +
    scale_y_continuous(limits = c(0, y_max), 
                       expand = c(0, 0)) +
    labs(y = y_title, x = element_blank()) +
    theme_classic() +
    theme(axis.text.x = element_text(angle = 45, hjust = 1, vjust = 1))
  
  if (length(states) > 1) {
    plt <- plt + facet_wrap(~State)
  }
  
  return(plt)
}

# Fig 4E: Control sleep baseline vs. light
paired_plot(paired, "Control", "Sleep", "black")

# Fig 4G: Brn3b-DTA sleep baseline vs. light
paired_plot(paired, "Brn3b-DTA", "Sleep", "navyblue")


# stats
paired_stats <- function(paired_df, geno, state) {
  paired_df <- paired_df %>%
    filter(Genotype == geno) %>%
    filter(State == state)
  
  t.test(
    filter(paired_df, Treatment == "Baseline")$Percent,
    filter(paired_df, Treatment == "Light")$Percent,
    paired = TRUE
  )$p.value
}

# paired t-test of effect of light on sleep amount
paired_stats(paired, "Control", "Sleep")
paired_stats(paired, "Brn3b-DTA", "Sleep")


# NREM vs. REM
# Fig 4-S1B: Control NREM/REM baseline vs. light
paired_plot(paired, "Control", c("NREM", "REM"), "black")

# Fig 4-S1C: Brn3b-DTA NREM/REM baseline vs. light
paired_plot(paired, "Brn3b-DTA", c("NREM", "REM"), "navyblue")

# stats
paired_stats(paired, "Control", "NREM")
paired_stats(paired, "Brn3b-DTA", "NREM")
paired_stats(paired, "Control", "REM")
paired_stats(paired, "Brn3b-DTA", "REM")


# - Testing Control vs. Brn3b-DTA for time spent asleep with light pulse ------
paired_combined_plot <- function(paired_df, state) {
  paired_df <- paired_df %>%
    filter(State == state)
  
  mean_values <-
    paired_df %>%
    group_by(Treatment, Genotype) %>%
    summarize(mean = mean(Percent),
              sem = sd(Percent)/sqrt(n()))
  
  mean_values %>%
    ggplot(aes(x = Treatment, y = mean)) +
    geom_col(aes(fill = Treatment), color = "black", show.legend = FALSE) +
    geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.4) +
    geom_jitter(data = paired_df, aes(color = Genotype, y = Percent),
                show.legend = FALSE, width = 0.2) +
    scale_color_manual(values = c("black", "navyblue")) +
    scale_fill_manual(values = c("white", "gray")) +
    scale_y_continuous(limits = c(0, max(paired_df$Percent) * 1.1),
                       expand = c(0, 0)) +
    labs(y = "Time spent asleep (%)", x = element_blank()) +
    facet_wrap(~Genotype) +
    theme_classic()
}

paired_combined_plot(paired, "Sleep")

# linear mixed model of effect of light by genotype
paired_model <- lmer(Percent ~ Genotype * Treatment + (1|Mouse),
                    data = filter(paired, State == "Sleep"))
pairs(emmeans(paired_model, ~ Treatment | Genotype))
pairs(emmeans(paired_model, ~ Genotype | Treatment))
