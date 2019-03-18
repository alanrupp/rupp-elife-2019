library(tidyverse)
library(lubridate)
library(emmeans); library(lme4)

# - Load in data --------------------------------------------------------------
genotypes <- read_csv("Genotypes.csv")
injections <- read_csv("Injections.csv")
injections <- filter(injections, !is.na(Time))
temps <- read_csv("Data.csv", col_types = c("Tcid"))

# - 48 hr plot ----------------------------------------------------------------
two_day_plot <- function(df, geno, color) {
  # add trial number to injection info
  injx <- 
    injections %>%
    arrange(Time) %>%
    group_by(Mouse, Injx) %>%
    mutate(Trial = seq(n()))
  
  # add trial info to data
  find_match <- function(injection) {
    which(df$Mouse == injection$Mouse &
            month(df$Date) == month(injection$Time) &
            day(df$Date) - day(injection$Time) <= 2 &
            day(df$Date) - day(injection$Time) >= 0)
  }
  
  df$Trial <- NA
  injx2 <- filter(injx, Injx == "PBS")
  for (i in 1:nrow(injx2)) {
    df[find_match(injx2[i,]), "Trial"] <- injx2[i, "Trial"]
  }
  df <- filter(df, !is.na(Trial))
  
  # plot data
  df %>%
    filter(!is.na(Temp)) %>%
    inner_join(., filter(genotypes, Genotype == geno), by = "Mouse") %>%
    mutate(Time = hour(Date)/24 + minute(Date)/24/60) %>%
    group_by(Mouse, Trial) %>%
    mutate(Day = day(Date) - min(day(Date))) %>%
    ungroup() %>%
    mutate(Time = Day + Time) %>%
    group_by(Mouse, Time) %>%
    summarize(Temp = mean(Temp)) %>%
    group_by(Time) %>%
    summarize(mean = mean(Temp),
              sd = sd(Temp)) %>%
    ggplot(aes(x = Time, y = mean)) +
    geom_line(color = color) +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                alpha = 0.4, fill = color) +
    geom_vline(aes(xintercept = 20/24), linetype = "dotted") +
    geom_vline(aes(xintercept = 1+20/24), linetype = "dotted") +
    geom_rect(aes(xmin = 0, xmax = 2.25, ymin = 35, ymax = 35.2),
              fill = "black", color = "black") +
    geom_rect(aes(xmin = 0.25, xmax = 0.75, ymin = 35, ymax = 35.2),
              fill = "white", color = "black") +
    geom_rect(aes(xmin = 1.25, xmax = 1.75, ymin = 35, ymax = 35.2),
              fill = "white", color = "black") +
    scale_x_continuous(limits = c(0, 2.25), expand = c(0, 0),
                       breaks = seq(0, 54/24, 6/24),
                       labels = rep(c(18, 0, 6, 12), length.out = 10)) +
    scale_y_continuous(limits = c(35, 38.9), expand = c(0, 0)) +
    labs(x= "Zeitgeber time (hr)", y = "Body temperature (°C)") +
    theme_classic()
}

# Fig 3B: Brn3b-Cre 54 hours
two_day_plot(temps, "Brn3b-Cre", "darkred")

# Fig 3-S2A: WT 54 hours
two_day_plot(temps, "Control", "black")


# - Relative change -----------------------------------------------------------
make_relative <- function(df, injx) {
  # add trial number to injection info
  injx <- 
    injx %>%
    arrange(Time) %>%
    group_by(Mouse, Injx) %>%
    mutate(Trial = seq(n()))
  
  # find matching rows
  find_match <- function(df, injx) {
    which(df$Mouse == injx$Mouse &
            (df$Date - injx$Time) <= 60*60*6 &
            (df$Date - injx$Time) >= -60*60)
  }
  
  each_relative <- function(df, injx) {
    df2 <- df[find_match(df, injx), ]
    df2 <- mutate(df2, Temp = Temp - Temp[which(df2$Date == injx$Time)])
    df2 <- mutate(df2, Date = Date - Date[which(df2$Date == injx$Time)])
    df2 <- mutate(df2, Time = as.numeric(Date)/60/60)
    df2 <- select(df2, Time, Mouse, Temp, -Date)
    df2 <- mutate(df2, "Injx" = injx$Injx, "Trial" = injx$Trial)
  }
  
  relative <-
    map(1:nrow(injx), ~ each_relative(df, injx[.x, ])) %>%
    bind_rows()
  
  return(relative)
}

relative <- make_relative(temps, injections)

# plot 
relative_plot <- function(relative_df, geno, injx, color) {
  relative_df %>%
    filter(Injx == injx) %>%
    inner_join(., filter(genotypes, Genotype == geno), by = "Mouse") %>%
    group_by(Mouse, Time) %>%
    summarize(Temp = mean(Temp)) %>%
    group_by(Time) %>%
    summarize(mean = mean(Temp),
              sd = sd(Temp)) %>%
    ggplot(aes(x = Time, y = mean)) +
    geom_line(color = color) +
    geom_rect(aes(xmin = -1, xmax = 6, ymin = -2, ymax = -1.8),
              fill = "black", color = "black") +
    geom_ribbon(aes(ymin = mean - sd, ymax = mean + sd),
                fill = color, alpha = 0.4) +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    geom_vline(aes(xintercept = 0), linetype = "dotted") +
    scale_x_continuous(breaks = seq(-1, 6, 1), expand = c(0, 0)) +
    scale_y_continuous(limits = c(-2, 1), expand = c(0, 0)) +
    labs(y = "Change in body temperature (°C)", 
         x = "Time since injection (hr)") +
    theme_classic()
}

# Fig 3C
relative_plot(relative, "Brn3b-Cre", "PBS", "gray")

# Fig 3D
relative_plot(relative, "Brn3b-Cre", "CNO", "darkred")

# Fig3-S2B
relative_plot(relative, "Control", "PBS", "gray")

# Fig3-S2C
relative_plot(relative, "Control", "CNO", "black")


# - Paired data ---------------------------------------------------------------
# only keep data that has PBS and CNO for each trial
make_relative_pairs <- function(relative_df) {
  find_pairless <- function(relative_df) {
    relative_df %>%
      group_by(Mouse, Trial, Injx) %>%
      count() %>%
      spread(key = Injx, value = n) %>%
      filter(is.na(CNO) | is.na(PBS))
  }
  
  no_pair <- find_pairless(relative_df)
  
  relative_df <- 
    relative_df %>%
    filter(!(paste(Mouse, Trial) %in% paste(no_pair$Mouse, no_pair$Trial)))
  
  return(relative_df)
}

relative_pairs <- make_relative_pairs(relative)

paired_plot <- function(relative_df, geno) {
  relative_df %>%
    inner_join(., filter(genotypes, Genotype == geno), by = "Mouse") %>%
    mutate(Injx = factor(Injx, levels = c("PBS", "CNO"))) %>%
    filter(Time > 0 & Time <= 6) %>%
    group_by(Mouse, Trial, Injx) %>%
    summarize(Temp = mean(Temp)) %>%
    group_by(Mouse, Injx) %>%
    summarize(Temp = mean(Temp, na.rm = TRUE)) %>%
    ggplot(aes(x = Injx, y = Temp, group = Mouse)) +
    geom_line() +
    geom_point() +
    geom_hline(aes(yintercept = 0), linetype = "dashed") +
    scale_y_continuous(expand = c(0, 0), limits = c(-1.5, 1)) +
    labs(x = element_blank(), y = "Change in body temperature (°C)") +
    theme_classic()
}

# Fig 3E
paired_plot(relative_pairs, "Brn3b-Cre")

# Fig 3-S2D
paired_plot(relative_pairs, "Control")

# paired t-tests for average body temperature change
paired_stats <- function(relative_df, geno) {
  mean_relative <-
    relative_df %>%
    inner_join(., filter(genotypes, Genotype == geno), by = "Mouse") %>%
    mutate(Injx = factor(Injx, levels = c("PBS", "CNO"))) %>%
    filter(Time > 0 & Time <= 6) %>%
    group_by(Mouse, Trial, Injx) %>%
    summarize(Temp = mean(Temp)) %>%
    group_by(Mouse, Injx) %>%
    summarize(Temp = mean(Temp))
  
  t.test(
    filter(mean_relative, Injx == "PBS")$Temp,
    filter(mean_relative, Injx == "CNO")$Temp,
    paired = TRUE
  )$p.value
}

paired_stats(relative_pairs, "Control")
paired_stats(relative_pairs, "Brn3b-Cre")


# - Mean relative change compared between genotypes ---------------------------
# Find relative change for each mouse
find_relative_mean <- function(relative_df) {
  relative_df %>%
    filter(!is.na(Temp)) %>%
    left_join(., genotypes, by = "Mouse") %>%
    mutate(Injx = factor(Injx, levels = c("PBS", "CNO"))) %>%
    mutate(Genotype = factor(Genotype, levels = c("Control", "Brn3b-Cre"))) %>%
    filter(Time > 0 & Time <= 6) %>%
    group_by(Genotype, Mouse, Trial, Injx) %>%
    summarize(Temp = mean(Temp)) %>%
    group_by(Genotype, Mouse, Injx) %>%
    summarize(Temp = mean(Temp))
}

relative_mean <- find_relative_mean(relative)

# plot relative change by genotype and injection
relative_mean_plot <- function(relative_df) {
  mean_values <- 
    relative_df %>%
    group_by(Genotype, Injx) %>%
    summarize(mean = mean(Temp),
              sd = sd(Temp))
    
  mean_values %>%
    ggplot(aes(x = Injx, y = mean, fill = Injx)) +
    geom_col(show.legend = FALSE, color = "black") +
    geom_errorbar(aes(ymin = mean - sd, ymax = mean + sd), width = 0.4) +
    geom_jitter(data = relative_df, aes(y = Temp, color = Genotype),
                show.legend = FALSE, width = 0.2) +
    geom_hline(aes(yintercept = 0)) +
    scale_fill_manual(values = c("white", "gray")) +
    scale_color_manual(values = c("black", "darkred")) +
    labs(y = "Mean change from baseline (°C)", x = element_blank()) +
    facet_wrap(~Genotype) +
    theme_classic() +
    theme(axis.text = element_text(color = "black"))
}

# Fig 3-S1D
relative_mean_plot(relative_mean)

# model of changes
relative_model <- lmer(Temp ~ Genotype * Injx + (1|Mouse), data = relative_mean)
pairs(emmeans(relative_model, ~ Injx | Genotype))
pairs(emmeans(relative_model, ~ Genotype | Injx))


# - Quantify dreadd infection ------------------------------------------------
counts <- read_csv("AAV_counts.csv")

counts <-
  counts %>%
  separate(Image, into = c("Mouse", "Mag", "Image")) %>%
  gather(`Red-not-green`, `Green-not-red`, Both, 
         key = "Cells", value = "Counts") 

# average of each type
counts %>%
  group_by(Mouse, Cells) %>%
  summarize(Counts = mean(Counts)) %>%
  group_by(Cells) %>%
  summarize(mean = mean(Counts),
            sem = sd(Counts)/sqrt(n())) %>%
  mutate(Cells = factor(Cells, levels = c("Red-not-green", "Green-not-red", 
                                          "Both"))) %>%
  arrange(desc(Cells)) %>%
  mutate("cumsum" = cumsum(mean)) %>%
  ggplot(aes(x = "", y = mean, fill = Cells)) +
  geom_col(color = "black", position = "stack") +
  geom_errorbar(aes(ymin = cumsum - sem, ymax = cumsum + sem), 
                width = 0.4) +
  scale_fill_manual(values = c("magenta", "green", "white"),
                    labels = c(expression("Brn3b"^"+"*"\nOpn4"^"-"),
                               expression("Brn3b"^"-"*"\nOpn4"^"+"), 
                               expression("Brn3b"^"+"*"\nOpn4"^"+"))) +
  theme_classic() +
  labs(y = "Cells / image", x = element_blank()) +
  scale_y_continuous(expand = c(0,0), limits = c(0, NA)) +
  theme(axis.text.x = element_text(color = "black", size = 8, angle = 45,
                                   hjust = 1, vjust = 1),
        axis.ticks.x = element_blank())

# proportion of ipRGCs infected by virus
counts %>%
  filter(Cells %in% c("Green-not-red", "Both")) %>%
  group_by(Mouse, Image) %>%
  mutate(Counts = Counts / sum(Counts)) %>%
  group_by(Mouse, Cells) %>%
  summarize(Counts = mean(Counts)) %>%
  group_by(Cells) %>%
  summarize(mean = mean(Counts),
            sem = sd(Counts)/sqrt(n())) %>%
  filter(Cells == "Both") %>%
  ggplot(aes(x = Cells, y = mean)) +
  geom_col(color = "black") +
  geom_errorbar(aes(ymin = mean - sem, ymax = mean + sem), width = 0.4) +
  theme_classic() +
  scale_x_discrete(label = "Infected") +
  scale_y_continuous(expand = c(0,0), limits = c(0,1),
                     breaks = seq(0, 1, 0.2)) +
  labs(y = expression("Opn4"^"+"*" mCherry"^"+"*" / Opn4"^"+"), x = element_blank()) +
  theme(axis.text.x = element_blank(),
        axis.ticks.x = element_blank())
