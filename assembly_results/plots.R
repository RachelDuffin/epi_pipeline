# Install and load packages ---------------------------------------------------

install.packages("ggrepel")
install.packages("lubridate")
install.packages("tidyverse")
install.packages("ggpubr")

library(lubridate)
library(tidyverse)
library(ggrepel)
library(grid)
library(ggpubr)

# Load in data ----------------------------------------------------------------

# read in dataset and remove any rows and columns where all values are NA
assembler_results <- 
  read.csv('/home/rachel/msc_project/assembly_results/assembler_testing_results.csv') %>% 
  mutate(Wall_clock_time_h=(period_to_seconds(lubridate::hms(Wall_clock_time_h_m_s)))/60/60, 
         CPU_time_h = CPU_time_s/60/60) %>%
  select(-c(Wall_clock_time_h_m_s, CPU_time_s, Dataset)) %>% 
  drop_na(Threads) %>% discard(~all(is.na(.) | . ==""))

# Functions -------------------------------------------------------------------

lower_ci <- function(mean, se, n, conf_level = 0.95){
  lower_ci <- mean - qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

upper_ci <- function(mean, se, n, conf_level = 0.95){
  upper_ci <- mean + qt(1 - ((1 - conf_level) / 2), n - 1) * se
}

summarise_data <- function(data, col_name){
  data %>% drop_na({{col_name}}) %>% group_by(Tool, Threads, Sample_type) %>%
    summarise(Mean = mean({{col_name}}, na.rm=TRUE), 
              SD = sd({{col_name}}, na.rm=TRUE),  
              Count = n()) %>% 
    mutate(SE = SD/sqrt(Count), 
           Lower_ci = lower_ci(Mean, SE, Count), 
           Upper_ci = upper_ci(Mean, SE, Count))
}

make_plot <- function(data, y_label, x_label, plot_name, width_scale, title) {
  axis_text_size <- 15
  facet_title_size = 20
  
  plot <- data %>% ggplot(aes(Threads, Mean, colour=Sample_type, label=Count)) +
    geom_point() +
    geom_text_repel(color='black') +
    geom_smooth(method = "loess", se = FALSE) +
    labs(y = y_label, x = x_label, title=title) +
    theme_bw() + 
    geom_errorbar(aes(ymin = Lower_ci, ymax = Upper_ci)) +
    theme(legend.position="bottom", 
          legend.text=element_text(size=width_scale),
          legend.title = element_blank(), 
          strip.text.x = element_text(size = facet_title_size), 
          axis.text = element_text(size=12), 
          axis.title = element_text(size=axis_text_size), 
          plot.title = element_text(size=18, face="bold", 
                                    margin = margin(10, 0, 10, 0))) +
    facet_wrap(vars(Tool))
  
  save_plot(plot, plot_name, width_scale)
  return(plot)
}

save_plot <- function(plt, name, width_scale) {
  ggsave(filename = name, 
         plot = plt, width = width_scale, height = 5, units = "in", dpi = 600)
} 

# Data wranging ---------------------------------------------------------------
# Create a tibble per measurement type to show mean, SD, count, SE, lower and
# upper confidence interval

wc_summarised <- assembler_results %>% summarise_data(Wall_clock_time_h)
cpu_summarised <- assembler_results %>% summarise_data(CPU_time_h)
ps_summarised <- assembler_results %>% summarise_data(Parallel_speedup_ratio)
pmu_summarised <- assembler_results %>% summarise_data(Peak_memory_usage_kbytes)

# Plot graphs -----------------------------------------------------------------
# Wall clock time = Time from start to end of program

# CPU time = Total time CPU spends in user and kernel mode (running program 
# code, and executing system calls respectively) (user time + system time)

# Parallel speedup = Ratio of time to run a program with one thread to time to 
# run it with N threads. Manually calculate for each tool using the ratio of 
# wall clock time for the process with 1 thread to the process with multiple 
# threads. Speedup is the ratio between sequential execution time and parallel 
# execution time (i.e. ratio between wall clock time when using one thread to 
# wall clock time when using multiple threads)

# Peak memory usage = Max memory used by a program during its lifetime
# (Maximum resident set size (kbytes))

width_scale <- 12

WC_plot <- make_plot(wc_summarised, 'Wall clock time (h)', 'Threads (n)', 
                     "Wall_clock_plot.png", width_scale, 'Wall clock time')

CPU_plot <- make_plot(cpu_summarised, 'CPU time (h)', 'Threads (n)', 
                      "CPU_time_plot.png", width_scale, 'Total CPU time')

PS_plot <- make_plot(ps_summarised, 'Parallel speedup (ratio)', 'Threads (n)', 
                     "Parallel_speedup_plot.png", width_scale, 
                     'Parallel speedup time')

PMU_plot <- make_plot(pmu_summarised, 'Peak memory usage (kbytes)', 
                      'Threads (n)', "Peak_memory_usage_plot.png", width_scale, 
                      'Peak memory usage')

ggarrange(WC_plot + rremove("xlab"), CPU_plot + rremove("xlab"), 
          PS_plot + rremove("xlab"), PMU_plot, 
          labels = c("A", "B", "C", "D"),
          font.label = list(size = 18, face="bold"),
          common.legend = TRUE, legend = "bottom", ncol = 1)

# Save plots ------------------------------------------------------------------

ggsave(filename = "Assembly_test_plots.png", width = width_scale, 
       height = 15, units = "in", dpi = 600, bg = "white")

