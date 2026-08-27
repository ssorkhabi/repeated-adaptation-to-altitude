og_tau <- read.csv("results/tau/orthogroup_tau_Zscore.csv")
picmin <- read.csv("results/picmin_categorical_results.csv")

sig_z <- og_tau %>% 
  filter(Orthogroup %in% picmin$Orthogroup[picmin$q < 0.01 & picmin$n_est == 11]) %>%
  pull(Z)

all_z <- og_tau$Z

message("Significant OGs mean tau Z: ", round(mean(sig_z, na.rm=TRUE), 3))
message("All OGs mean tau Z: ", round(mean(all_z, na.rm=TRUE), 3))
