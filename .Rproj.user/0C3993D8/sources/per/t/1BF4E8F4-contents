
# required packages
# install.packages("readxl")
# install.packages("writexl")
# install.packages("dplyr")
# install.packages("reshape")
# install.packages("ggpubr")
# install.packages("viridis")
# install.packages("hrbrthemes")
# install.packages("Cairo")
# install.packages("formattable")
# install.packages("plotly")
# install.packages("ggiraph")
# install.packages("Hmisc")


library(readxl)
library(writexl)
library(dplyr)
library(reshape)
library(ggpubr)
library(viridis)
library(hrbrthemes)
library(ggrepel)
# library(Cairo)
library(formattable)
library(plotly)
library(ggiraph)
library(Hmisc)

# read in user-defined functions
source("helper_functions.R")

# read in data
master_csf_df <- 
  readxl::read_excel("./Data/Master CSF Metabolites 2020.xlsx"
                     , sheet = "Master with Survival") %>%
  clean_names() %>%  # convert to lower case letters with underscore
  dplyr::filter(sample_id %in% c("Control", "Immuno", "NoIm")) %>% # get rid of invalid rows below the data in Excel
  na_if("N/A") %>% # replace character "N/A" with NA value
  mutate_at(vars(-c("mass_spect_id", "sample_id"))
            , as.numeric) %>% # convert measurement fields imported as character to 
  as.data.frame() # convert tibble to data.frame

# vector of metabolite variables
metab_vars <- 
  master_csf_df %>% 
  select(-c("mass_spect_id", "sample_id", "os_months", "b2d_months", "censor")) %>% 
  colnames()  


# melt for easy group_by operations on metabolities
master_csf_melted_df <- 
  melt(data = master_csf_df %>% dplyr::select(-c("os_months", "b2d_months", "censor"))
       , id.vars = c("mass_spect_id"
                     ,"sample_id"))


# rebase to control mean value for each metabolite
master_csf_rebased_df <-
  master_csf_df %>%
  mutate_at(metab_vars
            , list(~  ./mean(ifelse(sample_id == "Control", ., NA)
                             , na.rm = TRUE)))

# calculate summary statistics for original values
summary_stats <- 
  master_csf_melted_df %>%
  group_by(variable, sample_id) %>%
  summarise(
    rows = n(),
    valid_values = n() - sum(is.na(value)),
    nas = sum(is.na(value)),
    min = min(value, na.rm = TRUE),
    median = median(value, na.rm = TRUE),
    max = min(value, na.rm = TRUE),
    mean = mean(value, na.rm = TRUE),
    IQR = IQR(value, na.rm = TRUE)
  ) %>% 
  mutate_if(is.numeric, list(~na_if(., Inf))) %>%
  as.data.frame()

# use a consistent NA when a summary stat can't be calculated
is.nan.data.frame <- function(x)
  do.call(cbind, lapply(x, is.nan))

summary_stats[is.nan(summary_stats)] <- NA

write_xlsx(list(original_vals = summary_stats %>% filter(!grepl("_rebase",variable))
                , rebase_values = summary_stats %>% filter(grepl("_rebase",variable)))
           , "./Output/summary_stats.xlsx")

master_csf_melted_df <- master_csf_melted_df %>%
  left_join(summary_stats) 

wilcox_df <-
  merge(data.frame(metabolite = metab_vars)
        , data.frame(group1 = c("Control", "NoIm", "Immuno"))
        , by=NULL) %>%
  merge(data.frame(group2 = c("Control", "NoIm", "Immuno"))
        , by=NULL) %>%
  filter(group1 != group2) %>%
  mutate(p_value = NA, chg_mean = NA, chg_median = NA)


for(i in seq(nrow(wilcox_df))){
  
  if(summary_stats %>% filter(variable == wilcox_df[i, 1] & sample_id %in% c(wilcox_df[i, 2], wilcox_df[i, 3])) %>%
     summarise(min(valid_values)) > 0){
    wilcox_df[i, 4] = wilcox_lookup(a = wilcox_df[i, 2]
                                    , b = wilcox_df[i, 3]
                                    , var = wilcox_df[i, 1]
                                    , df = master_csf_melted_df)
    res <- metab_change(metab = wilcox_df[i, 1]
                        , df = summary_stats
                        , a = wilcox_df[i, 2]
                        , b = wilcox_df[i, 3])
    wilcox_df[i, 5] = res[1]
    wilcox_df[i, 6] = res[2]
    
    
    
  }else
  {wilcox_df[i, 4] = NA
  wilcox_df[i, 5] = NA
  wilcox_df[i, 6] = NA}
  
}


wilcox_df <- 
  wilcox_df %>% 
  filter(group1 == "Control" & group2 == "NoIm" |
           group1 == "Control" & group2 == "Immuno" |
           group1 == "NoIm" & group2 == "Immuno") %>%
  mutate(compare = paste(group2, group1, sep="/")) %>% 
  select(compare, everything(.), -c(group1, group2))



# all
metab_corr <- rcorr(x = master_csf_df %>% 
                      dplyr::select(-c(mass_spect_id:censor)) %>% 
                      as.matrix()
                    , type="pearson")

metab_corr_r_melt_df <- metab_corr$r %>%  
  melt() %>% 
  rename("r" = "value") %>%  
  mutate(X1 = as.character(X1), X2 = as.character(X2))

metab_corr_n_melt_df <- metab_corr$n %>% 
  melt() %>% 
  transmute(n = value)

metab_corr_p_melt_df <- metab_corr$P %>% 
  melt() %>% 
  transmute(p = value) 

metab_corr_melt_df_top <-
  cbind(
    metab_corr_r_melt_df,
    metab_corr_n_melt_df,
    metab_corr_p_melt_df
  ) %>% 
  filter(!is.na(p) & X1 < X2 & n > 5 & p < 0.05) %>% 
  arrange(desc(abs(r))) %>% 
  head(100) %>% 
  mutate(correlation_pair = paste(X1, X2, sep = " & "))

metab_corr_melt_df <- cbind(
  metab_corr_r_melt_df,
  metab_corr_n_melt_df,
  metab_corr_p_melt_df
) %>% 
  mutate(correlation_pair = paste(X1, X2, sep = " & "))

metab_signif_vars <- 
  wilcox_df %>% 
  filter(p_value < 0.05) %>% 
  arrange(desc(abs(chg_median))) %>% 
  distinct(metabolite)


# control

metab_corr_control <- rcorr(x = master_csf_df %>% 
                              dplyr::filter(sample_id == "Control") %>% 
                              dplyr::select(-c(mass_spect_id:censor)) %>% 
                              as.matrix()
                            , type="pearson")

metab_corr_control_r_melt_df <- metab_corr_control$r %>%  
  melt() %>% 
  rename("r" = "value") %>%  
  mutate(X1 = as.character(X1), X2 = as.character(X2))

metab_corr_control_n_melt_df <- metab_corr_control$n %>% 
  melt() %>% 
  transmute(n = value)

metab_corr_control_p_melt_df <- metab_corr_control$P %>% 
  melt() %>% 
  transmute(p = value) 

metab_corr_control_melt_df_top <-
  cbind(
    metab_corr_control_r_melt_df,
    metab_corr_control_n_melt_df,
    metab_corr_control_p_melt_df
  ) %>% 
  filter(!is.na(p) & X1 < X2 & n > 5 & p < 0.05) %>% 
  arrange(desc(abs(r))) %>% 
  head(100) %>% 
  mutate(correlation_pair = paste(X1, X2, sep = " & "))

metab_corr_control_melt_df <- cbind(
  metab_corr_control_r_melt_df,
  metab_corr_control_n_melt_df,
  metab_corr_control_p_melt_df
) %>% 
  mutate(correlation_pair = paste(X1, X2, sep = " & "))



# NoIm

metab_corr_noim <- rcorr(x = master_csf_df %>% 
                              dplyr::filter(sample_id == "NoIm") %>% 
                              dplyr::select(-c(mass_spect_id:censor)) %>% 
                              as.matrix()
                            , type="pearson")

metab_corr_noim_r_melt_df <- metab_corr_noim$r %>%  
  melt() %>% 
  rename("r" = "value") %>%  
  mutate(X1 = as.character(X1), X2 = as.character(X2))

metab_corr_noim_n_melt_df <- metab_corr_noim$n %>% 
  melt() %>% 
  transmute(n = value)

metab_corr_noim_p_melt_df <- metab_corr_noim$P %>% 
  melt() %>% 
  transmute(p = value) 

metab_corr_noim_melt_df_top <-
  cbind(
    metab_corr_noim_r_melt_df,
    metab_corr_noim_n_melt_df,
    metab_corr_noim_p_melt_df
  ) %>% 
  filter(!is.na(p) & X1 < X2 & n > 5 & p < 0.05) %>% 
  arrange(desc(abs(r))) %>% 
  head(100) %>% 
  mutate(correlation_pair = paste(X1, X2, sep = " & "))

metab_corr_noim_melt_df <- cbind(
  metab_corr_noim_r_melt_df,
  metab_corr_noim_n_melt_df,
  metab_corr_noim_p_melt_df
) %>% 
  mutate(correlation_pair = paste(X1, X2, sep = " & "))


# Immuno

metab_corr_immuno <- rcorr(x = master_csf_df %>% 
                           dplyr::filter(sample_id == "Immuno") %>% 
                           dplyr::select(-c(mass_spect_id:censor)) %>% 
                           as.matrix()
                         , type="pearson")

metab_corr_immuno_r_melt_df <- metab_corr_immuno$r %>%  
  melt() %>% 
  rename("r" = "value") %>%  
  mutate(X1 = as.character(X1), X2 = as.character(X2))

metab_corr_immuno_n_melt_df <- metab_corr_immuno$n %>% 
  melt() %>% 
  transmute(n = value)

metab_corr_immuno_p_melt_df <- metab_corr_immuno$P %>% 
  melt() %>% 
  transmute(p = value) 

metab_corr_immuno_melt_df_top <-
  cbind(
    metab_corr_immuno_r_melt_df,
    metab_corr_immuno_n_melt_df,
    metab_corr_immuno_p_melt_df
  ) %>% 
  filter(!is.na(p) & X1 < X2 & n > 5 & p < 0.05) %>% 
  arrange(desc(abs(r))) %>% 
  head(100) %>% 
  mutate(correlation_pair = paste(X1, X2, sep = " & "))

metab_corr_immuno_melt_df <- cbind(
  metab_corr_immuno_r_melt_df,
  metab_corr_immuno_n_melt_df,
  metab_corr_immuno_p_melt_df
) %>% 
  mutate(correlation_pair = paste(X1, X2, sep = " & "))


# combine

metab_corr_melt_df <- metab_corr_melt_df %>% mutate(sample_id = "All")
metab_corr_control_melt_df <- metab_corr_control_melt_df %>% mutate(sample_id = "Control")
metab_corr_noim_melt_df <- metab_corr_noim_melt_df %>% mutate(sample_id = "NoIm")
metab_corr_immuno_melt_df <- metab_corr_immuno_melt_df %>% mutate(sample_id = "Immuno")

metab_corr_combine_melt_df <- rbind(
  metab_corr_melt_df,
  metab_corr_control_melt_df,
  metab_corr_noim_melt_df,
  metab_corr_immuno_melt_df
) %>% arrange(correlation_pair)



metab_corr_melt_df_top <- metab_corr_melt_df_top %>% mutate(sample_id = "All")
metab_corr_control_melt_df_top <- metab_corr_control_melt_df_top %>% mutate(sample_id = "Control")
metab_corr_noim_melt_df_top <- metab_corr_noim_melt_df_top %>% mutate(sample_id = "NoIm")
metab_corr_immuno_melt_df_top <- metab_corr_immuno_melt_df_top %>% mutate(sample_id = "Immuno")

metab_corr_combine_melt_df_top <- rbind(
  metab_corr_melt_df_top,
  metab_corr_control_melt_df_top,
  metab_corr_noim_melt_df_top,
  metab_corr_immuno_melt_df_top
) %>% arrange(correlation_pair)


ggplotly(ggplot(master_csf_df, aes(x=adenine
                                   , y=methylcysteine
                                   , color=sample_id
                                   , text=paste("mass_spect_id:", mass_spect_id)
                                   )) +
  geom_point())




volcano_df <- wilcox_df %>% filter(!is.na(p_value) & compare == "Immuno/Control")

volcano_df$change <- "NO"
volcano_df$change[log2(volcano_df$chg_median) > 0.5 & volcano_df$p_value < 0.05] <- "UP"
volcano_df$change[log2(volcano_df$chg_median) < -0.5 & volcano_df$p_value < 0.05] <- "DOWN"

volcano_df$delabel <- NA
volcano_df$delabel[volcano_df$change != "NO"] <- volcano_df$metabolite[volcano_df$change != "NO"]


# p_volcano <- 
#   ggplot(data=volcano_df %>% arrange(metabolite), aes(x=log2(chg_median), y=-log10(p_value), col=change, label=delabel, tooltip=metabolite)) +
#   geom_point_interactive() +
#   theme_minimal() +
#   geom_text_repel() +
#   scale_color_manual(values=c("blue", "black", "red")) +
#   geom_vline(xintercept=c(-0.5, 0.5), col="gray") +
#   geom_hline(yintercept=-log10(0.05), col="darkgray")
# 
# p_volcano



master_csf_df <- master_csf_df %>% mutate(surv_status = ifelse(censor, 0, 1))

Surv(master_csf_df$os_months, master_csf_df$surv_status)

plot(survfit(Surv(os_months, surv_status) ~ 0, data = master_csf_df), 
     xlab = "Months", 
     ylab = "Overall survival probability")


Surv(master_csf_df$b2d_months, master_csf_df$surv_status)

plot(survfit(Surv(b2d_months, surv_status) ~ 1, data = master_csf_df), 
     xlab = "Months", 
     ylab = "Overall survival probability")



save.image(file = "workspace_for_shiny.RData")
