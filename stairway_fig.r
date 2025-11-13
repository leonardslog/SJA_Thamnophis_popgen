rm(list=ls())
setwd("~/Documents/sji_thamnophis/stairwayplot/")

# load data files
library(dplyr)
library(purrr)
library(ggplot2)
library(ggh4x)

# sum_files <- list.files("stairway",pattern = ".*_East.*\\.*final.summary*")
ordinoides_files <- list.files("ordinoides/", pattern = ".*final.summary*")
elegans_files <- list.files("elegans/",pattern = ".*final.summary*")
sirtalis_files <- list.files("sirtalis/",pattern = ".*final.summary*")

# name individual results by taxa
ordinoides_results <- lapply(paste0("ordinoides/",ordinoides_files[1]), read.table, header=TRUE)
# ordinoides_taxa <- gsub("_Ne_change.final.summary", "", ordinoides_files)
names(ordinoides_results) <- gsub("_Ne_change.final.summary", "", ordinoides_results) # give population name as name of each element of the list

elegans_results <- lapply(paste0("elegans/",elegans_files[1]), read.table, header=TRUE)
elegans_taxa <- gsub("_Ne_change.final.summary", "", elegans_files)
names(elegans_results) <- gsub("_Ne_change.final.summary", "", elegans_results) # give population name as name of each element of the list

sirtalis_results <- lapply(paste0("sirtalis/",sirtalis_files[1]), read.table, header=TRUE)
sirtalis_taxa <- gsub("_Ne_change.final.summary", "", sirtalis_files)
names(sirtalis_results) <- gsub("_Ne_change.final.summary", "", sirtalis_results) # give population name as name of each element of the list

# test: combine species data
ordinoides_named_results <-  imap(ordinoides_results, ~mutate(.x, Taxon = "ordinoides"))
elegans_named_results <-  imap(elegans_results, ~mutate(.x, Taxon = "elegans"))
sirtalis_named_results <-  imap(sirtalis_results, ~mutate(.x, Taxon = "sirtalis"))

ordinoides_result_df_test <- bind_rows(ordinoides_named_results)
elegans_result_df_test <- bind_rows(elegans_named_results)
sirtalis_result_df_test <- bind_rows(sirtalis_named_results)

# transform axes
ordinoides_result_df_test <- mutate(ordinoides_result_df_test, kya = year/1000)
ordinoides_result_df_test <- mutate(ordinoides_result_df_test, logNe = log(Ne_median))

elegans_result_df_test <- mutate(elegans_result_df_test, kya = year/1000)
elegans_result_df_test <- mutate(elegans_result_df_test, logNe = log(Ne_median))

sirtalis_result_df_test <- mutate(sirtalis_result_df_test, kya = year/1000)
sirtalis_result_df_test <- mutate(sirtalis_result_df_test, logNe = log(Ne_median))

total_result_df_test <- bind_rows(ordinoides_result_df_test,elegans_result_df_test,
                                  sirtalis_result_df_test)
total_result_df_test$Taxon <- factor(total_result_df_test$Taxon) # make the Taxon a factor

# plot population level data for each taxon
cols = c("#d73027", "#1f78b4", "gold")


plot(sirtalis_result_df_test$year/1000,sirtalis_result_df_test$Ne_median/1000, 
     log=c("xy"), type="n", xlim=c(1,50),ylim=c(1,50),
     xlab="Time (1k years ago)", ylab="Effective population size (1k individuals)"
     )

lines(ordinoides_result_df_test$year/1000,ordinoides_result_df_test$Ne_median/1000,type="s",col="#d73027",lwd = 2)
# lines(ordinoides_result_df_test$year/1000,ordinoides_result_df_test$Ne_2.5./1000,type="s",col="darkgrey",lty=3)
# lines(ordinoides_result_df_test$year/1000,ordinoides_result_df_test$Ne_97.5./1000,type="s",col="darkgrey",lty=3)

lines(elegans_result_df_test$year/1000,elegans_result_df_test$Ne_median/1000,type="s",col="gold",lwd = 2)
# lines(elegans_result_df_test$year/1000,elegans_result_df_test$Ne_2.5./1000,type="s",col="darkgrey",lty=3)
# lines(elegans_result_df_test$year/1000,elegans_result_df_test$Ne_97.5./1000,type="s",col="darkgrey",lty=3)

lines(sirtalis_result_df_test$year/1000,sirtalis_result_df_test$Ne_median/1000,type="s",col="#1f78b4",lwd = 2)
# lines(sirtalis_result_df_test$year/1000,sirtalis_result_df_test$Ne_2.5./1000,type="s",col="darkgrey",lty=3)
# lines(sirtalis_result_df_test$year/1000,sirtalis_result_df_test$Ne_97.5./1000,type="s",col="darkgrey",lty=3)

legend(1,3,legend = c("T. ordinoides","T. elegans","T. sirtalis"), col=c("#d73027","gold","#1f78b4"),lty=c(1,1,1,3),lwd=c(4,4,4,1), cex=0.8)

##################################
# ggplot w/ confidence intervals #
##################################

ordinoides_results <- read.table("ordinoides/Tordinoides_Ne_change.final.summary", header=TRUE)
ordinoides_results$Taxon <- "T. ordinoides"
ordinoides_results <- mutate(ordinoides_results, kya = year/1000)

elegans_results <- read.table("elegans/Telegans_Ne_change.final.summary", header=TRUE)
elegans_results$Taxon <- "T. elegans"
elegans_results <- mutate(elegans_results, kya = year/1000)

sirtalis_results <- read.table("sirtalis/Tsirtalis_Ne_change.final.summary", header=TRUE)
sirtalis_results$Taxon <- "T. sirtalis"
sirtalis_results <- mutate(sirtalis_results, kya = year/1000)

A <- rbind(ordinoides_results, elegans_results, sirtalis_results)
cols = c("#d73027", "#1f78b4", "gold")

# Make the plot
ggplot(data=A, aes(x=year, y=Ne_median/1000, ymin=Ne_2.5./1000, ymax=Ne_97.5./1000, 
                   fill=Taxon, linetype=Taxon)) + 
  geom_line() + 
  scale_color_manual(values = c("T. ordinoides" = "#d73027", "T. sirtalis" = "#1f78b4", "T. elegans" = "gold")) +
  scale_fill_manual(values = c("T. ordinoides" = "#d73027", "T. sirtalis" = "#1f78b4", "T. elegans" = "gold")) +
  geom_ribbon(alpha=0.3) + 
  scale_x_log10(expand = c(0, 0)) + 
  scale_y_log10(expand = c(0, 0)) + 
  xlab("Time ago (years)") + 
  ylab("Effective population size Ne (1K individuals)") +
  theme_bw()
