######################### Undirected Unequal Variance ######################### 
rm(list = ls())
library(BayesFactor)
library(rethinking)
library(progress)

# Notes: 
# Cohens's D = 2 * (t)/sqrt(N)
# pvalue from t-stats: 2*pt(t_stats, N, lower=FALSE) # two-tailed
# one-tailed 2*pt(-2,500, lower=TRUE)
set.seed(123)
n_sim <- 1000
n1 <- c(500, 1000)
n2 <- c(500, 1000)
diff_years <- seq(1, 5, length.out = 17)
twoSided <- TRUE

df <- expand.grid(n1, n2, diff_years)
colnames(df) <- c("N1", "N2", "diff_years")
df$median_t_stats <- NA
df$frequency_detected_t <- NA
df$median_cohens_d <- NA
df$median_BF <- NA
df$frequency_detected_BF <- NA

pb <- progress::progress_bar$new(total = nrow(df))
for (i in 1:nrow(df)){
  pop_sd1 <- 10
  pop_sd2 <- 15
  se <- sqrt(((pop_sd1^2) / df[i, "N1"]) + ((pop_sd2^2) / df[i, "N2"]))
  # Here we test for a undirected difference
  peak_chess_mean <- rnorm(n_sim, 30, se) # sample mean of individual age at peak skill in chess
  peak_go_mean <- rnorm(n_sim, 30 + df[i, "diff_years"], se) # sample mean of individual age at peak skill in go
  t_statistics <- (peak_go_mean-peak_chess_mean) / se
  df$median_t_stats[i] <- median(t_statistics)
  cohens_d <- (peak_go_mean-peak_chess_mean) / sqrt((pop_sd1^2 + pop_sd2^2)/2)
  df$median_cohens_d[i] <- median(cohens_d)
  
  p <- 0.001
  if(twoSided){
    p <- p/2
  }
  t_cutoff <- qt(p, df = df[i, "N1"] + df[i, "N2"] - 2, lower.tail = FALSE)
  
  # df$frequency_detected_t[i] <- ifelse(test = is.na(table(t_statistics < t_cutoff)["TRUE"]),
  # yes = 0,
  # no = table(t_statistics < t_cutoff)["TRUE"]/ length(t_statistics)) # power to detect difference at 0.001 level (approx)
  if(twoSided) {
    df$frequency_detected_t[i] <- mean(t_statistics > t_cutoff | t_statistics < -t_cutoff, na.rm = TRUE)
  } else {
    df$frequency_detected_t[i] <- mean(t_statistics > t_cutoff, na.rm = TRUE)
  }
  
  if(twoSided) {
    bf <- sapply(t_statistics, function(t) exp(BayesFactor::ttest.tstat(t, df[i, "N1"], df[i, "N2"],
                                                                        rscale = 0.1)$bf))
  } else {
    bf <- sapply(t_statistics, function(t) exp(BayesFactor::ttest.tstat(t, df[i, "N1"], df[i, "N2"],
                                                                        rscale = 0.1,
                                                                        nullInterval = c(0, Inf))$bf))
  }
  df$median_BF[i] <- median(bf)
  # Only Bayes Factors above 5
  df$frequency_detected_BF[i] <- ifelse(test = is.na(table(bf > 10)["TRUE"]),
                                        yes = 0,
                                        no =  table(bf > 10)["TRUE"]/length(bf))
  
  pb$tick()
}
save(df, file = "output/BF_power_unequal_undirected.RData")

# png("output/bf_power_unequal_undirected.png", width = 10, height = 5, units = "in", res = 800)
# postscript("output/bf_power_unequal_undirected.eps", onefile = TRUE, horizontal = TRUE, width = 10, height = 5)
cairo_ps("output/bf_power_unequal_undirected.eps", width = 10, height = 5)
plot_power <- function(age_diff, power_freq, power_bayes) {
  #browser()
  plot(age_diff, power_freq, type = "n", xaxt='n',
       ylim =c(0,1), xlab = "", ylab = "", cex.axis = 1.2,
       bty = "l")
  lines(x = age_diff, y = power_freq, col = "#BD0026", lwd = 2)
  lines(x = age_diff, y = power_bayes, col = "#4575B4", lwd = 2)
  points(x = age_diff, y = power_freq, col = "#BD0026", lwd = 2, pch = 16)
  points(x = age_diff, y = power_bayes, col = "#4575B4", lwd = 2, pch = 16)
  mtext(expression(delta~bar(phi)), side = 1, line = 3)
  axis(side = 1, at = pretty(age_diff), labels = pretty(age_diff), cex.axis = 1.2)
  mtext("Probability of Detection", side = 2, line = 2.2)
  #axis(side = 2, at = pretty(seq(0, 1, by = 0.2)),
  #     labels = seq(0, 1, by = 0.2), cex.axis = 1.2)
  mtext("Cohen's d", side = 3, line = 2.2)
  axis(side = 3, at = pretty(1:5), labels = round(pretty(age_diff)/sqrt((pop_sd1^2 + pop_sd2^2)/2), 2), cex.axis = 1.2)
  rect(xleft = 0, ybottom = 0.80,
       xright = 6, ytop = 1,
       col = col.alpha("#D95F02", 0.2), border = NA)
}

par(mfrow = c(1, 3), mar=c(4,4,6,1), oma=c(0,0,0,0))
with(subset(df, N1 == 500  & N2 == 500 ), plot_power(diff_years, frequency_detected_t, frequency_detected_BF))
mtext("A", side = 3, line = 4, adj = 0.5, font = 2)
with(subset(df, N1 == 1000 & N2 == 500 ), plot_power(diff_years, frequency_detected_t, frequency_detected_BF))
mtext("B", side = 3, line = 4, adj = 0.5, font = 2)
with(subset(df, N1 == 1000 & N2 == 1000), plot_power(diff_years, frequency_detected_t, frequency_detected_BF))
mtext("C", side = 3, line = 4, adj = 0.5, font = 2)
dev.off()

df[which(df$N1 == 500 & df$N2 == 500), c("diff_years", "frequency_detected_t", "frequency_detected_BF")]
df[which(df$N1 == 500 & df$N2 == 1000), c("diff_years", "frequency_detected_t", "frequency_detected_BF")]
df[which(df$N1 == 1000 & df$N2 == 1000), c("diff_years", "frequency_detected_t", "frequency_detected_BF")]
