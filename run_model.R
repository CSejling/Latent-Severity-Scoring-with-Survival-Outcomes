setwd("H:/PhD CGM/r_shiny_head_to_head")

library(rstan) # observe startup messages
library(dplyr)

options(mc.cores = parallel::detectCores())

rstan_options(auto_write = TRUE)

vary_rep <- function(vec, each){
  
  out <- numeric(0)
  
  k <- 1
  
  for (i in vec){
    
    out <- c(out, rep(i, each[k]))
    
    k <- k + 1
    
  }
  
  return(out)
  
}

# Data preparation:
pl <- read.csv2("PL_2023_24.csv", header = TRUE)
pl$W <- as.numeric(pl$Home.Goal > pl$Away.Goal)
pl$W[pl$Home.Goal == pl$Away.Goal] <- 3

pl <- subset(pl, Home.Team != "")

#translater_teams <- data.frame("Team"=unique(pl$Home.Team), "Index"=order(unique(pl$Home.Team)))
#translater_officials <- data.frame("Team"=unique(pl$Official), "Index"=order(unique(pl$Official)))

pl$og_hometeams <- pl$Home.Team
pl$og_awayteams <- pl$Away.Team

M <- n_distinct(pl$Home.Team)
pl$Home.Team[order(pl$Home.Team)] <- rep(1:M, each = 19)
pl$Home.Team <- as.numeric(pl$Home.Team)
pl$Away.Team[order(pl$Away.Team)] <- rep(1:M, each = 19)
pl$Away.Team <- as.numeric(pl$Away.Team)

K <- n_distinct(pl$Official)
pl$og_official <- pl$Official
pl$Official[order(pl$Official)] <- vary_rep(1:K, each = table(pl$Official[order(pl$Official)]))
pl$Offical <- as.numeric(pl$Official)

D <- sum(pl$W==1) + sum(pl$W==0)  
U <- sum(pl$W==3)
yd <- pl$W[pl$W!=3]
rd <- pl$Official[pl$W!=3]
ru <- pl$Official[pl$W==3]
ad <- pl$Home.Team[pl$W!=3]
au <- pl$Home.Team[pl$W==3]
bd <- pl$Away.Team[pl$W!=3]
bu <- pl$Away.Team[pl$W==3]
futd <- rep(1,D)
futu <- rep(1,U)

pl_dat <- list(M = M, K = K, D = D, U = U, yd = yd, rd = rd, ru = ru, ad = ad, au = au, bd = bd, bu = bu, futd = futd, futu = futu)

for (i in 1:length(pl_dat)){
  
  pl_dat[[i]] <- as.numeric(pl_dat[[i]])
  
}

set.seed(2101671)

# Model fit: 
fit <- stan(file = 'head_to_head.stan', data = pl_dat, chains = 4, iter = 5000)

library(stringr)

print(fit)
plot(fit, pars = c(str_c("lambda[",1:M,"]"), str_c("eta[",1:K,"]")))
traceplot(fit, pars = c(str_c("lambda[",1:M,"]"), str_c("eta[",1:K,"]")))

summary(summary(fit)$summary[, "n_eff"])

placings <- order(summary(fit)$summary[1:20,1], decreasing = FALSE)
placings_ref <- order(summary(fit)$summary[21:50,1], decreasing = FALSE)

library(bayesplot)
library(ggplot2)

x <- as.matrix(fit, pars = str_c("lambda[",1:M,"]"))  
#mcmc_intervals(x, pars = str_c("lambda[",placings,"]")) + 
#  scale_y_discrete(labels = left_join(data.frame("Index"=placings), unique(data.frame("Index"=pl$Home.Team,"Team" = pl$og_hometeams)), by = "Index")$Team) # my_labels is character vector 
mcmc_areas(x, pars = str_c("lambda[",placings,"]")) + scale_y_discrete(labels = left_join(data.frame("Index"=placings), unique(data.frame("Index"=pl$Home.Team,"Team" = pl$og_hometeams)), by = "Index")$Team) # my_labels is character vector 

x_ref <- as.matrix(fit, pars = str_c("eta[",1:K,"]")) 
mcmc_areas(x_ref, pars = str_c("eta[",placings_ref,"]")) + scale_y_discrete(labels = left_join(data.frame("Index"=placings_ref), unique(data.frame("Index"=as.numeric(pl$Official),"Official" = pl$og_official)), by = "Index")$Official) # my_labels is character vector 

#P(lambda_arsenal > lambda_mancity)
mean(x[,1]>x[,13])

# Posterior predictive sampling of point tallies

point_tally <- matrix(nrow=nrow(x),ncol=ncol(x))

for (i in 1:nrow(x)){
  
  log_power_levels <- x[i,]
  ref_effects <- x_ref[i,]
  
  pl$power_level_home <- exp(log_power_levels[pl$Home.Team] + ref_effects[pl$Offical])
  pl$power_level_away <- exp(log_power_levels[pl$Away.Team] + ref_effects[pl$Offical])
  
  sprob <- exp(-pl$power_level_home-pl$power_level_away)
  S <- rbinom(n = nrow(pl), size = 1, prob = sprob)
  
  pl$result <- rep(NA, nrow(pl))
  pl$result[S==1] <- 3
  
  wprob <- pl$power_level_home / (pl$power_level_home+pl$power_level_away)
  pl$result[S==0] <- rbinom(n = sum(1-S), size = 1, prob = wprob[S==0])
  
  pl$home_points <- S + 3*(pl$result==1)
  pl$away_points <- S + 3*(pl$result==0)
  
  for (j in 1:ncol(x)){
  
  point_tally[i,j] <- sum(pl$home_points[pl$Home.Team==j]) + sum(pl$away_points[pl$Away.Team==j])
    
  }
  
}

colnames(point_tally) <- str_c("Team",as.character(1:20))

mcmc_areas(point_tally, pars = str_c("Team",placings)) + scale_y_discrete(labels = left_join(data.frame("Index"=placings), unique(data.frame("Index"=pl$Home.Team,"Team" = pl$og_hometeams)), by = "Index")$Team) # my_labels is character vector 

summary(point_tally)

## Coverage of lambda parameters:

set.seed(5106912)

L <- 100

cov_frame <- data.frame("2.5" = rep(NA,M*S), "25" = rep(NA,M*S), "75" = rep(NA,M*S), "97.5" = rep(NA,M*S), "cov_50" = rep(NA,M*S), "cov_95" = rep(NA,M*S))

si <- sample(1:nrow(x), size = L, replace = FALSE)

k <- 1

for (i in si){
  
  log_power_levels <- x[i,]
  ref_effects <- x_ref[i,]
  
  pl$power_level_home <- exp(log_power_levels[pl$Home.Team] + ref_effects[pl$Offical])
  pl$power_level_away <- exp(log_power_levels[pl$Away.Team] + ref_effects[pl$Offical])
  
  sprob <- exp(-pl$power_level_home-pl$power_level_away)
  S <- rbinom(n = nrow(pl), size = 1, prob = sprob)
  
  pl$result <- rep(NA, nrow(pl))
  pl$result[S==1] <- 3
  
  wprob <- pl$power_level_home / (pl$power_level_home+pl$power_level_away)
  pl$result[S==0] <- rbinom(n = sum(1-S), size = 1, prob = wprob[S==0])
  
  # Artifical data:
  
  D <- sum(pl$result==1) + sum(pl$result==0)  
  U <- sum(pl$result==3)
  yd <- pl$result[pl$result!=3]
  rd <- pl$Official[pl$result!=3]
  ru <- pl$Official[pl$result==3]
  ad <- pl$Home.Team[pl$result!=3]
  au <- pl$Home.Team[pl$result==3]
  bd <- pl$Away.Team[pl$result!=3]
  bu <- pl$Away.Team[pl$result==3]
  futd <- rep(1,D)
  futu <- rep(1,U)
  
  pl_dat_test <- list(M = M, K = K, D = D, U = U, yd = yd, rd = rd, ru = ru, ad = ad, au = au, bd = bd, bu = bu, futd = futd, futu = futu)
  
  for (i in 1:length(pl_dat_test)){
    
    pl_dat_test[[i]] <- as.numeric(pl_dat_test[[i]])
    
  }
  
  #
  
  fit_test <- stan(file = 'head_to_head.stan', data = pl_dat_test, chains = 4, iter = 2000)
  
  cov_frame[((k-1)*M + 1):(k*M),1:4] <- summary(fit_test)$summary[1:M, c("2.5%", "25%", "75%", "97.5%")]
  
  cov_frame$cov_50[((k-1)*M + 1):(k*M)] <- as.numeric((log_power_levels > cov_frame$X25[((k-1)*M + 1):(k*M)]) & (log_power_levels < cov_frame$X75[((k-1)*M + 1):(k*M)]))
  
  cov_frame$cov_95[((k-1)*M + 1):(k*M)] <- as.numeric((log_power_levels > cov_frame$X2.5[((k-1)*M + 1):(k*M)]) & (log_power_levels < cov_frame$X97.5[((k-1)*M + 1):(k*M)]))
  
  k <- k + 1
  
}

# Gathering up results:

team_coverage_50 <- numeric(0)
team_coverage_95 <- numeric(0)

for (team in 1:20){
  
  team_coverage_50[team] <- mean(cov_frame$cov_50[team + M*(0:(L-1))])
  
  team_coverage_95[team] <- mean(cov_frame$cov_95[team + M*(0:(L-1))])
  
}

summary(team_coverage_50)
summary(team_coverage_95)

### return an array of three dimensions: iterations, chains, parameters 
a <- extract(fit, permuted = FALSE) 

summary(c(a[,,1]) - c(a[,,13]))

### use S3 functions on stanfit objects
a2 <- as.array(fit)
m <- as.matrix(fit)
d <- as.data.frame(fit)