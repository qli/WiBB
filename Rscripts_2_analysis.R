### R scripts for supporting the findings of article:
### Title: WiBB: an integrated method for quantifying the relative importance of predictive variables
### Authors: Qin Li, Xiaojun Kou
### Ecography；DOI: 10.1111/ecog.05651

# content
# 4. fit LM/GLM models and calculate variable importance by four indices: SW, SWi, beta*, WiB, all with bootstrap replicates
# 4.1 LM
# 4.2 GLM
# 5. variable ranking by 4 indices and plot figures


#---------------------------------------------------------------------------------------------
# 4. fit LM/GLM models
#---------------------------------------------------------------------------------------------

library(MuMIn)
library(mvtnorm)
library(parallel)
# library(qpcR) # akaike.weights
# library(fitdistrplus)

### parameters and functions before cl

slev <- c(25,50,100,200,500,1000) # 5 levels
n_set = 200

#---------------------------------------------------------------------------------------------
# 4.1 fit LM models and calculate variable importance by 4 indices
#---------------------------------------------------------------------------------------------
### --- run for lm on cluster with 200 datasets --- ###
# and with 6 sample size (increased to 1000) & bootstrap resamples

# Functions to simulate data sets
# source("~/Documents/PhD_UBC/Research/Rcodes_PhD/Codes_WiBB/01_Func_data_simulatuion.R")

# arrays to save each importance index
# first two dimensions = n_boot * WiB of 4 predictor;
# the third one = independent dataset; the fourth one = sample size;

# LOOP order:
# 1. n_rlev = 3 # r structure
# 2. n_slev = 6 # sample size
# 3. set.index = 200  # datasets
# 4. n_boot = boot.max # bootstrap number

# set.index * n_slev = 200 * 6

# function with sample size as the input
lm.index <- function(n_slev){
    res_list <- list() # result list
    
    sam.size <- slev[n_slev]
    boot.max = 320
    
    ### load data
    sub.sets <- readRDS(paste0("sub_datasets/dataset_lm/dt_r",n_rlev,"_s",slev[i],"_sub.rds"))

    # arrays to save each importance index
    # first two dimensions = n_boot * WiB of 4 predictor;
    # the third one = independent dataset;
    # create list by "parLapply" for 5 sample size
    
    # index WiB
    wi <- array(data = 0.0, dim = c(boot.max,4,n_set),
                dimnames=list(NULL, c("x1","x2","x3","x4"),NULL))
    
    # index beta_star
    b_star = wi
    
    # index standardized SW
    s_wi = wi
    
    # index SW: sum of weight
    sw = wi
    
    # save model parameters & statistics
    ms.rank <- array(data = 0.0, dim = c(16,11,boot.max,n_set)) # 16,11 = dim(data.frame(ms1))
    
    # LOOP order: set.index = 200 * n_boot = 320
    nb.simul = boot.max
    
    # intermediate store of sum of weights
    wi_x <- array(data = 0.0, dim = c(boot.max,4))
    
    for(set.index in 1:n_set){ # 1:200
        # cat("dataset",set.index,"| s_size:",slev[n_slev],"| n_rlev: ", n_rlev, "\n") 
        if(n_slev == 6){
            sub.dt <- sub.sets[[set.index]]
        }else{
            sub.dt <- sub.sets[[set.index]]$dat0
        }
        
        # loop.4 overall boostrap: boot.max in total
        for(m in 1:nb.simul){
            set.seed(813 + m + set.index + n_slev)
            
            # sample from dat0 further by 75%
            dat1 <-  sub.dt[sample(1:nrow(sub.dt),size=as.integer(slev[n_slev]*0.75),replace=F),]
            dat2 <- as.data.frame(scale(dat1)) # scale all variables
            
            ###--------- model averaging ---------### 
            reg <- lm(y ~  x1 + x2 + x3 + x4, data = dat2) # linear regression
            
            # extract coefficients of the full model and save in b_star array
            # if there is negative values, should convert to absolute values
            b_star[m, , set.index] <- abs(round(coef(reg)[2:5],4))
            
            # model fitting all combinations of predictors, rank models by AICc
            ms1 <- dredge(reg, rank = "AICc", extra = "adjR^2")
            
            # save data.frame(ms1)
            ms.rank[, , m, set.index] <- as.matrix(as.data.frame(ms1))
            # (Intercept),x1,x2,x3,x4,adjR^2,df,logLik,AICc,delta,weight
            
            var_ind <- array(data = 0, dim =c(nrow(ms1), 11))
            
            # extract coefficients from models
            coef.raw <- abs(ms1[,c("x1","x2","x3","x4")])
            
            # standardize coefficients (sum=1 for each model)
            # if the sum is too small (< 0.00001), set coefs of x1 ~ x4 = 0
            sumbyrow <- rowSums(coef.raw, na.rm=T)
            sumbyrow <- ifelse(sumbyrow < 0.00001, 0, sumbyrow)
            coef.w <- as.matrix(coef.raw/sumbyrow)
            coef.w[!is.finite(coef.w)] <- 0 # convert NaN/NA/Inf to 0
            var_ind[,1:4] <- coef.w
            
            var_ind[,5] <- ms1$AICc
            min_AICc <- min(var_ind[,5])
            
            # calculate relative likelihood and Akaike weights
            var_ind[,6] <- exp(-0.5*ms1$delta) # = akaike.weights(var_ind[,5])$rel.LL
            var_ind[,7] <- ms1$weight # = akaike.weights(var_ind[,5])$weights
            
            # weighted coefficients
            var_ind[,8:11] <- var_ind[,7] * var_ind[,1:4]
            var_ind <- round(var_ind, digits = 4)
            
            ###
            # estimate sum of weights
            ###
            wi_x[m,] <- colSums(var_ind[, 8:11])
            
            ###
            # SWi: sum of weights (SW != SWi)
            # SWi = SW/sum(SW)
            ###
            
            # matrix of presence or absence of each predictor in each model
            s_sig <- ifelse(is.na(ms1[,c("x1","x2","x3","x4")]), 0, 1)
            
            # for one particular model, all coefs share the same weights (the weight of the model)
            sw.temp <- colSums(var_ind[,7] * s_sig) # calculate SW = (P/A * weight)
            sw[m, , set.index] <- sw.temp
            s_wi[m, , set.index] <- round(sw.temp/sum(sw.temp), 4)
            
        } # finish bootstrap (n_boot)
        
        wi[,, set.index] <- wi_x # wi_x: boot.max rows
    }
    
    # return result list
    res.list$n_slev <- n_slev
    res.list$wi <- wi
    res.list$b_star <- b_star
    res.list$s_wi <- s_wi
    res.list$sw <- sw
    res.list$ms.rank <- ms.rank
    
    return(res.list)
    #saveRDS(res.list, file = paste0("fit_result/r",n_rlev,"_s",n_slev,"_lm_res.rds")) 
}

## ## ## ## ## ## ##
## start parallel ##
## ## ## ## ## ## ##

cat("Start parallel running for lm \n")

### r = 1
n_rlev = 1 # if separate for 3 r, define n_rlev first
print("n_rlev:")
print(n_rlev)

cl <- makeCluster(detectCores() - 3, type="FORK") # on cluster; 6 cores for 6 n_slev
lm_res <- parLapply(cl, 1:length(slev), lm.index)
saveRDS(lm_res, file = paste0("fit_result/r", n_rlev, "_lm.rds"))
print("FINISHED r1!")
stopCluster(cl)

### r = 2
n_rlev = 2 # if separate for 3 r, define n_rlev first
print("n_rlev:")
print(n_rlev)

cl <- makeCluster(detectCores() - 3, type="FORK") # on cluster; need 5 cores for 5 n_slev
lm.re <- parLapply(cl, 1:length(slev), lm.index)
saveRDS(lm.re, file = paste0("fit_result/r", n_rlev, "_lm.rds"))
print("FINISHED r2!")
stopCluster(cl)

### r = 3
n_rlev = 3 # if separate for 3 r, define n_rlev first
print("n_rlev:")
print(n_rlev)

cl <- makeCluster(detectCores() - 3, type="FORK") # on cluster; need 5 cores for 5 n_slev
lm.re <- parLapply(cl, 1:length(slev), lm.index)
saveRDS(lm.re, file = paste0("fit_result/r", n_rlev, "_lm.rds"))
print("FINISHED r3!")
stopCluster(cl)
# end run
cat("Finished! \n")


### one example of model results
n_rlev = 1
n_slev = 6
ex_dt = readRDS(paste0("fit_result/r", n_rlev, "_lm.rds"))
ex_dt = ex_dt[[n_slev]]
saveRDS(ex_dt, file = paste0("fit_result/r",n_rlev,"_s",n_slev,"_lm_res.rds")) 


#---------------------------------------------------------------------------------------------
# 4.2 fit GLM models and calculate variable importance by 4 indices
#---------------------------------------------------------------------------------------------


glm.index <- function(n_slev){
    # result list
    res.list <- list()
    
    sam.size <- slev[n_slev]
    
    ### load data from sub datasets ###
    sub.sets <- readRDS(paste0("sub_datasets/dataset_glm/dt_logit_r",n_rlev,"_s",slev[i],"_sub.rds"))
    
    # index WiB
    wi <- array(data = 0.0, dim = c(400,4,n_set),
                dimnames=list(NULL, c("x1","x2","x3","x4"),NULL))
    
    # index beta_star
    b_star = wi
    
    # index standardized SW
    s_wi = wi
    
    # index SW: sum of weight
    sw = wi
    
    # save model parameters & statistics
    ms.rank <- array(data = NA, dim = c(16,21,boot.max,n_set)) # 16,11 = dim(data.frame(ms1))
    
    # LOOP order: 2 left
    boot.max = 320
    nb.simul = boot.max
    
    # intermediate store of sum of weights
    wi_x <- array(data = 0.0, dim = c(boot.max,4))
    
    for(set.index in 1:n_set){
        # cat("dataset",set.index,"| s_size:",slev[n_slev],"| n_rlev: ", n_rlev, "\n") 
        sub.dt <- sub.sets[[set.index]]
        
        dat1 <- cbind(y=sub.dt$occ, sub.dt[,c("x1","x2","x3","x4")])
        
        # loop.4 overall boostrap: 400 in total, 
        # and then sample from the result based on n_boot?
        for(m in 1:nb.simul){
            set.seed(813 + m + set.index + n_slev)
            
            # sample from dat0 further
            logistic_dat <- dat1[sample(1:nrow(dat1),size=as.integer(slev[n_slev]*0.75),replace=F),]
            logistic_dat[2:5] <- scale(logistic_dat[2:5]) # scale all variables
            
            ### ------------ Start fitting glm model ------------ ###
            
            # to save intedmediate results of Wi
            var_ind <- data.frame(matrix(data = 0, nrow = 16, ncol = 21))
            names(var_ind)[5:21] <- c("adjR^2","df","AICc","rel.LL", "rWi",
                                      "coef_x1","coef_x2","coef_x3","coef_x4",
                                      "rB_x1","rB_x2","rB_x3","rB_x4",
                                      "WiB_x1","WiB_x2","WiB_x3","WiB_x4")
            
            glm_fit <- glm(y~., family=binomial(link='logit'), data=logistic_dat)
            
            # b_star: two options?
            # full model: the last one
            b_star[m, , set.index] <- abs(round(coef(glm_fit)[2:5],4))
            
            ###
            ### model rank
            ms1 <- dredge(glm_fit, rank = "AICc", extra = "adjR^2")
            
            # matrix of presence or absence of each predictor in each model
            s_sig <- ifelse(is.na(ms1[,c("x1","x2","x3","x4")]), 0, 1)
            
            # calculate deltaAIC, rel.LL, and relative weight for each model
            var_ind[,c("adjR^2","df","AICc","rWi")] <- as.data.frame(ms1)[,c("adjR^2","df", "AICc","weight")]
            var_ind[,"rel.LL"] <- exp(-0.5 * as.data.frame(ms1)[,"delta"])
            
            # save raw coefficients
            var_ind[,1:4] <- s_sig
            coef.raw <- as.matrix(abs(as.data.frame(ms1)[,c("x1","x2","x3","x4")]))
            coef.raw[!is.finite(coef.raw)] <- 0
            var_ind[,c("coef_x1","coef_x2","coef_x3","coef_x4")] <- coef.raw
            
            # standardize coefficients (sum=1 for each model)
            # calculate the relative contribution of correlation for each variable
            # (by dividing the sum)
            # (before accounting for the weight of the corresponding model)
            # if the sum is too small (< 0.00001), set coefs of x1 ~ x4 = 0
            sumbyrow <- rowSums(coef.raw, na.rm=T)
            sumbyrow <- ifelse(sumbyrow < 0.00001, 0, sumbyrow)
            coef.w <- as.matrix(coef.raw/sumbyrow)
            coef.w[!is.finite(coef.w)] <- 0
            var_ind[,c("rB_x1","rB_x2","rB_x3","rB_x4")] <- coef.w
            
            # calculate Wi * B for each predictor: coef.standardized * weight
            var_ind[,c("WiB_x1","WiB_x2","WiB_x3","WiB_x4")] <- var_ind[,"rWi"]*var_ind[,c("rB_x1","rB_x2","rB_x3","rB_x4")]
            
            var_ind <- round(var_ind, digits = 6)
            
            ### save calculations
            ms.rank[, , m, set.index] <- as.matrix(var_ind)
            
            ### --- summarize results --- ###
            # sum Wi*B based on all models
            wi_x[m,] <- colSums(var_ind[,c("WiB_x1","WiB_x2","WiB_x3","WiB_x4")])
            
            # SWi: for one particular model, all coefs share the same weights (the weight of the model)
            # var_ind[,1:4] = present/absent of the predictor in the model
            sw.temp <- colSums(var_ind[,"rWi"] * var_ind[,1:4]) # * var_ind[,1:4]
            sw[m, , set.index] <- sw.temp
            s_wi[m, , set.index] <- round(sw.temp/sum(sw.temp), 4)
            
            # SW results with model.ave()
            # confset <- get.models(ms1, subset=TRUE) # 
            # avgmod <- model.avg(confset)
            # summary(avgmod)$importance # parameters' sum of weigth
            
        } # finish 400 simulation (n_boot)
        wi[,, set.index] <- wi_x # wi_x: boot_max rows
    }
    
    # return result list
    res.list$n_slev <- n_slev
    res.list$wi <- wi
    res.list$b_star <- b_star
    res.list$s_wi <- s_wi
    res.list$sw <- sw
    res.list$ms.rank <- ms.rank
    
    return(res.list)
    #saveRDS(res.list, file = paste0("fit_result/r",n_rlev,"_s",n_slev,"_glm_res.rds"))
}

## ## ## ## ## ## ##
## start parallel ##
## ## ## ## ## ## ##

cat("Start parallel running for glm \n")

### r = 1
n_rlev = 1 # if separate for 3 r, define n_rlev first
cl <- makeCluster(detectCores() - 3, type="FORK") # on cluster; 6 cores for 6 n_slev
logit.re <- parLapply(cl, 1:length(slev), glm.index)
saveRDS(logit.re, file = paste0("fit_result/r", n_rlev, ".glm.rds"))
print("FINISHED r1!")
stopCluster(cl)

### r = 2
n_rlev = 2 # if separate for 3 r, define n_rlev first
cl <- makeCluster(detectCores() - 3, type="FORK") # on cluster
logit.re <- parLapply(cl, 1:length(slev), imp.index)
saveRDS(logit.re, file = paste0("fit_result/r", n_rlev, ".glm.rds"))
print("FINISHED r2!")
stopCluster(cl)

### r = 3
n_rlev = 3 # if separate for 3 r, define n_rlev first
cl <- makeCluster(detectCores() - 3, type="FORK") # on cluster
logit.re <- parLapply(cl, 1:length(slev), imp.index)
saveRDS(logit.re, file = paste0("fit_result/r", n_rlev, ".glm.rds"))
print("FINISHED r3!")
stopCluster(cl)

# end run
cat("Finished! \n")

### one example of model results
n_rlev = 1
n_slev = 6
ex_dt = readRDS(paste0("fit_result/r", n_rlev, "_glm.rds"))
ex_dt = ex_dt[[n_slev]]
saveRDS(ex_dt, file = paste0("fit_result/r",n_rlev,"_s",n_slev,"_glm_res.rds")) 

#---------------------------------------------------------------------------------------------
# 5. variable ranking by 4 indices and plot figures
#---------------------------------------------------------------------------------------------

name_ind <- c("wi","b_star","s_wi","sw") # wi = WiBB
mod_2 <- c("lm","glm")

# calculate bootstrap medians for each index
rank_mod = list()

for(mod_id in 1:length(mod_2)){
    name_mod = mod_2[mod_id]
    rank_rlev = list()
    
    for(n_rlev in 1:3){
        res_dt = readRDS(paste0("fit_result/r", n_rlev, "_", name_mod,".rds"))
        
        rank_index = list()
        
        for(ind_id in 1:length(name_ind)){
            index_name = name_ind[ind_id]
            print(paste(name_mod, "r =", n_rlev, index_name))
            
            # calculate median of variable importance for each index with bootstrap replicates
            # include values without bootstrap
            ind_med = median.cal(res_dt=res_dt,index.name=index_name, n_rlev=n_rlev)
            
            # rank variable importance
            rank_res = rank.f(m.array=ind_med$boot_med, n_slev=6, n_boot=6, n_set=200)
            
            # rank variable importance for each index without bootstrap replicates
            rank_b0 = rank.f(m.array=ind_med$boot0, n_slev=6, n_boot=1, n_set=200)
            
            # combine ranking results: freq.tb (x1 > x2 > x3 > x4)
            rank_res$freq.tb$n_boot = row.names(rank_res$freq.tb)
            rank_b0$freq.tb$n_boot = "b0"
            rank_index[[ind_id]] = rbind(rank_b0$freq.tb, rank_res$freq.tb)
            rank_index[[ind_id]]$index_name = index_name
        }
        
        rank_rlev[[n_rlev]] = do.call(rbind, rank_index)
        rank_rlev[[n_rlev]]$rlev = paste0("r",n_rlev)
    }
    
    rank_mod[[mod_id]] = do.call(rbind, rank_rlev)
    rank_mod[[mod_id]]$model = name_mod
}

rank_mod_df = do.call(rbind, rank_mod)
rank_mod_df = rank_mod_df %>% 
    gather(key = "s_size", value = "freq", s1:s6) %>% 
    arrange(model, index_name, rlev, s_size, n_boot)

saveRDS(rank_mod_df, file = "fit_result/variable_importance_rank_res.rds")


#---------------------------------------------------------------------------------------------
### visualization
#---------------------------------------------------------------------------------------------
library(tidyverse)
library(gridExtra)

theme_overall <- theme_bw() +
    theme(plot.title = element_text(hjust = 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          #strip.background = element_blank(),
          strip.text = element_text(size=rel(1.1)),
          axis.title = element_text(size = rel(1.1)),
          axis.text.x = element_text(size=rel(1)),
          axis.text.y = element_text(size=rel(1)),
          legend.title = element_text(size=rel(1.0)),
          legend.key = element_blank(),
          title = element_text(vjust = 2))


### Figure 1. frequency of correct ranking with WiBB
# select three bootstrap replicates for better visualization (to avoid too much overlap)
boot_sel <- c("b0","b1","b3","b5")
rank_wibb = rank_mod_df %>% filter(n_boot %in% boot_sel & index_name =="wi")

p_lm <- ggplot(rank_wibb %>% filter(model == "lm"), aes(x=s_size, y=freq, group=n_boot)) + 
    facet_wrap(~ rlve, labeller=label_bquote(cols = Delta~r:~.(rlve)), nrow=1) + 
    geom_line(linetype = 3,position = position_dodge(width = 0.5)) +
    geom_point(aes(shape = n_boot),size = 2, position = position_dodge(width = 0.5)) +
    scale_shape_manual(name="bootstrap \nreplicate", values = c(1,3,0,4)) +
    xlab("Sample size") +
    ylab("Relative frequency") +
    scale_y_continuous(expand = c(0,0), limits=c(0,1.05)) +
    ggtitle("(a) LM") +
    theme_overall

p_glm <- p_lm %+% (rank_wibb %>% filter(model == "glm")) + ggtitle("(b) GLM")

#ggsave(filename = "Figure1.pdf", grid.arrange(p_lm, p_glm, ncol=1), width = 10, height = 8)

### Figure 2. index comparison with n_boot = 320
rank_index = rank_mod_df %>% filter(n_boot == "b6")

rank_index$index_name <- factor(rank_index$index_name, levels = c("sw","s_wi","b_star","wi"))
labels_index <- list(bquote(italic("SW")), bquote(italic("SWi")), 
                     bquote(paste(italic(beta), "*")), bquote(italic("WiBB")))

p_ind_lm <- ggplot(rank_index %>% filter(model == "lm"), 
                   aes(x = s_size, y = freq, group = index_name)) + 
    facet_wrap(~rlve,labeller=label_bquote(cols = Delta~r:~.(rlve)), nrow=1) + 
    geom_line(linetype = 3,position = position_dodge(width = 0.5)) +
    geom_point(aes(shape = index_name),size = 1.5,
               position = position_dodge(width = 0.5)) +
    scale_shape_manual(name="importance \nindex",
                       labels = labels_index,
                       values = c(1,13,17,19)) +
    xlab("Sample size") +
    ylab("Relative frequency") +
    scale_y_continuous(expand = c(0, 0), limits=c(-0.05,1.05)) +
    ggtitle("(a) LM") +
    theme_overall

p_ind_glm <- p_ind_lm %+% (rank_index %>% filter(model == "glm")) + ggtitle("(b) GLM")

#ggsave(filename = "Figure2.pdf", grid.arrange(p_ind_lm, p_ind_glm, ncol=1), width = 10, height = 8)
