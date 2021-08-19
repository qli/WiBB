### R scripts for supporting the findings of article:
### Title: WiBB: an integrated method for quantifying the relative importance of predictive variables
### Authors: Qin Li, Xiaojun Kou
### Ecography；DOI: 10.1111/ecog.05651

# content
# 1. LM dataset simulation: 200 datasets
# 2. data resampling (by svarying ample size)
# 3. GLM dataset simulation


#---------------------------------------------------------------------------------------------
# 1. LM dataset simulation with three levels of delta_r
#---------------------------------------------------------------------------------------------
# three levels of differences of correlation coefficients of consecutive predictors
# delta_r = c(0.1, 0.2, 0.3), meaning three sets of coefficients:

library(mvtnorm)
library(parallel)
#----------------------------
# parameter settings
r1 = c(0.3, 0.2, 0.1, 0.0)
r2 = c(0.6, 0.4, 0.2, 0.0)
r3 = c(0.8, 0.6, 0.3, 0.0) # replace c(0.9, 0.6, 0.3, 0.0)
delta_r = list(r1,r2,r3)

for(i in 1:length(delta_r)){
    r_x1 <- delta_r[[i]][1]
    r_x2 <- delta_r[[i]][2]
    r_x3 <- delta_r[[i]][3]
    
    raw_dat <- list()
    simu.size <- rep(1000, times = 200) # 200 datasets, each with 1000 observations
    
    simu.set <- function(s_size){
        # use data.simulation() to generate dataset with preset sample size & delta_r level
        data.simulation(sample.size = s_size, corr_y_x1 = r_x1, corr_y_x2 = r_x2, corr_y_x3 = r_x3)
    }
    
    cl <- makeCluster(detectCores()-2, type="FORK")
    tryCatch(
        raw_dat <- parLapply(cl, simu.size, simu.set),
        error = function(e) print(e)
    )
    save(raw_dat, file = paste0("dt_lm_r",i,"rds"))
    stopCluster(cl)
}


# assessment of the collinearity between predictors
library(HH)
vif_df = list()
for(n_rlev in 1:3){
    raw_dt = readRDS(paste0("data/data_manu/full_datasets/dt_lm_r",n_rlev,".rds"))
    vif_df[[n_rlev]] = data.frame(n_rlev = rep(1,200), x1=NA,x2=NA,x3=NA,x4=NA)
        
    for(dt_id in 1:200){
        #mat_collinearity = model.matrix(~ x1 + x2 + x3 + x4, data = raw_dt[[dt_id]])
        # Correlation values (off-diagonal elements) of at least .4 are sometimes interpreted as indicating a multicollinearity problem.
        #cov(mat_collinearity[,-1])
        
        # Using the fact that the determinant of a diagonal matrix is the product of the eigenvalues => The presence of one or more small eigenvalues indicates collinearity
        #eigen(cov(mat_collinearity[,-1]))$values
        
        # variance inflation factor : a VIF of 5 or 10 and above indicates a multicollinearity problem (but see O'Brien 2007).
        # install.packages('HH')
        
        vif_df[[n_rlev]][dt_id,2:5] = vif(raw_dt[[dt_id]][,-1])
    }
}    

vif_df[[1]][,2:5] %>% range() # 1.000026 1.023857
vif_df[[2]][,2:5] %>% range() # 1.000024 1.018540
vif_df[[3]][,2:5] %>% range() # 1.000026 1.023857



#---------------------------------------------------------------------------------------------
# 2. data subsampling (by varying sample sizes) for LM dataset
#---------------------------------------------------------------------------------------------
slev <- c(25,50,100,200,500,1000)

# start parallel #
set.seq <- c(1:200) # 200 datasets
for(n_rlev in 1:3){
    raw_dt = readRDS(paste0("full_datasets/dt_lm_r",n_rlev,".rds"))
    for(i in 1:length(slev)){
        if(slev[i] == 1000){
            saveRDS(raw_dt, file = paste0("sub_datasets/dt_r",n_rlev,"_s", slev[i], "_sub.rds"))
        }else{ # apply function f.subsam() for subsample
            sub.sample.f <- function(set.id){
                dat <- raw_dt[[set.id]]
                sub_dt <- f.subsam(dat = dat, sample.size = slev[i])
                return(sub_dt)
            }
            cl <- makeCluster(detectCores()-2, type="FORK") # or other proper number of cores
            tryCatch(
                sub_set <- parLapply(cl, set.seq, sub.sample.f),
                error = function(e) print(e)
            )
            saveRDS(sub_set,file=paste0("sub_datasets/dataset_lm/dt_r",n_rlev,"_s",slev[i],"_sub.rds"))
            stopCluster(cl)
        }
    }
}


#---------------------------------------------------------------------------------------------
# 3. generate 1/0 binary datasets based on LM datasets
#---------------------------------------------------------------------------------------------
# generate GLM dataset from LM dataset
# 1) within each original LM dataset (n=1000), generate y3 & probability & generate binary response value (occ)
# 2) for subset dataset, left_join with above occ as resampling is done

for(n_rlev in 1:3){
    # 1) within each raw lm dataset (n=1000), generate y3 & p & occ
    raw_dt = readRDS(paste0("full_datasets/dt_lm_r",i,"rds"))
    logit_dt <- list() # to save logit dataset
    
    for(set.index in 1:200){
        #cat("dataset",set.index,"| s_size:",slev[n_slev],"| n_rlev: ", n_rlev, "\n") 
        lm_dt <- raw_dt[[set.index]]
        y <- scale(lm_dt$y)
        y3 <- y * 100
        p3 <- exp(y3) / (1 + exp(y3))
        
        # generate an uniform distribution and set all p3 > occ_thred to be 1, else be 0;
        occ_thred <- runif(1000, min=0, max=1)
        occ <- (p3 > occ_thred) 
        occ <- as.integer(occ)
        logit_dt[[set.index]] <- cbind(occ, raw.dt)
    }
    saveRDS(logit_dt, paste0("full_datasets/dt_logit_r",n_rlev,".rds"))
}

slev <- c(25,50,100,200,500,1000)

for(n_rlev in 1:3){
    logit_dt <- readRDS(paste0("full_datasets/dt_logit_r",n_rlev,".rds"))
    
    for(n_slev in 1:length(slev)){
        # 2) for subset dataset, left_join with above occ
        sub_sets <- readRDS(paste0("sub_datasets/dt_r",n_rlev,"_s", slev[i], "_sub.rds"))
        logit_dt_sub <- list()
        for(set.index in 1:1000){
            sub_dt <- sub_sets[[set.index]]$dat0
            logit_dt_sub[[set.index]] <- sub_dt %>% left_join(logit_dt[[set.index]])
        }
        saveRDS(logit_dt_sub, paste0("sub_datasets/dt_logit_r",n_rlev,"_s",slev[n_slev],"_sub.rds"))
    }
}

