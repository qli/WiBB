### custom functions of R scripts for supporting the findings of article:
### Title: WiBB: an integrated method for quantifying the relative importance of predictive variables
### Authors: Qin Li, Xiaojun Kou
### Ecography；DOI: 10.1111/ecog.05651

#---------------------------------------------------------------------------------------------
# function to simulate LM dataset
#---------------------------------------------------------------------------------------------
data.simulation <- function(sample.size=1000,corr_y_x1,corr_y_x2,corr_y_x3,precision=0.01){	
    ### correlation matrix
    # corr_y_x1 = expected correlation coefficient between y and x1
    # corr_y_x2 = expected correlation coefficient between y and x2
    # corr_y_x3 = expected correlation coefficient between y and x3
    corr_xi_xj <- 0  # = expected  correlation coefficient between predictor variables
    
    mat <- matrix(cbind(1, corr_y_x1, corr_y_x2, corr_y_x3,
                        corr_y_x1, 1, corr_xi_xj,  corr_xi_xj,  
                        corr_y_x2, corr_xi_xj,  1, corr_xi_xj, 
                        corr_y_x3 ,corr_xi_xj, corr_xi_xj, 1), nrow=4)
    cor1 <- 999
    cor2 <- 999
    cor3 <- 999
    
    # loop repeating data simulations until the requested correlations are reached
    while(abs(cor1 - corr_y_x1) > precision |  abs(cor2 - corr_y_x2) > precision | abs(cor3 - corr_y_x3) > precision ){
        dat <- rmvnorm(n=sample.size, mean=c(0,0,0,0), sigma=mat, method="chol")
        dat <- as.data.frame(dat); names(dat) <- c("y", "x1", "x2", "x3")
        
        y <- dat$y + 5
        x1 <- dat$x1 + 10
        x2 <- dat$x2 + 10
        x3 <- dat$x3 + 10
        cor1 <- cor.test(y,x1)$estimate
        cor2 <- cor.test(y,x2)$estimate
        cor3 <- cor.test(y,x3)$estimate
    }
    x4 <- rnorm(sample.size , 10, 1)	
    dat <- data.frame(y, x1, x2, x3, x4)
    return(dat)
}


#---------------------------------------------------------------------------------------------
# function of sub-sampling
#---------------------------------------------------------------------------------------------
f.subsam <- function(dat, sample.size){
    sub.list <- list()
    min_err = 100.0 # preset a large error/distance
    for(i_loop in 1:100000){
        dat10 <- dat[sample(1:nrow(dat), size=sample.size, replace=F),]
        cor_mat <- cor(dat10, y = NULL)
        the_err <- (cor_mat[1,2]-r_x1)**2 + (cor_mat[1,3]-r_x2)**2 + (cor_mat[1,4]-r_x3)**2 + (cor_mat[1,5])**2
        # bubble sort
        if(the_err < min_err){
            min_err <- the_err
            dat0 <- dat10
            min_mat <- cor_mat
        }
    }
    sub.list$min_err <- min_err # normally it's quite small
    sub.list$dat0 <- dat0
    sub.list$min_mat <- min_mat
    return(sub.list)
}



#---------------------------------------------------------------------------------------------
# function to return bootstrap median
#---------------------------------------------------------------------------------------------
# function to extract values of each index over sample sizes

index.arr <- function(model.res, index.name, s_slev=6, boot_max=320, n_set=200){
    res_array = array(data = NA, dim = c(boot_max, 4, n_set, s_slev))
    for(i in 1:s_slev){
        ind_id = which(names(res_dt[[s_slev]]) == index.name)
        res_array[,,,i] = res_dt[[s_slev]][[ind_id]][1:boot_max,,1:n_set]
    }
    return(res_array)
}

median.cal <- function(res_dt, index.name, n_rlev, n_slev=6, n_boot=6, n_set=200){
    boot_max = 320
    
    result <- list()
    
    # 0. basic parameters
    boot_no <- c(10, 20, 40, 80, 160, 320)
    slev <- c(25,50,100,200,500,1000)
    
    n_boot = length(boot_no)
    
    # 1. load file, depending on defined index
    #res_dt = readRDS(paste0("fit_result/r", n_rlev, "_", name_mod,".rds"))
    
    # ind_arr[boot_n 1:320, x1:x4, set.index 1:200, n_slev 1:6]
    ind_arr <- index.arr(model.res=res_dt, index.name = index.name)
    
    # 2. set empty array to save result
    # 1) median; 2) se;
    
    #boot.mu <- array(data = 0.0, dim =c(n_set,4,n_boot,n_slev)) # set.index, x1:4, n_boot, s_size
    boot_med <- array(data = 0.0, dim =c(n_set,4,n_boot,n_slev)) # set.index, x1:4, n_boot, s_size
    boot_se <- array(data = 0.0, dim =c(n_set,4,n_boot,n_slev))
    
    boot0 <- array(data = 0.0, dim =c(n_set,4,1,n_slev)) # set.index, x1:4, 1, s_size
    
    # 3. loop to calculate mean and se based on boot numbers
    
    for(set.index in 1:n_set){ # no. of data sets
        for(slev.id in 1:n_slev){
            for(boot.id in 1:n_boot){
                
                set.seed(813 + slev.id + boot.id + set.index)
                s_boot <- sample(1:boot_max, size = boot_no[boot.id], replace=F)
                
                # calculate median values over bootstrap replicates
                df.boot <- ind_arr[s_boot, , set.index, slev.id] # WiB with N_boot for x1:x4
                boot_med[set.index, , boot.id, slev.id] <- apply(df.boot, 2, median)
                boot_se[set.index, , boot.id, slev.id] <- apply(df.boot, 2, sd)/boot_no[boot.id]
                
                # extract one row as values without bootstrap
                boot0.id <- sample(1:boot_max, size = 1)
                boot0[set.index, , 1, slev.id] <- ind_arr[boot0.id, , set.index, slev.id]
            }
        }
    }
    
    result$boot_med  <- boot_med 
    result$boot_se <- boot_se
    result$boot0 <- boot0
    
    return(result)
}


#---------------------------------------------------------------------------------------------
#function to count the correct ranking frequency
#---------------------------------------------------------------------------------------------

rank.f <- function(m.array, n_slev=6, n_boot=6, n_set=200){
    
    result <- list()
    
    freq.tb <- as.data.frame(matrix(NA, nrow = n_boot, ncol = n_slev))
    names(freq.tb) <- paste0("s",1:n_slev)
    row.names(freq.tb) <- paste0("b",1:n_boot)
    
    freq.x13 <- freq.tb # with the correct order for x1-x3
    freq.x4 <- freq.tb # x4 is the last one
    
    for(N_slev in 1:n_slev){
        for(boot.id in 1:n_boot){
            
            df.set <- m.array[, , boot.id, N_slev]
            # set.no <- dim(df.set)[1] # default = n_set
            
            rank.T <- 0 # x1 > x2 > x3 > x4
            rank.x13 <- 0 # x1 > x2 > x3, & x4 might be anywhere (include rank.T)
            rank.x4c <- 0 # (x1-x3) > x4 (include rank.T)
            
            for(set.index in 1:n_set){
                
                x.vec <- df.set[set.index,]
                
                if(x.vec[1]>x.vec[2] & x.vec[2]>x.vec[3] & x.vec[3]>x.vec[4]){
                    rank.T = rank.T + 1
                }
                
                if(x.vec[1]>x.vec[2] & x.vec[2]>x.vec[3]){
                    rank.x13 = rank.x13 + 1
                }
                
                if(x.vec[1]>x.vec[4] & x.vec[2]>x.vec[4] & x.vec[3]>x.vec[4]){
                    rank.x4c = rank.x4c + 1
                }
                
            }
            
            # summary freq. over 200 datasets
            freq.tb[boot.id, N_slev] <- rank.T/200
            freq.x13[boot.id, N_slev] <- rank.x13/200
            freq.x4[boot.id, N_slev] <- rank.x4c/200
            
        }
    }
    
    result$freq.tb <- freq.tb
    result$freq.x13 <- freq.x13
    result$freq.x4 <- freq.x4
    
    return(result)
}

