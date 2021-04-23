### R scripts for supporting the findings of article:
### Title: WiBB: an integrated method for quantifying the relative importance of predictive variables
### Authors: Qin Li, Xiaojun Kou
### DOI: (Ecography under review with revision)

# content:
# apply WiBB method to an empirical dataset of 71 Mimulus species
# the goal is to identify important climatic variables in explaining species distributions

library(tidyverse)
library(MuMIn)
library(timeSeries)

#---------------------------------------------------------------------------------------------
# load Mimulus dataset
#---------------------------------------------------------------------------------------------
# cleaned occurrence data and climate variables for 71 Mimulus species
mim_dt = read.csv("empirical_dataset/mimulus_occ_var.csv")

# random sampled background pseudo-absence points with climatic data
backgroup_var = read.csv("empirical_dataset/background_pts_var.csv")

#---------------------------------------------------------------------------------------------
# carry out GLM and calculate WiBB for each species
#---------------------------------------------------------------------------------------------

# set parameters
sp_list = unique(mim_dt$species)
N_sp = length(sp_list) # species number
N_b <- 20 # bootstrap replicates

var_names <- c("T_cold","GDD0","P_season","TP_syn","Aridity","ISO")
pred_x <- paste0("x",1:6)
N_v <- length(var_names)

# create result matrix
wi <- array(data = 0.0, dim = c(N_b,N_v,N_sp), dimnames=list(NULL, pred.x,NULL))

ms.rank <- array(data = NA, dim = c(2^6,29,N_b,N_sp)) # 16,11 = dim(data.frame(ms1))

# start loop
for(i in 1:N_sp){
    
    df_sp <- mim_dt %>% filter(species == sp_list[i])
    df_sp <- df_sp %>% select(one_of(var.names))
    df_sp$y <- 1.0
    
    # intermediate results of sum of weights
    wi_x <- array(data = 0.0, dim = c(N_b,N_v))
    
    for(b.id in 1:N_b){
        
        # ----- read data and prepare for analaysis -----
        
        # two different sub-sampling for species with different sample sizes
        if(nrow(df_sp) <= 800){
            df_sp.sub <- df_sp[sample(1:nrow(df_sp), size=as.integer(nrow(df_sp)*0.75),replace=F),]
        }
        
        if(nrow(df_sp) > 800){ # max number = 640
            df_sp.sub <- df_sp[sample(1:nrow(df_sp), size=as.integer(800*0.75),replace=F),]
        }
        
        # sample 2000 background pseudo-absence points
        df_ev <- backgroup_var[sample(1:nrow(backgroup_var), size=2000,replace=F),]
        df_ev <- df_ev %>% select(one_of(var.names))
        df_ev$y <- 0.0
        
        df_all <- rbind(df_sp.sub, df_ev)
        names(df_all) <- c("x1", "x2", "x3", "x4","x5", "x6", "y")
        for (k in 1:6){df_all[,k] <- scale(df_all[,k])} # scale predictors
        test <- sweep(df_all[1:6],2,colMins(df_all[1:6])-0.01) # non-negative values by minus min
        test$y = df_all$y
        df_all <- test
        
        # --- end data preparition block ----
        
        # --- Start fit glm model ---
        names4col <- c(pred.x, # present/absent coef in each model/row
                       c("adjR^2","df","AICc","rel.LL", "rWi"),
                       paste0("coef_",pred.x), # raw abs(coef) values
                       paste0("rB_",pred.x), # relative/standarized ceof for each model
                       paste0("WiB_",pred.x)) # 
        
        var_ind <- data.frame(matrix(data = 0, nrow = 2^6, ncol = length(names4col)))
        names(var_ind) <- names4col
        
        ### model fitting and model ranking
        glm_fit <- glm(y~., data=df_all, family=binomial(link='logit'), na.action = "na.fail")
        #display(glm_fit)
        ms1 <- dredge(glm_fit, rank = "AICc", extra = "adjR^2")
        
        # matrix of presence or absence of each predictor in each model
        s_sig <- ifelse(is.na(ms1[,pred.x]), 0, 1)
        
        # calculate deltaAIC, rel.LL, and relative weight for each model
        var_ind[,c("adjR^2","df","AICc","rWi")] <- as.data.frame(ms1)[,c("adjR^2","df", "AICc","weight")]
        var_ind[,"rel.LL"] <- exp(-0.5 * as.data.frame(ms1)[,"delta"])
        
        # save raw coefficients
        var_ind[,pred.x] <- s_sig
        coef.raw <- as.matrix(abs(as.data.frame(ms1)[,pred.x]))
        coef.raw[!is.finite(coef.raw)] <- 0
        var_ind[,paste0("coef_",pred.x)] <- coef.raw
        
        # standardize coefficients (sum=1 for each model)
        # calculate the relative contribution of correlation for each variable
        # (by dividing the sum)
        # (before accounting for the weight of the corresponding model)
        # if the sum is too small (< 0.00001), set coefs of x1 ~ x6 = 0
        sumbyrow <- rowSums(coef.raw, na.rm=T)
        sumbyrow <- ifelse(sumbyrow < 0.00001, 0, sumbyrow)
        coef.w <- as.matrix(coef.raw/sumbyrow)
        coef.w[!is.finite(coef.w)] <- 0
        var_ind[,paste0("rB_",pred.x)] <- coef.w
        
        # calculate Wi * B for each predictor: coef.standardized * weight
        var_ind[,paste0("WiB_",pred.x)] <- var_ind[,"rWi"]*var_ind[,paste0("rB_",pred.x)]
        
        ### save calculations
        ms.rank[, , b.id, i] <- as.matrix(var_ind)
        
        ### --- summarize results --- ###
        # sum Wi*B based on all models
        wi_x[b.id, ] <- colSums(var_ind[,paste0("WiB_",pred.x)])
        
        #cat("N_sp =", i, "& boot.id =", b.id ,"is done! \n")
    }
    
    wi[, , i] <- wi_x # wi_x: 400 rows
    
    write.csv(wi_x, file = paste("empirical_dataset/wibb_sp",i,".csv",sep=''),row.names = FALSE)
    saveRDS(wi, file = paste("empirical_dataset/wibb_71sp_all.rds",sep=''))
}


#---------------------------------------------------------------------------------------------
# summary results
#---------------------------------------------------------------------------------------------
# calculate WiBB medians of each variables for each species
mim_wibb_mean <- as.data.frame(matrix(NA, ncol = N_v, nrow = N_sp))
#mim_wibb_med <- as.data.frame(matrix(NA, ncol = N_v, nrow = N_sp))
names(mim_wibb_mean) = var_names
    
for(i in 1:N_sp){
    wi_x <- wi[, , i]
    #mim_wibb_med[i,] <- apply(wi_x, 2, median)
    mim_wibb_mean[i,] <- apply(wi_x, 2, mean)
}

# Table S3
#write.csv(round(mim_wibb_med, 3), file = "result_mim_glm/WiBB_Mim_med.csv",row.names = FALSE)

# summary WiBB mean over species
mim.wi <- mim_wibb_mean %>% mutate(sp=1:nrow(mim_wibb_mean)) %>% 
    melt(id = "sp") %>% 
    group_by(sp, variable) %>% 
    summarise(wibb_mean = mean(value), 
              wibb_sd = sd(value)) %>% 
    mutate(cv = wibb_sd/wibb_mean)

#---------------------------------------------------------------------------------------------
# PLOT
#---------------------------------------------------------------------------------------------
theme_overall <- theme_bw() +
    theme(plot.title = element_text(hjust = 0),
          panel.grid.major = element_blank(),
          panel.grid.minor = element_blank(),
          strip.text = element_text(size=rel(1.1),face = "italic"),
          axis.title = element_text(size = rel(1.2)),
          axis.text.x = element_text(size=rel(1.2)),
          axis.text.y = element_text(size=rel(1.2)),
          legend.title = element_text(size=rel(1.0)),
          legend.key = element_blank())

wibb_ldf <- mim_wibb_mean %>% mutate(sp=1:71) %>% melt(id = "sp")
wibb_ldf$variable = factor(wibb_ldf$variable, levels = var_names,
                           labels = c("T[cold]","GDD0", "P[season]", "TP[syn]", "Aridity", "ISO"))
wibb_var_mean = wibb_ldf %>% group_by(variable) %>% 
    summarise(wibb_mean = mean(value))

p_mim_wibb <- ggplot(wibb_ldf, aes(x = value)) + 
    facet_wrap(~ variable, nrow = 2, labeller = label_parsed) + 
    geom_histogram(alpha=0.4, bins = 20) +
    geom_vline(data = wibb_var_mean, aes(xintercept = wibb_mean), linetype = 2) +
    xlab(expression(paste("Estimates of ", italic("WiBB")))) +
    ylab("Frequency") +
    scale_y_continuous(expand = c(0, 0)) +
    theme_overall

ggsave(filename = "result_mim_glm_summary/mimulus_wibb.pdf",p_mim_wibb,width=12,height=6)

