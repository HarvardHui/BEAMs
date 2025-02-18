#### Assessing the Impact of Batch Effect Associated Missing Values on Functional Analysis in High-Throughput Biomedical Data ####
# Author: Harvard Wai Hann Hui
# Date: 2024-10-18
# Script: Functions


####### Batch effect correction by ComBat----
# pdat should include class and batch factors
do.combat.sim <- function(df, pdat, par=T, plot=F, cov=T){
  df = as.data.frame(df)
  batch <- pdat$batch
  
  if("class" %in% colnames(pdat)){
    if(cov == TRUE){
      class <- pdat$class
      mod = model.matrix(~as.factor(class), data=df)
    } else{
      mod = NULL
    }
  } else{
    mod = NULL
  }
  
  (combat_edata = ComBat(dat=df, batch=batch, mod=mod, par.prior = par, prior.plots = plot))
  
  return(combat_edata)
}


####### Perform t-test based on class factor----
# output either "features" for significant features or "pvals"
# for p-values of all features
t.test.sim <- function(df, class_factor, output="features", fc.thres=0.5){
  unique_classes <- unique(class_factor)
  class1 = df[,which(class_factor==unique_classes[1])]
  class2 = df[,which(class_factor==unique_classes[2])]
    result=c()
    for(i in 1:nrow(df))
    {
      if(sd(class1[i,])<1e-6 && sd(class2[i,])<1e-6)
      {
        class1[i,1]=jitter(class1[i,1])
      }
      a=t.test(class1[i,],class2[i,], var.equal =TRUE)
      
      logfc <- a$estimate[2] - a$estimate[1]
      
      result=rbind(result,c(a$p.value,logfc))

    }
    result[,1]=p.adjust(result[,1],method = "BH")
    result=as.data.frame(result)
    rownames(result)=rownames(class1)
    result = cbind(rownames(result), result)
    colnames(result)=c("feature","pval","logfc")

    
    sig.feats = result$feature[which(result[,'pval']<0.05 & abs(result[,'logfc']) > fc.thres)]
    
    if (output == "features"){
      return(sig.feats)
    }
    if (output == "pvals"){
      out <- result[,'pval']
      names(out) <- result[,'feature']
      return(out)
    }
    if (output == "logfc"){
      out <- result[,'logfc']
      names(out) <- result[,'feature']
      return(out)
    }
}


####### Imputation by KNN----
# Samples in columns, features in rows
do.knn<-function(df, k, weighted =TRUE){
  mv_ind<-which(is.na(df), arr.ind = TRUE)
  mv_rows<-unique(mv_ind[,1])
  
  klist<-list()
  
  dist_mat = as.matrix(dist(t(df), method="euclidean"))
  imp.knn <- df
  imp.knn[is.finite(df) == FALSE] <- NA
  
  for (i in mv_rows){
    klist[[i]]<-list()
    feature_mvs<-which(is.na(df[i,]))
    feature_obs<-which(!is.na(df[i,]))
    # cand_vectors<-df[,feature_obs]
    for (j in 1:length(feature_mvs)){
      feature_mvs_ind<-feature_mvs[j]
      dist <- dist_mat[feature_mvs_ind, feature_obs] # Get distances for every observed sample in feature
      
      if (length(dist) < k) {
        stop(message = "Fewer than K finite distances found")
      } else {
        k_sample_ind <- order(dist)[1:k]
        k_samples <- feature_obs[k_sample_ind]

        wghts <- 1/dist[k_sample_ind]/sum(1/dist[k_sample_ind])
        
        if (weighted == TRUE){
          imp.knn[i, feature_mvs_ind] <- wghts %*% df[i, k_samples]
        } else if (weighted == FALSE){
          imp.knn[i, feature_mvs_ind] <- mean(df[i, k_samples])
        }
      }
      
      klist[[i]][[feature_mvs_ind]]<-k_samples
    }
  }

  return(imp.knn)
}


####### Check for BEAMs in the data----
# if by single batch is set to TRUE, then it will look for severe BEAms
check_beams <- function(xmis, class_factor, batch_factor, by_batch=TRUE, by_single_batch=TRUE){
  #### We require a minimum of 2 batches to define BEAMs
  if (length(unique(batch_factor)) < 2){
    stop("Number of batches less than the minimum required 2!")
  }
  
  xmis = as.matrix(xmis)
  
  if(by_batch == TRUE){
    mv_ind = xmis
    mv_ind[is.finite(mv_ind)] = 1
    mv_ind[is.na(mv_ind)] = 0
    mv_ind_t = data.frame(t(mv_ind), batch=batch_factor)
    sums = aggregate(. ~batch, data=mv_ind_t, FUN=sum)
    sums = sums[,2:ncol(sums)]
    
    # no. of samples in each batch
    batch_lengths = as.numeric(table(batch_factor))
    sums2 = cbind(sums, batch_lengths)
    sums2 = apply(sums2, 1, function(x) x[length(x)] - x[1:(length(x)-1)])
    colnames(sums2) = paste0("Batch_",unique(batch_factor))
    beam_features = rbind(batch_lengths, sums2)
    beam_features = apply(beam_features, 2, function(x) x[2:length(x)] == x[1])
    
    beam_features_inds = apply(beam_features,1,function(x) sum(x) > 0)
    
    if (by_single_batch == TRUE){
      beam_features_inds = apply(beam_features,1,function(x) sum(x) == (length(x) - 1))
    }
    
    return(beam_features_inds)
  }
  
  beam_features = c()
  for (h in 0:(length(unique(class_factor))-1)){
    b = as.factor(batch_factor[class_factor == h])
    
    current_class = xmis[,class_factor == h]
    beam_features_ind = c()
    
    for (j in 1:nrow(current_class)){
      # iterate by each feature
      vec = current_class[j,]
      
      # split feature vector by batches
      c <- vector("list", length(unique(b)))
      
      # split the vector by factor
      for (i in 1:length(unique(b))) {
        c[[i]] <- vec[b == unique(b)[i]]
      }
      
      # count number of MVs in each batch
      MV_count = lapply(c,is.na) %>% lapply(sum)
      beams = which(MV_count == length(c[[1]]))
      if (length(beams) > 0){
        beam_features_ind = c(beam_features_ind, j)
      }
    }
    beam_features = c(beam_features,(rownames(current_class))[beam_features_ind])
  }
  beam_inds = unique(beam_features)

  return(beam_inds)
}


####### Get TPR/FPR values of imputed data----
# true and false cases are determined using the supplied true differential features
tprfpr <- function(imp, class_factor, truediff, output="tpr", fc.thres=0.5){
  imp_true_cases = intersect(rownames(imp), truediff)
  imp_true_cases_inds = match(imp_true_cases, rownames(imp))
  
  # Set positive cases
  tps <- imp_true_cases
  # Set negative cases
  tns <- rownames(imp)[-imp_true_cases_inds]
  
  # Get positive cases from imputed data
  imp.ps <- t.test.sim(imp, class_factor, output="features", fc.thres=fc.thres)
  imp.ps_inds <- match(imp.ps, rownames(imp))
  # Get negative cases from imputed data
  imp.ns <- rownames(imp)[-imp.ps_inds]
  
  
  # True positives
  imp.tps <- intersect(tps, imp.ps)
  # False positives
  imp.fps <- intersect(tns, imp.ps)
  
  # False negatives (non-significant but true positives)
  imp.fn <- intersect(imp.ns, tps)
  # True negatives
  imp.tn <- intersect(imp.ns, tns)
  
  # TPR (aka Recall)
  tpr <- (length(imp.tps)/(length(imp.tps) + (length(imp.fn)))) * 100
  
  # FPR
  fpr <- (length(imp.fps)/(length(imp.fps) + length(imp.tn))) * 100
  
  # Precision
  precision <- (length(imp.tps)/(length(imp.tps) + length(imp.fps))) * 100
  
  # False Discovery Rate
  fdr <- (length(imp.fps)/(length(imp.fps) + length(imp.tps))) * 100
  
  # F-score
  f1 <- 2*((precision * tpr)/(precision + tpr))
  
  if (output == "tpr"){
    return(tpr)
  }
  if (output == "fpr"){
    return(fpr)
  }
  if (output == "precision"){
    return(precision)
  }
  if (output == "fdr"){
    return(fdr)
  }
  if (output == "fscore"){
    return(f1)
  }
  if (output == "confusion"){
    return(list(TP = length(imp.tps),
                FP = length(imp.fps),
                TN = length(imp.tn),
                FN = length(imp.fn)))
  }
}


####### Data simulation function---- 
simulate_microarray <- function(
    m,
    crosstab,
    delta = 1, # Batch Mean
    gamma = 0.5, # Batch SD
    phi = 0.2, # Differential feature proportion
    c = 10,
    d = 6,
    epsilon = 0.5, # limit = (, 1)
    kappa = 0.2, # limit = (, 0.3)
    a = 40,
    b = 5,
    dropout = FALSE,
    r = 2,
    s = -6,
    seed = NA
) {
  # Record parameters
  params <- c(
    crosstab = crosstab,
    delta = delta, gamma = gamma,
    phi = phi, c = c, d = d,
    epsilon = epsilon, kappa = kappa,
    a = a, b = b,
    dropout = dropout, r = r, s = s,
    seed = seed
  )
  
  if (!is.na(seed))
    set.seed(seed)
  
  n <- sum(crosstab)
  n_class <- nrow(crosstab)
  n_batch <- ncol(crosstab)
  class_ids <- 0:(n_class-1)
  gs <- rep(rep(seq_len(n_class), n_batch), crosstab) # class encoding
  ks <- rep(rep(seq_len(n_batch), each = n_class), crosstab) # batch encoding
  
  # Metadata
  gs_alphabet <- class_ids[gs]
  sid <- paste(paste0("ID", seq_len(n)), gs_alphabet, ks, sep = "_")
  metadata <- data.frame(
    class = gs_alphabet,
    batch = as.factor(ks),
    row.names = sid
  )
  
  log_psi <- rgamma(m, a, rate = b)
  # Log fold-change factors for each class
  # Class A has zero log fold change w.r.t. itself
  log_rho <- matrix(0, m, n_class)
  colnames(log_rho) <- LETTERS[seq_len(n_class)]
  diff.features <- NULL
  if (n_class > 1) {
    n_upreg <- n_downreg <- round(phi * m / 2, 0)
    n_diffexpr <- n_upreg + n_downreg
    for (g in seq(2, n_class)) {
      diff.features <- sort(sample(seq_len(m), n_diffexpr))
      upreg.features <- sort(sample(diff.features, n_upreg))
      downreg.features <- setdiff(diff.features, upreg.features)
      for (i in upreg.features) {
        log_rho[i, g] <- rgamma(1, c, rate = d)
      }
      for (i in downreg.features) {
        log_rho[i, g] <- -rgamma(1, c, rate = d)
      }
    }
  }
  
  # Base expression values with class effects
  Z <- matrix(0, m, n)
  colnames(Z) <- sid
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      g <- gs[j]
      Z[i, j] <- rnorm(1, log_psi[i] + log_rho[i, g], epsilon)
    }
  }
  
  # Sample specific scaling term (in log space)
  log_alpha <- rnorm(n, 0, kappa)
  Z <- sweep(Z, 2, log_alpha, `+`)
  
  # Batch effects
  log_beta <- matrix(rnorm(m * n_batch, 0, delta), m, n_batch)
  omega <- matrix(0, m, n)
  for (i in seq_len(m)) {
    for (j in seq_len(n)) {
      k <- ks[j]
      omega[i, j] <- rnorm(1, log_beta[i, k], gamma)
    }
  }
  
  X <- Z + omega
  X[X < 0] <- 0 # set negative values to zero
  
  if (dropout) {
    P <- sigmoid(X, c, d)
    indicator <- matrix(0, m, n)
    for (i in seq_len(m)) {
      for (j in seq_len(n)) {
        indicator[i, j] <- rbinom(1, 1, P[i, j])
      }
    }
    X <- X * indicator
  }
  
  list(
    X = X, metadata = metadata,
    diff.features = diff.features,
    Z = Z, batch.terms = omega,
    class.logfc = log_rho, batch.logfc = log_beta,
    params = params
  )
}


####### Filter features with missing proportions over the set threshold---- 
filter.mvs = function(df, thres){
  mv_amt = apply(df, 1, function(x) length(which(is.na(x))))
  mv_perc = (mv_amt/ncol(df))
  drop_inds = which(mv_perc > thres)
  
  if (length(drop_inds) > 0){
    filtered_df = df[-drop_inds,]
  } else {
    filtered_df = df
  }
  return(filtered_df)
}


####### Class-specific quantile normalization---- 
quantileNormByClass <- function(expr_matrix, class_factor) {
  
  # Check for input validity
  if (!is.matrix(expr_matrix) || !is.factor(class_factor)) {
    stop("expr_matrix must be a matrix and class_factor must be a factor")
  }
  
  if (ncol(expr_matrix) != length(class_factor)) {
    stop("Number of rows in expr_matrix must equal length of class_factor")
  }
  
  # Split the matrix by class
  classes = unique(class_factor)
  split_data = list()
  for (i in classes){
    split_data[[i]] = expr_matrix[,class_factor == i]
  }
  
  # Quantile normalize each split
  norm_data <- lapply(split_data, function(x) normalizeQuantiles(x))
  
  # Combine normalized data back into a matrix, preserving original order
  combined_data <- do.call(cbind, norm_data[names(split_data)])
  
  # Reorder columns to match original order
  combined_data <- combined_data[, colnames(expr_matrix)]
  
  return(combined_data)
}


####### Mean imputation function---- 
mean.imp = function(df){
  global.imp<-df
  for (i in 1:nrow(global.imp)){
    global.imp[i,which(is.na(global.imp[i,]))]<-mean(global.imp[i,],na.rm=TRUE)
  }
  out_mat<-global.imp
  return(out_mat)
}


####### Aggregates peptide expression to protein level---- 
aggregate_protein_expression <- function(peptide_data, feature_data, log=TRUE) {
  # Ensure row names in feature_data match column names in peptide_data
  if (!all(rownames(peptide_data) %in% feature_data$peptide)) {
    stop("Peptide IDs in feature_data do not match column names in peptide_data.")
  }
  peptide_data <- apply(peptide_data, 2, as.numeric)
  
  if (log == FALSE){
    peptide_data = 2^peptide_data
  }
  
  # Merge peptide expression data with feature data to get protein-level expression
  protein_data <- matrix(nrow=length(unique(feature_data[,2])), ncol=0)
  for (i in 1:ncol(peptide_data)){
    helper <- cbind(feature_data, intensity=peptide_data[,i])
    pep_to_prot <- helper %>%
      group_by(protein) %>%
      summarise(total_intensity = sum(intensity, na.rm=TRUE))
    protein_data <- cbind(protein_data, pep_to_prot$total_intensity)
  }
  rownames(protein_data) <- pep_to_prot$protein
  protein_data = log2(protein_data)
  protein_data[protein_data==0] <- NA
  protein_data[is.infinite(protein_data)] <- NA
  return(protein_data)
}


####### Missing value simulation---- 
jin.MNAR <- function(edata, mv_prop, mnar_ratio, beams = FALSE, batch_factor = NULL, seed=NULL, sample.MNAR = F){
  if(!is.null(seed)){
    set.seed(seed)
  }
  Actual_edata = edata
  if(beams == TRUE){
    batches = unique(batch_factor)
    beams_data = c()
    for(i in batches){
      batch_data = apply(edata[,which(batch_factor == i)], 1, mean, na.rm=TRUE)
      beams_data = cbind(beams_data, batch_data)
    }
    colnames(beams_data) = seq(length(batches))
    edata = beams_data
  }
  D = as.matrix(edata)
  n = length(D)
  alpha = mv_prop
  beta = mnar_ratio #ratio relative to mv_prop
  if (sample.MNAR == F){
    tmat = matrix(rnorm(n, mean = quantile(D, alpha, na.rm=TRUE), sd = 0.3), ncol=ncol(D), nrow=nrow(D))
  } else if (sample.MNAR == T){
    tmat <- matrix(NA, ncol=ncol(D), nrow=nrow(D))
    for (col in 1:ncol(tmat)){
      tmat[,col] = rnorm(nrow(tmat), mean = quantile(edata[,col], alpha, na.rm=TRUE), sd = 0.3)
    }
  }
  pmat = matrix(rbinom(n, 1, beta), ncol=ncol(D), nrow=nrow(D))
  
  
  # inject MNAR
  diff_mat = D < tmat
  diff_mat[diff_mat == FALSE] = 0
  diff_mat[diff_mat == TRUE] = 1
  pmat[pmat == FALSE] = 0
  pmat[pmat == TRUE] = 1
  
  indicator = diff_mat + pmat
  
  msdata = D
  # Inject MNAR
  msdata[indicator == 2] = NA
  
  
  if (beams == TRUE){
    msdata = Actual_edata
    for (j in batches){
      msdata[which(indicator[,j] == 2),which(batch_factor == j)] = NA
    }
  }
  
  # Determine MCAR/MAR percentage
  total_mcar = n * alpha * (1-beta)
  remaining_observed = which(is.finite(msdata))
  mcar_inds = sample(remaining_observed, round(total_mcar), replace = FALSE)
  
  # Inject MCAR
  if (beams == TRUE){
    remaining_observed = which(is.finite(indicator))
    mcar_inds = sample(remaining_observed, round(total_mcar), replace = FALSE)
    indicator[mcar_inds] = 3
    for (j in batches){
      msdata[which(indicator[,j] == 3),which(batch_factor == j)] = NA
    }
  } else{
    msdata[mcar_inds] = NA
  }
  
  return(list(mcar_ind = mcar_inds,
              mnar_ind = indicator,
              msdata = msdata))
}

####### Convert DEA outputs into plot data---- 
DEA_to_plot <- function(input, output){
  output_option = c("TPR", "FPR", "Precision", "Fscore","FDR")
  metric = output_option[match(output, output_option)]
  save_res <- lapply(input, function(x) x[[metric]]) %>% do.call(rbind, .)
  plot_dat <- aggregate(save_res$value, by=list(group=save_res$Group, imp=save_res$variable), FUN=mean)
  TPR_sd <- aggregate(save_res$value, by=list(group=save_res$Group, imp=save_res$variable), FUN=sd)
  plot_dat$sd <- TPR_sd$x
  colnames(plot_dat) <- c("group", "x", "mean", "sd")
  labels <- as.character(plot_dat$x)
  labels[labels == "Ground.truth"] <- "Ground truth"
  plot_dat$x <- factor(labels, levels=c("Ground truth","KNN","Mean","MinProb","SVD","MICE","RF"))
  
  plot_dat$data <- case_when(grepl("Control",plot_dat$group) == T ~ "Simulated (control)",
                             grepl("BEAM",plot_dat$group) == T ~ "Simulated (BEAMs)") %>% factor(., levels = c("Simulated (control)","Simulated (BEAMs)"))
  plot_dat$bec <- case_when(grepl("Uncorrected",plot_dat$group) == T ~ "Uncorrected",
                            grepl("Corrected",plot_dat$group) == T ~ "Corrected") %>% factor(., levels = c("Uncorrected","Corrected"))
  
  return(plot_dat)
}


####### MV simulation using Kong et al method---- 
# df = matrix without MVs ; total = total desired MV% ; 
# mcar = % of MVs as MCAR (remaining proportion is simulated as MNAR) 
mv.sim<-function(df,total,mcar, batch_factor = NULL, beams = F, seed=1234){
  set.seed(seed)
  # browser()
  df2<-2^df
  # Indicator matrix
  subdfhole=c()
  
  if(beams == TRUE){
    batches = unique(batch_factor)
    beams_data = c()
    for(i in batches){
      batch_data = apply(df[,which(batch_factor == i)], 1, mean, na.rm=TRUE)
      beams_data = cbind(beams_data, batch_data)
    }
    colnames(beams_data) = seq(length(batches))
  }
  
  
  if (beams == F){
    for(number in 1:ncol(df2)){
      #here generating MNAR
      mis_prop=total-mcar
      data_res = df2[,number]
      a=data_res
      if(mis_prop!=0){
        cutoff = quantile(a, mis_prop,na.rm=T)
        a[a< cutoff] = NA
        data_res=a
      }
      #here generating MCAR
      mis_prop=mcar
      q=round(length(data_res)*mis_prop)
      tmplis=which(!is.na(data_res))
      mi = sample(tmplis, q)
      data_res[mi] = NA
      
      subdfhole=cbind(subdfhole,data_res)
    }}else if (beams == T){
        indicator_mat <- matrix(NA, nrow=nrow(beams_data), ncol=ncol(beams_data))
        for(number in 1:ncol(beams_data)){
          #here generating MNAR
          mis_prop=total-mcar
          data_res = beams_data[,number]
          a=data_res
          if(mis_prop!=0){
            cutoff = quantile(a, mis_prop,na.rm=T)
            a[a< cutoff] = NA
            data_res=a
          }
          #here generating MCAR
          mis_prop=mcar
          q=round(length(data_res)*mis_prop)
          tmplis=which(!is.na(data_res))
          mi = sample(tmplis, q)
          data_res[mi] = NA
          
          # Update indicator matrix
          subdfhole=cbind(subdfhole,data_res)
        }
          # Insert BEAMs
          for (i in 1:length(batches)){
            df[is.na(subdfhole[,i]),which(batch_factor == batches[i])] <- NA
          }

          return(df)
    }
  
  colnames(subdfhole)=colnames(df2)
  
  # jid=c()
  # ### Which rows are completely NA
  # for(j in 1:nrow(subdfhole)){
  #   tmp=subdfhole[j,]
  #   ge=tmp[is.na(tmp)]
  #   if(length(ge)==ncol(subdfhole))
  #     jid=append(jid,j)
  # }
  # if (length(jid > 0)){
  #   subdfhole <- subdfhole[-jid,]
  #   df2 <- df2[-jid,]
  # }
  rownames(subdfhole)<-rownames(df2)
  
  df2<-log(df2,base=2)
  subdfhole<-log(subdfhole,base=2)
  return(subdfhole)
}
mv.sim<-function(df,total,mcar, batch_factor, beams = F, seed=1234){
  set.seed(seed)
  # browser()
  # Indicator matrix
  subdfhole=c()
  

  batches = unique(batch_factor)
  beams_data = c()
  for(i in batches){
    batch_data = apply(df[,which(batch_factor == i)], 1, mean, na.rm=TRUE)
    beams_data = cbind(beams_data, batch_data)
  }
  colnames(beams_data) = seq(length(batches))
  
  
  
  if (beams == F){

    for(number in 1:ncol(df)){
      batch_means <- beams_data[,which(batches == batch_factor[number])]
      #here generating MNAR
      mis_prop=total-mcar
      data_res = df[,number]
      a=data_res
      if(mis_prop!=0){
        cutoff = quantile(batch_means, mis_prop,na.rm=T)
        a[a< cutoff] = NA
        data_res=a
      }
      #here generating MCAR
      mis_prop=mcar
      q=round(length(data_res)*mis_prop)
      tmplis=which(!is.na(data_res))
      mi = sample(tmplis, q)
      data_res[mi] = NA
      
      subdfhole=cbind(subdfhole,data_res)
    }
    colnames(subdfhole)=colnames(df)
    rownames(subdfhole)<-rownames(df)
      
    subdfhole<-log(subdfhole,base=2)
      
    return(subdfhole)
      
    }else if (beams == T){
      indicator_mat <- matrix(NA, nrow=nrow(beams_data), ncol=ncol(beams_data))
      for(number in 1:ncol(beams_data)){
        #here generating MNAR
        mis_prop=total-mcar
        data_res = beams_data[,number]
        a=data_res
        if(mis_prop!=0){
          cutoff = quantile(a, mis_prop,na.rm=T)
          a[a< cutoff] = NA
          data_res=a
        }
        #here generating MCAR
        mis_prop=mcar
        q=round(length(data_res)*mis_prop)
        tmplis=which(!is.na(data_res))
        mi = sample(tmplis, q)
        data_res[mi] = NA
        
        # Update indicator matrix
        subdfhole=cbind(subdfhole,data_res)
      }
      # Insert BEAMs
      for (i in 1:length(batches)){
        df[is.na(subdfhole[,i]),which(batch_factor == batches[i])] <- NA
      }
      
      return(df)
    }

}
