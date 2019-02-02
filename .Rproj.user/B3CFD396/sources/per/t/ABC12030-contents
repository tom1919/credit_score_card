

LoadPackages <- function(packages) {
  # Load or install packages if they aren't already loaded.
  #
  # Args:
  #   packages: a vector of package names
  #
  for (package in packages) {
    if (!require(package, character.only=T, quietly=T)) {
      if (!package %in% installed.packages()) install.packages(package)
      library(package, character.only=T)
    }
  }
}

summarize_df <- function(df, r = 4) {
  # Create summary information about a data set
  #
  # Args:
  #   df = data frame 
  #   r = number of decimal places to round
  #
  # Return: A data frame containing various summary info about each column
  require(dplyr)
  df <- as.data.frame(df)
  summary_names <- c("col_name", 
                     "type", 
                     "num_unq",
                     "mode", 
                     "mode_ratio",                     
                     "num_missing", 
                     "num_na", 
                     "num_inf", 
                     "num_nan",
                     "min", 
                     "q1", 
                     "median", 
                     "mean", 
                     "q3", 
                     "max", 
                     "std_dev")
  col_summary <- data.frame(matrix(ncol = length(summary_names), nrow = ncol(df)))
  names(col_summary) <- summary_names
  
  for(i in 1:ncol(df)) {
    col <- df[,i]
    not_inf <- df[!is.infinite(col), i]
    freq_table <- sort(table(not_inf), decreasing=TRUE)[1]
    
    col_name <- names(df)[i]
    type <- class(col)
    num_inf <- length(df[is.infinite(col),i])
    num_nan <- length(df[is.nan(col),i])
    num_na <- length(df[is.na(col),i]) - num_nan
    num_missing <- num_na + num_inf + num_nan
    num_unq <- length(unique(not_inf[!is.na(not_inf)])) # NAs and INF values not included
    mode <- names(freq_table) # NAs and INF values not included
    mode <- ifelse(is.null(mode), NA, mode)
    mode_ratio <- unname(freq_table) / (length(col) - num_missing)
    
    
    if(is.numeric(col) == TRUE) {
      min <- min(col, na.rm = TRUE)
      q1 <- quantile(not_inf, .25, na.rm = TRUE) %>% unname()
      median <- median(not_inf, na.rm = TRUE)
      mean <- mean(not_inf, na.rm = TRUE)
      q3 <- quantile(not_inf, .75, na.rm = TRUE) %>% unname()
      max <- max(not_inf, na.rm = TRUE)
      std_dev <- sd(not_inf, na.rm = TRUE)
    } else {
      min <- NA
      q1 <- NA
      median <- NA
      mean <- NA
      q3 <- NA
      max <- NA
      std_dev <- NA
    }
    
    col_summary[i,] <- c(col_name, 
                         type, 
                         num_unq,
                         mode,
                         mode_ratio,
                         num_missing, 
                         num_na,
                         num_inf, 
                         num_nan,
                         min, 
                         q1, 
                         median, 
                         mean, 
                         q3, 
                         max, 
                         std_dev)
  }
  numerics <- dplyr::setdiff(names(col_summary),c("col_name", "type", "mode")) 
  col_summary <- dplyr::mutate_at(col_summary, numerics, funs(as.double(.)))
  col_summary <- dplyr::mutate_at(col_summary, numerics, funs(round(., r)))
  
  return(col_summary)
}

sig_bins <- function(df, col_names, numeric = T, good_col ){
  
  # loop through each col and use smbinning to bin them into categories
  # only create bins if unique vales > 5, there are significant splits and 
  # iv value is greater than 0.1
  # smbinning expcts "good variable"
  result_all_sig <- list() # Creating empty list to store all results #
  
  for(i in 1:length(col_names)){
    if(numeric == T){
      check_res <- smbinning(df = df, y = "good", x = col_names[i])
    }
    else{
      check_res <- smbinning.factor(df = df, y = "good", x = col_names[i])  
    }
    
    if(check_res == "Uniques values < 5") {
      print(paste0(col_names[i], " has less than 5 unq values"))
      next
    }
    else if(check_res == "No significant splits") {
      print(paste0(col_names[i], " has no sig splits"))
      next
    }
    else if(check_res$iv < 0.1) {
      print(paste0(col_names[i], " has iv less than .1"))
      next
    }
    else {
      result_all_sig[[col_names[i]]] <- check_res
      print(paste0(col_names[i], " has splits created"))
    }
  }
  return(result_all_sig)
}

bin_cols <- function(smbinning_result, df, numeric = T) {
  # takes result from smbinning and creates binned variables
  
  if(numeric == T){
    for(i in 1:length(smbinning_result)) {
      df <- smbinning.gen(df = df, ivout = smbinning_result[[i]], 
                          chrname = paste(smbinning_result[[i]]$x, 
                                          "_bin", sep = ""))
    }} else {
      for(i in 1:length(smbinning_result)) {
        df <- smbinning.factor.gen(df = df, ivout = smbinning_result[[i]], 
                                   chrname = paste(smbinning_result[[i]]$x, 
                                                   "_bin", sep = ""))
      }
    }
  return(df)
}

gen_woe <- function(sm_result, df){
  # take result from smbinning and generate WOE values
  for (j in 1:length(sm_result)) {
    for (i in 1:nrow(df)) {
      bin_name <- paste(sm_result[[j]]$x, "_bin", sep = "")
      bin <- substr(df[[bin_name]][i], 2, 2)
      
      woe_name <- paste(sm_result[[j]]$x, "_WOE", sep = "")
      
      if(bin == 0) {
        bin <- dim(sm_result[[j]]$ivtable)[1] - 1
        df[[woe_name]][i] <- sm_result[[j]]$ivtable[bin, "WoE"]
      } else {
        df[[woe_name]][i] <- sm_result[[j]]$ivtable[bin, "WoE"]
      }
    }
  }
  return(df)
}

allocate_points <- function(mod, df, pdo, score, odds) {
  # allocates points for score card
  fact <- pdo/log(2)
  os <- score - fact*log(odds)
  var_names <- names(mod$coefficients[-1])
  
  for(i in var_names) {
    beta <- mod$coefficients[i]
    beta0 <- mod$coefficients["(Intercept)"]
    nvar <- length(var_names)
    WOE_var <- df[[i]]
    points_name <- paste(str_sub(i, end = -4), "points", sep="")
    
    df[[points_name]] <- -(WOE_var*(beta) + (beta0/nvar))*fact + os/nvar
  }
  
  colini <- (ncol(df)-nvar + 1)
  colend <- ncol(df)
  df$Score <- rowSums(df[, colini:colend])
  return(df)
}

draw_confusion_matrix <- function(cm) {
  
  layout(matrix(c(1,1,2)))
  par(mar=c(2,2,2,2))
  plot(c(100, 345), c(300, 450), type = "n", xlab="", ylab="", xaxt='n', yaxt='n')
  title('CLASSIFICATION TABLE', cex.main=2)
  
  # create the matrix 
  rect(150, 430, 240, 370, col='#3F97D0')
  text(195, 440, 'No Default', cex=1.4, font=1) # level 0
  rect(250, 430, 340, 370, col='#F7AD50')
  text(295, 440, 'Defaulted', cex=1.4, font=1) # level 1
  text(125, 370, 'Predicted', cex=1.8, srt=90, font = 2)
  text(245, 450, 'Actual', cex=1.8, font = 2)
  rect(150, 305, 240, 365, col='#F7AD50')
  rect(250, 305, 340, 365, col='#3F97D0')
  text(140, 400, 'No Default', cex=1.4, srt=90, font=1) # level 0
  text(140, 335, 'Defaulted', cex=1.4, srt=90, font=1) # level 1
  
  # add in the cm results 
  res <- as.numeric(cm$table)
  text(195, 400, res[1], cex=1.9, font=2, col='white')
  text(195, 335, res[2], cex=1.9, font=2, col='white')
  text(295, 400, res[3], cex=1.9, font=2, col='white')
  text(295, 335, res[4], cex=1.9, font=2, col='white')
  
  # add in the specifics 
  plot(c(100, 0), c(100, 0), type = "n", xlab="", ylab="", main = 
         "ASSOCIATED STATISTICS", xaxt='n', yaxt='n', cex=1.8, font = 2)
  text(10, 85, names(cm$byClass[1]), cex=1.2, font=3)
  text(10, 70, round(as.numeric(cm$byClass[1]), 3), cex=1.2)
  text(30, 85, names(cm$byClass[2]), cex=1.2, font=3)
  text(30, 70, round(as.numeric(cm$byClass[2]), 3), cex=1.2)
  text(50, 87, names(cm$byClass[5]), cex=1.2, font=3)
  text(50, 71, round(as.numeric(cm$byClass[5]), 3), cex=1.2)
  text(70, 87, names(cm$byClass[6]), cex=1.2, font=3)
  text(70, 71, round(as.numeric(cm$byClass[6]), 3), cex=1.2)
  text(90, 87, names(cm$byClass[7]), cex=1.2, font=3)
  text(90, 71, round(as.numeric(cm$byClass[7]), 3), cex=1.2)
  
  # add in the accuracy information 
  text(30, 35, names(cm$overall[1]), cex=1.5, font=3)
  text(30, 19, round(as.numeric(cm$overall[1]), 3), cex=1.4)
  text(70, 35, names(cm$overall[2]), cex=1.5, font=3)
  text(70, 19, round(as.numeric(cm$overall[2]), 3), cex=1.4)
}  
