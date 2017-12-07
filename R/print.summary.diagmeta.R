print.summary.diagmeta <- function(x, ...) {


  cat("\nModel: ", model,"\n\n")
  
  
  if(!x$logistic) {
    if (x$log)
      cat("The optimal cut-off value is:",
          round(cutlog,3), "[",
          exp(cut - 1.96 * sqrt(varcut)), ",",
          exp(cut + 1.96 * sqrt(varcut)),"]", "\n\n")
    else
      {cat("The optimal cut-off value is:", round(cut,3), "[", cut-1.96 * sqrt(varcut), ",", cut+1.96 * sqrt(varcut),"]", "\n\n")}
    } else {
      if(log) {cat("The optimal cut-off value is:", round(cutlog,3),"\n\n")} else
      {cat("The optimal cut-off value is:", round(cut,3), "\n\n")}
    }
    
    cat("REML criterion: ", REMLcrit(lmeModel),"\n\n")
    
    if(log)
      cat("Pooled Sensitivity and Specificity at cutoff =", round(exp(evaluateCutoffHere),3), ":", "\n", SESP$Sens, "\n", SESP$Spec,"\n\n")
    else
      cat("Pooled Sensitivity and Specificity at cutoff =", round(evaluateCutoffHere,3), ":", "\n", SESP$Sens, "\n", SESP$Spec,"\n\n")
    
    cat("--------------------------------", "\n\n")
  }
