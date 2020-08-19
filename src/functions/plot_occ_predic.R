plot_occ_predic <- function(cov,covRange,coeffname1,coeffname2,intname,xlab,ylab,ylim) {
  
  vec <- which(!is.na(cov))
  vec <- cov[vec]
  
  mean <- attr(scale(vec), "scaled:center") # mean 
  sd <- attr(scale(vec), "scaled:scale")  # standard deviation
  orig.pred <- seq(from = covRange[1], to = covRange[2], by = covRange[3])# Predictor (x-values) which are unscaled 
  sc.pred <- (orig.pred - mean)/sd 
  
  coef1samples <- c(ttd_modfit$samples[[1]][,coeffname1],
                    ttd_modfit$samples[[2]][,coeffname1],
                    ttd_modfit$samples[[3]][,coeffname1])
  
  intsamples <- c(ttd_modfit$samples[[1]][,intname],
                  ttd_modfit$samples[[2]][,intname],
                  ttd_modfit$samples[[3]][,intname])
  
  predictions <- array(dim = c(length(sc.pred), length(coef1samples)))
  dim(predictions)
  
  if (is.na(coeffname2)) {
    
    for(i in 1:length(sc.pred)){
      
      predictions[i,] <- plogis(logit(intsamples) +
                                  coef1samples * sc.pred[i])
      
    }  
  } else {
    
    coef2samples <- c(ttd_modfit$samples[[1]][,coeffname2],
                      ttd_modfit$samples[[2]][,coeffname2],
                      ttd_modfit$samples[[3]][,coeffname2])
    
    for(i in 1:length(sc.pred)){
      predictions[i,] <- plogis(logit(intsamples) +
                                  coef1samples * sc.pred[i] + 
                                  coef2samples * sc.pred[i]^2)
    }  
    
  }
  
  LPB <-  apply(predictions, 1, quantile, probs = 0.025) # Lower bound
  UPB <-  apply(predictions, 1, quantile, probs = 0.975) # Upper bound
  y <- apply(predictions, 1, mean) 
  ylim <- ylim
  
  plotdf <- data_frame(xcol = orig.pred, ycol = y, lower = LPB, upper = UPB)
  
  finalplot <- ggplot(plotdf, aes(x = xcol, y = ycol)) +
    geom_line(colour="black", size = 2) +
    geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
    xlab(xlab) +
    ylab(ylab) +
    #scale_x_continuous(breaks = c(8:14))+
    scale_y_continuous(limits = ylim)+
    theme(axis.text.x=element_text(size=16, colour = "black"),
          axis.text.y=element_text(size=16, colour = "black"),
          axis.title = element_text(size = 16,margin = margin(t = 0, r = 20, b = 100, l = 20)),
          panel.grid = element_blank(),panel.background = element_blank(),
          panel.border = element_blank(),
          axis.line = element_line(size = 1),
          axis.ticks = element_line(color = "black",size = 1.2),
          axis.ticks.length = unit(0.2,"cm"))
  
  return(finalplot)
  
}