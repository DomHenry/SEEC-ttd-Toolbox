plot_detec_probability <- function(spp, tmax) {
  
  ref <- which(lam_est_spp$param %in% spp)
  mean.lambda <- as.numeric(lam_est_spp[ref,2])
  upper <- as.numeric(lam_est_spp[ref,4])
  lower <- as.numeric(lam_est_spp[ref,3])
  
  duration <- 1:tmax 
  p.pred <- 1 - exp(-mean.lambda*duration)
  upp <- 1-exp(-upper*duration)
  low <- 1-exp(-lower*duration)
    
  plotdf <- data_frame(xcol = duration, ycol = p.pred, lower = low, upper = upp)
    
  ggplot(plotdf, aes(x = xcol, y = ycol)) +
      geom_line(colour="black", size = 1) +
      geom_ribbon(aes(ymin=lower, ymax=upper), alpha=0.2) +
      labs(title = spp) +
      xlab("\nSurvey duration (seconds)") +
      ylab("Detection probability\n") +
      scale_y_continuous(limits = c(0,0.5))+
      theme(axis.text.x=element_text(size=12, colour = "black"),
            axis.text.y=element_text(size=12, colour = "black"),
            axis.title = element_text(size = 12,margin = margin(t = 0, r = 20, b = 100, l = 20)),
            panel.grid = element_blank(),panel.background = element_blank(),
            panel.border = element_blank(),
            axis.line = element_line(size = 1),
            axis.ticks = element_line(color = "black",size = 1.2),
            axis.ticks.length = unit(0.2,"cm"),
            plot.title = element_text(size = 15))
  
 
}

