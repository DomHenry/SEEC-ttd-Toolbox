plot_random_effects <- function(pm,coeff_name,main_lab,xlims,n){
  
    plot(pm[["mean"]] ,1:n, xlim = xlims, 
       xlab = "Parameter estimate", ylab = "", 
       main = main_lab, pch = 16)
  abline(v = 0, lwd = 2, col = "green")
  
  segments(pm[1:n,"lower"], 1:n, pm[1:n,"upper"], 1:n, col = "grey", lwd = 1.5)
  sig1 <- (pm[1:n,"lower"] * pm[1:n,"upper"]) > 0
  
  segments(pm[1:n,"lower"][sig1 == 1], (1:n)[sig1 == 1], pm[1:n,"upper"][sig1 == 1], (1:n)[sig1 == 1], col = "blue", lwd = 1.5)
  points(pm[["mean"]], 1:n, pch = 16)
  
  abline(v = plot_data[which(plot_data$param == paste0("mu.",coeff_name)),2], col = "red", lwd = 3)
  abline(v = plot_data[which(plot_data$param == paste0("mu.",coeff_name)),3:4], col = "red", lwd = 2, lty = 2)
  
  }