zero_downx <- function(gauge){
  library(data.table)
  SignalLagged <- gauge[1:nrow(gauge)-1,] 
  Signal <- gauge[2:nrow(gauge),]
  zero_logical <- SignalLagged[,2] >= 0 & Signal[,2] <=0
  zero_ind <- which(zero_logical)
  zero_crossings <- Signal[zero_ind,] 
  
  setDT(zero_crossings, keep.rownames = "Index")
  Tstart <- zeros(nrow(zero_crossings)-1,1)
  Tend <- zeros(nrow(zero_crossings)-1,1)
  height <- zeros(nrow(zero_crossings)-1,1)
  surf_min <- zeros(nrow(zero_crossings)-1,1)
  surf_min_ind <- zeros(nrow(zero_crossings)-1,1)
  surf_max <- zeros(nrow(zero_crossings)-1,1)
  surf_max_ind <- zeros(nrow(zero_crossings)-1,1)
  for (i in 2:nrow(zero_crossings)){
    low <- as.integer(zero_crossings[(i-1),1])
    high <- as.integer(zero_crossings[i,1])
    Tstart[i-1] <- gauge[low,1]
    Tend[i-1] <- gauge[high,1]
    temp <- gauge[low:high,2]
    surf_min_ind[i-1] <- low + which.min(temp)
    surf_min[i-1] <- min(temp)
    surf_max_ind[i-1] <- low + which.max(temp)
    surf_max[i-1] <- max(temp)
    height[i-1] <- max(temp) - min(temp)
  }
  surf_min_time <- gauge[surf_min_ind,1]
  surf_max_time <- gauge[surf_max_ind,1]
  waves <- data.frame(height)
  waves$Tstart <- Tstart
  waves$Tend <- Tend
  waves$surf_min <- surf_min
  waves$surf_min_time <- surf_min_time
  waves$surf_max <- surf_max
  waves$surf_max_time <- surf_max_time
  plot <- ggplot() + geom_line(data = gauge, aes(x=Time, y=Surf, color = "Gauge")) + geom_point(data = zero_crossings, aes(x=Time, y=Surf, color = "Crossings"))+
    geom_point(data = waves, aes(x=surf_max_time, y=surf_max, color = "Peak")) +  geom_point(data = waves, aes(x=surf_min_time, y=surf_min, color = "Trough")) +
    scale_color_manual(values = c('Gauge' = 'black', "Crossings" = "orange", "Peak" = "red", "Trough" = "blue")) + ggtitle("Zero Down Crossings")
  
  return(list(waves,plot))
}