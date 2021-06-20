#Load libraries----------------------------------------------------------------------------------------------------------------------
library(phonTools)
library(zoo)
library(ggplot2)
library(ggpubr)
library(copula)
library(POT)
library(parallel)
#Load user defined functions---------------------------------------------------------------------------------------------------------
cwd <- getwd()
source(file.path(cwd,"zero_downx.R"))
source(file.path(cwd,"calc_overtopping.R"))
#Load experimental data--------------------------------------------------------------------------------------------------------------
benchfiles = Sys.glob("/home/george/Thesis/tsosMi_Experiments/*.ASC")
#numCores <- detectCores()

for (benchfile in benchfiles){
#main_function <- function(benchfile){
  k <- match(benchfile, benchfiles)
  BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
  WG7 <- BenchRaw[1:(nrow(BenchRaw)-1000),c(1,8)]
  names(WG7)[1] <- "Time"
  names(WG7)[2] <- "WG7"
  WG7$WG7 <- WG7$WG7/100
  WG6 <- BenchRaw[,c(1,7)]
  names(WG6)[1] <- "Time"
  names(WG6)[2] <- "Surf"
  WG6$Surf <- WG6$Surf*2.27/100
  correction <- mean(WG6[1:2000,2])
  WG6$Surf <- WG6$Surf - correction
  WG5 <- BenchRaw[,c(1,6)]
  names(WG5)[1] <- "Time"
  names(WG5)[2] <- "Surf"
  WG5$Surf <- WG5$Surf*2.38/100
  correction <- mean(WG5[1:2000,2])
  WG5$Surf <- WG5$Surf - correction
  
  #Calculate, filter and plot the cumulative overtopping in m^3/m---------------------------------------------------------------------
  Signal <- BenchRaw[,c(10)]
  box_type = 3
  V <- calc_overtopping(Signal, box_type)
  
  #Calculate the measurement time step
  step = V$Time[2] - V$Time[1] 
  Vm <- data.frame(rollmean(V, 1000))
  Vmm <- data.frame(rollmedian(V, 1000))
  der <- data.frame(diff(Vm$V, lag =1)/step)
  der$Time <- Vm[2:nrow(Vm),1]
  names(der)[1]<- "der"
  
  # ggplot() +
  #   geom_line(data = V, aes(x=Time, y=V, color = "Raw")) + labs(color = 'Overtopping') + ggtitle("Cumulative Overtopping")
  # ggplot() +
  #   geom_line( data = Vmm, aes(x=Time, y=V, color = "Median")) + labs(color = 'Overtopping') + ggtitle("Cumulative Overtopping")
  # ggplot(data = V) + geom_line(aes(x=Time, y=V, color = "Raw")) + geom_line( data = Vm, aes(x=Time, y=V, color = "Mean")) +
  #   geom_line( data = Vmm, aes(x=Time, y=V, color = "Median")) + labs(color = 'Overtopping') + ggtitle("Cumulative Overtopping") +
  #   scale_color_manual(values = c('Raw' = 'black','Mean' = 'red', "Median" = "blue"))

  # Isolate the individual overtopping events-----------------------------------------------------------------------------------------
  # Find the time delay between WG7 signal and the overtopping volume signal. It employs the derivative of the 
  # overtopping signal
  nlag <- 1500
  correls = zeros(nlag,1)
  for (i in 1:nlag){
    correls[i] <- cor(der[1:(nrow(der)+1-i),1], WG7[(1+i-1):(nrow(WG7)),1], method = "pearson")
  }
  
   WG7$Time <- WG7$Time + which.max(correls)*step 
  # ggplot() + geom_line(data = V,aes(x = Time, y=q, color = "Volume")) +
  #   geom_line(data = der, aes(x=Time, y=der, color = "Derivative")) + geom_line(data = WG7, aes(x= Time, y = WG7, color = "WG7"))
  
  # Find the time of the overtopping events 
  dat <- WG7
  names(dat)[1] <- "time"
  names(dat)[2] <- "obs"
  a <- clust(dat, u = 0.003, tim.cond = 1, clust.max = TRUE)
  peaks <- data.frame(WG7[a[,3],2])
  names(peaks)[1] <- "peaks"
  peaks$Time <- WG7[a[,3],1]
  peaks$inds <- a[,3] 
  peaks <- peaks[order(peaks$Time),]
  # ggplot() + geom_point(data = peaks, aes(x=Time, y=peaks, color = "Peaks")) +
  #    geom_line(data = WG7, aes(x= Time, y = WG7, color = "WG7"))

  #Calculate the mean and median individual overtopping volume. 
  #Median appears to provide better results 
  final_obs <- data.frame(0,V[nrow(V),1],0)
  names(final_obs) <- c("peaks","Time","inds")
  peaks <- rbind(peaks, final_obs)
  peaks$Time <- round(peaks$Time, digits = 4 )
  ind <- match(peaks$Time, V$Time)
  Vmean <- zeros((nrow(peaks)),1)
  Vmedian <- zeros((nrow(peaks)),1)
  for (i in 2:nrow(peaks)){
    Vmean[i] <- mean(V[ind[i-1]:ind[i],2])  
    Vmedian[i] <- median(V[ind[i-1]:ind[i],2])  
  }
  
  Vevent <- data.frame(Vmean[2:nrow(Vmean)])
  Vevent$Median <- Vmedian[2:nrow(Vmedian)]
  Vevent$Time <- V[ind[(1:nrow(peaks)-1)],1]
  names(Vevent)[1] <- "Mean"
  # ggplot(data = V) + geom_line( data = Vmm, aes(x=Time, y=V, color = "Median")) +
  #   labs(color = 'Overtopping') + ggtitle("Cumulative Overtopping") +
  #   geom_point(data = Vevent, aes(x=Time, y=Median, color = "event")) +
  #   scale_color_manual(values = c("Median" = "black", "event" = "orange"))
  # ggplot(data = V) + geom_line(aes(x=Time, y=V, color = "Raw")) + geom_line( data = Vm, aes(x=Time, y=V, color = "Mean")) +
  #   geom_line( data = Vmm, aes(x=Time, y=V, color = "Median")) + labs(color = 'Overtopping', y = expression(V (m^3/s))) + ggtitle("Cumulative Overtopping") +
  #   geom_point(data = Vevent, aes(x=Time, y=Median, color = "event")) +
  #   scale_color_manual(values = c('Raw' = 'black','Mean' = 'red', "Median" = "blue", "event" = "orange"))

  Vindiv <- as.data.frame(diff(Vmedian, lag = 1)*1000)
  names(Vindiv)[1] <- "V"
  Vindiv$Time <- Vevent$Time
  #Find the zero-down crossings--------------------------------------------------------------------------------------------------------
  out <- zero_downx(WG6)
  waves <- out[[1]]
  #out[[2]] + xlim(100,150) + labs(color = "Legend", x = "Time (sec.)", y = "Surface Elevation (m.)")
  # Clean the record from the noise
  waves <- waves[!(waves$height<0.01),]
  
  #Associate waves and volumes-------------------------------------------------------------------------------------------------------
  wave_index <- zeros(nrow(Vindiv),1)
  for (i in 1:nrow(wave_index)){
    wave_index[i] <- which((waves$Tstart <= Vindiv[i,2] & waves$Tend >= Vindiv[i,2]))
  }
  waves_paired <- waves[wave_index-1,]
  waves_paired$V <- Vindiv$V
  waves_paired$Case <- substr(benchfile,55,59)
  pso_paired <- as.data.frame(pobs(waves_paired))
  pso_paired$Case <- substr(benchfile,55,59)
  assign(sprintf("waves_paired%i",k), waves_paired)
  assign(sprintf("pso_paired%i",k), pso_paired)
}

#results = mclapply(benchfiles, main_function, mc.cores = numCores)

  
#Plot the correlations---------------------------------------------------------------------------------------------------------------
#IR34JS - PM
rho1 <- round(cor(pso_paired5$V, pso_paired5$height, method = "spearman"), digits = 2)
rho2 <- round(cor(pso_paired6$V, pso_paired6$height, method = "spearman"), digits = 2)
set <- rbind(waves_paired5,waves_paired6)
set1 <- rbind(pso_paired5,pso_paired6)
plot11 <- ggplot(data = set, aes(x=V, y=height, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
  scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
  annotate("label", x = 17.5, y = 0.1, label = sprintf("atop(rho[IR34J]==%1.2f ,rho[IR34P]==%1.2f)",rho1,rho2), parse = TRUE) + 
  labs(color = 'Case', y = "Wave Height (m.)", x = "V (lts)") + ggtitle("Volume-Height") + ylim(0.08,0.35)
plot21 <- ggplot(data = set1, aes(x=V, y=height, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
  scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
  labs(color = 'Case', y = "Wave Height (m.)", x = "V (lts)") + ggtitle("Volume-Height") + xlim(0, 1) 
 #+ annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")

rho1 <- round(cor(pso_paired5$V, pso_paired5$surf_min, method = "spearman"), digits = 2)
rho2 <- round(cor(pso_paired6$V, pso_paired6$surf_min, method = "spearman"), digits = 2)

plot12 <- ggplot(data = set, aes(x=V, y=surf_min, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
  scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
  annotate("label", x = 17, y = -0.105, label = sprintf("atop(rho[IR34J]==%1.2f ,rho[IR34P]==%1.2f)",rho1,rho2), parse = TRUE) + 
  labs(color = 'Case', y = "Trough Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Trough")
plot22 <- ggplot(data = set1, aes(x=V, y=surf_min, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
  scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
  labs(color = 'Case', y = "Trough Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Trough") + xlim(0, 1) 
  # annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")


rho1 <- round(cor(pso_paired5$V, pso_paired5$surf_max, method = "spearman"), digits = 2)
rho2 <- round(cor(pso_paired6$V, pso_paired6$surf_max, method = "spearman"), digits = 2)

plot13 <- ggplot(data = set, aes(x=V, y=surf_max, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
  scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
  annotate("label", x = 17.8, y = 0.1, label = sprintf("atop(rho[IR34J]==%1.2f ,rho[IR34P]==%1.2f)",rho1,rho2), parse = TRUE) + 
  labs(color = 'Case', y = "Crest Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Crest")
plot23 <- ggplot(data = set1, aes(x=V, y=surf_max, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
  scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
  labs(color = 'Case', y = "Crest Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Crest") + xlim(0, 1) 
  # annotate(geom="text",x=0.8, y=0.92, label=sprintf("Rho = %.2f", rho), color = "blue")

figure <- ggarrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3 , nrow = 2)

annotate_figure(figure,top = text_grob("Individual Overtopping Scatter Plots (Hs=0.11 m., Ts=3.8 sec., s=0.030)", color = "black", face = "bold", size = 12),
                left = text_grob("         Copula                                                                   Joint Distribution",  face = "bold", rot = 90))

# #IR14JS - PM
# rho1 <- round(cor(pso_paired1$V, pso_paired1$height, method = "spearman"), digits = 2)
# rho2 <- round(cor(pso_paired3$V, pso_paired3$height, method = "spearman"), digits = 2)
# set <- rbind(waves_paired1,waves_paired3)
# set1 <- rbind(pso_paired1,pso_paired3)
# plot11 <- ggplot(data = set, aes(x=V, y=height, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) +
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   annotate("label", x = 17, y = 0.26, label = sprintf("atop(rho[IR14J]==%1.2f ,rho[IR14P]==%1.2f)",rho1,rho2), parse = TRUE) +
#   labs(color = 'Case', y = "Wave Height (m.)", x = "V (lts)") + ggtitle("Volume-Height")
# plot21 <- ggplot(data = set1, aes(x=V, y=height, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   labs(color = 'Case', y = "Wave Height (m.)", x = "V (lts)") + ggtitle("Volume-Height") + xlim(0, 1) 
# 
# #+ annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")
# 
# rho1 <- round(cor(pso_paired1$V, pso_paired1$surf_min, method = "spearman"), digits = 2)
# rho2 <- round(cor(pso_paired3$V, pso_paired3$surf_min, method = "spearman"), digits = 2)
# 
# plot12 <- ggplot(data = set, aes(x=V, y=surf_min, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   annotate("label", x = 17, y = -0.008, label = sprintf("atop(rho[IR14J]==%1.2f ,rho[IR14P]==%1.2f)",rho1,rho2), parse = TRUE) + 
#   labs(color = 'Case', y = "Trough Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Trough")
# plot22 <- ggplot(data = set1, aes(x=V, y=surf_min, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   labs(color = 'Case', y = "Trough Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Trough") + xlim(0, 1) 
# # annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")
# 
# 
# rho1 <- round(cor(pso_paired1$V, pso_paired1$surf_max, method = "spearman"), digits = 2)
# rho2 <- round(cor(pso_paired3$V, pso_paired3$surf_max, method = "spearman"), digits = 2)
# 
# plot13 <- ggplot(data = set, aes(x=V, y=surf_max, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   annotate("label", x = 17, y = 0.175, label = sprintf("atop(rho[IR14J]==%1.2f ,rho[IR14P]==%1.2f)",rho1,rho2), parse = TRUE) + 
#   labs(color = 'Case', y = "Crest Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Crest")
# plot23 <- ggplot(data = set1, aes(x=V, y=surf_max, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   labs(color = 'Case', y = "Crest Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Crest") + xlim(0, 1) 
# # annotate(geom="text",x=0.8, y=0.92, label=sprintf("Rho = %.2f", rho), color = "blue")
# 
# figure <- ggarrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3 , nrow = 2)
# 
# 
# annotate_figure(figure,top = text_grob("Individual Overtopping Scatter Plots (Hs=0.11 m., Ts=1.9 sec., s=0.040)", color = "black", face = "bold", size = 12),
#                 left = text_grob("         Copula                                                                   Joint Distribution",  face = "bold", rot = 90))
# 
# #IR24JS - PM
# rho1 <- round(cor(pso_paired2$V, pso_paired2$height, method = "spearman"), digits = 2)
# rho2 <- round(cor(pso_paired4$V, pso_paired4$height, method = "spearman"), digits = 2)
# set <- rbind(waves_paired2,waves_paired4)
# set1 <- rbind(pso_paired2,pso_paired4)
# plot11 <- ggplot(data = set, aes(x=V, y=height, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   annotate("label", x = 17, y = 0.26, label = sprintf("atop(rho[IR24J]==%1.2f ,rho[IR24P]==%1.2f)",rho1,rho2), parse = TRUE) + 
#   labs(color = 'Case', y = "Wave Height (m.)", x = "V (lts)") + ggtitle("Volume-Height")
# plot21 <- ggplot(data = set1, aes(x=V, y=height, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   labs(color = 'Case', y = "Wave Height (m.)", x = "V (lts)") + ggtitle("Volume-Height") + xlim(0, 1) 
# 
# #+ annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")
# 
# rho1 <- round(cor(pso_paired2$V, pso_paired2$surf_min, method = "spearman"), digits = 2)
# rho2 <- round(cor(pso_paired4$V, pso_paired4$surf_min, method = "spearman"), digits = 2)
# 
# plot12 <- ggplot(data = set, aes(x=V, y=surf_min, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   annotate("label", x = 17, y = -0.008, label = sprintf("atop(rho[IR24J]==%1.2f ,rho[IR24P]==%1.2f)",rho1,rho2), parse = TRUE) + 
#   labs(color = 'Case', y = "Trough Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Trough")
# plot22 <- ggplot(data = set1, aes(x=V, y=surf_min, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   labs(color = 'Case', y = "Trough Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Trough") + xlim(0, 1) 
# # annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")
# 
# 
# rho1 <- round(cor(pso_paired2$V, pso_paired2$surf_max, method = "spearman"), digits = 2)
# rho2 <- round(cor(pso_paired4$V, pso_paired4$surf_max, method = "spearman"), digits = 2)
# 
# plot13 <- ggplot(data = set, aes(x=V, y=surf_max, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   annotate("label", x = 17, y = 0.175, label = sprintf("atop(rho[IR24J]==%1.2f ,rho[IR24P]==%1.2f)",rho1,rho2), parse = TRUE) + 
#   labs(color = 'Case', y = "Crest Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Crest")
# plot23 <- ggplot(data = set1, aes(x=V, y=surf_max, group = Case)) + geom_point(aes(shape= Case, color = Case), size = 2) + 
#   scale_shape_manual(values = c(15,16)) + scale_color_manual(values = c("black","red")) +
#   labs(color = 'Case', y = "Crest Elevation (m.)", x = "V (lts)") + ggtitle("Volume-Crest") + xlim(0, 1) 
# # annotate(geom="text",x=0.8, y=0.92, label=sprintf("Rho = %.2f", rho), color = "blue")
# 
# figure <- ggarrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3 , nrow = 2)
# 
# 
# annotate_figure(figure,top = text_grob("Individual Overtopping Scatter Plots (Hs=0.11 m., Ts=3 sec., s=0.035)", color = "black", face = "bold", size = 12),
#                 left = text_grob("         Copula                                                                   Joint Distribution",  face = "bold", rot = 90))
