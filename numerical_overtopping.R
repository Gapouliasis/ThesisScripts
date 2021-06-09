#Load libraries----------------------------------------------------------------------------------------------------------------------
library(phonTools)
library(zoo)
library(ggplot2)
library(ggpubr)
library(copula)
library(POT)
#Load user defined functions---------------------------------------------------------------------------------------------------------
cwd <- getwd()
source(file.path(cwd,"zero_downx.R"))
source(file.path(cwd,"calc_overtopping.R"))
#Load numerical data--------------------------------------------------------------------------------------------------------------
case_folder <- "IRCase_Test"
filename = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"postProcessing/surfaceElevation/0/surfaceElevation.dat")
Meas = read.table(filename,header = TRUE)
gauge_pos <- Meas[c(1,2,3),]
Meas <- Meas[-c(1,2,3,4),]
rownames(Meas) <- NULL

# Isolate the individual overtopping events-----------------------------------------------------------------------------------------
#Find the zero-down crossings--------------------------------------------------------------------------------------------------------
gauge <- Meas[,c(1,251)]
names(gauge)[2] <- "Surf"
#out <- zero_downx(dat) #Requires the first column to be Time and the second surface elevation
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

ggplot() + geom_line(data = gauge, aes(x=Time, y=Surf, color = "Gauge")) + geom_point(data = zero_crossings, aes(x=Time, y=Surf, color = "Crossings"))+
  geom_point(data = waves, aes(x=surf_max_time, y=surf_max, color = "Peak")) +  geom_point(data = waves, aes(x=surf_min_time, y=surf_min, color = "Trough")) +
  scale_color_manual(values = c('Gauge' = 'black', "Crossings" = "orange", "Peak" = "red", "Trough" = "blue")) + ggtitle("Zero Down Crossings")

waves <- out[[1]]
#out[[2]]
# Clean the record from the noise
waves <- waves[!(waves$height<0.01),]


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
#   geom_line( data = Vmm, aes(x=Time, y=V, color = "Median")) + labs(color = 'Overtopping') + ggtitle("Cumulative Overtopping") +
#   geom_point(data = Vevent, aes(x=Time, y=Median, color = "event")) +
#   scale_color_manual(values = c('Raw' = 'black','Mean' = 'red', "Median" = "blue", "event" = "orange"))

Vindiv <- as.data.frame(diff(Vmedian, lag = 1)*1000)
names(Vindiv)[1] <- "V"
Vindiv$Time <- Vevent$Time
#Find the zero-down crossings--------------------------------------------------------------------------------------------------------
out <- zero_downx(WG6)
waves <- out[[1]]
#out[[2]]
# Clean the record from the noise
waves <- waves[!(waves$height<0.01),]

#Associate waves and volumes--------------------------------------------------------------------------------------------------------
wave_index <- zeros(nrow(Vindiv),1)
for (i in 1:nrow(wave_index)){
  wave_index[i] <- which((waves$Tstart <= Vindiv[i,2] & waves$Tend >= Vindiv[i,2]))
}
waves_paired <- waves[wave_index-2,]
waves_paired$V <- Vindiv$V
pso_paired <- as.data.frame(pobs(waves_paired))
plot11 <- ggplot(data = waves_paired) + geom_point(aes(x=V, y=height)) + ggtitle("Volume-Height")
rho <- round(cor(pso_paired$V, pso_paired$height, method = "spearman"), digits = 2)
plot21 <- ggplot(data = pso_paired) + geom_point(aes(x=V, y=height)) + xlim(0, 1) +
  annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")

plot12 <- ggplot(data = waves_paired) + geom_point(aes(x=V, y=surf_min)) + ggtitle("Volume-Trough")
rho <- round(cor(pso_paired$V, pso_paired$surf_min, method = "spearman"), digits = 2)
plot22 <- ggplot(data = pso_paired) + geom_point(aes(x=V, y=surf_min)) + xlim(0, 1) +
  annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")

plot13 <- ggplot(data = waves_paired) + geom_point(aes(x=V, y=surf_max)) + ggtitle("Volume-Peak")
rho <- round(cor(pso_paired$V, pso_paired$surf_max, method = "spearman"), digits = 2)
plot23 <- ggplot(data = pso_paired) + geom_point(aes(x=V, y=surf_max)) + xlim(0, 1) +
  annotate(geom="text",x=0.8, y=0.95, label=sprintf("Rho = %.2f", rho), color = "blue")

figure <- ggarrange(plot11, plot12, plot13, plot21, plot22, plot23, ncol = 3 , nrow = 2)

annotate_figure(figure,top = text_grob("Individual Overtopping Scatter Plots (Hs=0.11 m., Ts=3 sec., PM)", color = "black", face = "bold", size = 12),
                left = text_grob("          Copula                                        Joint Distribution",  face = "bold", rot = 90))
