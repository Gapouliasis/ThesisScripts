#overtopping <- function(Signal, box_type)
# Function translate WG9 signal to overtopping, for box 2 & 3. WG8 and
# results and box=1 results should not be used here. 
# Function calculates volume of water in the tank/meter of mouth width. 
# Number of box (G column of data code)
#{
library(phonTools)
library(zoo)
library(ggplot2)
library(ggpubr)
library(copula)
library(POT)
source("/home/george/OpenFOAM/george-v1912/run/RScripts/zero_downx.R")
#benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200902_S1_WG_R16NA.ASC"
benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200825_S1_WG_R14NA.ASC"
#benchfile = Sys.glob("/home/george/Thesis/tsosMi_Experiments/*.ASC")
#BenchRaw = read.table(benchfile[4],header = TRUE, sep = ";", skip = 6)
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

Signal <- BenchRaw[,c(10)]
box_type = 2

if (box_type==2) {
  mouth = 0.28 
} else if (box_type==3) {
  mouth = 0.112}

#Calibration function of WG9 
cald = - 4.38/100

#Calculate water depth in box
d = 4.4/100 + (Signal - 9.86172602739726) * cald

#Apply equations for different  water depths
f1 = which(d<0.47)
f2 = which(d>=0.47 & d<0.48)
f3 = which(d>=0.48)
v = zeros(length(d),1)
v[f1] = 0.216*0.128*d[f1]
v[f2] = 0.296*0.128*(d[f2]-0.47)+0.47*0.216*0.128
v[f3] = 0.47*0.216*0.128+0.01*0.296*0.128+(0.19097*(d[f3]-0.48)^2+(d[f3]-0.48)*0.52)*0.28
v=v-mean(v[1:2000])
q=v/mouth

V <- data.frame(BenchRaw$Measurement.time.s.)
names(V)[1] <- "Time"
V$V <- q #Cumulative overtopping volume, in m^3/m

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

# Isolate the individual overtopping events 
# Find the time delay between WG7 signal and the overtopping volume signal. It employs the derivative of the 
# overtopping signal
nlag <- 1500
correls = zeros(nlag,1)
for (i in 1:nlag){
  correls[i] <- cor(der[1:(nrow(der)+1-i),1], WG7[(1+i-1):(nrow(WG7)),1], method = "pearson")
}

# WG7$Time <- WG7$Time + which.max(correls)*step 
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
ggplot() + geom_point(data = peaks, aes(x=Time, y=peaks, color = "Peaks")) +
   geom_line(data = WG7, aes(x= Time, y = WG7, color = "WG7"))

#Calculate the mean and median individual overtopping volume. Median appears to provide better results 
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
ggplot(data = V) + geom_line( data = Vmm, aes(x=Time, y=V, color = "Median")) +
  labs(color = 'Overtopping') + ggtitle("Cumulative Overtopping") +
  geom_point(data = Vevent, aes(x=Time, y=Median, color = "event")) +
  scale_color_manual(values = c("Median" = "black", "event" = "orange"))
# ggplot(data = V) + geom_line(aes(x=Time, y=V, color = "Raw")) + geom_line( data = Vm, aes(x=Time, y=V, color = "Mean")) +
#   geom_line( data = Vmm, aes(x=Time, y=V, color = "Median")) + labs(color = 'Overtopping') + ggtitle("Cumulative Overtopping") +
#   geom_point(data = Vevent, aes(x=Time, y=Median, color = "event")) +
#   scale_color_manual(values = c('Raw' = 'black','Mean' = 'red', "Median" = "blue", "event" = "orange"))

Vindiv <- as.data.frame(diff(Vmedian, lag = 1)*1000)
names(Vindiv)[1] <- "V"
Vindiv$Time <- Vevent$Time
#Find the zero-down crossings 
out <- zero_downx(WG6)
waves <- out[[1]]
#out[[2]]
# Clean the record from the noise
waves <- waves[!(waves$height<0.01),]

#Associate waves and volumes
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
