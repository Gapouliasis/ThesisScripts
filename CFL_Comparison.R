#Load libraries---------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(seewave)
library(phonTools)
library(zoo)
library(ggpubr)
#Load user defined functions---------------------------------------------------------------------------------------------------------
cwd <- getwd()
source(file.path(cwd,"calc_overtopping.R"))
#Load cases for comparison-----------------------------------------------------------------------------------------------------------
case_folders <- c("R24NA_MESH15_ThC010","R24NA_MESH15_ThC025","R24NA_MESH15_Th","R24NA_MESH15_ThC075","R24NA_MESH15_ThC1")
for (case_folder in case_folders){
  filenames = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"postProcessing/overtopping/0/overtopping.dat")
  tx  <- readLines(filenames)
  tx2  <- gsub(pattern = '\\(', replace = " ", x = tx)
  tx3  <- gsub(pattern = '\\)', replace = " ", x = tx2)
  writeLines(tx3, con=filenames)
  assign(case_folder,read.table(filenames[1],header=FALSE, skip = 1))
}

vars <- colnames(Meas)

benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200825_S1_WG_R24NA.ASC"
BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
expmeas <-BenchRaw[,c(1,2,3,4,5,6,7,10)]
rm(BenchRaw)

calibfile = "/home/george/Thesis/tsosMi/Calibration.txt"
calibration = read.table(calibfile,header = FALSE)
calibration <- as.vector(t(calibration))
wgun <- expmeas[,c(2:8)]
wgun <- as.matrix(wgun)
temp <- wgun%*%diag(calibration)
temp <- as.data.frame(temp)
BenchCalib <- temp
colnames(BenchCalib) <- c("WG1","WG2","WG3","WG4","WG5","WG6","WG10")

mu <- colMeans(BenchCalib[1:5000, ])
BenchCalib$WG1 <- BenchCalib$WG1 - mu[1]
BenchCalib$WG2 <- BenchCalib$WG2 - mu[2]
BenchCalib$WG3 <- BenchCalib$WG3 - mu[3]
BenchCalib$WG4 <- BenchCalib$WG4 - mu[4]
BenchCalib$WG5 <- BenchCalib$WG5 - mu[5]
BenchCalib$WG6 <- BenchCalib$WG6 - mu[6]

BenchCalib <- BenchCalib/100

BenchCalib$Time <- expmeas[,c(1)] 
rm(expmeas)

#Compare the experimental and numerical overtopping volumes----------------------------------------------------------------------
#Transform the instantaneous overtopping discharge to cumulative overtopping volume 
temp <- R24NA_MESH15_Th[,c(2,3)]
steps <- diff(R24NA_MESH15_Th$V1, lag = 1)
temp1 <- temp[1:(nrow(temp)-1),1]
temp2 <- temp[2:nrow(temp),1]
temp <- cumsum(0.5*(temp1 + temp2)*steps)
V1 <- data.frame(temp)
#Vnum <- data.frame(apply(temp, MARGIN = 1, sum))/0.8
colnames(V1) <- "V"
V1$Time <- R24NA_MESH15_Th[2:nrow(R24NA_MESH15_Th),1] + 35

temp <- R24NA_MESH15_ThC010[,c(2,3)]
steps <- diff(R24NA_MESH15_ThC010$V1, lag = 1)
temp1 <- temp[1:(nrow(temp)-1),1]
temp2 <- temp[2:nrow(temp),1]
temp <- cumsum(0.5*(temp1 + temp2)*steps)
V2 <- data.frame(temp)
#Vnum <- data.frame(apply(temp, MARGIN = 1, sum))/0.8
colnames(V2) <- "V"
V2$Time <- R24NA_MESH15_ThC010[2:nrow(R24NA_MESH15_ThC010),1] + 35

temp <- R24NA_MESH15_ThC025[,c(2,3)]
steps <- diff(R24NA_MESH15_ThC025$V1, lag = 1)
temp1 <- temp[1:(nrow(temp)-1),1]
temp2 <- temp[2:nrow(temp),1]
temp <- cumsum(0.5*(temp1 + temp2)*steps)
V3 <- data.frame(temp)
#Vnum <- data.frame(apply(temp, MARGIN = 1, sum))/0.8
colnames(V3) <- "V"
V3$Time <- R24NA_MESH15_ThC025[2:nrow(R24NA_MESH15_ThC025),1] + 35

temp <- R24NA_MESH15_ThC075[,c(2,3)]
steps <- diff(R24NA_MESH15_ThC075$V1, lag = 1)
temp1 <- temp[1:(nrow(temp)-1),1]
temp2 <- temp[2:nrow(temp),1]
temp <- cumsum(0.5*(temp1 + temp2)*steps)
V4 <- data.frame(temp)
#Vnum <- data.frame(apply(temp, MARGIN = 1, sum))/0.8
colnames(V4) <- "V"
V4$Time <- R24NA_MESH15_ThC075[2:nrow(R24NA_MESH15_ThC075),1] + 35

temp <- R24NA_MESH15_ThC1[,c(2,3)]
steps <- diff(R24NA_MESH15_ThC1$V1, lag = 1)
temp1 <- temp[1:(nrow(temp)-1),1]
temp2 <- temp[2:nrow(temp),1]
temp <- cumsum(0.5*(temp1 + temp2)*steps)
V5 <- data.frame(temp)
#Vnum <- data.frame(apply(temp, MARGIN = 1, sum))/0.8
colnames(V5) <- "V"
V5$Time <- R24NA_MESH15_ThC1[2:nrow(R24NA_MESH15_ThC1),1] + 35

#Load the experimental wave gauge data in the overtopping box 
#benchfile = file.path("/home/george/Thesis/tsosMi/Shape 1", bench_case)
BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
Signal <- BenchRaw[,c(10)]
#Transform the water level to cumulative overtopping volume
box_type = 2

V <- calc_overtopping(Signal,box_type)

#Calculate the measurement time step
step = V$Time[2] - V$Time[1] 
Vm <- data.frame(rollmean(V, 1000))
Vmm <- data.frame(rollmedian(V, 1000))
der <- data.frame(diff(Vm$V, lag =1)/step)
der$Time <- Vm[2:nrow(Vm),1] 
names(der)[1]<- "der"
imark <- floor(seq(from=1,to=dim(V1)[1],length=100))
ggplot() + geom_line( data = V, aes(x=Time, y=V, color = "Experimental")) + 
  #geom_point(data = V1[imark,], aes(x=Time, y=V, color = "CFL=0.50")) +
  geom_line(data = V2, aes(x=Time, y=V, color = "CFL=0.10")) + geom_line(data = V3, aes(x=Time, y=V, color = "CFL=0.25")) + 
  geom_line(data = V4, aes(x=Time, y=V, color = "CFL=0.75")) + geom_line(data = V5, aes(x=Time, y=V, color = "CFL=1.00")) +
  geom_line(data = V1, aes(x=Time, y=V, color = "CFL=0.50")) +
  scale_color_manual(values = c('CFL=0.10' = 'deeppink','CFL=0.25' = 'darkorchid1','CFL=0.50' = 'red',
                                'CFL=0.75' = 'darkorange1','CFL=1.00' = 'darkgreen','Experimental' = 'black')) + 
  ggtitle(sprintf("Cumulative Overtopping %s", case_folder)) + ylim(0,0.012) +xlim(50,100)

#Load the wave gauge data---------------------------------------------------------------------------------------------------------
case_folders <- c("R24NA_MESH15_ThC010","R24NA_MESH15_ThC025","R24NA_MESH15_Th","R24NA_MESH15_ThC075","R24NA_MESH15_ThC1")
for (case_folder in case_folders){
  filenames = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"postProcessing/surfaceElevation/0/surfaceElevation.dat")
  tx  <- readLines(filenames)
  tx2  <- gsub(pattern = '\\(', replace = " ", x = tx)
  tx3  <- gsub(pattern = '\\)', replace = " ", x = tx2)
  writeLines(tx3, con=filenames)
  assign(case_folder,read.table(filenames[1],header=FALSE, skip = 1))
}
                        
var <- BenchCalib$WG6
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.0005/Fr
temp$Fr <- Fr*2000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_exp <- temp

var <- R24NA_MESH15_Th$gauge_48
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.01/Fr
temp$Fr <- Fr*100
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_V1 <- temp

var <- R24NA_MESH15_ThC010$gauge_48
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.01/Fr
temp$Fr <- Fr*100
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_V2 <- temp

var <- R24NA_MESH15_ThC025$gauge_48
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.01/Fr
temp$Fr <- Fr*100
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_V3 <- temp

var <- R24NA_MESH15_ThC075$gauge_48
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.01/Fr
temp$Fr <- Fr*100
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_V4 <- temp

var <- R24NA_MESH15_ThC1$gauge_48
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.01/Fr
temp$Fr <- Fr*100
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_V5 <- temp
