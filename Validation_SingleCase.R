#Load libraries---------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(seewave)
library(phonTools)
library(zoo)
#Load user defined functions---------------------------------------------------------------------------------------------------------
cwd <- getwd()
source(file.path(cwd,"calc_overtopping.R"))
#Load data--------------------------------------------------------------------------------------------------------------------------
case_folder <- "R35NA_MESH14_2"
bench_case <- "20200902_S1_WG_R35NA.ASC"
# case_folder <- "R14NA_MESH11_4"
# bench_case <- "20200825_S1_WG_R14NA.ASC"
#Load OpenFoam data
filenames = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"postProcessing/surfaceElevation/0/surfaceElevation.dat")
Meas = read.table(filenames[1],header = TRUE)
gauge_pos <- Meas[c(1,2,3),]
Meas <- Meas[-c(1,2,3,4),]
rownames(Meas) <- NULL

#Load OceanWaves3D data
filenames = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"waveGauges.dat")
Meas_OCW = read.table(filenames[1],header=TRUE, skip = 8)
colnames(Meas_OCW) <- c("Time","WG1","WG2","WG3","WG4","WG5","WG6")

#Load and calibrate experimental data
benchfile = file.path("/home/george/Thesis/tsosMi/Shape 1", bench_case)
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

#Plot experimental and numerical Water Elevations-----------------------------------------------------------------------------------
limit = 150000
a <- BenchCalib[1:limit, ]
b <- Meas[1:round(limit)/20, ]

b$Time2 <- b$Time + 16.3

Meas_OCW$Time2 <- Meas_OCW$Time + 16.55
# Meas_1$Time <- Meas_1$Time + 19.5
ggplot() + geom_line(data=a,aes(x=Time, y=WG2, color="Experimental")) + geom_line(data=Meas_OCW,aes(x=Time2, y=WG5, color="Numerical (OCW3D)")) +
  ggtitle(sprintf("WG5 (17.6 m.) %s", case_folder), case_folder) + xlim(25,60) + 
  scale_color_manual(values = c('Experimental' = 'black','Numerical (OCW3D)' = 'red')) + labs(color = 'Legend')

ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + geom_line(data=b,aes(x=Time2, y=gauge_38, color="Numerical (OF)")) + 
  ggtitle(sprintf("WG6-Gauge39 (33.38 m.) %s", case_folder)) +  
  scale_color_manual(values = c('Experimental' = 'black','Numerical (OF)' = 'red')) + xlim(40,80)

#OpenFoam spectra in the foot of the breakwater----------------------------------------------------------------------------------
var <- Meas$gauge_38
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.01/Fr
temp$Fr <- Fr*100
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_num <- temp
var <- BenchCalib$WG6
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.0005/Fr
temp$Fr <- Fr*2000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_exp <- temp

ggplot() + geom_line(data = spectra_exp, aes(x=T, y=Pscaled, color = "Experimental")) + 
  geom_line(data = spectra_num, aes(x=T, y=Pscaled, color = "Numerical (OF)")) + labs(x = "T (sec.)", y = "Pscaled") + 
  scale_color_manual(values = c('Numerical (OF)' = 'red','Experimental' = 'black')) + labs(color = 'Legend') + 
  xlim(0,5) + ggtitle(sprintf("Water Elevation Spectrum vs Wave Period (33.33 m.) %s", case_folder))

ggplot() + geom_line(data = spectra_exp, aes(x=Fr, y=Pscaled, color = "Experimental")) + 
  geom_line(data = spectra_num, aes(x=Fr, y=Pscaled, color = "Numerical (OF)")) + labs(x = "Fr (Hz)", y = "Pscaled") + 
  scale_color_manual(values = c('Numerical (OF)' = 'red','Experimental' = 'black')) + labs(color = 'Legend') + 
  xlim(0,3) + ggtitle(sprintf("Water Elevation Spectrum vs Frequency (33.33 m.) %s", case_folder))

#OCW3D spectra in the middle of the flume---------------------------------------------------------------------------------------- 
equ_time <- seq(from = 1, to = 70 , by = 0.001)
trans_OCW <- as.data.frame(approx(Meas_OCW$Time,Meas_OCW$WG5, xout = equ_time)) #Linearly interpolate data in equidistant axis

var <- trans_OCW$y
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.001/Fr
temp$Fr <- Fr*1000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_num <- temp
var <- BenchCalib$WG5
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.0005/Fr
temp$Fr <- Fr*2000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_exp <- temp

ggplot() + geom_line(data = spectra_exp, aes(x=T, y=Pscaled, color = "Experimental")) + 
  geom_line(data = spectra_num, aes(x=T, y=Pscaled, color = "Numerical (OCW3D)")) + labs(x = "T (sec.)", y = "Pscaled") + 
  scale_color_manual(values = c('Numerical (OCW3D)' = 'red','Experimental' = 'black')) + labs(color = 'Legend') + 
  xlim(0,5) + ggtitle(sprintf("Water Elevation Spectrum vs Wave Period (17.6 m.) %s", case_folder))

ggplot() + geom_line(data = spectra_exp, aes(x=Fr, y=Pscaled, color = "Experimental")) + 
  geom_line(data = spectra_num, aes(x=Fr, y=Pscaled, color = "Numerical (OCW3D)")) + labs(x = "Fr (Hz)", y = "Pscaled") + 
  scale_color_manual(values = c('Numerical (OCW3D)' = 'red','Experimental' = 'black')) + labs(color = 'Legend') + 
  xlim(0,3) + ggtitle(sprintf("Water Elevation Spectrum vs Frequency (17.6 m.) %s", case_folder))

#Compare the experimental and numerical overtopping volumes----------------------------------------------------------------------
#Load the OpenFoam data 
filenames =  file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"postProcessing/overtopping/0/overtopping.dat")
tx  <- readLines(filenames)
tx2  <- gsub(pattern = '\\(', replace = " ", x = tx)
tx3  <- gsub(pattern = '\\)', replace = " ", x = tx2)
writeLines(tx3, con=filenames)
QMeas = read.table(filenames[1],header=FALSE, skip = 1)
colnames(QMeas) <- c("Time","Qx","Qy","Qz")
#Transform the instantaneous overtopping discharge to cumulative overtopping volume 
temp <- QMeas[,c(2,3)]
# temp <- QMeas[,c(2,3)]**2
# temp <- apply(temp, MARGIN = 1, sum)
# temp <- sqrt(temp)/0.8
steps <- diff(QMeas$Time, lag = 1)
temp1 <- temp[1:(nrow(temp)-1),1]
temp2 <- temp[2:nrow(temp),1]
temp <- cumsum(0.5*(temp1 + temp2)*steps)
Vnum <- data.frame(temp)/0.8
#Vnum <- data.frame(apply(temp, MARGIN = 1, sum))/0.8
colnames(Vnum) <- "V"
Vnum$Time <- QMeas[2:nrow(QMeas),1] + 16

#Load the experimental wave gauge data in the overtopping box 
benchfile = file.path("/home/george/Thesis/tsosMi/Shape 1", bench_case)
BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
Signal <- BenchRaw[,c(10)]
#Transform the water level to cumulative overtopping volume
box_type = 3

V <- calc_overtopping(Signal,box_type)
# if (box_type==2) {
#   mouth = 0.28 
# } else if (box_type==3) {
#   mouth = 0.112}
# 
# #Calibration function of WG9 
# cald = - 4.38/100
# 
# #Calculate water depth in box
# d = 4.4/100 + (Signal - 9.86172602739726) * cald
# 
# #Apply equations for different  water depths
# f1 = which(d<0.47)
# f2 = which(d>=0.47 & d<0.48)
# f3 = which(d>=0.48)
# v = zeros(length(d),1)
# v[f1] = 0.216*0.128*d[f1]
# v[f2] = 0.296*0.128*(d[f2]-0.47)+0.47*0.216*0.128
# v[f3] = 0.47*0.216*0.128+0.01*0.296*0.128+(0.19097*(d[f3]-0.48)^2+(d[f3]-0.48)*0.52)*0.28
# v=v-mean(v[1:2000])
# q=v/mouth
# 
# V <- data.frame(BenchRaw$Measurement.time.s.)
# names(V)[1] <- "Time"
# V$V <- q #Cumulative overtopping volume, in m^3/m

#Calculate the measurement time step
step = V$Time[2] - V$Time[1] 
Vm <- data.frame(rollmean(V, 1000))
Vmm <- data.frame(rollmedian(V, 1000))
der <- data.frame(diff(Vm$V, lag =1)/step)
der$Time <- Vm[2:nrow(Vm),1] 
names(der)[1]<- "der"
ggplot() + geom_line( data = Vmm, aes(x=Time, y=V, color = "Experimental")) + geom_line(data = Vnum, aes(x=Time, y=V, color = "Numerical")) +
  scale_color_manual(values = c('Numerical' = 'black','Experimental' = 'red')) + 
  ggtitle(sprintf("Cumulative Overtopping %s", case_folder))



                        

