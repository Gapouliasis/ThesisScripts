library(ggplot2)
library(reshape2)
library(seewave)
library(phonTools)
library(zoo)
filenames = "/home/george/OpenFOAM/george-v1912/run/R14NA_MESH11/postProcessing/surfaceElevation/0/surfaceElevation.dat"
Meas = read.table(filenames[1],header = TRUE)
gauge_pos <- Meas[c(1,2,3),]
Meas <- Meas[-c(1,2,3,4),]
rownames(Meas) <- NULL

# filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH9/postProcessing/surfaceElevation/0/surfaceElevation.dat"
# Meas1 = read.table(filenames[1],header=TRUE)
# gauge_pos <- Meas1[c(1,2,3),]
# Meas1 <- Meas1[-c(1,2,3,4),]
# rownames(Meas1) <- NULL
# filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH8/waveGauges.dat"
# Meas_OCW = read.table(filenames[1],header=TRUE, skip = 8)
# names(Meas_OCW)[1] <- "Time"
# names(Meas_OCW)[2] <- "WG1"
# names(Meas_OCW)[3] <- "WG2"
# names(Meas_OCW)[4] <- "WG3"
# names(Meas_OCW)[5] <- "WG4"
# names(Meas_OCW)[6] <- "WG5"
# names(Meas_OCW)[7] <- "WG6"
# 
# filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH10/waveGauges.dat"
# Meas_1 = read.table(filenames[1],header=TRUE, skip = 8)
# names(Meas_1)[1] <- "Time"
# names(Meas_1)[2] <- "WG1"
# names(Meas_1)[3] <- "WG2"
# names(Meas_1)[4] <- "WG3"
# names(Meas_1)[5] <- "WG4"
# names(Meas_1)[6] <- "WG5"
# names(Meas_1)[7] <- "WG6"

#longmeas <- melt(Meas, id.vars="Time")
#ggplot(longmeas, aes(Time,value, col=variable)) + geom_line()
   
#ggplot(longmeas, aes(Time,value)) + geom_line() + facet_wrap(~variable)
  
#rm(longmeas)

vars <- colnames(Meas)

benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200825_S1_WG_R14NA.ASC"
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
names(BenchCalib)[1] <- "WG1"
names(BenchCalib)[2] <- "WG2"
names(BenchCalib)[3] <- "WG3"
names(BenchCalib)[4] <- "WG4"
names(BenchCalib)[5] <- "WG5"
names(BenchCalib)[6] <- "WG6"
names(BenchCalib)[7] <- "WG10"

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

#longbench <- melt(BenchCalib[1:1000000,], id.vars="Time")
#ggplot(longbench, aes(Time,value, col=variable)) + geom_line() 
  
#ggplot(longbench, aes(Time,value)) + geom_line() + facet_wrap(~variable)
limit = 150000
a <- BenchCalib[1:limit, ]
b <- Meas[1:round(limit)/20, ]

b$Time2 <- b$Time + 16

# Meas_OCW$Time <- Meas_OCW$Time + 20.5
# Meas_1$Time <- Meas_1$Time + 19.5
# ggplot() + geom_line(data=a,aes(x=Time, y=WG2, color="Measured")) + geom_line(data=Meas_OCW,aes(x=Time, y=WG2, color="Mesh 5")) +
#   ggtitle("WG2 (16.2 m.)") + geom_line(data=Meas_1,aes(x=Time, y=WG2, color="Mesh 10")) + 
#   scale_color_manual(values = c('Measured' = 'red','Mesh 5' = 'black', 'Mesh 10' = 'blue')) + labs(color = 'Legend')

ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + geom_line(data=b,aes(x=Time2, y=gauge_38, color="Numerical")) + 
  ggtitle("WG6-Gauge39 (33.38 m.)") +  scale_color_manual(values = c('Experimental' = 'red','Numerical' = 'black')) + xlim(40,80)

filenames = "/home/george/OpenFOAM/george-v1912/run/R14NA_MESH12/postProcessing/overtopping/0/overtopping.dat"
tx  <- readLines(filenames)
tx2  <- gsub(pattern = '\\(', replace = " ", x = tx)
tx3  <- gsub(pattern = '\\)', replace = " ", x = tx2)
writeLines(tx3, con=filenames)
QMeas = read.table(filenames[1],header=FALSE, skip = 1)
colnames(QMeas) <- c("Time","Qx","Qy","Qz")

temp <- QMeas[,c(2,3)]**2
temp <- apply(temp, MARGIN = 1, sum)
temp <- sqrt(temp)/0.8
steps <- diff(QMeas$Time, lag = 1)
temp <- cumsum(temp[2:length(temp)]*steps)
Vnum <- data.frame(temp)
colnames(Vnum) <- "V"
Vnum$Time <- QMeas[2:nrow(QMeas),1] + 16
benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200825_S1_WG_R14NA.ASC"
BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
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
ggplot() + geom_line( data = Vmm, aes(x=Time, y=V, color = "Experimental")) + geom_line(data = Vnum, aes(x=Time, y=V, color = "Numerical")) +
  scale_color_manual(values = c('Numerical' = 'black','Experimental' = 'red')) + ggtitle("Cumulative Overtopping")
# a <- BenchCalib[80001:limit, ]
# b <- Meas[1700:nrow(Meas), ]
# exp_spec <- meanspec(a, f = 2000, wl = 1024, PSD = TRUE, plot = FALSE, FUN = "mean")
# spectra <- as.data.frame(exp_spec)
# names(spectra)[2] <- "Exp"
# num_spec <- meanspec(b$gauge_38, f = 100, wl = 1024, PSD = TRUE, plot = FALSE, FUN = "mean")
# temp <- as.data.frame(num_spec)
# spectra$Num <- temp$y
# 
# longmeas <- melt(spectra, id.vars="x")
# ggplot(longmeas, aes(x,value, col=variable)) + geom_line()


                        

