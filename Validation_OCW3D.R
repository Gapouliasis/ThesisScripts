#Load libraries---------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(seewave)
library(ggpubr)
#Read the wave gauge data ----------------------------------------------------------------------------------------------------------
case_folder <- "OCW3D_11_4INT"
filenames = file.path("/home/george/OpenFOAM/george-v1912/run/OCW3D",case_folder,"waveGauges.dat")
# case_folder <- "RUN35NA_MESH11_6"
# filenames = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"waveGauges.dat")
Meas = read.table(filenames[1],header=TRUE, skip = 8)
colnames(Meas) <- c("Time","WG1","WG2","WG3","WG4","WG5","WG6")
vars <- colnames(Meas)

# case_folder <- "OCW3D_6"
# filenames = file.path("/home/george/OpenFOAM/george-v1912/run/OCW3D",case_folder,"waveGauges.dat")
# Meas2 = read.table(filenames[1],header=TRUE, skip = 8)
# colnames(Meas1) <- c("Time","WG1","WG2","WG3","WG4","WG5","WG6")
# 
# case_folder <- "OCW3D_7"
# filenames = file.path("/home/george/OpenFOAM/george-v1912/run/OCW3D",case_folder,"waveGauges.dat")
# Meas3 = read.table(filenames[1],header=TRUE, skip = 8)
# colnames(Meas3) <- c("Time","WG1","WG2","WG3","WG4","WG5","WG6")

#Load and calibrate experimental data
benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200902_S1_WG_R35NA.ASC"
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
colnames(BenchCalib) <- c("WG1","WG2","WG3","WG4","WG5","WG6")

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

#Plot experimental and OCW3D data--------------------------------------------------------------------------------------------------
a <- BenchCalib[1:120000, ]

Meas$Time2 <- Meas$Time + 17
# Meas2$Time2 <- Meas2$Time + 17
# Meas3$Time2 <- Meas3$Time + 17

# ggplot() + geom_line(data = a,aes(x=Time, y=WG4, color = "Measured")) + geom_line(data = Meas,aes(x=Time2, y=WG4, color = "OCW3D")) +
#   geom_point(data = Meas2,aes(x=Time2, y=WG4, color = "OCW3D2")) + geom_line(data = Meas3,aes(x=Time2, y=WG4, color = "OCW3D3")) +
#   scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black','OCW3D2' = 'blue','OCW3D3' = 'green')) + 
#   labs(color = 'OCW3D Mesh') + ggtitle("WG1 (17.3 m.)")

plot1 <- ggplot() + geom_line(data = a,aes(x=Time, y=WG1, color = "Measured")) + geom_line(data = Meas,aes(x=Time2, y=WG1, color = "OCW3D")) +
  scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("WG1 (3.35 m.)")

plot2 <- ggplot() + geom_line(data = a,aes(x=Time, y=WG2, color = "Measured")) + geom_line(data = Meas,aes(x=Time2, y=WG2, color = "OCW3D")) +
  scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("WG2 (16.2 m.)")

plot3 <- ggplot() + geom_line(data = a,aes(x=Time, y=WG3, color = "Measured")) + geom_line(data = Meas,aes(x=Time2, y=WG3, color = "OCW3D")) +
  scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("WG3 (16.9 m.)")

plot4 <- ggplot() + geom_line(data = a,aes(x=Time, y=WG4, color = "Measured")) + geom_line(data = Meas,aes(x=Time2, y=WG4, color = "OCW3D")) +
  scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("WG1 (17.3 m.)")

figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)

annotate_figure(figure,top = text_grob("OCW3D", color = "black", face = "bold", size = 12))

#Read the spatial data from the simulation and plot-------------------------------------------------------------------------------
# filenames = Sys.glob("/home/george/OpenFOAM/george-v1912/run/OCW3D/OCW3D_9/fort.*")
# Meas21 = read.table(filenames[450],header=FALSE)
# colnames(Meas21) <- c("X","Y","WL","R2")
# ggplot() + geom_line(data = Meas21, aes(x=X, y=WL))
# 
# filenames = Sys.glob("/home/george/OpenFOAM/george-v1912/run/OCW3D/OCW3D12_3/fort.*")
# Meas21 = read.table(filenames[1116],header=FALSE)
# colnames(Meas21) <- c("X","Y","WL","R2")
# filenames = Sys.glob("/home/george/OpenFOAM/george-v1912/run/OCW3D/OCW3D_11/fort.*")
# Meas22 = read.table(filenames[1000],header=FALSE)
# colnames(Meas22) <- c("X","Y","WL","R2")
# 
# ggplot() + geom_line(data=Meas21, aes(x=X, y=WL, color = "Mesh12")) +geom_line(data=Meas22, aes(x=X, y=WL, color = "Mesh1")) +
#   scale_color_manual(values = c("Mesh12" = "black", "Mesh1" = "red")) +
#   labs(x = "x (m.)", y = "Water Elevation (m.)",color = 'OCW3D Mesh') + ggtitle("Flume Water Elevation at t=28 sec.")

#Load the wave gauge data from all the OCW3D cases. Keep WG2 data for each case---------------------------------------------------
# filenames = Sys.glob("/home/george/OpenFOAM/george-v1912/run/OCW3D/OCW3D_*/waveGauges.dat")
# temp = read.table(filenames[1],header=TRUE, skip = 8)
# equ_time <- seq(from = 1, to = 70 , by = 0.01)
# temp2 <- approx(temp[,1],temp[,3], xout = equ_time)
# Meas <- as.data.frame(temp2)
# colnames(Meas) <- c("Time","OCW3D_1")
# for (i in 2:length(filenames)){
#   temp = read.table(filenames[i],header=TRUE, skip = 8)
#   temp2 <- approx(temp[,1],temp[,3], xout = equ_time)
#   Meas$temp <- temp2$y
#   names(Meas)[i+1] <- sprintf("OCW3D_%s",i)
# }

# longmeas <- melt(Meas, id.vars="Time")
# ggplot(longmeas, aes(Time,value, col=variable)) + geom_line()

#ggplot(longmeas, aes(Time,value)) + geom_line() + facet_wrap(~variable)

#rm(longmeas)



