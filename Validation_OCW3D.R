library(ggplot2)
library(reshape2)
library(seewave)
library(ggpubr)
filenames = "/home/george/OpenFOAM/george-v1912/run/OCW3D/OCW3D_7/waveGauges.dat"
Meas = read.table(filenames[1],header=TRUE, skip = 8)
names(Meas)[1] <- "Time"
names(Meas)[2] <- "WG1"
names(Meas)[3] <- "WG2"
names(Meas)[4] <- "WG3"
names(Meas)[5] <- "WG4"
names(Meas)[6] <- "WG5"
names(Meas)[7] <- "WG6"

filenames = Sys.glob("/home/george/OpenFOAM/george-v1912/run/OCW3D/OCW3D_1/fort.*")
Meas2 = read.table(filenames[400],header=FALSE)
colnames(Meas2) <- c("X","Y","WL","R2")

ggplot(data = Meas2) + geom_line(aes(x=X, y=WL))

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
a <- BenchCalib[1:120000, ]

Meas$Time <- Meas$Time +10.5
Meas2$Time <- Meas2$Time +10.5

plot1 <- ggplot() + geom_line(data = a,aes(x=Time, y=WG1, color = "Measured")) + geom_line(data = Meas,aes(x=Time, y=WG1, color = "OCW3D")) +
  geom_line(data = Meas2,aes(x=Time, y=WG1, color = "OCW3D 2")) +
  scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black', 'OCW3D 2' = 'blue')) + labs(color = 'OCW3D Mesh') + ggtitle("WG1 (3.35 m.)")

plot2 <- ggplot() + geom_line(data = a,aes(x=Time, y=WG2, color = "Measured")) + geom_line(data = Meas,aes(x=Time, y=WG2, color = "OCW3D")) +
  geom_line(data = Meas2,aes(x=Time, y=WG2, color = "OCW3D 2")) +
  scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black', 'OCW3D 2' = 'blue')) + labs(color = 'OCW3D Mesh') + ggtitle("WG2 (16.2 m.)")

plot3 <- ggplot() + geom_line(data = a,aes(x=Time, y=WG3, color = "Measured")) + geom_line(data = Meas,aes(x=Time, y=WG3, color = "OCW3D")) +
  geom_line(data = Meas2,aes(x=Time, y=WG3, color = "OCW3D 2")) +
  scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black', 'OCW3D 2' = 'blue')) + labs(color = 'OCW3D Mesh') + ggtitle("WG3 (16.9 m.)")

plot4 <- ggplot() + geom_line(data = a,aes(x=Time, y=WG4, color = "Measured")) + geom_line(data = Meas,aes(x=Time, y=WG4, color = "OCW3D")) +
  geom_line(data = Meas2,aes(x=Time, y=WG4, color = "OCW3D 2")) +
  scale_color_manual(values = c('Measured' = 'red','OCW3D' = 'black', 'OCW3D 2' = 'blue')) + labs(color = 'OCW3D Mesh') + ggtitle("WG1 (17.3 m.)")

figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)

annotate_figure(figure,top = text_grob("OCW3D", color = "black", face = "bold", size = 12))

# 
# spectra <- as.data.frame(exp_spec)
# names(spectra)[2] <- "Exp" 
# num_spec <- meanspec(Meas$gauge_39, f = 100, wl = 1024, PSD = TRUE, plot = FALSE, FUN = "mean")
# temp <- as.data.frame(num_spec)
# spectra$Num <- temp$y






