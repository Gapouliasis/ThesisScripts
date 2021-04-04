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
#Load cases for comparison
case_folders <- c("R14NA_MESH11_2","R14NA_MESH11_3","R14NA_MESH11_4","R14NA_MESH11_5","R14NA_MESH11_6")
for (case_folder in case_folders){
  filenames = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"postProcessing/surfaceElevation/0/surfaceElevation.dat")
  assign(case_folder,read.table(filenames[1],header = TRUE))
}

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

#longbench <- melt(BenchCalib[1:1000000,], id.vars="Time")
#ggplot(longbench, aes(Time,value, col=variable)) + geom_line() 
  
#ggplot(longbench, aes(Time,value)) + geom_line() + facet_wrap(~variable)
limit = 150000
a <- BenchCalib[1:limit, ]

MESH1$Time2 <- MESH1$Time + 16
MESH2$Time2 <- MESH2$Time + 16
MESH3$Time2 <- MESH3$Time + 16
MESH4$Time2 <- MESH4$Time + 16
MESH5$Time2 <- MESH5$Time + 16
MESH6$Time2 <- MESH6$Time + 16
MESH7$Time2 <- MESH7$Time + 16
MESH8$Time2 <- MESH8$Time + 16

ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + geom_line(data=MESH1,aes(x=Time2, y=gauge_38, color="MESH1")) + 
  geom_line(data=MESH2,aes(x=Time2, y=gauge_38, color="MESH2")) + geom_line(data=MESH3,aes(x=Time2, y=gauge_38, color="MESH3")) +
  geom_line(data=MESH4,aes(x=Time2, y=gauge_38, color="MESH4")) + ggtitle("Water level comparison on WG6 (33.38 m.) MESH1-4") +
  scale_color_manual(values = c('Experimental' = 'black','MESH1' = 'red','MESH2' = 'blue', 'MESH3' = 'darkorange1', 'MESH4' = 'darkgreen')) + xlim(40,80)

plot1 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + 
  geom_line(data=MESH1,aes(x=Time2, y=gauge_38, color="MESH1")) + scale_color_manual(values = c('Experimental' = 'black','MESH1' = 'cornflowerblue')) + 
   xlim(40,80)
plot2 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + 
  geom_line(data=MESH2,aes(x=Time2, y=gauge_38, color="MESH2")) + scale_color_manual(values = c('Experimental' = 'black','MESH2' = 'cornflowerblue')) + 
   xlim(40,80)
plot3 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + 
  geom_line(data=MESH3,aes(x=Time2, y=gauge_38, color="MESH3")) + scale_color_manual(values = c('Experimental' = 'black','MESH3' = 'cornflowerblue')) + 
  xlim(40,80)
plot4 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + 
  geom_line(data=MESH4,aes(x=Time2, y=gauge_38, color="MESH4")) + scale_color_manual(values = c('Experimental' = 'black','MESH4' = 'cornflowerblue')) + 
  xlim(40,80)

figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)

annotate_figure(figure,top = text_grob("Water level comparison on WG6 (33.38 m.) MESH1-4", color = "black", face = "bold", size = 12))

ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + geom_line(data=MESH5,aes(x=Time2, y=gauge_38, color="MESH5")) + 
  geom_line(data=MESH6,aes(x=Time2, y=gauge_38, color="MESH6")) + geom_line(data=MESH7,aes(x=Time2, y=gauge_38, color="MESH7")) +
  geom_line(data=MESH8,aes(x=Time2, y=gauge_38, color="MESH8")) + ggtitle("Water level comparison on WG6 (33.38 m.) MESH5-8") +
  scale_color_manual(values = c('Experimental' = 'black','MESH5' = 'red','MESH6' = 'blue', 'MESH7' = 'darkorange1', 'MESH8' = 'darkgreen')) + xlim(40,80)

plot1 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + 
  geom_line(data=MESH5,aes(x=Time2, y=gauge_38, color="MESH5")) + scale_color_manual(values = c('Experimental' = 'black','MESH5' = 'cornflowerblue')) + 
  xlim(40,80)
plot2 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + 
  geom_line(data=MESH6,aes(x=Time2, y=gauge_38, color="MESH6")) + scale_color_manual(values = c('Experimental' = 'black','MESH6' = 'cornflowerblue')) + 
  xlim(40,80)
plot3 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + 
  geom_line(data=MESH7,aes(x=Time2, y=gauge_38, color="MESH7")) + scale_color_manual(values = c('Experimental' = 'black','MESH7' = 'cornflowerblue')) + 
  xlim(40,80)
plot4 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + 
  geom_line(data=MESH8,aes(x=Time2, y=gauge_38, color="MESH8")) + scale_color_manual(values = c('Experimental' = 'black','MESH8' = 'cornflowerblue')) + 
  xlim(40,80)

figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)

annotate_figure(figure,top = text_grob("Water level comparison on WG6 (33.38 m.) MESH5-8", color = "black", face = "bold", size = 12))

plot1 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + geom_line(data=MESH1,aes(x=Time2, y=gauge_38, color="MESH1")) +
  geom_line(data=MESH5,aes(x=Time2, y=gauge_38, color="MESH5")) + scale_color_manual(values = c('Experimental' = 'black','MESH1' = 'red', 'MESH5' = 'cornflowerblue')) + 
  xlim(40,80)
plot2 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + geom_line(data=MESH2,aes(x=Time2, y=gauge_38, color="MESH2")) +
  geom_line(data=MESH6,aes(x=Time2, y=gauge_38, color="MESH6")) + scale_color_manual(values = c('Experimental' = 'black','MESH2' = 'red','MESH6' = 'cornflowerblue')) + 
  xlim(40,80)
plot3 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + geom_line(data=MESH3,aes(x=Time2, y=gauge_38, color="MESH3")) +
  geom_line(data=MESH7,aes(x=Time2, y=gauge_38, color="MESH7")) + scale_color_manual(values = c('Experimental' = 'black','MESH3' = 'red','MESH7' = 'cornflowerblue')) + 
  xlim(40,80)
plot4 <- ggplot() + geom_line(data=a,aes(x=Time, y=WG6, color="Experimental")) + geom_line(data=MESH4,aes(x=Time2, y=gauge_38, color="MESH4")) +
  geom_line(data=MESH8,aes(x=Time2, y=gauge_38, color="MESH8")) + scale_color_manual(values = c('Experimental' = 'black','MESH4' = 'red','MESH8' = 'cornflowerblue')) + 
  xlim(40,80)

figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)

annotate_figure(figure,top = text_grob("Water level comparison on WG6 (33.38 m.) Pairs", color = "black", face = "bold", size = 12))

filenames = "/home/george/OpenFOAM/george-v1912/run/R14NA_MESH1/postProcessing/overtopping/0/overtopping.dat"
tx  <- readLines(filenames)
tx2  <- gsub(pattern = '\\(', replace = " ", x = tx)
tx3  <- gsub(pattern = '\\)', replace = " ", x = tx2)
writeLines(tx3, con=filenames)
QMESH1 = read.table(filenames[1],header=FALSE, skip = 1)
colnames(QMESH1) <- c("Time","Qx","Qy","Qz")

temp <- QMESH1[,c(2,3)]**2
temp <- apply(temp, MARGIN = 1, sum)
temp <- sqrt(temp)/0.8
steps <- diff(QMESH1$Time, lag = 1)
temp <- cumsum(temp[2:length(temp)]*steps)
VMESH1 <- data.frame(temp)
colnames(VMESH1) <- "V"
VMESH1$Time <- QMESH1[2:nrow(QMESH1),1] + 16

filenames = "/home/george/OpenFOAM/george-v1912/run/R14NA_MESH4/postProcessing/overtopping/0/overtopping.dat"
tx  <- readLines(filenames)
tx2  <- gsub(pattern = '\\(', replace = " ", x = tx)
tx3  <- gsub(pattern = '\\)', replace = " ", x = tx2)
writeLines(tx3, con=filenames)
QMESH2 = read.table(filenames[1],header=FALSE, skip = 1)
colnames(QMESH2) <- c("Time","Qx","Qy","Qz")

temp <- QMESH2[,c(2,3)]**2
temp <- apply(temp, MARGIN = 1, sum)
temp <- sqrt(temp)/0.8
steps <- diff(QMESH2$Time, lag = 1)
temp <- cumsum(temp[2:length(temp)]*steps)
VMESH4 <- data.frame(temp)
colnames(VMESH4) <- "V"
VMESH4$Time <- QMESH2[2:nrow(QMESH2),1] + 16

benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200825_S1_WG_R14NA.ASC"
BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
Signal <- BenchRaw[,c(10)]
box_type = 2
V <- calc_overtopping(Signal, box_type)

#Calculate the measurement time step
step = V$Time[2] - V$Time[1] 
Vm <- data.frame(rollmean(V, 1000))
Vmm <- data.frame(rollmedian(V, 1000))
der <- data.frame(diff(Vm$V, lag =1)/step)
der$Time <- Vm[2:nrow(Vm),1] 
names(der)[1]<- "der"
ggplot() + geom_line( data = Vmm, aes(x=Time, y=V, color = "Experimental")) + geom_line(data = VMESH1, aes(x=Time, y=V, color = "MESH1")) +
  geom_line(data = VMESH4, aes(x=Time, y=V, color = "MESH4")) + ggtitle("Cumulative Overtopping") +
  scale_color_manual(values = c('Experimental' = 'black','MESH1' = 'red', 'MESH4' = 'cornflowerblue'))
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


                        

