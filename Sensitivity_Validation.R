#Load libraries---------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(reshape2)
library(seewave)
library(phonTools)
library(zoo)
library(ggpubr)
library(POT)
#Load user defined functions---------------------------------------------------------------------------------------------------------
cwd <- getwd()
source(file.path(cwd,"calc_overtopping.R"))
#Load numerical overtopping data for comparison--------------------------------------------------------------------------------------
case_folders <- c("R24NA_MESH15_Th","R24NA_MESH15_Th_Dt","R24NA_MESH15_Th_Coarse2","R24NA_MESH15_Th_Fine")
i <- 1
for (case_folder in case_folders){
  filenames = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"postProcessing/overtopping/0/overtopping.dat")
  tx  <- readLines(filenames)
  tx2  <- gsub(pattern = '\\(', replace = " ", x = tx)
  tx3  <- gsub(pattern = '\\)', replace = " ", x = tx2)
  writeLines(tx3, con=filenames)
  assign(case_folder,read.table(filenames[1],header=FALSE, skip = 1))
  i <- i + 1 
}

#Transform the instantaneous overtopping discharge to cumulative overtopping volume 
i <- 1
for (case_folder in case_folders){
  dat <- eval(parse(text = case_folder))
  temp <- dat[,c(2,3)]
  steps <- diff(dat$V1, lag = 1)
  temp1 <- temp[1:(nrow(temp)-1),1]
  temp2 <- temp[2:nrow(temp),1]
  V <- cumsum(0.5*(temp1 + temp2)*steps)
  temp1 <- temp[1:(nrow(temp)-1),2]
  temp2 <- temp[2:nrow(temp),2]
  V <- abs(cumsum(0.5*(temp1 + temp2)*steps)) + V
  Time <- dat[2:nrow(dat),1] + 36.5
  assign(sprintf("V%i",i), as.data.frame(cbind(V,Time)))
  i <- i + 1
}

#Load the experimental wave gauge data---------------------------------------------------------------------------------------------
benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200825_S1_WG_R24NA.ASC"
BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
expmeas <-BenchRaw[,c(1,2,3,4,5,6,7,10)]
rm(BenchRaw)

#Calibrate the experimental data
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

#Load and transform the water level to cumulative overtopping volume
BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
Signal <- BenchRaw[,c(10)]
box_type = 2
V <- calc_overtopping(Signal,box_type)

#Smooth the raw cumulative signal
step = V$Time[2] - V$Time[1] 
Vm <- data.frame(rollmean(V, 1000))
Vmm <- data.frame(rollmedian(V, 1000))
der <- data.frame(diff(Vm$V, lag =1)/step)
der$Time <- Vm[2:nrow(Vm),1] 
names(der)[1]<- "der"

#Plot the comparison of experimental and numerical overtopping volumes--------------------------------------------------------------
imark <- floor(seq(from=1,to=dim(V1)[1],length=100))
ggplot() + geom_line( data = V, aes(x=Time, y=V, color = "Experimental")) + 
  #geom_point(data = V1[imark,], aes(x=Time, y=V, color = "CFL=0.50")) +
  geom_line(data = V2, aes(x=Time, y=V, color = "Coarse")) + geom_line(data = V3, aes(x=Time, y=V, color = "Coarse2")) + 
  geom_line(data = V4, aes(x=Time, y=V, color = "Fine")) + geom_line(data = V1, aes(x=Time, y=V, color = "Base")) +
  scale_color_manual(values = c('Base' = 'red','Coarse' = 'blue','Coarse2' = 'chartreuse4',
                                'Fine' = 'darkorange1','Experimental' = 'black')) +
   ggtitle("Cumulative Overtopping Comparison") + ylim(0,0.012) +xlim(50,100)
  

#Calculate and plot relative error of individual overtopping volume----------------------------------------------------------------
i <- 1
for (case_folder in case_folders){
  temp <- eval(parse(text = case_folder))
  dat <- temp[,c(1,2)]
  names(dat)[1] <- "time"
  names(dat)[2] <- "obs"
  b <- clust(dat, u = 0.001, tim.cond = 1, clust.max = FALSE)
  
  temp_clust <- 0
  nclust <- length(b)
  
  for (j in 1:nclust){
    temp <- b[[j]][2,]
    steps <- diff(b[[j]][1,], lag = 1)
    temp1 <- temp[1:(length(temp)-1)]
    temp2 <- temp[2:length(temp)]
    temp <- sum(0.5*(temp1 + temp2)*steps)
    temp_clust[j] <- temp
  }
  
  if (i == 1){
    temp <- zeros(23,1)
    VClust <- data.frame(temp)
    colnames(VClust) <- "V1"
  }
  
  VClust[sprintf("V%i",i)] <- as.data.frame(temp_clust)
  i <- i + 1
}

V12<- abs((VClust$V1 - VClust$V2)/VClust$V1)*100
V13 <- abs((VClust$V1 - VClust$V3)/VClust$V1)*100
V14 <- abs((VClust$V1 - VClust$V4)/VClust$V1)*100

Ver <- as.data.frame(cbind(V12,V13,V14)) 

mer <- melt(Ver)
names(mer)[1] <- "Case"

ggplot(mer, aes(Case,value, fill = Case)) +  geom_boxplot() + geom_jitter(position=position_jitter(0.2)) + 
  labs(x = "Case", y = "Relative Error %") + ggtitle("Relative Error of individual overtopping volume compared to CFL010")
#+ scale_fill_brewer(palette="Dark2") 

#Load the wave gauge data---------------------------------------------------------------------------------------------------------
for (case_folder in case_folders){
  filenames = file.path("/home/george/OpenFOAM/george-v1912/run",case_folder,"postProcessing/surfaceElevation/0/surfaceElevation.dat")
  tx  <- readLines(filenames)
  tx2  <- gsub(pattern = '\\(', replace = " ", x = tx)
  tx3  <- gsub(pattern = '\\)', replace = " ", x = tx2)
  writeLines(tx3, con=filenames)
  assign(case_folder,read.table(filenames[1],header=TRUE))
}

R24NA_MESH15_Th <- R24NA_MESH15_Th[-c(1,2,3,4),]    
R24NA_MESH15_Th_2 <- R24NA_MESH15_Th_2[-c(1,2,3,4),] 
R24NA_MESH15_Th_3 <- R24NA_MESH15_Th_3[-c(1,2,3,4),] 
R24NA_MESH15_Th_4 <- R24NA_MESH15_Th_4[-c(1,2,3,4),] 

#Calculate experimental and numerical spectra--------------------------------------------------------------------------------------
var <- BenchCalib$WG6
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.0005/Fr
temp$Fr <- Fr*2000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_exp <- temp

i <- 1
for (case_folder in case_folders){
  case_folder <- case_folders[1]
  dat <- eval(parse(text = case_folder))
  var <- dat$gauge_48
  Pscaled <- Mod(4*fft(var)/length(var))
  Fr <- 0:(length(var)-1)/length(var)
  temp <- as.data.frame(Pscaled)
  temp$T <- 0.01/Fr
  temp$Fr <- Fr*100
  temp <- temp[-1,]
  colnames(temp) <- c("Pscaled","T","Fr")
  assign(sprintf("spectra_V%i",i),temp)
  i <- i + 1
}

#Plot spectra and water surface elevation------------------------------------------------------------------------------------------
ggplot() + geom_line(data = spectra_exp, aes(x=Fr, y=Pscaled, color = "Experimental")) + 
  geom_line(data = spectra_V1, aes(x=Fr, y=Pscaled, color = "V1")) + 
  geom_line(data = spectra_V2, aes(x=Fr, y=Pscaled, color = "V2")) + geom_line(data = spectra_V3, aes(x=Fr, y=Pscaled, color = "V3")) +
  geom_line(data = spectra_V4, aes(x=Fr, y=Pscaled, color = "V4")) +
  labs(x = "Fr (Hz)", y = "Pscaled") + labs(color = 'Legend') + scale_color_brewer(palette="RdGy")  +  
  xlim(0,3) + ggtitle("Water Elevation Spectrum vs Frequency Comparison (33.33 m.)") + xlim(0,3) 

ggplot() + geom_line(data = BenchCalib, aes(x=Time, y=WG6, color = "Experimental")) + 
  geom_line(data = R24NA_MESH15_Th, aes(x=Time, y=gauge_48, color = "V1")) + 
  geom_line(data = R24NA_MESH15_Th_2, aes(x=Time, y=gauge_48, color = "V2")) + geom_line(data = R24NA_MESH15_Th_3, aes(x=Time, y=gauge_48, color = "V3")) +
  geom_line(data = R24NA_MESH15_Th_4, aes(x=Time, y=gauge_48, color = "V4")) +
  labs(x = "Fr (Hz)", y = "Pscaled") + labs(color = 'Legend') +
  scale_color_manual(values = c('V1' = 'deeppink','V2' = 'darkorchid1','V3' = 'red',
                                'V4' = 'darkorange1','Experimental' = 'black'))  +  
  xlim(40,80) + ggtitle("Water Surface Elevation Comparison (33.33 m.)")

