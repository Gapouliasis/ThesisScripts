#Load libraries----------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(seewave)
library(reshape2)
library(signal)
#Load and calibrate the wave pedal data----------------------------------------------------------------------------------------------
filenames = Sys.glob("/home/george/Thesis/tsosMi/Shape 1/*FP*NA.ASC")
filename <- filenames[1]
# filename <- filenames[26]
Meas = read.table(filename,header=TRUE, sep = ";",skip = 6)
print(filename)
vars <- colnames(Meas)
wavepedal <-Meas[,c(1,10)]
rm(Meas)
colnames(wavepedal) <- c("Time","Pedal")
wavepedal$Pedal <- wavepedal$Pedal*10.16/100
wavepedal$Pedal <- wavepedal$Pedal - mean(wavepedal[1:5000,2])
#Filter data and calculate pedal velocity--------------------------------------------------------------------------------------------
temp <- bwfilter(wavepedal$Pedal, 2000, n =2, to = 2) #Apply Butterworth low pass filter
wavepedal$FPedal1 <- temp

#temp <- diff(temp, lag =1)/0.0005 #Forward difference
temp2 <- temp[5:nrow(temp)]
temp1 <- temp[4:(nrow(temp)-1)]
tempb1 <- temp[2:(nrow(temp)-3)]
tempb2 <- temp[1:(nrow(temp)-4)]
temp <- (-temp2 + 8*temp1 - 8*tempb1 + tempb2)/(12*0.0005) #Fourth order central difference
pedalflux <- wavepedal[3:(nrow(wavepedal)-2),1]
pedalflux <- as.data.frame(pedalflux)
names(pedalflux)[1]<-"Time"
equ_time <- seq(from = 1, to = 100 , by = 0.01)
temp_equ <- approx(pedalflux$Time,temp, xout = equ_time) #Linearly interpolate data in equidistant axis
#temp_equ <- temp_equ[30000,200000]
pedals <- as.data.frame(temp_equ)
names(pedals)[1] <- "Time"
names(pedals)[2] <- substr(filename,41,55)
i=3
for (filename in filenames[2:28]){
  Meas = read.table(filename,header=TRUE, sep = ";",skip = 6)
  print(filename)
  vars <- colnames(Meas)
  wavepedal <-Meas[,c(1,10)]
  rm(Meas)
  colnames(wavepedal) <- c("Time","Pedal")
  wavepedal$Pedal <- wavepedal$Pedal*10.16/100
  wavepedal$Pedal <- wavepedal$Pedal - mean(wavepedal[1:5000,2])
  #Filter data and calculate pedal velocity--------------------------------------------------------------------------------------------
  temp <- bwfilter(wavepedal$Pedal, 2000, n =2, to = 2) #Apply Butterworth low pass filter
  wavepedal$FPedal1 <- temp
  
  #temp <- diff(temp, lag =1)/0.0005 #Forward difference
  temp2 <- temp[5:nrow(temp)]
  temp1 <- temp[4:(nrow(temp)-1)]
  tempb1 <- temp[2:(nrow(temp)-3)]
  tempb2 <- temp[1:(nrow(temp)-4)]
  temp <- (-temp2 + 8*temp1 - 8*tempb1 + tempb2)/(12*0.0005) #Fourth order central difference
  pedalflux <- wavepedal[3:(nrow(wavepedal)-2),1]
  pedalflux <- as.data.frame(pedalflux)
  names(pedalflux)[1]<-"Time"
  equ_time <- seq(from = 1, to = 100 , by = 0.01)
  temp_equ <- approx(pedalflux$Time,temp, xout = equ_time) #Linearly interpolate data in equidistant axis
  #temp_equ <- temp_equ[30000,200000]
  pedals$temp <- temp_equ$y
  names(pedals)[i] <- substr(filename,41,55)
  i <- i + 1
}

longmeas <- melt(pedals[,1:10], id.vars="Time")
ggplot(longmeas, aes(Time,value, col=variable)) + geom_line() + xlim(20,60)

plot(pedals[,1], pedals[,1])

write.table(out,file = "/home/george/OpenFOAM/george-v1912/run/R24NA_PaddleSignal_FP1inp",
            row.names = FALSE, col.names = FALSE)