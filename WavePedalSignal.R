#Load libraries----------------------------------------------------------------------------------------------------------------------
library(ggplot2)
library(seewave)
library(reshape2)
library(signal)
#Load and calibrate the wave pedal data----------------------------------------------------------------------------------------------
filenames = Sys.glob("/home/george/Thesis/tsosMi/Shape 1/*FP*NA.ASC")
filename <- filenames[10]
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
pedalflux$FPedal1 <- temp

temp <- bwfilter(wavepedal$Pedal, 2000, n =2, to = 10) #Apply Butterworth low pass filter
wavepedal$FPedal2 <- temp

#temp <- diff(temp, lag =1)/0.0005
temp2 <- temp[5:nrow(temp)]
temp1 <- temp[4:(nrow(temp)-1)]
tempb1 <- temp[2:(nrow(temp)-3)]
tempb2 <- temp[1:(nrow(temp)-4)]
temp <- (-temp2 + 8*temp1 - 8*tempb1 + tempb2)/(12*0.0005) #Fourth order central difference
pedalflux$FPedal2 <- temp

temp <- sgolayfilt(wavepedal$Pedal) #Apply Savitzky-Golay smoothing filter
wavepedal$FPedal3 <- temp

#temp <- diff(temp, lag =1)/0.0005
temp2 <- temp[5:length(temp)]
temp1 <- temp[4:(length(temp)-1)]
tempb1 <- temp[2:(length(temp)-3)]
tempb2 <- temp[1:(length(temp)-4)]
temp <- (-temp2 + 8*temp1 - 8*tempb1 + tempb2)/(12*0.0005) #Fourth order central difference
pedalflux$FPedal3 <- temp

pedalflux$FPedal4 <- sgolayfilt(pedalflux$FPedal1)

temp <- wavepedal$Pedal
temp2 <- temp[5:length(temp)]
temp1 <- temp[4:(length(temp)-1)]
tempb1 <- temp[2:(length(temp)-3)]
tempb2 <- temp[1:(length(temp)-4)]
temp <- (-temp2 + 8*temp1 - 8*tempb1 + tempb2)/(12*0.0005) #Fourth order central difference
pedalflux$Pedal <- temp

#Plot the raw and filtered data
# longmeas <- melt(wavepedal[,], id.vars="Time")
# ggplot(longmeas, aes(Time,value, col=variable)) + geom_line()
ggplot(data = wavepedal) + geom_line(aes(x = Time, y = Pedal, color = "Experimental")) + 
  geom_line(aes(x = Time, y = FPedal1, color = "Filtered")) + xlim(20,75) + 
  geom_line(data = pedalflux, aes(x = Time, y = FPedal4, color = "Velocity")) + 
  scale_color_manual(values = c("Experimental" = "black", "Filtered" = "darkorange1", "Velocity" = "deepskyblue"))


#Plot the paddle velocity from the filtered data
zero <- 30000
ggplot() + geom_line(data=pedalflux[zero:160000, ],aes(x=Time,y=FPedal2, color = "Alt")) +
  geom_line(data=pedalflux[zero:160000, ],aes(x=Time,y=FPedal4, color = "Current"))
pedalflux$Time <- pedalflux$Time - pedalflux$Time[zero]
pedalflux <- pedalflux[zero:nrow(pedalflux),]
longmeas <- melt(pedalflux, id.vars="Time")
#ggplot(longmeas, aes(Time,value, col=variable)) + geom_line()
out <- pedalflux[c("Time","FPedal4")]
equ <- seq(from = 0, to = 100 , by = 0.00005)
out <- approx(pedalflux$Time,pedalflux$FPedal4, xout = equ)

# #Write pedal velocity to file-------------------------------------------------------------------------------------------------------
# write.table(out,file = "/home/george/OpenFOAM/george-v1912/run/R14NA_PaddleSignal_SG2INT.inp",
#             row.names = FALSE, col.names = FALSE)

#Calculate and plot wave pedal displacement spectra---------------------------------------------------------------------------------
var <- wavepedal$FPedal1
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.0005/Fr
temp$Fr <- Fr*2000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_filt <- temp
var <- wavepedal$Pedal
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.0005/Fr
temp$Fr <- Fr*2000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_raw <- temp
dis_spectra <- temp
# ggplot() + geom_point(data = spectra_raw, aes(x=T, y=Pscaled, color = "Raw")) + 
#   geom_line(data = spectra_filt, aes(x=T, y=Pscaled, color = "Filtered")) + labs(x = "T (sec.)", y = "Pscaled") + 
#   scale_color_manual(values = c('Filtered' = 'red','Raw' = 'black')) + xlim(0,10) + ggtitle("Pedal displacement Spectrum vs Period")

ggplot() + geom_line(data = spectra_raw, aes(x=Fr, y=Pscaled, color = "Raw")) + 
  geom_line(data = spectra_filt, aes(x=Fr, y=Pscaled, color = "Filtered")) + labs(x = "Fr (Hz)", y = "Pscaled") + 
  scale_color_manual(values = c('Filtered' = 'red','Raw' = 'black')) + labs(color = 'Legend') + 
  xlim(0,5) + ggtitle("Pedal displacement Spectrum vs Frequency")

#Calculate and plot wave pedal velocity spectra---------------------------------------------------------------------------------
var <- pedalflux$FPedal4
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.0005/Fr
temp$Fr <- Fr*2000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_filt <- temp
vel_spectra <- temp
var <- pedalflux$Pedal
Pscaled <- Mod(4*fft(var)/length(var))
Fr <- 0:(length(var)-1)/length(var)
temp <- as.data.frame(Pscaled)
temp$T <- 0.0005/Fr
temp$Fr <- Fr*2000
temp <- temp[-1,]
colnames(temp) <- c("Pscaled","T","Fr")
spectra_raw <- temp
# ggplot() + geom_point(data = spectra_raw, aes(x=T, y=Pscaled, color = "Raw")) + 
#   geom_line(data = spectra_filt, aes(x=T, y=Pscaled, color = "Filtered")) + labs(x = "T (sec.)", y = "Pscaled") + 
#   scale_color_manual(values = c('Filtered' = 'red','Raw' = 'black')) + xlim(0,10) + ggtitle("Pedal displacement Spectrum vs Period")

ggplot() + geom_line(data = spectra_raw, aes(x=Fr, y=Pscaled, color = "Raw")) + 
  geom_line(data = spectra_filt, aes(x=Fr, y=Pscaled, color = "Filtered")) + labs(x = "Fr (Hz)", y = "Pscaled") + 
  scale_color_manual(values = c('Filtered' = 'red','Raw' = 'black')) + labs(color = 'Legend') + 
  xlim(0,5) + ggtitle("Pedal velocity Spectrum vs Frequency")

#Plot the paddle displacement and velocity spectra--------------------------------------------------------------------------------- 
ggplot() + geom_line(data = dis_spectra, aes(x=Fr, y=Pscaled, color = "Displacement")) + 
  geom_line(data = vel_spectra, aes(x=Fr, y=Pscaled, color = "Velocity")) + labs(x = "Fr (Hz)", y = "Pscaled") + 
  scale_color_manual(values = c('Velocity' = 'red','Displacement' = 'black')) + labs(color = 'Legend') + 
  xlim(0,5) + ggtitle("Pedal displacement velocity Spectrum vs Frequency")
