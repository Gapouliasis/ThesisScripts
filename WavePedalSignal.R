library(ggplot2)
library(seewave)
library(reshape2)
require(signal)
filenames = Sys.glob("/home/george/Thesis/tsosMi/Shape 1/*FP*.ASC")
filename <- filenames[4]
Meas = read.table(filename,header=TRUE, sep = ";",skip = 6)
print(filename)
vars <- colnames(Meas)

wavepedal <-Meas[,c(1,10)]
rm(Meas)
names(wavepedal)[1]<-"Time"
names(wavepedal)[2]<-"Pedal"
wavepedal$Pedal <- wavepedal$Pedal*10.16/100
#ggplot() + geom_line(data=wavepedal,aes(x=Time,y=Pedal))
temp <- bwfilter(wavepedal$Pedal, 2000, n =2, to = 5)
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
pedalflux$FPedal09 <- temp

temp <- bwfilter(wavepedal$Pedal, 2000, n =2, to = 10)
wavepedal$FPedal05 <- temp

#temp <- diff(temp, lag =1)/0.0005
temp2 <- temp[5:nrow(temp)]
temp1 <- temp[4:(nrow(temp)-1)]
tempb1 <- temp[2:(nrow(temp)-3)]
tempb2 <- temp[1:(nrow(temp)-4)]
temp <- (-temp2 + 8*temp1 - 8*tempb1 + tempb2)/(12*0.0005) #Fourth order central difference
pedalflux$FPedal05 <- temp

temp <- diff(wavepedal$Pedal, lag = 1)/0.0005
#pedalflux$Pedal <- temp

longmeas <- melt(wavepedal[,], id.vars="Time")
ggplot(longmeas, aes(Time,value, col=variable)) + geom_line()

#longmeas <- melt(pedalflux, id.vars="Time")
#ggplot(longmeas, aes(Time,value, col=variable)) + geom_line()
zero <- 30000
ggplot() + geom_line(data=pedalflux[zero:160000, ],aes(x=Time,y=FPedal05, color = "Alt")) +
  geom_line(data=pedalflux[zero:160000, ],aes(x=Time,y=FPedal09, color = "Current"))
pedalflux$Time <- pedalflux$Time - pedalflux$Time[zero]
pedalflux <- pedalflux[zero:nrow(pedalflux),]
longmeas <- melt(pedalflux, id.vars="Time")
#ggplot(longmeas, aes(Time,value, col=variable)) + geom_line()
out <- pedalflux[c("Time","FPedal09")]

write.table(out,file = "/home/george/OpenFOAM/george-v1912/run/R14NA_MESH11/R14NA_PaddleSignal_5.inp",
            row.names = FALSE, col.names = FALSE)

# a <- BenchCalib[80001:limit, ]
# b <- Meas[1700:nrow(Meas), ]
# exp_spec <- meanspec(wavepedal$Pedal, f = 2000, wl = 1024, PSD = TRUE, plot = FALSE, FUN = "mean")
# spectra <- as.data.frame(exp_spec)
# ggplot() + geom_line(data = spectra, aes(x=x ,y=y))
# names(spectra)[2] <- "Exp"
# num_spec <- meanspec(b$gauge_38, f = 100, wl = 1024, PSD = TRUE, plot = FALSE, FUN = "mean")
# temp <- as.data.frame(num_spec)
# spectra$Num <- temp$y
# 
# longmeas <- melt(spectra, id.vars="x")
# ggplot(longmeas, aes(x,value, col=variable)) + geom_line()
