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
temp <- bwfilter(wavepedal$Pedal, 2000, n =2, to = 10)
wavepedal$FPedal1 <- temp

temp <- diff(temp, lag =1)/0.0005
pedalflux <- wavepedal[-1,1]
pedalflux <- as.data.frame(pedalflux)
names(pedalflux)[1]<-"Time"
pedalflux$FPedal09 <- temp

temp <- bwfilter(wavepedal$Pedal, 2000, n =2, to = 25)
wavepedal$FPedal05 <- temp

temp <- diff(temp, lag =1)/0.0005
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

write.table(out,file = "/home/george/OpenFOAM/george-v1912/run/R14NA_MESH1/R14NA_PaddleSignal.inp",
            row.names = FALSE, col.names = FALSE)
