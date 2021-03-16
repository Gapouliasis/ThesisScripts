library(ggplot2)

Q = read.table("/home/george/OpenFOAM/george-v1912/run/Model2/postProcessing/overtopping/0/overtopping.dat",header=TRUE, skip = 1)
colnames(Q)
names(Q)[1]<- "Time"
names(Q)[2]<- "X"
names(Q)[3]<- "Y"
names(Q)[4]<- "Z"

ggplot() + geom_point(data=Q,aes(x=Time,y=Y))
