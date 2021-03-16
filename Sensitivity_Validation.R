library(ggplot2)
library(reshape2)
library(seewave)

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH1/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH1 = read.table(filenames[1],header=TRUE)
gauge_pos1 <- MESH1[c(1,2,3),]
MESH1 <- MESH1[-c(1,2,3,4),]
rownames(MESH1) <- NULL

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH1/postProcessing/overtopping/0/overtopping.dat"
temp = read.table(filenames[1],header=FALSE, skip = 1)
Q1 <- data.frame(temp$V1) 
names(Q1)[1] <- "Time"
temp$V2 <- '^'(temp$V2,2)
temp$V3 <- '^'(temp$V3,2)
temp$V4 <- '^'(temp$V4,2)
sum <- temp$V2 + temp$V3 + temp$V3
Q1$Q <- '^'(sum,0.5) 

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH2/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH2 = read.table(filenames[1],header=TRUE)
gauge_pos2 <- MESH2[c(1,2,3),]
MESH2 <- MESH2[-c(1,2,3,4),]
rownames(MESH2) <- NULL

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH2/postProcessing/overtopping/0/overtopping.dat"
temp = read.table(filenames[1],header=FALSE, skip = 1)
Q2 <- data.frame(temp$V1) 
names(Q2)[1] <- "Time"
temp$V2 <- '^'(temp$V2,2)
temp$V3 <- '^'(temp$V3,2)
temp$V4 <- '^'(temp$V4,2)
sum <- temp$V2 + temp$V3 + temp$V3
Q2$Q <- '^'(sum,0.5) 

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH3/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH3 = read.table(filenames[1],header=TRUE)
gauge_pos3 <- MESH3[c(1,2,3),]
MESH3 <- MESH3[-c(1,2,3,4),]
rownames(MESH3) <- NULL

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH3/postProcessing/overtopping/0/overtopping.dat"
temp = read.table(filenames[1],header=FALSE, skip = 1)
Q3 <- data.frame(temp$V1) 
names(Q3)[1] <- "Time"
temp$V2 <- '^'(temp$V2,2)
temp$V3 <- '^'(temp$V3,2)
temp$V4 <- '^'(temp$V4,2)
sum <- temp$V2 + temp$V3 + temp$V3
Q3$Q <- '^'(sum,0.5) 

# filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH4/postProcessing/surfaceElevation/0/surfaceElevation.dat"
# MESH4 = read.table(filenames[1],header=TRUE)
# gauge_pos4 <- MESH4[c(1,2,3),]
# MESH4 <- MESH4[-c(1,2,3,4),]
# rownames(MESH4) <- NULL
# filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH4/postProcessing/overtopping/0/overtopping.dat"
# temp = read.table(filenames[1],header=FALSE, skip = 1)
# Q4 <- data.frame(temp$V1) 
# names(Q4)[1] <- "Time"
# temp$V2 <- '^'(temp$V2,2)
# temp$V3 <- '^'(temp$V3,2)
# temp$V4 <- '^'(temp$V4,2)
# sum <- temp$V2 + temp$V3 + temp$V3
# Q4$Q <- '^'(sum,0.5) 

#longmeas <- melt(Meas, id.vars="Time")
#ggplot(longmeas, aes(Time,value, col=variable)) + geom_line()

#ggplot(longmeas, aes(Time,value)) + geom_line() + facet_wrap(~variable)

#rm(longmeas)


ggplot() + geom_line(data=MESH1,aes(x=Time, y=gauge_10), color="red") + geom_line(data=MESH2,aes(x=Time, y=gauge_10), color="black") +
  geom_line(data=MESH3,aes(x=Time, y=gauge_10), color="blue") + ggtitle("gauge_10 (26.405 m.)")

ggplot() + geom_line(data=MESH1,aes(x=Time, y=gauge_20), color="red") + geom_line(data=MESH2,aes(x=Time, y=gauge_20), color="black") +
  geom_line(data=MESH3,aes(x=Time, y=gauge_20), color="blue") + ggtitle("gauge_20 (28.81 m.)")

ggplot() + geom_line(data=MESH1,aes(x=Time, y=gauge_30), color="red") + geom_line(data=MESH2,aes(x=Time, y=gauge_30), color="black") +
  geom_line(data=MESH3,aes(x=Time, y=gauge_30), color="blue") + ggtitle("gauge_30 (31.22 m.)")

ggplot() + geom_line(data=MESH1,aes(x=Time, y=gauge_39), color="red") + geom_line(data=MESH2,aes(x=Time, y=gauge_39), color="black") +
  geom_line(data=MESH3,aes(x=Time, y=gauge_39), color="blue") + ggtitle("gauge_39 (33.38 m.)")


ggplot() + geom_line(data=Q1,aes(x=Time, y=Q), color="red") + geom_line(data=Q2,aes(x=Time, y=Q), color="black")
 + geom_line(data=Q3,aes(x=Time, y=Q), color="blue") + geom_line(data=Q4,aes(x=Time, y=Q), color="green")






