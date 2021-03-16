library(ggplot2)
library(reshape2)
library(seewave)
library(ggpubr)
filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH1/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH1 = read.table(filenames[1],header=TRUE)
gauge_pos1 <- MESH1[c(1,2,3),]
MESH1 <- MESH1[-c(1,2,3,4),]
rownames(MESH1) <- NULL

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH2/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH2 = read.table(filenames[1],header=TRUE)
gauge_pos2 <- MESH2[c(1,2,3),]
MESH2 <- MESH2[-c(1,2,3,4),]
rownames(MESH2) <- NULL

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH3/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH3 = read.table(filenames[1],header=TRUE)
gauge_pos3 <- MESH3[c(1,2,3),]
MESH3 <- MESH3[-c(1,2,3,4),]
rownames(MESH3) <- NULL

# filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH4/postProcessing/surfaceElevation/0/surfaceElevation.dat"
# MESH4 = read.table(filenames[1],header=TRUE)
# gauge_pos4 <- MESH4[c(1,2,3),]
# MESH4 <- MESH4[-c(1,2,3,4),]
# rownames(MESH4) <- NULL

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH5/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH5 = read.table(filenames[1],header=TRUE)
gauge_pos1 <- MESH5[c(1,2,3),]
MESH5 <- MESH5[-c(1,2,3,4),]
rownames(MESH5) <- NULL

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH6/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH6 = read.table(filenames[1],header=TRUE)
gauge_pos2 <- MESH6[c(1,2,3),]
MESH6 <- MESH6[-c(1,2,3,4),]
rownames(MESH6) <- NULL

filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH7/postProcessing/surfaceElevation/0/surfaceElevation.dat"
MESH7 = read.table(filenames[1],header=TRUE)
gauge_pos3 <- MESH7[c(1,2,3),]
MESH7 <- MESH7[-c(1,2,3,4),]
rownames(MESH7) <- NULL

# filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH8/postProcessing/surfaceElevation/0/surfaceElevation.dat"
# MESH8 = read.table(filenames[1],header=TRUE)
# gauge_pos4 <- MESH8[c(1,2,3),]
# MESH8 <- MESH8[-c(1,2,3,4),]
# rownames(MESH8) <- NULL
# filenames = "/home/george/OpenFOAM/george-v1912/run/R17NA_MESH6/postProcessing/overtopping/0/overtopping.dat"
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

plot1 <- ggplot() + geom_line(data = MESH1,aes(x=Time, y=gauge_10, color = "Coarse")) + geom_line(data = MESH5,aes(x=Time, y=gauge_10, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_10 (26.405 m.)")

plot2 <- ggplot() + geom_line(data = MESH1,aes(x=Time, y=gauge_20, color = "Coarse")) + geom_line(data = MESH5,aes(x=Time, y=gauge_20, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_20 (28.81 m.)")

plot3 <- ggplot() + geom_line(data = MESH1,aes(x=Time, y=gauge_30, color = "Coarse")) + geom_line(data = MESH5,aes(x=Time, y=gauge_30, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_30 (31.22 m.)")

plot4 <- ggplot() + geom_line(data = MESH1,aes(x=Time, y=gauge_39, color = "Coarse")) + geom_line(data = MESH5,aes(x=Time, y=gauge_39, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_39 (33.38 m.)")

figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)

annotate_figure(figure,top = text_grob("MESH 1", color = "black", face = "bold", size = 12))

plot1 <- ggplot() + geom_line(data = MESH2,aes(x=Time, y=gauge_10, color = "Coarse")) + geom_line(data = MESH6,aes(x=Time, y=gauge_10, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_10 (26.405 m.)")

plot2 <- ggplot() + geom_line(data = MESH2,aes(x=Time, y=gauge_20, color = "Coarse")) + geom_line(data = MESH6,aes(x=Time, y=gauge_20, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_20 (28.81 m.)")

plot3 <- ggplot() + geom_line(data = MESH2,aes(x=Time, y=gauge_30, color = "Coarse")) + geom_line(data = MESH6,aes(x=Time, y=gauge_30, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_30 (31.22 m.)")

plot4 <- ggplot() + geom_line(data = MESH2,aes(x=Time, y=gauge_39, color = "Coarse")) + geom_line(data = MESH6,aes(x=Time, y=gauge_39, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_39 (33.38 m.)")

figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)

annotate_figure(figure,top = text_grob("MESH 2", color = "black", face = "bold", size = 12))

plot1 <- ggplot() + geom_line(data = MESH3,aes(x=Time, y=gauge_10, color = "Coarse")) + geom_line(data = MESH6,aes(x=Time, y=gauge_10, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_10 (26.405 m.)")

plot2 <- ggplot() + geom_line(data = MESH3,aes(x=Time, y=gauge_20, color = "Coarse")) + geom_line(data = MESH6,aes(x=Time, y=gauge_20, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_20 (28.81 m.)")

plot3 <- ggplot() + geom_line(data = MESH3,aes(x=Time, y=gauge_30, color = "Coarse")) + geom_line(data = MESH6,aes(x=Time, y=gauge_30, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_30 (31.22 m.)")

plot4 <- ggplot() + geom_line(data = MESH3,aes(x=Time, y=gauge_39, color = "Coarse")) + geom_line(data = MESH6,aes(x=Time, y=gauge_39, color = "Fine")) +
  scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_39 (33.38 m.)")

figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)

annotate_figure(figure,top = text_grob("MESH 3", color = "black", face = "bold", size = 12))

# plot1 <- ggplot() + geom_line(data = MESH4,aes(x=Time, y=gauge_10, color = "Coarse")) + geom_line(data = MESH8,aes(x=Time, y=gauge_10, color = "Fine")) +
#   scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_10 (26.405 m.)")
# 
# plot2 <- ggplot() + geom_line(data = MESH4,aes(x=Time, y=gauge_20, color = "Coarse")) + geom_line(data = MESH8,aes(x=Time, y=gauge_20, color = "Fine")) +
#   scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_20 (28.81 m.)")
# 
# plot3 <- ggplot() + geom_line(data = MESH4,aes(x=Time, y=gauge_30, color = "Coarse")) + geom_line(data = MESH8,aes(x=Time, y=gauge_30, color = "Fine")) +
#   scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_30 (31.22 m.)")
# 
# plot4 <- ggplot() + geom_line(data = MESH4,aes(x=Time, y=gauge_39, color = "Coarse")) + geom_line(data = MESH8,aes(x=Time, y=gauge_39, color = "Fine")) +
#   scale_color_manual(values = c('Coarse' = 'red','Fine' = 'black')) + labs(color = 'OCW3D Mesh') + ggtitle("gauge_39 (33.38 m.)")
# 
# figure <- ggarrange(plot1, plot2, plot3, plot4, ncol = 2 , nrow = 2)
# 
# annotate_figure(figure,top = text_grob("MESH 4", color = "black", face = "bold", size = 12))




