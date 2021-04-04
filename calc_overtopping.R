calc_overtopping <- function(Signal, box_type)
# Function translate WG9 signal to overtopping, for box 2 & 3. WG8 and
# results and box=1 results should not be used here. 
# Function calculates volume of water in the tank/meter of mouth width. 
# Number of box (G column of data code) 
#Modified from D. Dermentzoglou
{
if (box_type==2) {
  mouth = 0.28 
} else if (box_type==3) {
  mouth = 0.112}

#Calibration function of WG9 
cald = - 4.38/100

#Calculate water depth in box
d = 4.4/100 + (Signal - 9.86172602739726) * cald

#Apply equations for different  water depths
f1 = which(d<0.47)
f2 = which(d>=0.47 & d<0.48)
f3 = which(d>=0.48)
v = zeros(length(d),1)
v[f1] = 0.216*0.128*d[f1]
v[f2] = 0.296*0.128*(d[f2]-0.47)+0.47*0.216*0.128
v[f3] = 0.47*0.216*0.128+0.01*0.296*0.128+(0.19097*(d[f3]-0.48)^2+(d[f3]-0.48)*0.52)*0.28
v=v-mean(v[1:2000])
q=v/mouth

V <- data.frame(BenchRaw$Measurement.time.s.)
names(V)[1] <- "Time"
V$V <- q #Cumulative overtopping volume, in m^3/m
return(V)
}