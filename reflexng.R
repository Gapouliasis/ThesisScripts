#function [F,Ar,Ai,Etai,Ampli,Etar,Amplr] = reflexng(Eta,dt,wgi,wgx,wgd,g,fn,fx,x_Sign)
#    reflexion analysis from n wave gauges
#    Given
# 	Eta	= gauge elevation matrix, gauge are columns
# 	dt		= time resolution
# 	wgi	= columns where used gauges are stored
# 	wgx	= position of gauges, increasing in incident direction considering
# 	0 for the first one
# 	wgd	= depth at wave gauges
# 	g		= gravity acceleration
# 	fn		= minimum frequency of interest, default: frequency resolution
# 	fx		= maximum frequency of interest, default: Nyquist frequency
#    x_Sign  = coordinate for the output 
# 
#    Evaluates
# 	F		= frequency values
# 	Ar		= reflected complex amplitude, phase is referred to x=0
# 	Ai		= incident complex amplitude, phase is referred to x=0
# 	Etai	= incident_wave water_elevation history at x=0
# 	Ampli	= incident_wave envelope amplitude history at x=0
# 	Etar	= reflected_wave water_elevation history at x=0
# 	Amplr	= reflected_wave envelope amplitude history at x=0

benchfile = "/home/george/Thesis/tsosMi/Shape 1/20200825_S1_WG_R14NA.ASC"
BenchRaw = read.table(benchfile,header = TRUE, sep = ";", skip = 6)
Eta <- BenchRaw[,c(2:6)]
dt <- BenchRaw[2,1] - BenchRaw[1,1] 
ng <- ncol(Eta)
le <- nrow(Eta)
df <- 1/(dt*le)
fn <- 1
fx <- 4000
x_Sign <- 33.3
wgx <- 

jn=max(1,round(fn/df)) 
fmin=jn*df 
jx=min(round(fx/df),le/2) 
fmax=jx*df 

# if (isfinite(fn)){ 
#   jn=max(1,round(fn/df)) 
#   fmin=jn*df 
# } else { 
#   jn=1 
#   fmin=df 
# }

# if (isfinite(fx)){ 
#   jx=min(round(fx/df),le/2) 
#   fmax=jx*df 
# }else {
#   jx=le/2 
#   fmax=1/(2*dt) 
# }

print(jn, jx)
Fr <- seq(from = fmin, to = fmax, by = df)
#Fr <- [fmin:df:fmax] #Needs transposition
lf=length(Fr)

# Thp1=phas1(wgx,Fr,wgd,g,x_Sign)
# %PHAS1	evaluate PHASe differences between wave-gauges and the 1-st one
# % input variables:
#   %	X = wave gauge positions raw vector (1,l)
# %	F = frequency column vector (n,1)
# %	d = water depth
# %	g = gravity
# % output variables
# %	Thp1	phase differences
# %	K	wave numbers
poly_coef <- c(0.00011,0.00039,0.00171,0.00654,0.02174,0.06320,0.16084,0.35550,0.66667,1)
D <- 0.5
Kod=((pi+pi)*Fr)^2*(D/g)
K=sqrt(Kod.*(Kod+polyval(poly_coef,Kod)^(-1)))*diag(d.^(-1))
Thp1=K*(X-x_Sign)

W=zeros(lf,ng)
for (i in 1:ng){
	for (j in 1:ng){
		W[,i]=W[,i]+sin(Thp1[,i]-Thp1[,j])^2./(1+(pi)^(-2)*(Thp1[,i]-Thp1[,j])^2)
}
}

D=zeros(lf,1)
for (i in 1:ng){
	for (j in 1:i){
		D=D+4*W[,i].*W[,j]*(sin(Thp1[,i]-Thp1[,j])^2)
	}
}

im=sqrt(-1)
C=zeros(lf,ng)
for (i in 1:ng){
  for (j in 1:ng){
		C[,i]=C[,i]+W[,j]*sin(Thp1[,i]-Thp1[,j])*exp(im*Thp1[,j])
  }	
}	
C=(W*C)
for (j in 1:lf){
	if (D(j)>0){
	  C[j,]=C[j,]*(2*im/D[j]) 
  }else{ 
    C[j,]=0 
}
}

#A=(2/le)*fft(inpaint_nans(Eta(:,wgi),3))
A=(2/le)*fft(Eta[,wgi])
Ar=zeros(lf,1)
Ai=zeros(lf,1)
for (j in 1:lf){
	Ar[j]=A[jn+j,]*transpose(C[j,])
  Ai[j]=A[jn+j,]*conj[transpose(C[j,])]
}

remove(A, C, W, Thp1)
remove(df, fmax, fmin, im, i, j, lf, ng) 

if (nargout>3){
	An=zeros(le,1)
	An(jn+1:1:jx+1)=Ai
	An(le+1-jn:-1:le+1-jx)=conj(Ai)
	Etai=(le/2)*real(ifft(An))
	remove(An) 
}

if (nargout>3){
	An=zeros(le,1)
	An(jn+1:1:jx+1)=Ar
	An(le+1-jn:-1:le+1-jx)=conj(Ar)
	Etar=(le/2)*real(ifft(An))
	remove(An) 
}
remove(le,jn,jx)

if (nargout>4){
	Ampli=abs(hilbert(Etai))
}

if (nargout>4){
	Amplr=abs(hilbert(Etar))
}