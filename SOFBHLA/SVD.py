## Running from the same server to avoid EPICS ARRARY SIZE issus #
## export EPICS_CA_MAX_ARRAY_BYTES=9000000
# reduce(min, ormXXRV.flat) for min in the array

############  SVD Calculation ###############
import sys
import time 
from numpy import *
import numpy as np
import cothread
from cothread.catools import caput, caget
import pickle


#Adjust options:
H_plane_only = 0
H_plane_disable_one_family = 0

V_plane_only = 0
V_plane_disable_one_family = 0



alphaPV= 'SR-SOFB{}alpha'
eigenValPV= 'SR-SOFB{}eigenValues'
alpha=caget(alphaPV)
alphaX=caget(alphaPV)
alphaY=caget(alphaPV)

#alpha=.005  #Beamline see noises
#alpha=.5  

print "0000000000000000"
print "alpha:"
print alpha

N=359
rowNumber=N+1;
colNumber=N;
index=0;

### Functions ###
def sv(filename,x):
	ff=open(filename,'w')
	pickle.dump(x,ff)
	ff.close()

def rl(filename):
	ff=open(filename)
	xx=pickle.load(ff)
	ff.close()
	return xx

def reverseR(R,row,col,alpha):
	U,s,V = linalg.svd(R,full_matrices=True) # SVD decomposition of R
	L=len(s)
	print L
	S = floor(random.rand(row+1,col+1)*0)
	S1 = floor(random.rand(row+1,col+1)*0)
	for i in range(0,L-1):
		S[i,i]=s[i]
		S1[i,i]= s[i]/(s[i]*s[i]+alpha)	
	SV=transpose(S1)
	RV=dot(dot(transpose(V),SV),transpose(U))
	return RV

def svdS(R,row,col,alpha):
	U,s,V = linalg.svd(R,full_matrices=True) # SVD decomposition of R
	L=len(s)
	print L
	S = floor(random.rand(row+1,col+1)*0)
	S1 = floor(random.rand(row+1,col+1)*0)
	for i in range(0,L-1):
		S[i,i]=s[i]
		S1[i,i]= s[i]/(s[i]*s[i]+alpha)	
	SV=transpose(S1)
	return S

def svdSV(R,row,col,alpha):
	U,s,V = linalg.svd(R,full_matrices=True) # SVD decomposition of R
	L=len(s)
	print L
	S = floor(random.rand(row+1,col+1)*0)
	S1 = floor(random.rand(row+1,col+1)*0)
	for i in range(0,L-1):
		S[i,i]=s[i]
		S1[i,i]= s[i]/(s[i]*s[i]+alpha)	
	SV=transpose(S1)
	return SV


### Input: ORM, SOFB orbit ###
### Units: urad, um ###
ORMPV= 'SR-SOFB{BPM}sofb:R'
#ormPV1d = caget(ORMPV)

SPWeightPV= 'SR-SOFB{BPM}sofb:SP_Weight'
SPWeight = caget(SPWeightPV)


#SAVE R
#sv("R_matrix_1010Inj_1011",ormPV1d);

#RESTORE R
ormPV1d = rl("R_matrix_1010Inj_1011")

print "++++++++++++++"
print ormPV1d

#Disbale/enable H/V plane
if(H_plane_only==1):
	for i in range(0,129600):
		if(i>64799):
			ormPV1d[i]=0
		elif((i%360)>179):
			ormPV1d[i]=0
if(V_plane_only==1):
	for i in range(0,129600):
		if(i<64800):
			ormPV1d[i]=0
		elif((i%360)<180):
			ormPV1d[i]=0

#Disbale/enable H/V corrector family
ormPV1d= ormPV1d.reshape(360,360)
if(H_plane_disable_one_family==1):
	for i in range(0,180):
		if( ((i%12)==1) or ((i%12)==11)):
			ormPV1d[:,i] = 0

if(V_plane_disable_one_family==1):
	for i in range(180,360):
		if( ((i%12)==1) or ((i%12)==11)):
			ormPV1d[:,i] = 0
		if( ((i%12)==3) or ((i%12)==10)):
			ormPV1d[:,i] = 0


ormPV1d= ormPV1d.flatten()

ormPV1d_361 = np.zeros(129960)

for i in range(0,129960):
	if(i<129600):
		ormPV1d_361[i] = ormPV1d[i]
	else:
		ormPV1d_361[i] = SPWeight
print "++++++++++++++"
print ormPV1d_361


caput(ORMPV,ormPV1d_361) 
R = np.reshape(ormPV1d_361,(-1,360))

RV = floor(random.rand(colNumber,rowNumber)*0)
RV = reverseR(R,rowNumber,colNumber,alpha)

print "----------------- (RV*R)  ----------------"
print "XXXXXXXXXXXXXX  YYYYYYYYYYYYYYYYYYYYYYYYYYY"
print dot(RV,R)

RV1PV=RV.flatten()

RV_PV='SR-SOFB{BPM}sofb:RV'
caput(RV_PV, RV1PV,wait=False)






