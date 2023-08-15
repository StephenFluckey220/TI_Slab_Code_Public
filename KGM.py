import numpy as np
import matplotlib.pyplot as plt
import scipy.linalg as la
import cmath
from matplotlib import cm
from matplotlib.ticker import LinearLocator, FormatStrFormatter
import Bi2Se3Parameters as bp
from Bi2Se3Parameters import *



#path K-G-M for x and y 
#kx=np.linspace(-0.1,0.1,2*N)*10**10/meter
#ky=np.linspace(-0.1,0.1,N)*10**10/meter
N=4500	  		
kxKG=np.linspace(0.0,0,1)*10**10/meter
kyKG=np.linspace(0.0,0,1)*10**10/meter
#kyKG=np.linspace(1/3,0,N)*10**10/meter

#kzKG=np.linspace(0,0,N)*10**10/meter
kxGM=np.linspace(0,1,N)*10**10/meter
kyGM=np.linspace(0.00,0,N)*10**10/meter
#kzGM=np.linspace(0,0,N)*10**10/meter
kx=kxGM
ky=kyGM
#kz=np.append(kzKG,kzGM)
kz=0
nseg=2

def get_kx():
	return kx
def get_ky():
	return ky
def get_kz():
	return kz
def get_N():
	return N
def get_nseg():
	return nseg
