from TI_slab_code.TI_slab_tunnel import *

from Bi2Se3Parameters import *
import matplotlib.pyplot
from KGM import *
import math
import argparse as ag

parser=ag.ArgumentParser()   
parser.add_argument("--test",nargs='?',default=False, const=True)
parser.add_argument("--save",nargs='?',default=False, const=True)
parser.add_argument("--band",nargs='?',default=False, const=True)
parser.add_argument("--wave",nargs='?',default=False, const=True)
parser.add_argument("--log",nargs='?',default=False, const=True)
parser.add_argument("--current",nargs='?',default=False, const=True)
parser.add_argument("--both",nargs='?',default=False, const=True)	
parser.add_argument("--nbands",default=12,type=int)
parser.add_argument("--nstates",default=4,type=int)
parser.add_argument("--ndiv",default=120,type=int)#increasing this decreases the total magnitude, at least up to a point, and decreases the size of the flat region seems to have leveled out by about 600 divisions
parser.add_argument("--outputt",default=1,type=int)
parser.add_argument("--outputb",default=1,type=int)
parser.add_argument("--wavestate",default=0,type=int)
parser.add_argument("--vacwidtha",default=6,type=float)
parser.add_argument("--width",default=12,type=float)
parser.add_argument("--dx",default=0,type=float)#could be changed to xwidth without changing units
parser.add_argument("--kx",default=0.048952909788227444,type=float) 
parser.add_argument("--ky",default=0,type=float)
parser.add_argument("--vacpotential",default=5,type=float)
parser.add_argument("--dflt",default=0.1,type=float)#change in the fermi level of top and bottom slab for vacuum 
parser.add_argument("--dflb",default=0,type=float)
parser.add_argument("--FermiChangei",default=0,type=float)
parser.add_argument("--FermiChangef",default=0.2,type=float)
parser.add_argument("--points",default=20,type=int)
parser.add_argument("--BandGap",default='0.28',type=float)
args=parser.parse_args()

#centre_vacuum,constants,kpoints,total_thickness,vacuum_thickness,nz,sigma)
kx=get_kx()
ky=get_ky()
k=[]

total_size=args.width#12.0

vacuum_size=args.vacwidtha#6.0

nbands=args.nbands

take_bands=args.nstates

#matparam['M']=args.BandGap*eV
vacparam['C']=args.vacpotential*eV
s=(1.00+0.16*abs(1-matparam['M']/(0.28*eV)))*0.034*matparam['M']/(0.28*eV)+matparam['delta']/2
nz=int(args.ndiv*total_size//12.0)
pi=math.pi
k.append(kx/6)
k.append(ky/6)
#print(k)
nm=1*10**-9
constants={'kT':0.025,'eV':eV}
print('nz:',nz)

print('self_dz'+str(total_size/nz))
T=TI_tunnel(matparam,matparam2,vacparam,cent_vacc,constants,k,total_size*nm,vacuum_size*nm,nz,s*eV,nbands,take_bands,args.both,args.test)
if args.band:
	E=T.Bandstructure()
	#print(np.shape(E))
	for i in range(len(E[0,:])):
		plt.plot(E[:,i].real/eV)
	plt.show()
	plt.close()

if args.wave:

	if args.both:
		zNT=T.Wave_func('top',args.kx*10**10,args.ky*10**10,args.outputt).reshape(-1,4)#reshape(int(nz*4*1.5+4*(nz-40)/20+4)//4,4)
		zNB=T.Wave_func('bottom',0.048952909788227444*10**10,args.ky*10**10,args.outputb).reshape(-1,4)#.reshape(int(nz*4*1.5+4*(nz-40)/20+4)//4,4)
	else:
		zNT=T.Wave_func('top',args.kx*10**10,args.ky*10**10,args.outputt).reshape(int(nz),4)
		zNT=(zNT.imag)**2+(zNT.real)**2
		zNB=T.Wave_func('bottom',args.kx*10**10,args.ky*10**10,args.outputb).reshape(int(nz),4)
		zNB=(zNB.imag)**2+(zNB.real)**2

		zNT2=T.Wave_func('top',args.kx*10**10,args.ky*10**10,args.outputt+1).reshape(int(nz),4)
		zNT2=(zNT2.imag)**2+(zNT2.real)**2
		zNB2=T.Wave_func('bottom',args.kx*10**10,args.ky*10**10,args.outputb+1).reshape(int(nz),4)
		zNB2=(zNB2.imag)**2+(zNB2.real)**2

	#plt.spy(T.H_slab('bottom',143595130.8987827,0))
	#plt.show()
	#plt.close()
	print('Overlap:',np.sum(np.multiply(zNT,zNB)))
	plt.figure(figsize=(3,2))
	if args.outputt!=0:
		for l in range(len(zNT)):
			if zNT[l,0]<10**-10:
				zNT[l,:]=0
		for l in range(len(zNT)):
			if zNT[l,1]<10**-10:
				zNT[l,:]=0
		for l in range(len(zNT)):
			if zNT[l,2]<10**-10:
				zNT[l,:]=0
		for l in range(len(zNT)):
			if zNT[l,3]<10**-10:
				zNT[l,:]=0
		for m in range(len(zNT[0,:])):
			if args.log:
				plt.semilogy(np.linspace(0,total_size,len(zNT[:,0])),zNT[:,m],color='g',linestyle='dashed')
				plt.semilogy(np.linspace(0,total_size,len(zNT[:,0])),zNT2[:,m],color='r',linestyle='dashed')
			else:
				plt.plot(np.linspace(0,total_size,len(zNT[:,0])),zNT[:,m],color='g',linestyle='dashed')
			#plt.plot(zNT[:,1],color='r')
			#plt.plot(zNT[:,2],color='g')
			#plt.plot(zNT[:,3],color='y')
	if args.outputb!=0:
		for l in range(len(zNB)):
			if zNB[l,0]<10**-10:
				zNB[l,:]=0
		for l in range(len(zNB)):
			if zNB[l,1]<10**-10:
				zNB[l,:]=0
		for l in range(len(zNB)):
			if zNB[l,2]<10**-10:
				zNB[l,:]=0
		for l in range(len(zNB)):
			if zNB[l,3]<10**-10:
				zNB[l,:]=0
		for m in range(len(zNB[0,:])):
			if args.log:
				plt.semilogy(np.linspace(0,total_size,len(zNB[:,0])),zNB[:,m],color='g')
				plt.semilogy(np.linspace(0,total_size,len(zNB2[:,0])),zNB[:,m],color='r')
			else:
			#plt.plot(zNB[:,0],color='b')
				plt.plot(np.linspace(0,total_size,len(zNB[:,0])),zNB[:,m],color='g')
			#plt.plot(zNB[:,2],color='g')
			#plt.plot(zNB[:,3],color='y')
	plt.ylabel('$|\Psi|^2$')
	plt.xlabel('X (nm)')
	#plt.legend(['Right','Left'])
	name=str('WaveFunction'+str(args.ndiv)+'['+str(args.outputt)+','+str(args.outputb)+'].pdf')
	plt.savefig(name,format='pdf',dpi=400,bbox_inches='tight',pad_inches=0.2)
	plt.show()
#add for loop to plot all 4 wave states
#plt.close()

#plt.plot(zNT[:,0].real*zNB[:,0].real)

iHf=0
if args.wave:
	for i in range(4):
		iHf=np.absolute((zNB[:,i].T.conj()).dot(zNT[:,i]))**2

		print(iHf)	


#plt.show()
#plt.close()

fig=plt.figure()

Vd=np.linspace(args.FermiChangei,args.FermiChangef,args.points)*eV
Vd1=np.linspace(0,args.dlft-0.01,60)*eV
Vd2=np.linspace(args.dlft-0.01,args.dflt+0.01,3)*eV
Vd3=np.linspace(args.dflt+0.01,0.20,20)*eV

Vd=np.append(Vd1,Vd2)
Vd=np.append(Vd,Vd3)
Vg=args.dflt*eV#0.1*eV
J=np.absolute(T.Current_calc(Vg,Vd,'top'))
#name=str('J'+args.dflt+args.BandGap)
np.savetxt('J.txt',J)
np.savetxt('V.txt',Vd)
if args.current:

	plt.plot(Vd/eV,(J/10**-6)/10000,'k')
	plt.xlabel('$V_{ds}$ (eV)')
	plt.ylabel('J ($\mu$A/$cm^2$)')
	#plt.ylim((0,10**-6))
	plt.savefig('J.png',format='png',dpi=200)

	plt.show()
	#plt.savefig('J.png',format='png',dpi=200)
