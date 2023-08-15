#from mpl_toolkits.mplot3d import Axes3D
import numpy as np
import scipy
import matplotlib.pyplot as plt
import logging
import os
from struct import *
import shutil
#import numpy as np
import matplotlib.mlab as mlab
from Bi2Se3Parameters import *
import matplotlib.gridspec as gs
import sys
from scipy.sparse.linalg import eigsh
from scipy import sparse
import cmath
import math
if os.path.exists('./results'):
  shutil.rmtree('./results')

os.mkdir('./results')
logging.basicConfig(filename='./results/out.log',filemode='w',level=logging.DEBUG)
logger=logging.getLogger()

class TI_tunnel:
	"""This code is to calculate TI_slab tunneling """
	def __init__(self,material_parameters_top,material_parameters_bottom,vacuum_parameters,centre_vacuum,constants,kpoints,total_thickness,vacuum_thickness,nz,sigma,nbands,take_bands,both,test):
	
		self.mat_param_top=material_parameters_top
		self.mat_param_bottom=material_parameters_bottom
		self.vac_param=vacuum_parameters
		self.cent_vac=centre_vacuum
		self.kpoints=kpoints#[0][0:-2],kpoints[1][0:-2]]
		self.nbands=nbands
		self.n_mat=int(np.ceil(nz*((total_thickness-vacuum_thickness)/total_thickness)))
		self.n_vac=int(nz-self.n_mat)
		self.dz=total_thickness/float(nz)
		print('self_dz'+str(self.dz))
		self.s=sigma
		self.constants=constants
		self.nz=nz
		self.take_bands=take_bands
		self.test=test
		self.both=both
		#if self.both==True:
		#	self.n_vac=self.n_vac//2
		print(self.n_mat,self.n_vac,self.mat_param_top,self.mat_param_bottom,self.vac_param,self.cent_vac)

	def H_mat(self,kx,ky,typ,side):
		if(typ=='TI'):	
			if side=='top':
				matparam=self.mat_param_top
			else:			
				matparam=self.mat_param_bottom
		elif(typ=='Vacuum'):
			matparam=self.vac_param

		elif(typ=='Vacc_mat'):
			matparam=self.cent_vac
#			print(matparam)
		kperp2=kx**2+ky**2
		H=np.zeros((4,4),dtype=np.complex)

		dz=self.dz#matparam['dz']
		A1=matparam['A1']
		A2=matparam['A2']
		Tz1=1j*matparam['A1']/(2*dz)
		M=matparam['M']
		B2=matparam['B2']
		B1=matparam['B1']
		D1=matparam['D1']
		D2=matparam['D2']
		C=matparam['C']

		Mtp=M-B2*kperp2-2*B1/(dz**2)+2*D1/(dz**2)+D2*kperp2+C#double deriveative B1*kz**2
		Mtm=-M+B2*kperp2+2*B1/(dz**2)+2*D1/(dz**2)+D2*kperp2+C
		delta=matparam['delta']

		kplus=kx+cmath.sqrt(-1)*ky
		kminus=kx-cmath.sqrt(-1)*ky	
		H=np.array([[Mtp+delta,0,0,A2*kminus],[0,Mtm+delta,A2*kminus,0],[0,A2*kplus,Mtp-delta,0],[A2*kplus,0,0,Mtm-delta]])

		Wt=np.array([[B1/(dz**2)-D1/(dz**2), -Tz1,0 ,0 ],[-Tz1,-B1/(dz**2)-D1/(dz**2),0 ,0 ],[0,0 ,B1/(dz**2)-D1/(dz**2),Tz1 ],[0,0 , Tz1,-B1/(dz**2)-D1/(dz**2)]])
		W=np.array([[B1/(dz**2)-D1/(dz**2),Tz1 ,0 ,0 ],[Tz1,-B1/(dz**2)-D1/(dz**2) ,0 ,0 ],[0,0 ,B1/(dz**2)-D1/(dz**2) ,-Tz1 ],[0,0 , -Tz1,-B1/(dz**2)-D1/(dz**2)]])

		return H,Wt,W
	
	def H_slab(self,side,kx,ky):
		
		
		H,Wt,W=self.H_mat(kx,ky,'TI','top')
		Hvac,Wtvac,Wvac=self.H_mat(kx,ky,'Vacuum','top')
		H_TI=np.kron(np.eye(self.n_mat),H)+np.kron(np.eye(self.n_mat,k=1),W)+np.kron(np.eye(self.n_mat,k=-1),Wt)
		H_Vac=np.kron(np.eye(self.n_vac),Hvac)+np.kron(np.eye(self.n_vac,k=1),Wtvac)+np.kron(np.eye(self.n_vac,k=-1),Wvac)
		HintR=np.zeros((len(H_TI[:,0]),len(H_Vac[0,:])),dtype=complex)
		HintL=np.zeros((len(H_Vac[:,0]),len(H_TI[0,:])),dtype=complex)
		HintLL=np.zeros((len(H_Vac[:,0]),len(H_TI[0,:])),dtype=complex)
		HintRR=np.zeros((len(H_TI[:,0]),len(H_Vac[0,:])),dtype=complex)

		H_B,Wt_B,W_B=self.H_mat(kx,ky,'TI','bottom')
		H_TI_B=np.kron(np.eye(self.n_mat),H_B)+np.kron(np.eye(self.n_mat,k=1),W_B)+np.kron(np.eye(self.n_mat,k=-1),Wt_B)
		HintR_B=np.zeros((len(H_TI_B[:,0]),len(H_Vac[0,:])),dtype=complex)
		HintL_B=np.zeros((len(H_Vac[:,0]),len(H_TI_B[0,:])),dtype=complex)
		HintLL_B=np.zeros((len(H_Vac[:,0]),len(H_TI_B[0,:])),dtype=complex)
		HintRR_B=np.zeros((len(H_TI_B[:,0]),len(H_Vac[0,:])),dtype=complex)

		HintR[len(HintR[:,0])-4:,0:4]=Wtvac
		HintL[0:4,len(HintL[0,:])-4:]=Wvac

		HintRR[0:4,len(HintRR[0,:])-4:]=Wvac
		HintLL[len(HintLL[:,0])-4:,0:4]=Wtvac
		
		HintR=sparse.csr_matrix(HintR)
		HintL=sparse.csr_matrix(HintL)

		HintRR=sparse.csr_matrix(HintRR)
		HintLL=sparse.csr_matrix(HintLL)

		HintR_B[len(HintR_B[:,0])-4:,0:4]=Wtvac
		HintL_B[0:4,len(HintL_B[0,:])-4:]=Wvac

		HintRR_B[0:4,len(HintRR_B[0,:])-4:]=Wvac
		HintLL_B[len(HintLL_B[:,0])-4:,0:4]=Wtvac
		
		HintR_B=sparse.csr_matrix(HintR_B)
		HintL_B=sparse.csr_matrix(HintL_B)

		HintRR_B=sparse.csr_matrix(HintRR_B)
		HintLL_B=sparse.csr_matrix(HintLL_B)


		if(side=='top'):
			
			H_big=sparse.bmat([[H_Vac,HintLL],[HintRR,H_TI]])
			if self.both==True:
				H_big=sparse.bmat([[H_Vac,HintLL,None],[HintRR,H_TI,HintRR],[None,HintLL,H_Vac]])
		elif(side=='bottom'):
			H_big=sparse.bmat([[H_TI_B,HintR_B],[HintL_B,H_Vac]])
			if self.both==True:
				H_big=sparse.bmat([[H_Vac,HintLL_B,None],[HintRR_B,H_TI_B,HintRR_B],[None,HintLL_B,H_Vac]])
		return H_big

	def Wave_func(self,side,kx,ky,state):
		if np.mean(self.H_slab('top',kx,ky)-self.H_slab('top',kx,ky).conj().T)!=0:
			print('Non-Hermitian')
			Evec=0
		else:
			val,vec = eigsh(self.H_slab(side,kx,ky),self.nbands,sigma=self.s)#expand around k or eV#eigsh(self.H_slab(side,kx,ky)/self.constants['eV'],self.nbands,sigma=self.s/self.constants['eV'])#expand around k or eV
			Ev=np.argsort(val)
			Evec=vec[:,Ev]
			Eval=val[Ev]
		if self.both==True:
			WF13=Evec[:len(Evec)//3,state]
			WF23=Evec[len(Evec)//3:int(2*len(Evec))//3,state]
			WF33=Evec[int(2*len(Evec))//3:len(Evec),state]
			if side=='top':
				Evec[:,state]=np.concatenate((WF23,WF13,WF33))
			if side=='bottom':
				Evec[:,state]=np.concatenate((WF33,WF13,WF23))
		#print('Energies:',Eval/eV)

		return(Evec[:,state])

	def Intersections(self,E1,E2):

#                BandsFermin=np.flip(BandsFermi)
		kx=self.kpoints[0]
		ky=self.kpoints[1]

		BandsNoFermi=E1
		BandsFermin=E2

		IntersectionEnergy=[]
		Intersectionkx=[]
		Intersectionky=[]
		IntersectionTheta=[]
		index=[]
		diffbands=BandsNoFermi-BandsFermin
		leng2=int(len(kx)/2)
		difff=[]
#		print(np.shape(BandsFermin))
		for i in range(len(BandsNoFermi[:,0])):
			for j in range(len(BandsFermin[:,0])):

				diff=(BandsNoFermi[i,:]-BandsFermin[j,:])
				cond=diff[1:]*diff[:1]<0
				diff2=(np.where(cond))
				diffE=(np.extract(cond,BandsNoFermi[i,:]))
				if(len(diff2[0])>0):
					
					difff_K=np.where(np.absolute(diff)<1e-4)
					difff.append([i,j])#,diff2[0])#indicates 2*10 intersections
		outF = open("Intersections.txt", "w")
		print(difff)
		#plt.figure(1)
		#for l in range(len(BandsNoFermi[:,0])):
		#	plt.plot(BandsNoFermi[l,:])
		#	plt.plot(BandsFermin[l,:])
		#plt.show()
		for i in range(len(difff)):

			df=abs(BandsNoFermi[difff[i][0]]-BandsFermin[difff[i][1]])#should be 0 or very close to it at intersection point
			index.append(np.where(df==min(df))[0][0])
			#print(np.real(kx[np.where(df==min(df))[0]]))
			Intersectionkx.append(float(np.real(kx[np.where(df==min(df))[0][0]])))
			Intersectionky.append(float(np.real(ky[np.where(df==min(df))[0][0]])))
			IntersectionEnergy.append(np.real(BandsNoFermi[difff[i][0],np.where(df==min(df))[0]]))
		
		if self.test==True:
			print('Number of Intersections:', len(Intersectionkx))

	#	Intersectionkx=np.asarray([Intersectionkx]).reshape(len(Intersectionkx),1)
	#	Intersectionky=np.asarray(Intersectionky)

		return Intersectionkx,Intersectionky,IntersectionEnergy,difff,index




	def Transition_matrix2(self,kx,ky,band1,band2):


		Wf1,E1=np.array(self.Wave_func('top',kx,ky,band1))
		Wf2,E1=np.array(self.Wave_func('bottom',kx,ky,band2))

		Hvac,Wtvac,Wvac=self.H_mat(kx,ky,'Vacc_mat','top')
		size=self.n_mat+self.n_vac

		#HV=np.kron(np.eye(size),Hvac)#+np.kron(np.eye(size,k=1),Wtvac)+np.kron(np.eye(size,k=-1),Wvac)	
		HV=np.eye(4*size)*self.cent_vac['C']#/self.constants['eV']

		plt.spy(HV)
		plt.show()
#		print(Wf1)
		D=HV.dot(Wf1)
		iHf=np.absolute((Wf2.conj().T).dot(D))**2

#		iHf=np.absolute(np.dot(Wf2.T.conj(),np.dot(HV,Wf1)))**2
		

		#iHf=np.absolute((Wf2.T.conj()).dot(HV.dot(Wf1)))**2
		return(iHf)

	def Transition_matrix(self,kx,ky,Wf1,Wf2):

		Hvac,Wtvac,Wvac=self.H_mat(kx,ky,'Vacc_mat','top')
		size=self.n_mat+self.n_vac
		HV=np.kron(np.eye(int(size)),Hvac)#/self.constants['eV']
		#if self.both==True:
		#	size=int(self.nz*4*1.5+4*(self.nz-40)/20+4)#self.n_mat+self.n_vac
		#	HV=np.kron(np.eye(int(size/4)),Hvac)#/self.constants['eV']
		
		#plt.spy(HV)
		#plt.show()
		#if(self.test==True):
		#	print(HV,Wf1)
		D=HV.dot(Wf1)
		M=((Wf2.T.conj()).dot(D))

		iHf=M.real**2+M.imag**2
		#print('Overlap:',np.real(np.sum(np.multiply(Wf1,Wf2) ) ) )
		#if(self.test==True):
		#zNT=Wf1.reshape(int(self.nz*4*1.5+4*(self.nz-40)/20+4)//4,4)
		#zNB=Wf2.reshape(int(self.nz*4*1.5+4*(self.nz-40)/20+4)//4,4)
		#plt.figure(1)			
		#plt.plot(zNT.real)
		#plt.plot(zNB.real)
		#plt.show()

		return iHf	

	def drawProgressBar(self,percent,barLen=20):
		barLen=int(barLen)
		sys.stdout.write("\r")
		progress = ""
		for i in range(barLen):
			if i < int(barLen * percent):
				progress += "="
			else:
				progress += " "
		sys.stdout.write("[ %s ] %.2f%%" % (progress, percent * 100))
		sys.stdout.flush()
	def eigen(self,kx,ky,typ='top',prev=0):
		if np.mean(self.H_slab('top',kx,ky)-self.H_slab('top',kx,ky).conj().T)!=0:
			print('Non-Hermitian')
			f=0
			Vec_sort=0
		else:
			#print(self.take_bands)
			Vec_sort=np.zeros((self.take_bands,self.take_bands),dtype=complex)
			Eigs,Vec=eigsh(self.H_slab(typ,kx,ky),self.nbands,sigma=self.s)#self.s)#eigsh(self.H_slab(typ,kx,ky)/self.constants['eV'],self.nbands,sigma=self.s/self.constants['eV'])#self.s)
			EV=np.argsort(Eigs)
			Eig=Eigs[EV]
			Vec=Vec[:,EV]
			#sortval=np.sort(val)
			#print((Eig))
			imiddle=np.where(Eig>=(self.s))[0][0]#np.where(Eig>=(self.s/self.constants['eV']))[0][0]
			#print(imiddle,self.take_bands,len(Eig),self.s,kx)
			f=Eig[imiddle-self.take_bands//2:imiddle-self.take_bands//2+self.take_bands]
			Vec_sort=Vec[:,imiddle-self.take_bands//2:imiddle-self.take_bands//2+self.take_bands]
			for i in range(len(Vec_sort[0,:])):
				C=np.sqrt(np.trapz(Vec_sort[:,i]*Vec_sort[:,i].conj(),None,self.dz))
				Vec_sort[:,i]=Vec_sort[:,i]/C
				Vec_Int=np.trapz(Vec_sort[:,i]*Vec_sort[:,i].conj(),None,self.dz)
				#print('Normalized Integral:',Vec_Int)
		#print(np.shape(f))
		return f,Vec_sort

	def Bandstructure(self,Vt=0,typ='top'):

		#E=[]
		E=np.zeros((len(self.kpoints[0]),self.take_bands),dtype=float)		
		for j in range(len(self.kpoints[0])):
			Eig,Vec=self.eigen(self.kpoints[0][j],self.kpoints[1][j],prev=np.mean(E[j-1,:]))
	#		print(len(Eig),self.take_bands)
			E[j,:]=Eig


		return(E)

	def dedk(self,band,index):
		
		dk=(self.kpoints[0][1]-self.kpoints[0][0])	
#		print(dk)
		if index+1!=len(band):
			derivative=(band[index-1]-band[index+1])/(2*dk)
		else:
			derivative=(band[index-1]-band[index])/(2*dk)
		return derivative
	def dJdV(self,Current,Vd):
		Jdiff=[]
		Jdiff=np.diff(Current)/np.diff(Vd)
		return Jdiff

	
	def Current_calc(self,Vg,Vd,side):
		Total_cur=[]
		#for i in range(len(Vd)):
		if 1:
			EE1=np.zeros((len(self.kpoints[0]),self.take_bands),dtype=float)
			EE2=np.zeros((len(self.kpoints[0]),self.take_bands),dtype=float)
			#print(-Vd[i]/self.constants['eV']+Vg/self.constants['eV'])
			if self.both==True:
				Vec1=np.zeros((len(self.kpoints[0]),int(736),self.take_bands),dtype=complex)
				Vec2=np.zeros((len(self.kpoints[0]),int(736),self.take_bands),dtype=complex)
			else:
				Vec1=np.zeros((len(self.kpoints[0]),int(self.nz*4),self.take_bands),dtype=complex)
				Vec2=np.zeros((len(self.kpoints[0]),int(self.nz*4),self.take_bands),dtype=complex)
			for j in range(len(self.kpoints[0])):

				if(side=='top'):

					Eig1,V1=self.eigen(self.kpoints[0][j],self.kpoints[1][j],typ='top')
					Eig2,V2=self.eigen(self.kpoints[0][j],self.kpoints[1][j],typ='bottom')
				#	print(Vd[i]/self.constants['eV'])
					Vec1[j,:,:]=V1
					Vec2[j,:,:]=V2
					EE2[j,:]=Eig2
					EE1[j,:]=Eig1#-Vd[i]/self.constants['eV']+Vg/self.constants['eV']
		for i in range(len(Vd)):
			
			E1=EE1-Vd[i]+(Vg)#-Vd[i]/(self.constants['eV'])+(Vg/self.constants['eV'])
			E2=EE2

			Intersection_kx,Intersection_ky,IntersectionEnergy,Intersection_band,index=self.Intersections(np.array(E1).T,np.array(E2).T)
		
			print('number of intersections:',len(IntersectionEnergy))
			#print('Number of intersections for',len(E1),'bands:',len(Intersection_kx), 'at Vd=',Vd[i]/eV)
			J=0
			EdCurrent=np.zeros(len(Intersection_kx))
			
			if(self.test==True):
				for j in range(len(E1[0,:])):
					#plt.plot(np.flip(self.kpoints[0]),np.flip(E1[:,j]),color='k')
					#plt.plot(np.flip(self.kpoints[0]),np.flip(E2[:,j]),color='k')
	
					plt.plot(self.kpoints[0]/10**10,E1[:,j]/eV,color='k')
					plt.plot(-self.kpoints[0]/10**10,-(E2[:,j])/eV,color='b')
					plt.scatter(Intersection_kx[j-1]/10**10,IntersectionEnergy[j-1]/eV,color='r')
					plt.vlines(0.0001,-0.6,0.6)
					plt.hlines(0,-0.1,0.1,linestyle='--',color='b')
					plt.hlines(-Vd[i]/eV,-0.1,0.1,linestyle='--',color='k')
				plt.ylim(-0.5,0.5)   
				plt.xlim(-0.1,0.1)
				plt.xlabel('$k_x$ $\dot A$',fontsize=16)
			
				plt.ylabel('$E$ (eV)',fontsize=16)
				plt.show()
				
				plt.close()

			for j in range(len(Intersection_kx)):
				if len(IntersectionEnergy[j])!=1:
					IntersectionEnergy[j]=IntersectionEnergy[j][0]
					
				
				Ieft=IntersectionEnergy[j]+(Vd[i])#+(Vd[i]/self.constants['eV'])#+Vg/self.constants['eV']#+Vt#need to add an option for top and bottom later
				Iefb=IntersectionEnergy[j]
				print('Ieft:',Ieft,'Iefb:',Iefb)
				iHf=self.Transition_matrix(Intersection_kx[j],Intersection_ky[j],Vec1[index[j],:,Intersection_band[j][0]],Vec2[index[j],:,Intersection_band[j][1]])#/(self.constants['eV']**2)
				if(self.test==True):
					print('Intersecting K-point:', Intersection_kx[j]/10**10,'Energy:',IntersectionEnergy[j]/eV,'Overlap:',np.real(np.sum(np.multiply(Vec1[index[j],:,Intersection_band[j][0]],Vec2[index[j],:,Intersection_band[j][1]]) ) ),'Intersecting Bands:', Intersection_band[j][0],Intersection_band[j][1])
				print('')
				print('Average Transition Element:',np.mean(iHf),'for Vd=',Vd[i]/eV,'Overlap:',np.real(np.sum(np.multiply(Vec1[index[j],:,Intersection_band[j][0]],Vec2[index[j],:,Intersection_band[j][1]]) ) ) )
				
				

				if(Intersection_band[j][0]==Intersection_band[j][1]):
					Ederi=np.absolute(self.dedk(E1[:,Intersection_band[j][0]],index[j]))
				else:
					Ederi=np.absolute(self.dedk(E1[:,Intersection_band[j][0]],index[j])-self.dedk(E2[:,Intersection_band[j][1]],index[j]))
				print('Average Derivative:',np.mean(Ederi),'for Vd=',Vd[i]/eV)
			#	if(self.test==True):
			#		plt.plot(self.kpoints[0],E1[:,Intersection_band[j][0]],color='b')
			#		plt.plot(self.kpoints[0],E2[:,Intersection_band[j][1]],color='r')
			#		plt.scatter(Intersection_kx[j],E1[index[j],Intersection_band[j][0]],color='black')
			#		plt.scatter(Intersection_kx[j],E2[index[j],Intersection_band[j][1]],color='black')

			#		plt.show()

				k1=np.absolute(Intersection_kx[j])#self.kpoints[0][index[j]]	
				fl=(np.exp(Ieft/(self.constants['kT']*self.constants['eV']))+1)**(-1)#(np.exp(Ieft/(self.constants['kT']))+1)**(-1)
				fr=(np.exp(Iefb/(self.constants['kT']*self.constants['eV']))+1)**(-1)#(np.exp(Iefb/(self.constants['kT']))+1)**(-1)
				print(np.absolute(fl-fr)/eV)
				if(np.absolute(fl-fr)/eV<1e-5):
					fl=0
					fr=0
	#			iHf=1
	#			Ederi=1
		#		k1=1
				#print(Ieft,Iefb,fr-fl,self.constants['kT'])
				#print(Intersection_kx[j],Intersection_ky[j],Intersection_band[j][0],Intersection_band[j][1],IntersectionEnergy[j],(fr-fl),self.kpoints[0][index[j]],Ederi,k1)
				J=-eV/(2*math.pi*hbar)*iHf*(fl-fr)*k1*(Ederi**(-1))*(self.constants['eV'])+J#J=-2*math.pi*iHf*(fl-fr)*k1*(Ederi**(-1))*(self.constants['eV'])*(10**10)+J # fix for hbar units later
			if(self.test==True):			
				print(J,Vd[i]/self.constants['eV'],len(Intersection_kx),Ieft,Iefb,iHf,Ederi,k1/(10**10),(fl-fr))
			#if np.mean(Ederi**-1)/Edlast>10:
			#	continue
			#if len(Intersection_kx)==1:
			#	Total_cur.append(2*J)
			#else:
			Total_cur.append(J)
			Edlast=np.mean(Ederi**-1)
			self.drawProgressBar(i/float(len(Vd)))
		#	print(Ieft)
		return(Total_cur)






