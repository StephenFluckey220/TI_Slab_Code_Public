# Bi2Se3Parameters
import cmath

pi=cmath.pi#unit predeclarations 
J=1
kg=1
s=1
meter=1#10**10

eV=1.602e-19*J
hbar=1.05*10**(-34)*J*s
hbar2 = 1#1.05 * 10 ** (-34)*J*s
m0= 9.31*10**(-31)*kg
kT = 0.026


V=1*eV
M=0.28*eV
A1=2.2*eV*(10**(-10)*meter)
A2=4.1*eV*(10**(-10)*meter)
B1=10*eV*(10**(-10)*meter)**2
B2=56.6*eV*(10**(-10)*meter)**2
C=-0.0068*eV
D1=1.3*eV*(10**(-10)*meter)**2
D2=19.6*eV*(10**(-10)*meter)**2
deltaT=0.0*eV
deltaB=0*eV
#print(matparam)


matparam={'A1':A1, 'A2':A2, 'M':M, 'B1':B1, 'B2':B2,'C':C, 'D1':D1, 'D2':D2,'delta':deltaT}
matparam2={'A1':A1, 'A2':A2, 'M':1*M, 'B1':B1, 'B2':B2,'C':C, 'D1':D1, 'D2':D2,'delta':deltaB}
print(matparam)

vacparam={'A1':0, 'A2':0, 'M':0, 'B1':0, 'B2':0,'C':V, 'D1':D1, 'D2':D2,'delta':deltaB}
cent_vacc={'A1':0, 'A2':0, 'M':0, 'B1':0, 'B2':0,'C':0, 'D1':D1, 'D2':D2,'delta':deltaB}
