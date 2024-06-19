#libraries
import numpy as np
import matplotlib.pyplot as plt
import sympy as sp
from sympy.physics.wigner import wigner_3j as Wig3
from sympy.physics.wigner import wigner_6j as Wig6

#variables
#variation variables
w = np.linspace(0,120,10)
#sat_z = [0.512, 0.576, 0.640, 0.704, 0.768, 0.832, 0.896, 0.960, 1.024]
#epilson =[0.012, 0.024, 0.036, 0.048, 0.060, 0.072, 0.084, 0.096, 0.108, 0.120]

#constant variables
M = 2.208518465e-25 #133#in au which is 2.208518465e-25 kg
phi = 0             #if phase is assumed
h = 1.055e-34       #Js
mu = 9.27E-24       #Joule/Tesla
S = 1/2             #spin
lg = 0              #angular orbital no, g - ground , e - excited
le = 1
jg = 1/2            #total angular momentum no
je = 3/2
I = 7/2             #total nuclear angular momentum
k = 11732.307104e2  #m**-1 wavenumber
Gamma = 2*np.pi*5.22227e6     # natural linwidth of ceasium
delta = 0.2 * 3.519654552e14  # here delta is detuning in the transition Fg = 4 to Fe = 5, here i assuming a 0.2 detuning , data from danielsteck  https://steck.us/alkalidata/cesiumnumbers.1.6.pdf
Delta = 151.21e6              # this is defined as hyperfine spacing of the excited states
delta_z = -2.49 * Gamma
time = 2/Gamma

fglow = 3
fgup = 4

s_z, s_t ,s_r, z, t, r = sp.symbols("s_z,s_t,s_r,z,t,r")                                                                       #saturation parameters of long, trans , repump
del_r, del_t, del_z = sp.symbols("del_r,del_t,del_z")                                                                          #detuning of long and repump
f1, f2, f3, f4, f5, f6, f7 = sp.symbols("f1, f2, f3, f4, f5, f6, f7")                                                          #fg=3 ground state populations
g1, g2, g3, g4, g5, g6, g7, g8, g9 = sp.symbols("g1, g2, g3, g4, g5, g6, g7, g8, g9")                                          #fg=4 ground level populations
p1, p2, p3, p4, p5 = sp.symbols("p1, p2, p3, p4,p5")                                                                           #fe=2 excited level populations
q1, q2, q3, q4, q5, q6, q7 = sp.symbols("q1, q2, q3, q4, q5, q6, q7")                                                          #fe=3 excited level populations
r1, r2, r3, r4 ,r5, r6, r7, r8, r9 = sp.symbols("r1, r2, r3, r4, r5, r6, r7, r8, r9")                                          #fe=4 excited level populations
t1, t2, t3, t4, t5, t6, t7, t8, t9, t10 ,t11 = sp.symbols("t1,t2,t3,t4,t5,t6,t7,t8,t9,t10,t11")                                #fe=5 excited level populations
lam, Ga, t = sp.symbols("lam,Ga,t")

#equations
def Amplitude(x,r,s_z,epsilon):
    return (F_0(s_z,epsilon)) / (M*sp.sqrt((w2(s_z)-x**2)**2+((G(s_z,r)**2)*x**2))) * sp.cos(x * time)

def w2(s_z):
    x = 0
    fg = 4
    for fe in range(fg-1,fg+2):
        for mg in range(-fg,fg+1):
           x = x + (w2Fg_Fe(fg+1,fe,s_z) * (((mg + 1)* Landeg(fe,le,je,I) - (mg * Landeg(fg,lg,jg,I))) * A(fe, mg, 1, fg, mg) - ((mg - 1) * Landeg(fe,le,je,I) - (mg * Landeg(fg,le,je,I))) * A(fe, mg, 1, fg, mg)) * kronckerdelta(fg,4))
    return x

def G(s_z,r):
    x = 0
    fg = 4
    for fe in range(fg-1,fg+2):
        for mg in range(-fg,fg+1):
            x = x + (gFg_Fe(fg+1,fe,s_z,r) * (A(fe, mg, 1, fg, mg) + A(fe, mg, -1, fg, mg)) * kronckerdelta(fg,4))
    return x

def F_0(s_z,epsilon):
    x = 0
    fg = 4
    for fe in range(fg-1,fg+2):
        for mg in range(-fg,fg+1):
            x = x + (FFg_Fe(fg+1,fe,s_z) * (A(fe, mg, +1, fg, mg) + A(fe, mg, -1, fg, mg)) * kronckerdelta(fg,4))
    return x

def Delta(a,b):                                     #hyperfine excited levels spacing
    if a == 4:
       return 352.45e6 #hz
    elif a == 5:
        return 452.24e6 #hz
    else :
        return 0

def gFg_Fe(a,b,s_z,r):
    return r * w2Fg_Fe(a,b,s_z)

def w2Fg_Fe(a,b,s_z):
    return ((-4 * mu * k * s_z *(delta_z + Delta(a,b))) / Gamma ) / (M * (1 + 4 * (delta_z + Delta(a,b))**2  / Gamma**2 )**2 )

def FFg_Fe(a,b,s_z):
    return (h * k * s_z * Gamma * epsilon * (1 + 4 *(delta_z + Delta(a,b))**2 / Gamma**2)) / (2 * (1 + 4 * (delta_z + Delta(a,b))**2 / Gamma**2 )**2 )

def A(fe,mg1,i,fg,mg):
    return R(fe,mg1,i,fg,mg) * (Ppop(fe,mg) - Qpop(fg,mg))

def R(fe,mg1,i,fg,mg):
    return 3 * (2 * je + 1) * (2 * jg + 1) * (2 * fe + 1) * (2 * fg + 1) * (Wig6(le,je,S,jg,lg,1) * Wig6(je,fe,I,fg,jg,1) * Wig3(fg,1,fe,mg,i,-(mg+i)))**2

def Landeg(f,l,j,i):
    return ((f*(f+1) + j*(j+1) - i*(i+1)) / (2*f*(f+1)) ) * ((3*j*(j+1) - S*(S+1) + l*(l+1)) / (2*j*(j+1)))

def kronckerdelta(a, b):
    if (a == b):
        return 1
    else:
        return 0

#new part of program
def Ppop(b,a):
    if (a == -4 and b == 4):
        return 0.068
    elif (a == -3 and b == 4):
        return 0.068
    elif (a == -2 and b == 4):
        return 0.068
    elif (a == -1 and b == 4):
        return 0.068
    elif (a == 0 and b == 4):
        return 0.068
    elif (a == 1 and b == 4):
        return 0.068
    elif (a == 2 and b == 4):
        return 0.068
    elif (a == 3 and b == 4):
        return 0.068
    elif (a == 4 and b == 4):
        return 0.068
    elif (a == -3 and b == 3):
        return 0.013
    elif (a == -2 and b == 3):
        return 0.013
    elif (a == -1 and b == 3):
        return 0.013
    elif (a == 0 and b == 3):
        return 0.013
    elif (a == 1 and b == 3):
        return 0.013
    elif (a == 2 and b == 3):
        return 0.013
    elif (a == 3 and b == 3):
        return 0.013
    else :
        return 0
def Qpop(b,a):
   if b == 2:
        if a == -2:
           return 0.002
        elif a == -1:
           return 0.002
        elif a == 0:
            return 0.002
        elif a == 1:
            return 0.002
        elif a == 2:
            return 0.002
        else:
            return 0
   elif b == 3:
        if a == -3:
            return 0.0014
        elif a == -2:
            return 0.0014
        elif a == -1:
            return 0.0014
        elif a == 0:
            return 0.0014
        elif a == 1:
            return 0.0014
        elif a == 2:
            return 0.0014
        elif a == 3:
            return 0.0014
        else:
            return 0
   elif b == 4:
        if a == -4:
            return 0.021
        elif a == -3:
            return 0.021
        elif a == -2:
            return 0.021
        elif a == -1:
            return 0.021
        elif a == 0:
            return 0.021
        elif a == 1:
            return 0.021
        elif a == 2:
            return 0.021
        elif a == 3:
            return 0.021
        elif a == 4:
            return 0.021
        else:
            return 0
   elif b == 5:
        if a == -5:
            return 0.017
        elif a == -4:
            return 0.017
        elif a == -3:
            return 0.017
        elif a == -2:
            return 0.017
        elif a == -1:
            return 0.017
        elif a == 0:
            return 0.017
        elif a == 1:
            return 0.017
        elif a == 2:
            return 0.017
        elif a == 3:
            return 0.017
        elif a == 4:
            return 0.017
        elif a == 5:
            return 0.017
        else:
            return 0
####################

#plotting
#plotting variation wrt b
s_z = 0.779
epsilon = 0.073
#b = [0.057,0.068,0.080,0.091,0.102,0.114,0.125,0.136,0.148]
b = [0.091]
rho = np.zeros(len(b))
for i in range(len(b)):
    rho[i] = (h * k) / (mu * b[i])
Amp = [Amplitude(w[j],b[0],s_z,epsilon) for j in range(len(w))]
Amp = np.delete(Amp,[0])
print(Amp)
w = np.delete(w,[0])
plt.plot(w,Amp,'o-')
