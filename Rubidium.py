#liibraries 
import numpy as np
import sympy as sp
import matplotlib.pyplot as plt
from sympy.physics.wigner import wigner_3j as Wig3
from sympy.physics.wigner import wigner_6j as Wig6
from sympy.functions.elementary.exponential import exp_polar

#constants
#Gamma =  36.10e6 #Hz
#delta = 1.5*Gamma #Hz
S = 1/2
lg = 0
le = 1
jg = 1/2
je = 3/2
I = 3/2

#variables 
f1, f2, f3 = sp.symbols("f1, f2, f3")                                                               #fg=1 ground state populations 
h1, h2, h3, h4, h5 = sp.symbols("h1, h2, h3, h4, h5")                                               #fg=2 ground level populations
q1 = sp.Symbol("q1")                                                                                #fe=0 excited level populations 
p1, p2 = sp.symbols("p1, p2")                                                                       #fe=1 excited level populations 
g1, g2, g3 = sp.symbols("g1, g2, g3")                                                               #fe=2 excited level populations
f1_dot, f2_dot, f3_dot = sp.symbols("f1_dot, f2_dot, f3_dot")                                       #derivative of the energy level populations 
h1_dot, h2_dot, h3_dot, h4_dot, h5_dot = sp.symbols("h1_dot, h2_dot, h3_dot, h4_dot, h5_dot")
q1_dot = sp.Symbol("q1_dot") 
p1_dot, p2_dot = sp.symbols("p1_dot, p2_dot") 
g1_dot, g2_dot, g3_dot = sp.symbols("g1_dot, g2_dot, g3_dot")
lam, lam1, lam2, lam3, G, t = sp.symbols("lam, lam1, lam2, lam3, G, t")                             # trial function  
s0, s1, s2 = sp.symbols("s0,s1,s2")                                                                 #saturation parameter

#Derivatives
def P_dot(a,b):
    if (a == -2 and b == 2) :
        h1_dot = sp.Derivative(h1, t)
        return h1_dot 
    elif (a == -1 and b == 2):
        h2_dot = sp.Derivative(h2, t)
        return h2_dot
    elif (a == 0 and b == 2):
        h3_dot = sp.Derivative(h3, t)
        return h3_dot
    elif (a == 1 and b == 2):
        h4_dot = sp.Derivative(h4, t)
        return h4_dot
    elif (a == 2 and b == 2):
        h5_dot = sp.Derivative(h5, t)
        return h5_dot
    elif (a == -1 and b == 1):
        f1_dot = sp.Derivative(f1, t)
        return f1_dot
    elif (a == 0 and b == 1):
        f2_dot = sp.Derivative(f2, t)
        return f2_dot
    elif (a == 1 and b == 1):
        f3_dot = sp.Derivative(f3, t)
        return f3_dot 
    else :
        return 0

def Q_dot(a,b):
    if (a == 0 and b == 2):
        g1_dot = sp.Derivative(g1, t)
        return g1_dot
    elif (a == 1 and b == 2):
        g2_dot = sp.Derivative(g2, t)
        return g2_dot
    elif (a == 2 and b == 2):
        g3_dot = sp.Derivative(g3, t)
        return g3_dot
    elif (a == 0 and b == 1):
        p1_dot = sp.Derivative(p1, t)
        return p1_dot
    elif (a == 1 and b == 1):
        p2_dot = sp.Derivative(p2, t)
        return p2_dot
    elif (a == 0 and b == 0):
        q1_dot = sp.Derivative(q1, t)
        return q1_dot
    else :
        return 0
    
#variable assignment
def P(a,b):
    if (a == -2 and b == 2): 
        return h1
    elif (a == -1 and b == 2):
        return h2
    elif (a == 0 and b == 2):
        return h3
    elif (a == 1 and b == 2):
        return h4
    elif (a == 2 and b == 2):
        return h5
    elif (a == -1 and b == 1):
        return f1
    elif (a == 0 and b == 1):
        return f2
    elif (a == 1 and b == 1):
        return f3 
    else :
        return 0
    
def Q(a,b):
    if (a == 0 and b == 2):
        return g1
    elif (a == 1 and b == 2):
        return g2
    elif (a == 2 and b == 2):
        return g3
    elif (a == 0 and b == 1):
        return p1
    elif (a == 1 and b == 1):
        return p2
    elif (a == 0 and b == 0):
        return q1
    else :
        return 0

#equations
def C(e,me,g,mg):
    return (2*le +1) * (2*je + 1) * (2*jg + 1) * (2*e + 1) * (2*g + 1) * (Wig6(le,je,S,jg,lg,1) * Wig6(je,e,I,g,jg,1) * Wig3(g,1,e,mg,me-mg,-me))**2

def sat(a):
    if a == 2:
        return s2
    elif a == 1 :
        return s1
    else :
        return s0

def RateEq1():
    #1st equation
    fg1eq = []
    for m in range(-1,2):
        ans = 0 
        for fe in range(0,3):
            lv = 0
            for me in range(m-1, m+2):
                lv = lv + (G * C(fe,me,1,m) * Q(me,fe))
            ans = ans + lv + ((-G/2) * C(fe,m+1,1,m) * sat(fe) * (P(m,1) - Q(m+1,fe)))
        ans = P_dot(m,1) - ans
        #print(ans)
        fg1eq.append(ans)
    return fg1eq
        
def RateEq2():
    fg2eq = []
    for m in range(-2,3):
        ans = 0
        for i in range(1,4):
            for j in range(m-1, m+2):
                ans = ans + (G * C(i,j,2,m) * Q(j,i))
        ans = P_dot(m, 2) - ans
        #print(ans)
        fg2eq.append(ans)
    return fg2eq
        
def RateEq3():
    fe_eq = []
    for fe in range(0,3):
        for m in range(-fe,fe+1):
            ans = 0   
            lv = 0
            for i in range(1,3):
                for j in range(m-1, m+2):
                    lv = lv + (G * C(fe,m,i,j) * Q(m,fe))
            ans = Q_dot(m,fe) - ((G/2) * C(fe,m,1,m-1) * sat(fe) * (P(m-1,1) - Q(m,fe))) + lv
            fe_eq.append(ans)
            #print(ans)
    return fe_eq       

fg1eq = RateEq1()
fg2eq = RateEq2()
fe_eq = RateEq3()

exp_h1 = fg2eq[0]
exp_h2 = fg2eq[1]
exp_h3 = fg2eq[2]
exp_h4 = fg2eq[3]
exp_h5 = fg2eq[4]

exp_f1 = fg1eq[0]
exp_f2 = fg1eq[1]
exp_f3 = fg1eq[2]

exp_q1 = fe_eq[0]
exp_p1 = fe_eq[2]
exp_p2 = fe_eq[3]
exp_g1 = fe_eq[6]
exp_g2 = fe_eq[7]
exp_g3 = fe_eq[8]

#trial function substitution
#first we need to find the equation for f1 in terms of either g1_dot, p1_dot or q1_dot, since g is in every transistion , g is used
tf = sp.exp(-lam * G * t)
tf_exp1 = exp_f1 + exp_g1 + exp_p1 + exp_q1
tf_exp2 = exp_f2 + exp_g2 + exp_p2
tf_exp3 = exp_f3 + exp_g3
#sp.pprint(tf_exp1)
tf_exp_g1 = exp_g1.subs([((f1 - g1), f1)])
tf_f1 = sp.solve(tf_exp_g1,f1,simplify = False)[0]
c_g1 = abs(exp_g1.coeff((s2 * G * (f1 - g1)) , 1)) * 2
tf_f1 = tf_f1.subs(g1,s2*tf*c_g1)
tf_f1 = tf_f1.doit()
tf_f1dot = sp.Derivative(tf_f1, t)
tf_f1dot = tf_f1dot.doit()

tf_d = sp.diff(tf,t)

c2 = abs(exp_p1.coeff((s1 * G * (f1 - p1)), 1)) * 2
c3 = abs(exp_q1.coeff((s0 * G * (f1 - q1)), 1)) * 2
tf_eq1 = tf_exp1.subs([(g1,s2*c_g1*tf),(p1,c2*s1*tf),(q1,c3*s0*tf),(f1,tf_f1),(sp.Derivative(g1,t),c_g1*tf_d),(sp.Derivative(p1,t),c2*s1*tf_d),(sp.Derivative(q1,t),c3*s0*tf_d),(sp.Derivative(f1,t),tf_f1dot)])
tf_eq1 = tf_eq1.doit()
tf_eq1 = tf_eq1 / (G * tf)
tf_eq1 = sp.simplify(tf_eq1)
roots1 = sp.solve(tf_eq1, lam)
#print(roots1)              #without approximation

lam_Eq0 = 1
lam_Eq1 = tf_eq1.coeff(lam,0)/2 #approximation

#Populations Calcution 
level0 = sp.exp(-G * t)
level1 = sp.exp(-lam1 * G * t)
level2 = sp.exp(-lam2 * G * t)
level3 = sp.exp(-lam3 * G * t)

f1_popeq = sp.Eq(f1,(1/8)*level1)
sp.pprint(f1_popeq)

g1_popeq = exp_g1
g1_popeq = g1_popeq.subs(sp.Derivative(g1,t),0)
g1_popeq = (1/8) * sp.solve(g1_popeq, G*g1)[0]
g1_popeq = g1_popeq.subs([(G,1),(f1,level1),(g1,level0)])
g1_popeq = sp.Eq(g1, g1_popeq)
sp.pprint(g1_popeq)

p1_popeq = exp_p1
p1_popeq = p1_popeq.subs(sp.Derivative(p1,t),0)
p1_popeq = (1/8) * sp.solve(p1_popeq, G*p1)[0]
p1_popeq = p1_popeq.subs([(G,1),(f1,level1),(p1,level0)])
p1_popeq = sp.Eq(p1, p1_popeq)
sp.pprint(p1_popeq)

q1_popeq = exp_q1
q1_popeq = q1_popeq.subs(sp.Derivative(q1,t),0)
q1_popeq = (1/8) * sp.solve(q1_popeq, G*q1)[0]
q1_popeq = q1_popeq.subs([(G,1),(f1,level1),(q1,level0)])
q1_popeq = sp.Eq(q1, q1_popeq)
sp.pprint(q1_popeq)

#for second level transitions
tf_exp_g2 = exp_g2.subs([((f2 - g2), f2)])
tf_f2 = sp.solve(tf_exp_g2,f2,simplify = False)[0]
c_g2 = abs(exp_g2.coeff(s2 * G * (f2 - g2) , 1))*2
tf_f2 = tf_f2.subs(g2,c_g2*s2*tf)
tf_f2 = tf_f2.doit()
tf_f2dot = sp.Derivative(tf_f2, t)
tf_f2dot = tf_f2dot.doit()

c2 = abs(exp_p2.coeff((s1 * G * (f2 - p2)), 1)) * 2
tf_eq2 = tf_exp2.subs([(g2,s2*c_g2*tf),(p2,c2*s1*tf),(f2,tf_f2),(sp.Derivative(g2,t),c_g2*tf_d),(sp.Derivative(p2,t),c2*s1*tf_d),(sp.Derivative(f2,t),tf_f2dot),(g1, (g1_popeq).rhs),(q1, (q1_popeq).rhs)])
tf_eq2 = tf_eq2.doit()
tf_eq2 = tf_eq2.subs([(G*tf, 1)])
tf_eq2 = sp.simplify(tf_eq2)
roots2 = sp.solve(tf_eq2, lam)
#sp.pprint(roots2)        #without aproximation

lam_eq3 = tf_eq2.coeff(lam,0)/2
lam_eq3 = lam_eq3.subs(G,0)
lam_eq3 = sp.simplify(lam_eq3)
#sp.pprint(lam_eq3)
