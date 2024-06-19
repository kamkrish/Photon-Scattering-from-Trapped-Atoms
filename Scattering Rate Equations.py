#libraries
import numpy as np
import matplotlib.pyplot as plt

#constants
#Na-3D
Na_diameter = 18e-3
Na_I_cool = 9.39
Na_I_re = 10e-3
Na_Gamma = 2*np.pi*9.795e6

#K-3D
K_diameter = 24e-3 
K_I_cool = 1.75
K_I_re = 13e-3
K_Gamma = 2*np.pi*6.035e6

delta_Na = Na_Gamma*6
delta_K = K_Gamma*6

s=np.linspace(0,100,10000)

#Rabi Frequency
def Rabi(s, Gamma):
    Omega = Gamma * np.sqrt(s/2)
    return Omega

#Scattering Rate equation
def Rate_equation(Gamma,Omega,delta):
    R = (Gamma/2)*((Omega**2/2)/(Gamma**2/4 + Omega**2/2 + delta**2))
    return R

#calculate
#Na
Na_Rabi_cool = Rabi(s, Na_Gamma)
Na_R_cool = Rate_equation(Na_Gamma, Na_Rabi_cool, delta_Na)

#K
K_Rabi_cool = Rabi(s, K_Gamma)
K_R_cool = Rate_equation(K_Gamma, K_Rabi_cool, delta_K)

#plot
fig,axs = plt.subplots(1,2,constrained_layout=True)
axs[0].plot(s, Na_R_cool)
axs[1].plot(s,K_R_cool)
fig.suptitle("Scattering Rate vs Intensity (δ = 6Γ)")
ax1=axs[0]
ax1.set_title("Sodium")
ax1.set_ylabel("Scattering Rate R")
ax1.set_xlabel(r"$\frac{I}{I_{sat}}$", size=15)
ax2=axs[1]
ax2.set_title("Potassium")
ax2.set_ylabel("Scattering Rate R")
ax2.set_xlabel(r"$\frac{I}{I_{sat}}$", size=15)
plt.show()
-------------------------------------------------------------------------------------------------------------------------------------------------------------

#libraries
import numpy as np
import matplotlib.pyplot as plt

#constants
#Na-3D
Na_diameter = 18e-3
Na_I_cool = 9.39e-3
Na_I_re = 10
Na_Gamma = 2*np.pi*9.795e6

#K-3D
K_diameter = 24e-3
K_I_cool = 17.5e-3
K_I_re = 13
K_Gamma = 2*np.pi*6.035e6

delta_K = np.linspace(-20*K_Gamma,20*K_Gamma,10000)
delta_K_plot = delta_K/K_Gamma

delta_Na = np.linspace(-20*Na_Gamma,20*Na_Gamma,10000)
delta_Na_plot = delta_Na/Na_Gamma

#Rabi Frequency
def Rabi(s, Gamma):
    Omega = Gamma * np.sqrt(s/(2))
    return Omega

#Scattering Rate equation
def Rate_equation(Gamma,Omega,delta):
    R = (Gamma/2)*((Omega**2/2)/(Gamma**2/4 + Omega**2/2 + delta**2))
    return R

#calculate
#Na
Na_Rabi_cool = Rabi(Na_I_cool, Na_Gamma)
Na_R_cool = Rate_equation(Na_Gamma, Na_Rabi_cool, delta_Na)

#K
K_Rabi_cool = Rabi(K_I_cool, K_Gamma)
K_R_cool = Rate_equation(K_Gamma, K_Rabi_cool, delta_K)

#plot
fig,axs = plt.subplots(1,2,constrained_layout=True)
axs[0].plot(delta_Na_plot, Na_R_cool)
axs[1].plot(delta_K_plot,K_R_cool)
fig.suptitle("Scattering Rate vs Detuning")
ax1=axs[0]
ax1.set_title("Sodium (s = 33e-3)")
ax1.set_ylabel("Scattering Rate R")
ax1.set_xlabel(r"$\delta/\Gamma$", size=15)
ax2=axs[1]
ax2.set_title("Potassium (s = 24e-3)")
ax2.set_ylabel("Scattering Rate R")
ax2.set_xlabel(r"$\delta/\Gamma$", size=15)
plt.show()

----------------------------------------------------------------------------------------------------------------------------------------------------------------------

#libraries
import numpy as np
import matplotlib.pyplot as plt
from mpl_toolkits import mplot3d

#constants
#Na-3D
Na_diameter = 18e-3
Na_I_cool = 33e-3
Na_I_re = 10e-3
Na_Gamma = 6e6

#K-3D
K_diameter = 24e-3
K_I_cool = 17e-3
K_I_re = 13e-3
K_Gamma = 23e6

delta = np.linspace(-500e6,500e6,1000)

s=np.linspace(0,100000,1000)

#Rabi Frequency
def Rabi(s, Gamma):
    Omega = Gamma * np.sqrt(s/(2))
    return Omega

#Scattering Rate equation
def Rate_equation(Gamma,Omega,delta):
    R = (Gamma/2)*((Omega**2/2)/(Gamma**2/4 + Omega**2/2 + delta**2))
    return R

#calculate
#Na
Na_Rabi_cool=np.zeros(1000)
Na_R_cool=np.zeros((1000,1000))
index1=0
for i in s:
  index2=0
  Na_Rabi_cool[index1] = Rabi(i, Na_Gamma)
  for j in delta:
    Na_R_cool[index1][index2] = Rate_equation(Na_Gamma, Na_Rabi_cool[index1], j)
    index2=index2+1
  index1=index1+1

#K
K_Rabi_cool = Rabi(s, K_Gamma)
K_R_cool = Rate_equation(K_Gamma, K_Rabi_cool, delta)

#plot
fig = plt.figure(figsize=(10,10), constrained_layout = True)
ax1 = plt.axes(projection ='3d')
result = ax1.contourf(delta, s, Na_R_cool)
ax1.set_title('Scattering Rate vs Detuning and s')
ax1.set_ylabel(r"$\frac{I}{I_{sat}}$",size=15)
ax1.set_xlabel(r"$\delta$",size=15)
ax1.set_zlabel("Scattering Rate R")

fig.colorbar(result)

plt.show()
