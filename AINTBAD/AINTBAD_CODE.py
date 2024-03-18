################################################################################################################
import numpy as np
import matplotlib.pyplot as plt
import matplotlib as mpl
import numpy as np
import scipy.optimize
import matplotlib.ticker 
import os
import warnings
from scipy.optimize import curve_fit
from scipy.optimize import minimize
from scipy import special
from scipy.integrate import quad
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter,AutoMinorLocator)
warnings.filterwarnings('ignore')
################################################################################################################
plt.rc('font', size = 30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)
Tm=273.15
minToSec=1/60
################################################################################################################
file=None
data=None
print("The AINTBAD code obtains A, B in Equation (7) that better fits the target data.")
dir_list = os.listdir('./input'); check=False
while(check is False):
  file_name=input("\nEnter the file's name (file format is .txt, the 1st column is temperature and the second is fice): ")
  check=file_name in dir_list
  if check is False:
     print("\nCould not find the file. Try again.")
     print("List of files:")
     for i in range(0,len(dir_list)):
        print(str(dir_list[i]))
file = str(file_name)
################################################################################################################
### User need to give cooling rate and droplet volume
coolingRate=input("\nWhat is the cooling rate in K/min? "); coolingRate=float(coolingRate)
dropletVolume=input("What is the mean droplet diameter in micrometers? "); dropletVolume=float(dropletVolume)
p0 = [3.45e40, 1.2]
################################################################################################################
initial=0
data = np.loadtxt("input/"+str(file),usecols=[0,1])
nsteps=data.shape[0]
temperature = np.copy(data[initial:nsteps,0]) 
SurvivalProb = np.copy(data[initial:nsteps,1])
dimensionlessTemperature=1-temperature/Tm
dimensionlessTemperature_fit=np.linspace(start=dimensionlessTemperature[0],stop=dimensionlessTemperature[-1],num=1000)
volumeFactor=(4.0*np.pi/3.0)*(dropletVolume*1e-4**3)
################################################################################################################
print("\nThis is the input data:\n")
fig1, ax = plt.subplots()
ax.plot(temperature,SurvivalProb,'o',color='red',ms=6,label='input data')
ax.set_ylabel('P')
ax.set_xlabel('$\delta_T$')
ax.xaxis.set_major_locator(MultipleLocator(1))
ax.legend(loc='best',ncol=1,fontsize=25)
fig1.set_size_inches(6.0, 6.0, forward=True)
fig1.tight_layout()
plt.show()
################################################################################################################
def exponent(dt,b):
    aux=b/((1-dt)*dt**2)
    if aux<1e-12 :
       return 1.0
    if aux>=1e-12 :
       return np.exp(-aux)
################################################################################################################
def integral(x, p):
    # p is array of A and B 
    x1=np.zeros(len(x))
    for i in range(0,len(x)):
        I1 = quad(exponent, 0.0001, x[i], args=(p[1])) # evaluating the integral 
        x1[i]=np.exp(-p[0]*(Tm/(coolingRate*minToSec))*volumeFactor*I1[0]) # eq.7
    return x1 
################################################################################################################
def objective_function(p, x, y): # p represent A and B, dT and Pexp respectively
     error = integral(x, p) - y
     score = np.sum(error**2)
     return score
################################################################################################################    
# Solving equation 7
################################################################################################################    
result = minimize(objective_function, p0, args=(dimensionlessTemperature, SurvivalProb), method='Nelder-Mead', options={'ftol':1e-4, 'maxiter':1000})
################################################################################################################
# Extract the optimized parameters
optimized_params = result.x
print("\nScore: ", result.fun)
fit=optimized_params
S=integral(dimensionlessTemperature_fit, fit)
A_param="{0:.4}".format(fit[0]); B_param="{0:.4}".format(fit[1])
################################################################################################################
fig2, ax = plt.subplots()
ax.plot(dimensionlessTemperature,SurvivalProb,'o',color='red',ms=6,label='input data')
ax.plot(dimensionlessTemperature_fit, S,'-',color='black',linewidth=2, label='fit to Eq. 7')
ax.set_ylabel('P')
ax.set_xlabel('$\delta_T$')
ax.legend(loc='best',ncol=1,fontsize=15)
fig2.set_size_inches(6.0, 6.0, forward=True)
fig2.tight_layout()
plt.show()
################################################################################################################
# Printing parameters
################################################################################################################
print("\n\nFrom the optimization\n")
print('A_param='+str(A_param)+', B_param='+str(B_param))
deltaT=dimensionlessTemperature_fit[(S>0.5-0.01) & (S<0.5+0.01)]
barrier=float(B_param)/float(np.mean(deltaT)**2)
barrier="{0:.4}".format(barrier)
T50=Tm*(1-np.mean(deltaT))
T50="{0:.4}".format(T50)
print("delta_G="+str(barrier)+" kT")
print("T50="+str(T50)+" K")
################################################################################################################
# Saving data file: temperature vs survival probability
################################################################################################################
np.savetxt("fit_"+str(file), np.column_stack((Tm*(1-dimensionlessTemperature_fit), S)), fmt='%10.8f', newline='\n') 
print("\n\nThe fitted data was saved in the same directory.\n")
################################################################################################################    
# Fitting the data to equation 8
################################################################################################################
check=True
while(check):
    decision=input("Do you want to fit the data to Equation 8? Yes or No? ")
    if decision == 'No':
        check=False; 
    if decision == 'Yes':
        check=False; 
################################################################################################################
if decision == 'Yes':
   def survivalProb(x, a, b): 
      return np.exp(-a*x*((np.exp(-b/(x*x))-(np.sqrt(np.pi*b)/x)*special.erfc(np.sqrt(b)/x))))
   fig3, ax = plt.subplots()
   ax.plot(dimensionlessTemperature,SurvivalProb,'o',color='red',ms=6,label='input data')
   fit, cov = scipy.optimize.curve_fit(survivalProb, dimensionlessTemperature, SurvivalProb, p0 = [3.45e20, 1.07],maxfev=10000)
   A_param="{0:.4}".format(fit[0]); B_param="{0:.4}".format(fit[1])
   S=survivalProb(dimensionlessTemperature_fit, *fit)
   ax.plot(dimensionlessTemperature_fit, S,'-',color='black',linewidth=2, label='fit to Eq. 8')
   ax.legend(loc='best',ncol=1,fontsize=15)
   ax.set_ylabel('P')
   ax.set_xlabel('$\delta_T$')
   fig3.set_size_inches(7.0, 7.0, forward=True)
   fig3.tight_layout()
   plt.show()
   print("\n\nFrom the optimization\n")
   print('\nA_param='+str(A_param)+', B_param='+str(B_param))
   deltaT=dimensionlessTemperature_fit[(S>0.5-0.01) & (S<0.5+0.01)]
   barrier=float(B_param)/float(np.mean(deltaT)**2)
   barrier="{0:.4}".format(barrier)
   T50=Tm*(1-np.mean(deltaT))
   T50="{0:.4}".format(T50)
   print("delta_G="+str(barrier)+" kT")
   print("T50="+str(T50)+" K")
else:
   print("\n\nYou are done!")
plt.show()
