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
from matplotlib.ticker import (MultipleLocator, FormatStrFormatter, AutoMinorLocator)
from scipy.special import gamma as gamma_function
warnings.filterwarnings('ignore')
################################################################################################################
plt.rc('font', size = 30)
plt.rc('xtick', labelsize=20) 
plt.rc('ytick', labelsize=20)
Tm=273.15
minToSec=1/60
prefactorVolume=(4.0*np.pi/3.0)

################################################################################################################
print("The IPA code uses any parametrization for the nucleation rate and a size distribution to generate survival probability curves.")
dir_list = os.listdir('./input'); check=False
print("\n Step 1 \n")
while(check is False):
  print("1 - CNT parametrization from Qiu et. al 2019")
  print("2 - CNT parametrization from Koop and Murray 2016")
  print("3 - A and B from the AINTBAD model")
  option=input("\nEnter the number for the nucleation rate parametrization option: "); option=int(option)
  if option < 1 or option > 3:
    check=False
  else:
    check=True
################################################################################################################
if option == 1:
  file = './input/CNT_parametrization_Qiu2019.txt'
  initial=0
  data = np.loadtxt(file,usecols=[0,1,2])
  nsteps=data.shape[0]-110#-170
  temperature = np.copy(data[initial:nsteps,0]) # in K
  nucleationBarriers = np.copy(data[initial:nsteps,1]) # in kT
  nucleationRate = np.copy(data[initial:nsteps,2]) # in (cm^(-3)s^(-1))
  dimensionlesstemperature=(1-temperature/Tm)
################################################################################################################
if option == 2:
  c0=-3020.684; c1=-425.921; c2=-25.9779; c3=-0.868451; c4=-1.66203*1e-2; c5=-1.71736*1e-4; c6=-7.46953*1e-7
  temperature=np.linspace(start=230,stop=250,num=100)
  dimensionlesstemperature=(1-temperature/Tm)
  nucleationRate=10**(c0+c1*(temperature-Tm)+c2*(temperature-Tm)**2+c3*(temperature-Tm)**3+c4*(temperature-Tm)**4+c5*(temperature-Tm)**5+c6*(temperature-Tm)**6)
################################################################################################################
if option == 3:
  A_param = input("\nWhat is the parameter A: ");
  A_param = float(A_param)
  B_param = input("\nWhat is the parameter B: ");
  B_param = float(B_param)
  temperature=np.linspace(start=230,stop=250,num=100)
  dimensionlesstemperature=(1-temperature/Tm)
  nucleationRate=A_param*np.exp(-B_param/(dimensionlesstemperature**2))
################################################################################################################
fig, ax = plt.subplots()
ax.plot(temperature,nucleationRate,'-',color='red',ms=5,label='input')
ax.set_yscale('log')
ax.grid(True)
ax.legend(loc='best',ncol=1,fontsize=15)
ax.set_ylabel('J ($cm^{-3} s^{-1}$)')
ax.set_xlabel('T (K)')
fig.set_size_inches(7.0, 7.0, forward=True)
fig.tight_layout()
fig.savefig('nucleationRate.pdf', dpi=300)
################################################################################################################
print("\n Step 2 \n")
# parameters of the distribution
micrometer_to_centimeter=1e-4
npoints=10000
print("These are the available droplet diameter distributions:")
print("1 - Uniform distribution")
print("2 - Gaussian distribution")
print("3 - Gamma distribution")
check=False
while(check is False):
   option = input("Choose the type of distribution: "); option = int(option)
   if option < 1 or option > 3:
      check=False
   else:
      check=True
if option==1:
    print("Give the upper and lower boundaries (in micrometers)")
    D_lower = input("Lower diameter: ");  D_lower=float(D_lower)*micrometer_to_centimeter
    D_upper = input("Upper diameter: ");  D_upper=float(D_upper)*micrometer_to_centimeter
if option==2:
   print("Give the parameters of the Gaussian distribution (in micrometers)")
   user_mode=input("Mean: "); user_mode=float(user_mode)*micrometer_to_centimeter
   user_scale=input("Standard deviation: "); user_scale=float(user_scale)*micrometer_to_centimeter
if option==3:
    print("Give the parameters of the Gamma distribution")
    shape_parameter = input("Shape parameter: "); shape_parameter=float(shape_parameter)
    scale_parameter = input("Scale parameter: "); scale_parameter=float(scale_parameter)
################################################################################################################
Dmin=0*micrometer_to_centimeter; Vmax=1000*micrometer_to_centimeter
################################################################################################################
def uniform_pdf(x, a, b):
    """
    Probability Density Function (PDF) of a continuous uniform distribution
    in the range [a, b].
    """
    pdf=np.zeros(len(x))
    for i in range(0,len(x)):
       if ((x[i]>=a) and (x[i]<b)):
        pdf[i]=1/(b - a)
    return pdf
################################################################################################################
def Gaussian_PDF(x,mode,scale):
    pdf=np.zeros(len(x))
    pdf=(1/(scale*np.sqrt(2*np.pi)))*np.exp(-0.5*((x-mode)/scale)**2)
    return pdf
################################################################################################################
def Gamma_PDF(x, a, b):
    numerator = b**a * (x)**(a-1) * np.exp(-b*x)
    denominator = gamma_function(a)
    y = numerator / denominator
    x=x[y>0]
    y=y[y>0]
    return x, y
################################################################################################################
fig, ax = plt.subplots()
if option==1:
   D_dist=np.linspace(start=0,stop=D_upper*50,num=npoints)
   freq_D_dist=uniform_pdf(D_dist,D_lower,D_upper)
   ax.plot(D_dist,freq_D_dist,'-',color='red')
   ax.set_xlim([0, D_upper+2*D_upper])
if option==2:
   D_dist = np.linspace(start=0, stop=user_mode * 50, num=npoints)
   freq_D_dist=Gaussian_PDF(D_dist,user_mode,user_scale); ii = np.isnan(freq_D_dist); freq_D_dist[ii] = 0
   ax.plot(D_dist*2,freq_D_dist,'-',color='red')
   ax.set_xlim([0, user_mode + 2 * user_mode])
if option==3:
   D_dist = np.linspace(start=0, stop=shape_parameter * 20, num=npoints)
   D_dist,freq_D_dist=Gamma_PDF(D_dist,shape_parameter,scale_parameter); ii = np.isnan(freq_D_dist); freq_D_dist[ii] = 0
   ax.plot(D_dist,freq_D_dist,'-',color='red')
################################################################################################################
ax.set_ylabel('PDF')
ax.set_xlabel('Droplet diameter ($cm$)')
fig.set_size_inches(7.0, 7.0, forward=True)
fig.tight_layout()
fig.savefig('diameterDistribution.pdf', dpi=300)
np.savetxt("Diameter_distribution.txt", np.column_stack((D_dist,freq_D_dist)), fmt='%10.8f', newline='\n')
################################################################################################################
print("\n Step 3 \n")
coolingRate=input("What is the cooling rate in K/min? "); coolingRate=float(coolingRate)
################################################################################################################
# OUTPUT
################################################################################################################
def trapezoidal_discrete_integral(x, y, a, b):
    """
    Approximate the integral of discrete data within the range [a, b] using the trapezoidal rule.

    Parameters:
    x (array): Array of x-values (data points).
    y (array): Array of corresponding y-values.
    a (float): Lower integration limit.
    b (float): Upper integration limit.

    Returns:
    integral (float): Approximated integral within [a, b].
    """
    mask = (x >= a) & (x <= b)
    x_interval = x[mask]
    y_interval = y[mask]
    integral = np.trapz(y_interval, x_interval)
    return integral
################################################################################################################
def func(deltaT,J):
    ycom=np.zeros(len(deltaT)) # initializing the output array to zeros
    for i in range(0,len(deltaT)): # loop over length of the data to find integral
        # Slice the data points up to the current deltaT[i] and use Trapezoidal's rule
        ycom[i]=trapezoidal_discrete_integral(deltaT, J, 0, deltaT[i])
    vfac=np.zeros(len(deltaT))
    vfac1=0
    for i in range(0,len(D_dist)):
        dR=D_dist[i]/2
        Pv=freq_D_dist[i]
        # sum of survival probability for all volumes
        vfac=vfac+Pv*np.exp(-prefactorVolume*(Tm/(coolingRate*minToSec))*((dR)**3)*ycom)
        vfac1=vfac1+Pv
    return vfac/vfac1 # normalization of probabilities
################################################################################################################
zipped = list(zip(dimensionlesstemperature, nucleationRate))
# Sort by the first array
zipped.sort(key=lambda x: x[0])
# Unzip the sorted zipped list
dimensionlesstemperature, nucleationRate = zip(*zipped)
# Convert the result back to NumPy arrays
dimensionlesstemperature = np.array(dimensionlesstemperature)
nucleationRate = np.array(nucleationRate)
SurvivalProbability=func(dimensionlesstemperature,nucleationRate)
################################################################################################################
fig, ax = plt.subplots()
ax.plot(dimensionlesstemperature,SurvivalProbability,'-o',color='blue')
def survivalProb(x, a, b):
    return np.exp(-a*x*((np.exp(-b/(x*x))-(np.sqrt(np.pi*b)/x)*special.erfc(np.sqrt(b)/x))))
fit, cov = scipy.optimize.curve_fit(survivalProb, dimensionlesstemperature, SurvivalProbability, p0 = [3.45e40, 1.07],maxfev=100000)
A_apparent="{0:.4}".format(fit[0]); B_apparent="{0:.4}".format(fit[1])
dimensionlessTemperature_fit=np.linspace(start=dimensionlesstemperature[0],stop=dimensionlesstemperature[-1],num=1000)
SurvivalProbability_fit=survivalProb(dimensionlessTemperature_fit, *fit)
deltaT=dimensionlessTemperature_fit[(SurvivalProbability_fit>0.5-0.01) & (SurvivalProbability_fit<0.5+0.01)]
T50=(1-np.mean(deltaT))*Tm
###############################################################################################################
T50="{0:.5}".format(T50)
print("\n OUTPUT")
print("\nThe plot and data of the generated survival probability were saved in the same directory.")
print("\nT50 = "+str(T50)+' K')
ax.grid(True)
ax.set_ylabel('Survival probability')
ax.set_xlabel('$\delta_T$')
fig.set_size_inches(8.0, 8.0, forward=True)
fig.tight_layout()
fig.savefig('survivalProbability.pdf', dpi=300)
np.savetxt("SurvivalProbability.txt", np.column_stack(((1-dimensionlesstemperature)*Tm,SurvivalProbability)), fmt='%10.8f', newline='\n')

