# -*- coding: utf-8 -*-
import sys
import math
import cmath


try:
	from matplotlib import pylab
	from pylab import *
	#import pylab
	pylab_available = "YES"
except:
	pylab_available = "NO"


####### PROGRAM DESCRIPTION ###########################
#
### CRest.py ###
# written by Robinson Peric 
# at Hamburg University of Technology (TUHH)
# version: 9 august 2019
#
#
### What the code does ###
# The code computes analytical predictions for the reflection coefficient
# for long-crested free-surface wave propagation when using 
# forcing zones (e.g. wave absorbing layers, damping zones, sponge layers, etc.) as described in [1]
#
### Please note ###
# Recommendations are given in [1] how to tune the forcing zone's parameters depending on the waves. See also the manual distributed along with the code.
# Although theory predictions are often quite accurate [1], please keep in mind that every theory has its limitations! Tuning the forcing zone parameters according to the present theory does NOT guarantee that the actual reflection coefficient in the simulation will equal the prediction. Other mechanisms can also lead to undesired wave reflections. Thus it is recommended to select zone thickness slightly larger than minimum necessary.
# The theory predictions were found to be of satisfactory accuracy for regular, irregular, linear, highly nonlinear, 2D-, and 3D-waves [1,2,3]. Indeed, the present code for 1D-wave-propagation (i.e. 2D-flow with wave propagation in one direction) was found to be of satisfactory accuracy for tuning forcing zones in typical 3D-flow simulations [2,3].
# Future research will help to better assess the accuracy of the theory's predictions.
# 
#
### Reference ###
# [1] Perić, R., & Abdel-Maksoud, M. (2018). Analytical prediction of reflection coefficients for wave absorbing layers in flow simulations of regular free-surface waves. Ocean Engineering, 147, 132-147.
# [2] Perić, R., & Abdel-Maksoud, M. (2019). Reducing Undesired Wave Reflection at Domain Boundaries in 3D Finite Volume–Based Flow Simulations via Forcing Zones. Journal of Ship Research.
# [3] Perić, R., Vukčević, V., Abdel-Maksoud, M., & Jasak, H. (2018). Tuning the Case-Dependent Parameters of Relaxation Zones for Flow Simulations With Strongly Reflecting Bodies in Free-Surface Waves. arXiv preprint arXiv:1806.10995.
#
### How to use the program ###
# Call program like this:
# python ./CRest.py
#
# to run the program in batch-mode (i.e. non-interactive), simply pass the corresponding arguments
#
### Requirements ###
# The programming language python version 2.7 or >=3.0  (https://www.python.org/downloads/) must be installed. 
# It is recommended to also have matplotlib (https://matplotlib.org/users/installing.html) installed. Then the code will plot results beautifully.
#
#
### Bug fixes ###
# The 2017 version of this code was inaccurate when forcing only the vertical velocities in shallow water. This has been amended, and now the code should work for all water depths.
#
# Please report errors to robinson.peric@tuhh.de
# 
#####################################################


# calculate the linear wave theory solution for the given period and water depth
def period(T, h):
	w = 2.0 * math.pi / T	# angular wave frequency
	# iterate to obtain corresponding intermediate depth wavelength
	tmp = (w**2)/9.81	# initial guess: deep water wave number
	for i in range(50):	# 50 iterations is usually enough...
		tmp2 = w*w/(9.81*math.tanh(tmp*h))
		tmp = (tmp+tmp2)/2.0	#bisection
	return 2.0 * math.pi / tmp	# return wavelength



# batch mode
if (len(sys.argv)-1 == 6):
	T = float(sys.argv[1])
	h = float(sys.argv[2])
	L = float(sys.argv[3])
	forcedEQs = str(sys.argv[4])
	zone_thickness_xd = float(sys.argv[5])
	blend = float(sys.argv[6])
	#pylab_available = "NO"
	
elif (len(sys.argv)-1 == 7):
	T = float(sys.argv[1])
	h = float(sys.argv[2])
	L = float(sys.argv[3])
	forcedEQs = str(sys.argv[4])
	zone_thickness_xd = float(sys.argv[5])
	blend = float(sys.argv[6])
	n = float(sys.argv[7])
	#pylab_available = "NO"
	
# interactive mode
else:
	T = float(input("\n\n\nPlease enter wave period (s):\n"))

	h = float(input("\nPlease enter water depth (m):\n"))

	L = float(input("\nPlease enter wavelength (m):\n"))

	w = 2 * math.pi / T			# angular wave frequency
	wavenumber_k = 2 * math.pi / L		# wave number
	c = L / T				# phase velocity
	
	# check wavelength
	Llin0L = L / period(T,h)
	if ( ( Llin0L < 0.99 ) | ( Llin0L > 1.2 ) ):
		print("\nThis wavelength is " + str(Llin0L*100) + "% percent of the prediction by linear wave theory.")
		print("Even for nonlinear waves, the wavelength rarely exceeds linear wave theory predictions by more than 20%. Please check input parameters.")
		try:
			input("Continue? Press any key:")
		except:
			print("")

	print("\nPlease indicate in which governing equations forcing source terms are applied, i.e. which flow quantities should be forced?")
	print("Enter 1 for horizontal velocities")
	print("Enter 2 for vertical velocity")
	print("Enter 3 for volume fraction")
	print("Enter 4 for horizontal + vertical velocity")
	print("Enter 5 for horizontal + vertical velocity + volume fraction")
	forcedEQs = str(input("\nPlease enter number:\n"))

	zone_thickness_xd = float(input("\nPlease enter forcing zone thickness (m):\n"))
	
	print("\nPlease select a blending function b(x)\n\n")
	print("--- Commonly used b(x) ---")
	print("Enter 1 for CONSTANT blending, i.e. b(x) = 1")
	print("Enter 2 for LINEAR blending, i.e. b(x) = x")
	print("Enter 3 for QUADRATIC blending, i.e. b(x) = x^2")
	print("Enter 4 for COSINE-SQUARED blending, i.e. b(x) = ( cos(x) )^2")
	print("Enter 5 for EXPONENTIAL blending, i.e. b(x)= ( exp(x^2) - 1 ) / ( exp(1) - 1 )")
	print("--- Expert b(x) ---")
	print("Enter 6 for POWER blending with exponent n, i.e. b(x) = x^n")
	print("Enter 7 for EXPONENTIAL blending with exponent n, i.e. b(x)= ( exp(x^n) - 1 ) / ( exp(1) - 1 )")
	print("Enter 8 for COSINE-SQUARED blending with exponent n, i.e. ( ( cos(x) )^2 )^n")
	print("Enter 9 for CUSTOM blending")
	blend = float(input("\nPlease enter number of blending function:\n"))
	if ( blend > 5.0):
		n = float(input("\nPlease enter blending exponent n:"))

w = 2 * math.pi / T			# angular wave frequency
wavenumber_k = 2 * math.pi / L		# wave number
c = L / T				# phase velocity

if ( h/L < 1.0 ):
	Ekinz0Ekinx = ( math.sinh( 2.0 * wavenumber_k * h ) - 2.0 * wavenumber_k * h ) / ( math.sinh( 2.0 * wavenumber_k * h ) + 2.0 * wavenumber_k * h )
else:
	Ekinz0Ekinx = 1.0
if (forcedEQs == "1"):
	forcedEQs = 1.0				# Ekinx / Ekinx
elif (forcedEQs == "2"):
	forcedEQs = Ekinz0Ekinx			# Ekinz / Ekinx
elif (forcedEQs == "3"):
	forcedEQs = 1.0 + Ekinz0Ekinx		# Epot/Ekinx =  Ekinx / Ekinx + Ekinz / Ekinx
elif (forcedEQs == "4"):
	forcedEQs = 1.0 + Ekinz0Ekinx		# ( Ekinx + Ekinz ) / Ekinx
elif (forcedEQs == "5"):
	forcedEQs = 2.0 * ( 1.0 + Ekinz0Ekinx )		# ( Ekinx + Ekinz + Epot ) / Ekinx
else:
	print("Wrong entry.")
	sys.exit(1)



### default settings
csv_file_separator = " "		# alt. "," or ";"
factor =1.05				# resolution of plot: very fine (1.01), very coarse (2.0)
gammaMin=10**-4/T*forcedEQs	# range of forcing strength gamma for which reflection coefficient is computed
gammaMax=10**7/T*forcedEQs	# -"-
dampres=200				# number of piece-wise constant blending zones, into which forcing zone is subdivided
Lx = 2.0*L 				# for plots: domain size outside the forcing zone







###### VARIABLE DECLARATIONS AND INITIALIZATIONS ###########################

# initialize vectors
xd = [0.0] * ( int(dampres) + 2 )		 	# initialize vector for thicknesses of zone segments, xd[1] is Lx
k = [0.0] * ( int(dampres) + 2 )		 	# initialize vector for wave number, k[1] is 2*pi/L
Ct = [0.0] * ( int(dampres) + 2 )		 	# initialize vector for transmission coefficient, Ct[1] is coefficient at the entrance to the forcing zone
beta = [0.0] * ( int(dampres) + 2 )			# initialize vector needed to compute Cr
Cr = [0.0] * ( int(dampres) + 2 )		 	# initialize vector for reflection coefficient, Cr[1] is the global reflection coefficient (i.e. for the entrance to the forcing zone)

# write initial conditions
dx = zone_thickness_xd / dampres	# thickness of single 'cell' within forcing zone
xd[1] = Lx
for i in range(2,dampres+2,1):
	xd[i] = dx
k[0] = wavenumber_k
k[1] = wavenumber_k
Ct[0] = 1.0
Cr[0] = 0.0
Ct[dampres+1] = 0.0
Cr[dampres+1] = 1.0


###### FUNCTION DECLARATIONS ###########################

# return absolute value of complex number x
def abs(x):
	return math.sqrt( (x.real)**2 + (x.imag)**2 )	

# evaluate blending function
def b(x):
	if (x >= Lx):

		# constant	
		if (blend == 1):
			return 1.0	# constant blending

		# linear
		elif (blend == 2):
			return  ( x - Lx ) / ( ( Lx + zone_thickness_xd ) - Lx )

		# quadratic
		elif (blend == 3):
			return  ( ( x - Lx ) / ( ( Lx + zone_thickness_xd ) - Lx ) )**2

		# cosine-squared
		elif (blend == 4):
			return math.cos( 0.5 * math.pi + 0.5 * math.pi * ( x - Lx ) / ( ( Lx + zone_thickness_xd ) - Lx ) )**2

		# exponential
		elif (blend == 5):
			return (math.exp( ( ( x - Lx ) / ( ( Lx + zone_thickness_xd ) - Lx )  )**2 ) -1) / ( math.exp(1) - 1 )
		# power with exponent n
		elif (blend == 6):
			return ( ( x - Lx ) / ( ( Lx + zone_thickness_xd ) - Lx ) )**n
		# exponential with exponent n
		elif (blend == 7):
			return (math.exp( ( ( x - Lx ) / ( ( Lx + zone_thickness_xd ) - Lx )  )**n ) -1) / ( math.exp(1) - 1 )
		# cosine-squared with exponent n
		elif (blend == 8):
			return ( math.cos( 0.5 * math.pi + 0.5 * math.pi * ( x - Lx ) / ( ( Lx + zone_thickness_xd ) - Lx ) )**2 )**n

		# custom
		elif (blend == 9):
			print("\nERROR: You have not yet implemented a custom blending function.")
			sys.exit(1)	# please replace this line and the line above by a return-statement for the desired custom blending function
		
		else:
			print("ERROR: Specified blending function not found.")
			sys.exit(1)
	
	else:
		return 0.0

# return \beta_{i}
def fbeta(i):
    return ( ( 1.0 + Cr[i] * cmath.exp( 1j * k[i] * xd[i] * 2.0 )  ) / ( 1.0 -  Cr[i] * cmath.exp( 1j * k[i] * xd[i] * 2.0 ) ) )

# return C_{T,i}
def fCt(i):
    return ( 1.0 - Cr[i] ) / (  1.0 - Cr[i+1] * cmath.exp( 1j * k[i+1] * xd[i+1] * 2.0 )  ) 

# return C_{R,i}
def fCr(i):
    return ( -k[i] + k[i+1] * beta[i+1] ) / ( k[i] + k[i+1] * beta[i+1] )





###### MAIN ##############################################
gamma = gammaMin
G=[]
CR=[]

gamma_opt = 0.0
C_R_opt = 1000.0

# create file
f = open("C_R.csv", 'w')
f.close()

# calculate C_R for different gamma values
while (gamma < gammaMax):
	tmp = Lx
	for i in range(2,dampres+2,1):
		k[i] = cmath.sqrt(  (w**2)/ ( c**2 ) + (1j * w * gamma * b( tmp + 0.5 * xd[i] ) ) / ( (c)**2 ))
		tmp += xd[i]

	# calculate transmission and reflection coefficients
	for i in range(dampres,0,-1):   # go cell-by-cell from the end of the forcing zone to its entrance. E.g. for 3 zones, start with zone 3, then 2, then 1.
		beta[i+1] = fbeta(i+1)
		Cr[i] = fCr(i)
		Ct[i] = fCt(i)

	# write reflection coefficients to file
	f = open("C_R.csv", 'a')
	f.write(str(gamma/forcedEQs) + csv_file_separator + str(abs(Cr[1])) +  "\n")
	G.append(gamma/forcedEQs)
	CR.append(abs(Cr[1]))
	f.close()

	if (abs(Cr[1]) < C_R_opt):
		gamma_opt = gamma/forcedEQs
		C_R_opt = abs(Cr[1])

	gamma = gamma * factor


# if pylab and matplotlib are installed, open window and plot C_R(gamma)
if (pylab_available == "YES"):
	pylab.figure(figsize=(20,8))
	pylab.xlabel('$\gamma\ (1/\mathrm{s})$ ', {'fontsize': 20})
	pylab.ylabel('$C_R$', {'fontsize': 20})
	pylab.title('Reflection coefficient $C_R$ versus forcing strength $\gamma$')
	pylab.semilogx()
	pylab.semilogy()
	pylab.plot(G,CR, linewidth=1.5)
	pylab.show()

print("\n\nProgram finished. \nReflection coefficients were written to file C_R.csv\n")

print("Global minimum for reflection coefficient C_R_opt = "+str(C_R_opt*100)+"%\noccurs at forcing strength gamma = "+str(gamma_opt)+" s^-1\n")

sys.exit(0)

