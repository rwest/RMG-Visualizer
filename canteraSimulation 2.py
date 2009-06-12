# this version updated by RHW in May 2008
"""
Batch reactor

        +------------------------+
        |                        |
        |                        |
        |  reactor      (T,V)    |
        |                        |
        |                        |
        +------------------------+

"""

######## SETTINGS

saveCSVever=False  # whether to save the [T,P,massfractions] profiles as CSV files (for MixMaster)
saveCSV=3 # saveCSV is False or an integer saying how many steps to do between saves (1=save every step)

if True:
    ctifile='chem' # the .cti extension is added later
    clearfigs=True
    format='-'


####### END OF SETTINGS


import sys
from Cantera import *
from Cantera.Reactor import *
from Cantera.Func import *
from math import *
import random
import pylab
import scipy
from pylab import *
from scipy import *
from scipy.optimize import leastsq

import time
def prettyTime(seconds):
	"""returns a string of the form HH:MM:SS when given a number of seconds"""
	return time.strftime("%H:%M:%S",time.gmtime(seconds))
	

# unit conversion factors
mm = 0.001
cm = 0.01
hours = 60*60
atm=101325

# reset()  # reset Cantera (clear all cached data)


# create a gas object
gas = importPhase(ctifile+'.cti', 'chem')

io2= gas.speciesIndex('O2(1)')
assert io2==0, 'Was expecting O2(1) to be first species'
maxtime=4*hours # seconds

iTarget=gas.speciesIndex('C16H26O2(119)')

fuelmassfracs="n-undecane(2):0.05, n-tridecane(3):0.19, SPC(4):0.11, n-hexadecane(5):0.25, SPC(6):0.12, n-nonadecane(7):0.18, n-heneicosane(8):0.10"
gas.set(Y=fuelmassfracs)

# Temperature
temperature= 135 +273 # Kelvin
#temperature=500
gas.set(T=temperature)


#The density of petroleum diesel is about 0.85 kg/l (7.09 lbs/gallon) en.wikipedia.org/wiki/Diesel
density= 0.85*1000 # [kg/m3]
gas.set(Rho=density)



#  Oxygen concentration:
# (doi: 10.1021/ie071481z ): At 25C and 1 atm of oxygen, O2 is sparsely 
# soluble in octane and decane and its liquid-phase mole fractions were 
# reported to be 2.05E-3 and 2.18E-3 respectively. Using this
# data and the partial pressure of O2 used under the given
# experimental conditions, the liquid-phase concentration of O2
# was calculated by assuming Henry's Law behavior. 
#  
#  using van't Hoff eqn with C=1700K gives k(140C) = 4.896*k(298K)
#  ie. concentration for same partial pressuer is ~5 times higher than at 298K
#  so for Po2 = 0.2 atm at 140C we have mole fraction of roughly 2E-3.
#
# not sure what to do for diesel, so use the same!
molefracO2=0.002



def setOneMoleFrac(gas,index,molefrac):
	molefracs=gas.moleFractions()
	molefracs*=(1-molefrac)/(1-molefracs[index])
	molefracs[index]=molefrac
	gas.set(X=molefracs)
	
setOneMoleFrac(gas,io2,molefracO2)

# create a reactor
reactor = Reactor(gas, volume = 1.0)
# create a network
sim = ReactorNet([reactor])  
           


if saveCSV:
    filehandle = open("MassFractionProfile.csv",'w')
    output=['Time (s)','T','P']
    output.extend(gas.speciesNames())
    output.extend(['foo'])
    writeCSV(filehandle, output)
    

if sim.time()>0.0:
    print "warning: t>0"
sim.setInitialTime(0.0)


def fuelmassfraction(gas):
	"""gives the mass fraction of unreacted fuel 
	   (components 1-7, because O2 is 0)"""
	fraction=0
	for i in range(1,8):
		fraction+=gas.massFraction(i)
	return fraction
	
	
initialFuelMassFraction=fuelmassfraction(gas)
def conversion(gas):
	return 1-fuelmassfraction(gas)/initialFuelMassFraction
	
steps=0

# observed deposit formation is 8E-4 grams per gram of fuel
finalConversion=8E-4

while conversion(gas)<finalConversion and sim.time()<maxtime:
    print "time: %.3gs (%s)\t frac used: %.2g"%(sim.time(),prettyTime(sim.time()),conversion(gas))
    
    sim.step(maxtime)
    steps+=1

    setOneMoleFrac(gas,io2,molefracO2)
    

    
    if saveCSV and steps%saveCSV:
        output=[sim.time(), gas.temperature(), gas.pressure()]
        #output.extend(gas.massFractions())
        output.extend([gas.massFraction(i) for i in range(gas.nSpecies())])
        output.extend([0.0])
        writeCSV(filehandle, output)

if saveCSV and not filehandle.closed: 
    filehandle.flush()
    filehandle.close()     
    
simtime=sim.time()
conversion=conversion(gas)

print "(%d steps)\t"%steps,

print "time (s)\t%.5g\t(%s)\tfrac used\t%.4g"%(simtime,prettyTime(simtime),conversion) 


