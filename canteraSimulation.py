# this version updated by RHW in May 2008
"""
Batch reactor
with peroxide!

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


molefracPeroxide=1e-6


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

reset()  # reset Cantera (clear all cached data)



for molefracPeroxide in [1]:# [10**i for i in range(-9,-1)]:
	# create a gas object
	gas = importPhase(ctifile+'.cti', 'chem')
	
	io2= gas.speciesIndex('O2(1)')
	assert io2==0, 'Was expecting O2(1) to be first species'
	maxtime=4*hours # seconds
	
#	iPeroxide=gas.speciesIndex('peroxydcbz(9)')
	
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
	
	
	
	#setOneMoleFrac(gas,iPeroxide,molefracPeroxide)
	#print "initial mole fraction of Peroxide:",molefracPeroxide
	
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
	    #print "time: %.3gs (%s)\t frac used: %.2g"%(sim.time(),prettyTime(sim.time()),conversion(gas))
	    
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
	#print "time (s)\t%.5g\t(%s)\tfrac used\t%.4g"%(simtime,prettyTime(simtime),conversion) 
	print "molefracPeroxide\t%g\ttime(s)\t%g"%(molefracPeroxide,simtime)



assert False, "stop"

# temperature=1500
# pressure=1.0e5

f=figure(1)
title("Effect of pressure")
Ts=range(901,2500,50)
Ps=[1e5*10**(i) for i in arange(-1,2.1,1)] # was ,0.25)]
Ps.insert(0,0) # Pratsinis
Ps.insert(0,-1) # Raghavan
if clearfigs: clf()
meshTemps=[]
meshPressures=[]
meshLogKs=[]
for P in Ps:
    temps=[]
    logKs=[]
    
    if P==-1:
        if not clearfigs: continue
        label="Raghavan (2001)"
        for T in [646, 2044]:
            temps.append(T)
            logKs.append(RaghavanLogK(T))    
    elif P==0:
        if not clearfigs: continue
        label="Pratsinis (1990)"
        for T in [973, 1273]:
            temps.append(T)
            logKs.append(PratsinisLogK(T))
    else:
        for T in Ts:
            label="%.2g bar"%(P/1e5)
            X="TiCl4:0.01, O2:0.05, Ar:0.94"
            s=Simulation(T,P,X)
            print "T (K)\t%d\tP (bar)\t%.2g\t "%(T,P/1e5),
            logK=s.getLogK()
            if logK>-Inf:
                temps.append(T)
                logKs.append(logK)
    
        meshTemps.append(temps)
        meshPressures.append([P for T in temps]) 
        meshLogKs.append(logKs)
    temps=array(temps)
    logKs=array(logKs)
    plot( 1000.0/temps, logKs, format, label=label, linewidth=2 )
    
ylabel("ln k")
xlabel("1000/T")
if clearfigs:
    legend(loc=0)
    drawPratsinisBox()
    drawRaghavanBox()


### 3D plot
import matplotlib.axes3d as p3
f=figure(2)
ax=p3.Axes3D(f)
show()


##RAGHAVAN'S
#temps=array([  1046.82525236,  1073.85699043,
#        1151.80924153,  1257.77744668,  1298.76346761,  1406.39240553,
#        1439.28607796,  1514.65377135,  1604.81167741,  1645.52873364,
#        1715.11423198,  1740.23047543,  1806.84217932,  2043.4652303 ])
#x=1000.0/temps
#y=log(array([   551205.35156326,
#         429845.99298875,   736792.54475485,   964331.41922609,
#         887997.53785585,  1198889.6613466 ,   942090.5684987 ,
#        1436794.95733164,  1706284.57763703,  1810558.77650401,
#        1458560.42080549,  1526135.31697436,  1695978.55319514,
#        2384852.8012917 ]))/log(10.0)
## z given by Raghavan's arrhenius rate       
#z=array([RaghavanLogK(T) for T in temps])


        
# z calculated from Raghavan's measurement (assuming lorentzian, with his Eact)
def RaghavanLogKexpt(Tp,NNo,tc):
    Eact=108681 # J/mol
   # Eact= 88800 # J/mol (pratsinis)
    molargasconstant=8.314472 # J/mol/K
    # NNo=NNo-0.1;
    f=0.2
    if NNo<f: return 15.0
    if NNo>=1.0: return -15.0
    numerator= -1.0*( log((NNo-f)/(1-f)) )
    denominator= tc * sqrt( pi * molargasconstant * Tp / Eact)
    #print "Tp=%g\tNNo=%g\tTc=%g"%(Tp,NNo,tc)
    k=numerator/denominator
    #print "k=",numerator,'/',denominator,"\t",k
    #print "log(k)=",log(k)
    return log(k)
    
temps=array([  646.,   757.,  1047.,  1074.,  1152.,  1258.,  1299.,  1407.,
        1440.,  1516.,  1607.,  1648.,  1718.,  1743.,  1809.,  2044.])       
measuredNNos=array([ 1.  ,  0.97,  1.02,  0.9 ,  0.96,  1.04,  0.94,  0.49,  0.44,
        0.6 ,  0.33,  0.35,  0.16,  0.24,  0.17,  0.17])
reactiontimes=array([ 0.467,  0.166,  0.058,  0.063,  0.04 ,  0.032,  0.032,  0.027,
        0.023,  0.024,  0.019,  0.019,  0.018,  0.018,  0.013,  0.01 ])
pressures=array([  165163.21252451,   252971.75364402,   551205.35156326,
         429845.99298875,   736792.54475485,   964331.41922609,
         887997.53785585,  1198889.6613466 ,   942090.5684987 ,
        1436794.95733164,  1706284.57763703,  1810558.77650401,
        1458560.42080549,  1526135.31697436,  1695978.55319514,
        2384852.8012917 ])
x=1000.0/temps 
y=log(pressures)/log(10.0)
z=array([RaghavanLogKexpt( temps[i], measuredNNos[i], reactiontimes[i]) for i in range(temps.size)])

RaghavanMeasuredInverseTs=x
RaghavanMeasuredPressures=y
RaghavanMeasuredLogKs=z

ax.plot3D(x,y,z,'o-',label="Raghavan (2001)", linewidth=2)
save('lineXrag.out',x)
save('lineYrag.out',y)
save('lineZrag.out',z)

temps=arange(973.0,1274,50)
x=1000.0/temps
y=log(array([OneAtm for T in temps]))/log(10.0)
z=array([PratsinisLogK(T) for T in temps])
ax.plot3D(x,y,z,label="Pratsinis (1990)", linewidth=2)

PratData=array([[ 0.883, 1.1112 ],[ 0.8981, 1.0121 ],[ 0.9149, 0.9306 ],[ 0.9317, 0.6962 ],[ 0.9499, 0.4882 ],[ 0.9677, 0.1884 ],[ 0.9832, -0.2687 ],[ 1.0073, -0.2841 ]])
x=PratData[:,0]
y=log(array([OneAtm for t in x]))/log(10.0)
z=PratData[:,1]
ax.plot3D(x,y,z,'o-',label="Pratsinis (1990) 1:5", linewidth=2)
save('lineXprat.out',x)
save('lineYprat.out',y)
save('lineZprat.out',z)

PratsinisMeasuredInverseTs=x
PratsinisMeasuredPressures=y
PratsinisMeasuredLogKs=z


def plotRaghavanMeasurements():
    temps=array([  646.,   757.,  1047.,  1074.,  1152.,  1258.,  1299.,  1407.,
            1440.,  1516.,  1607.,  1648.,  1718.,  1743.,  1809.,  2044.])       
    measuredNNos=array([ 1.  ,  0.97,  1.02,  0.9 ,  0.96,  1.04,  0.94,  0.49,  0.44,
            0.6 ,  0.33,  0.35,  0.16,  0.24,  0.17,  0.17])
    reactiontimes=array([ 0.467,  0.166,  0.058,  0.063,  0.04 ,  0.032,  0.032,  0.027,
            0.023,  0.024,  0.019,  0.019,  0.018,  0.018,  0.013,  0.01 ])
    pressures=array([  165163.21252451,   252971.75364402,   551205.35156326,
             429845.99298875,   736792.54475485,   964331.41922609,
             887997.53785585,  1198889.6613466 ,   942090.5684987 ,
            1436794.95733164,  1706284.57763703,  1810558.77650401,
            1458560.42080549,  1526135.31697436,  1695978.55319514,
            2384852.8012917 ])
            
    x=1000.0/temps 
    for i in range(temps.size):
        temp=temps[i]
        measuredNNo=measuredNNos[i]
        pressure=pressures[i]
        reactiontime=reactiontimes[i]
        errorInNNo=0.1
        x=1000.0/temp
        y=RaghavanLogKexpt( temp, measuredNNo, reactiontime )
        ymin =RaghavanLogKexpt( temp, measuredNNo+errorInNNo, reactiontime )
        ymax =RaghavanLogKexpt( temp, measuredNNo-errorInNNo, reactiontime )
        ybar=array([ymin,ymax])
        xbar=array([x, x])
        if y>9:
            plot( [x], [y], 'b^', label="Raghavan", linewidth=1 )
        elif y<-5:
            plot( [x], [y], 'bv', label="Raghavan", linewidth=1 )
        else:
            plot( [x], [y], 'ro', label="Raghavan", linewidth=1 )
        
        plot( xbar, ybar, 'b-_', label="Raghavan", linewidth=1 )
        print "Tp\t%dK\tNNo\t%g\tt\t%g\tk_eff\t%.4g\tmin\t%.4g\tmax\t%.4g"%(temp,measuredNNo, reactiontime, y,ymin,ymax)
    
def plotPratsinisMeasurements():
    plot( PratsinisMeasuredInverseTs, PratsinisMeasuredLogKs, 'go', label="Pratsinis", linewidth=2 )


X=1000.0/array(meshTemps)
Y=log(array(meshPressures))/log(10.0)
Z=array(meshLogKs)
surf1=ax.plot_surface(X,Y,Z)
frame1=ax.plot_wireframe(X,Y,Z)

Zprat=array([[PratsinisLogK(1000.0/i) for i in j] for j in X])
surf2=ax.plot_surface(X,Y,Zprat)
frame2=ax.plot_wireframe(X,Y,Zprat)
frame1._wrapped.set_color('r')

Zrag=array([[RaghavanLogK(1000.0/i) for i in j] for j in X])
surf3=ax.plot_surface(X,Y,Zrag)
frame3=ax.plot_wireframe(X,Y,Zrag)
frame2._wrapped.set_color('g')

save('meshX.out', X)
save('meshY.out', Y)
save('meshZme.out',Z)
save('meshZprat.out',Zprat)
save('meshZrag.out',Zrag)


for surf in [surf1, surf2, surf3]:
    surf._wrapped.set_alpha(0.2)
    surf._wrapped.set_linewidth(0)
    
for frame in [frame1, frame2, frame3]:
    frame._wrapped.set_alpha(0.5)
    
ax.set_xlabel('1000/T')
ax.set_ylabel('log_{10} P')
ax.set_zlabel('ln k')


## now add things
figure(1)
if clearfigs:
    plotRaghavanMeasurements()
    plot( PratsinisMeasuredInverseTs, PratsinisMeasuredLogKs, 'go', label="Pratsinis", linewidth=2 )

#### O2 
figure(3)
if clearfigs: clf()
P=1.0e5
title("Effect of excess O2 (no argon)")
o2xs=[-1,0,1,4,9]
for o2x in o2xs:
    temps=[]
    logKs=[]
    if o2x==-1:
        if not clearfigs: continue
        label="Raghavan (2001)"
        for T in [646, 2044]:
            temps.append(T)
            logKs.append(RaghavanLogK(T)) 
    elif o2x==0:
        if not clearfigs: continue
        label="Pratsinis (1990)"
        for T in [973, 1273]:
            temps.append(T)
            logKs.append(PratsinisLogK(T))
    else:
        for T in Ts:
            #label="O2/TiCl4=%.2g"%(o2x)
            ar=0.0
            ticl=(1.0-ar)/(1+o2x)
            o2=ticl*o2x
            X="TiCl4:%f, O2:%f, Ar:%f"%(ticl,o2,ar)
            label=X
            s=Simulation(T,P,X)
            print "T\t%dK\tO2/TiCl4\t%.2g\t"%(T,o2x),
            logK=s.getLogK()
            if logK>-Inf:
                temps.append(T)
                logKs.append(logK)
    temps=array(temps)
    logKs=array(logKs)
    plot( 1000.0/temps, logKs, format, label=label, linewidth=2 )
    


ylabel("ln k")
xlabel("1000/T")
if clearfigs:
    legend(loc=0)
    drawPratsinisBox()
    drawRaghavanBox()
    plotRaghavanMeasurements()
    plot( PratsinisMeasuredInverseTs, PratsinisMeasuredLogKs, 'go', label="Pratsinis", linewidth=2 )

figure(4)
if clearfigs: clf()
title("Effect of argon fraction (4:1 O2:TiCl4)")
P=1.0e5
argons=[-2,-1,0,0.9,0.99]
o2x=4
for argon in argons:
    temps=[]
    logKs=[]
    if argon==-2:
        if not clearfigs: continue
        label="Raghavan (2001)"
        for T in [646, 2044]:
            temps.append(T)
            logKs.append(RaghavanLogK(T)) 
    elif argon==-1:
        if not clearfigs: continue
        label="Pratsinis (1990)"
        for T in [973, 1273]:
            temps.append(T)
            logKs.append(PratsinisLogK(T))
    else:
        for T in Ts:
            #label="O2/TiCl4=%.2g"%(o2x)
            ar=argon
            ticl=(1.0-ar)/(1+o2x)
            o2=ticl*o2x
            X="TiCl4:%f, O2:%f, Ar:%f"%(ticl,o2,ar)
            label=X
            s=Simulation(T,P,X)
            print "T (K)\t%d\tArgon\t%.2g\t"%(T,argon),
            logK=s.getLogK()
            if logK>-Inf:
                temps.append(T)
                logKs.append(logK)
    temps=array(temps)
    logKs=array(logKs)
    plot( 1000.0/temps, logKs, format, label=label, linewidth=2 )



ylabel("ln k")
xlabel("1000/T")
if clearfigs:
    legend(loc=0)
    drawPratsinisBox()
    drawRaghavanBox()
    plotRaghavanMeasurements()
    plot( PratsinisMeasuredInverseTs, PratsinisMeasuredLogKs, 'go', label="Pratsinis", linewidth=2 )
