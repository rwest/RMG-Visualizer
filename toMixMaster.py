# this version updated by RHW in January 2009
"""
load the chem.inp file from RMG and get it ready for mixmaster
"""
import os, sys, shutil
#from Cantera import *
#from Cantera.Reactor import *o
#from Cantera.Func import *
#from math import *
#import random
#import pylab
#import scipy
#from pylab import *	
#from scipy import *
#from scipy.optimize import leastsq

# Get folder with RMG results in
from RMG_results_path import RMGworkingDir

# convert the chemkin file from RMG into a cantera file chem.cti
from Cantera import ck2cti
infile='chem.inp'
#oldpath=os.path.join(RMGworkingDir,'chemkin',infile)
oldpath=os.path.join(RMGworkingDir,infile)
newpath=os.path.join(os.getcwd(),infile)
print "copying %s to %s"%(oldpath,newpath)
shutil.copy2(oldpath, newpath) # copy it to the current folder

thermodb=''
trandb=''
nm='chem'
ck2cti.ck2cti(infile = infile, thermodb = thermodb,  trandb = trandb, idtag = nm, debug=0, validate=1)


# convert the Final_Model.txt into approprita CSV file

print "NB  ForMixMaster.csv is wrong. MixMaster wants MASS fractions and we are giving it MOLE fractions" 

temperature=273+150
pressure=208*101325
print " using these settings:\n Temperature: %f K \t Pressure: %f Pa\n"%(temperature,pressure)

# load file
resultFile=file(os.path.join(RMGworkingDir,'Final_Model.txt'))
# search for "Mole Fraction Profile Output"
line=resultFile.next()
while (line.find('Mole Fraction Profile Output')<0):
    line=resultFile.next()
# add "T \t P" to the  following line
titles=resultFile.next()
print "Species:",titles
output=titles.strip()+"\tT\tP\tnothing\n"
# add the temperature IN KELVIN and pressure IN PASCAL to all the following nonblank lines
line=resultFile.next()
while (line.strip()):
    output += line.strip() + "\t%f\t%f\t0\n"%(temperature,pressure)
    line=resultFile.next()
# turn whitespaces into commas
# save the output
outputFile=file('ForMixMaster.csv','w')
outputFile.write(output.replace('\t',','))
outputFile.close()


### make pretty table of species
import ctml_writer
from ctml_writer import *
# these lists store top-level entries
ctml_writer._elements = []
ctml_writer._species = []
ctml_writer._speciesnames = []
ctml_writer._phases = []
ctml_writer._reactions = []
ctml_writer._atw = {}
ctml_writer._enames = {}
ctml_writer._valsp = ''
ctml_writer._valrxn = ''
ctml_writer._valexport = ''
ctml_writer._valfmt = ''

execfile('chem.cti')

for rxn in ctml_writer._reactions:
	for side in [rxn._r,rxn._p]:
		for sp in side:
			num=side[sp]
			if num!=1: print num,
			print sp
		
import jinja2
env = jinja2.Environment(loader = jinja2.FileSystemLoader('templates'))
template = env.get_template('rxnlist.html')
outstring=template.render(title=infile, reactionList=ctml_writer._reactions)
outfile=file('ReactionList'+'.html','w')
outfile.write(outstring)
outfile.close()

print outstring


# load mixmaster
from MixMaster import MixMaster
#o=MixMaster()
#o.loadmech('','chem.cti')