# this version updated by RHW in November 2008
"""
load the chem.inp file from RMG and get it ready for mixmaster
"""
import os, sys, shutil, re
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

import openbabel, pybel
# please cite:
# Pybel: a Python wrapper for the OpenBabel cheminformatics toolkit
# Noel M O'Boyle, Chris Morley and Geoffrey R Hutchison
# Chemistry Central Journal 2008, 2:5
# doi:10.1186/1752-153X-2-5


# Get folder with RMG results in
from RMG_results_path import RMGworkingDir


picfolder='pics'
molfolder='mols'
for path in [picfolder,molfolder]:
    os.path.isdir(path) or os.mkdir(path)


periodicTableByNumber={ 1: 'H',  2: 'He',  3: 'Li',  4: 'Be',  5: 'B',  6: 'C',  7: 'N',  8: 'O',  9: 'F',  10: 'Ne',  11: 'Na',  12: 'Mg',  13: 'Al',  14: 'Si',  15: 'P',  16: 'S',  17: 'Cl',  18: 'Ar',  19: 'K',  20: 'Ca',  21: 'Sc',  22: 'Ti',  23: 'V',  24: 'Cr',  25: 'Mn',  26: 'Fe',  27: 'Co',  28: 'Ni',  29: 'Cu',  30: 'Zn',  31: 'Ga',  32: 'Ge',  33: 'As',  34: 'Se',  35: 'Br',  36: 'Kr',  37: 'Rb',  38: 'Sr',  39: 'Y',  40: 'Zr',  41: 'Nb',  42: 'Mo',  43: 'Tc',  44: 'Ru',  45: 'Rh',  46: 'Pd',  47: 'Ag',  48: 'Cd',  49: 'In',  50: 'Sn',  51: 'Sb',  52: 'Te',  53: 'I',  54: 'Xe',  55: 'Cs',  56: 'Ba',  57: 'La',  58: 'Ce',  59: 'Pr',  60: 'Nd',  61: 'Pm',  62: 'Sm',  63: 'Eu',  64: 'Gd',  65: 'Tb',  66: 'Dy',  67: 'Ho',  68: 'Er',  69: 'Tm',  70: 'Yb',  71: 'Lu',  72: 'Hf',  73: 'Ta',  74: 'W',  75: 'Re',  76: 'Os',  77: 'Ir',  78: 'Pt',  79: 'Au',  80: 'Hg',  81: 'Tl',  82: 'Pb',  83: 'Bi',  84: 'Po',  85: 'At',  86: 'Rn',  87: 'Fr',  88: 'Ra',  89: 'Ac',  90: 'Th',  91: 'Pa',  92: 'U',  93: 'Np',  94: 'Pu',  95: 'Am',  96: 'Cm',  97: 'Bk',  98: 'Cf',  99: 'Es',  100: 'Fm',  101: 'Md',  102: 'No',  103: 'Lr',  104: 'Rf',  105: 'Db',  106: 'Sg',  107: 'Bh',  108: 'Hs',  109: 'Mt',  110: 'Ds',  111: 'Rg',  112: 'Uub',  113: 'Uut',  114: 'Uuq',  115: 'Uup',  116: 'Uuh',  117: 'Uus',  118: 'Uuo'}
periodicTableBySymbol=dict([(val, key) for key, val in periodicTableByNumber.items()])   

OBMolBondTypes={'S':1, 'D':2, 'T':3, 'B':5 }



# copy the RMG dictionary file
infile='RMG_Dictionary.txt'
oldpath=os.path.join(RMGworkingDir,infile)
newpath=os.path.join(os.getcwd(),infile)
print "copying %s to %s"%(oldpath,newpath)
shutil.copy2(oldpath, newpath) # copy it to the current folder

# load file
RMGfile=file(newpath)


for i in range(1,30000):
    print 'Molecule', i,'\t',
    name=''
    try:
        while name=='':
            name=RMGfile.next().strip()
    except StopIteration:
        print 'No more molecules'
        break
    print name
    graph=[]
    line=RMGfile.next()
    while line.strip():
        graph.append(line)
        line=RMGfile.next()
    # now have 'name' and 'graph'
    
    mol = openbabel.OBMol()
    #print 'Should print 0 (atoms)'
    #print mol.NumAtoms()

    re_bond=re.compile('\{(?P<atomnum>\d+),(?P<bondtype>[SDTB])\}')
    for line in graph:
        #print 'line:',line.strip()
        if len(line.split())>3:
            (number, element, radical, bonds)=line.split(None,3)
        else:
            (number, element, radical )=line.split(None)
        a = mol.NewAtom()
        a.SetAtomicNum(periodicTableBySymbol[element])  # 6 for a carbon atom
        #a.SetVector(0.0, 1.0, 2.0) # coordinates
        if int(radical[0]): # the [0] is so we take the first character of the string, in case it's something like "2T"
            a.SetSpinMultiplicity(int(radical[0])+1)
            # note that for non-radicals it's 0, but single radicals are 2, double radicals are 3...
            # http://openbabel.org/wiki/Radicals_and_SMILES_extensions#How_OpenBabel_does_it
        for bond in bonds.split():
            #print bond
            matchobject=re_bond.match(bond)
            if matchobject:
                fromAtom=int(number)
                toAtom=int(matchobject.group('atomnum'))
                bondType=matchobject.group('bondtype')
                if toAtom>fromAtom:
                    continue # because toAtom hasn't been placed yet!
                # print "%s bond from %d to %d"%(bondType,fromAtom,toAtom)
                mol.AddBond(fromAtom,toAtom,OBMolBondTypes[bondType])
            else:
                raise "couldn't figure out this bond: %s"%bond
    pymol=pybel.Molecule(mol)
    print pymol.write().strip(), 
    chemkinformula=pymol.formula+'J'*(pymol.spin-1)
    print chemkinformula
    if pymol.OBMol.NumHvyAtoms()>1:
        pymol.removeh()
    pymol.draw(filename=os.path.join('pics',name+'.png'), update=True, show=False)
    pymol.write(format='mol',filename=os.path.join('mols',name+'.mol'),overwrite=True)
        
