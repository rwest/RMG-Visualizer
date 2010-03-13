# this version updated by RHW in January 2009
"""
load the chem.inp file from RMG and get it ready for mixmaster
"""
import os, sys, shutil

def convertChemkin2Cantera(RMG_results):
    """Convert the Chemkin file into a Cantera file.
    
    Does its work inside RMG_results/chemkin"""
    
    from Cantera import ck2cti
    starting_dir = os.path.getcwd()
    chemkin_dir = os.path.join(RMG_results,'chemkin')
    os.path.chdir(chemkin_dir)
    try:
        infile='chem.inp'
        thermodb=''
        trandb=''
        nm='chem'
        ck2cti.ck2cti(infile = infile, thermodb = thermodb,  trandb = trandb, idtag = nm, debug=0, validate=1)
    finally:
        os.path.chdir(starting_dir)

def convertFinalModel2MixMaster(RMG_results):
    """Convert the Final_Model.txt into appropriate CSV data file for mixmaster.
    
    Needs a MolarMasses.txt file, which is created in another function"""
    
    massesfilename=os.path.join(RMG_results,'MolarMasses.txt')
    print "Reading molar masses from",massesfilename
    massesfile=file(massesfilename)
    massesdict=dict()
    for line in massesfile:
        (species,mass)=line.split()
        massesdict[species]=mass
    massesfile.close()
    
    temperature=273+150
    pressure=208*101325
    print "Using these settings:\n Temperature: %f K \t Pressure: %f Pa\n"%(temperature,pressure)
    
    # load file
    filename ='Final_Model.txt'
    filepath = os.path.join(RMG_results,filename)
    resultFile=file(filepath)
    
    # search for "Mole Fraction Profile Output"
    line=resultFile.next()
    while (line.find('Mole Fraction Profile Output')<0):
        line=resultFile.next()
    # add "T \t P" to the  following line
    titles=resultFile.next()
    print "Species:",titles
    output=titles.strip()+"\tT\tP\tnothing\n"
    items=titles.split()
    assert items[0]=='Time'
    speciesnames=items[1:]
    masses=list()
    for species in speciesnames:
        masses.append(float(massesdict[species]))
    	
    # add the temperature IN KELVIN and pressure IN PASCAL to all the following nonblank lines
    line=resultFile.next()
    while (line.strip()):
        massfractions=[]
        massfractionsum = 0
        items = line.split()
        time = items[0]
        molefracs = items[1:]
        for i,molefrac in enumerate(molefracs):
            massfrac = float(molefrac)*masses[i]
            massfractions.append(massfrac)
            massfractionsum += massfrac
        massfractions = [str(m/massfractionsum) for m in massfractions]
        output += str(time)+'\t'
        output += '\t'.join(massfractions)
        output +=  "\t%f\t%f\t0\n"%(temperature,pressure)
        line=resultFile.next()
    # turn whitespaces into commas
    # save the output
    outputFile=file(os.path.join(RMG_results,'ForMixMaster.csv'),'w')
    outputFile.write(output.replace('\t',','))
    outputFile.close()
    print "ForMixMaster.csv now contains mass fractions, as required by MixMaster"
    

def makeTableOfSpecies(RMG_results):
    """Make a pretty table of species"""
    ### make pretty table of species
    import ctml_writer
    from ctml_writer import *
    # these lists store top-level entries. Empty them!
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
    
    import jinja2
    env = jinja2.Environment(loader = jinja2.FileSystemLoader('templates'))
    template = env.get_template('rxnlist.html')
    outstring=template.render(title=infile, reactionList=ctml_writer._reactions)
    outfile=file('ReactionList'+'.html','w')
    outfile.write(outstring)
    outfile.close()

def loadMixMaster(RMG_results):
    """Load MixMaster"""
    os.path.chdir(RMG_results)
    from MixMaster import MixMaster
    o=MixMaster()
    o.loadmech('','chem.cti')
    
if __name__ == "__main__":
    RMG_results = "RMG_result"
    
    convertChemkin2Cantera(RMG_results)
    convertFinalModel2MixMaster(RMG_results)
    makeTableOfSpecies(RMG_results)