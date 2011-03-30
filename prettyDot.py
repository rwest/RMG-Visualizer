# this version updated by RHW in November 2008
"""
make the dot file prettier
"""
import os, sys, re

StripLineLabels=False
if StripLineLabels:
    print "stripping edge (line) labels"

infile=file('rxnpath.dot')
outfile=file('rxnpath2.dot','w')

# replace this:
#  s10 [ fontname="Helvetica", label="C11H23J"];
# with this:
#  s10 [ shapefile="mols/C11H23J.png" label="" width="1" height="1" imagescale=true fixedsize=true color="white" ];
 
reSize=re.compile('size=\"5,6\"\;page=\"5,6\"')
reNode=re.compile('(?P<node>s\d+)\ \[\ fontname=\"Helvetica\",\ label=\"(?P<label>[^\"]*)\"\]\;')

for line in infile:
    (line,changedSize)=reSize.subn('size="12,12";page="12,12"',line)
    match=reNode.search(line)
    if match:
        if os.path.isfile("pics/%s.png"%match.group('label')):
            line='%s [ image="pics/%s.png" label="" width="0.5" height="0.5" imagescale=false fixedsize=false color="none" ];\n'%(match.group('node'),match.group('label'))

    # rankdir="LR" to make graph go left>right instead of top>bottom

    if StripLineLabels:
        line=re.sub('label\s*=\s*\"\s*[\d.]+\"','label=""',line)
        
    # change colours
    line=re.sub('color="0.7,\ (.*?),\ 0.9"',r'color="1.0, \1, 0.7*\1"',line)
        
    outfile.write(line)

outfile.close()
infile.close()
        
print "now try:\n dot -O -Tpdf rxnpath2.dot"