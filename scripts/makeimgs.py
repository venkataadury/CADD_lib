#!/usr/bin/env python3
import os,sys
from rdkit import Chem
from rdkit.Chem.Draw import rdMolDraw2D
import re

pattern = re.compile("<\?xml.*\?>")
def drawMol(mc, molSize=(600,250), kekulize=True):
    print("Type:",type(mc))
    print(mc)
    if kekulize:
        try:
            Chem.Kekulize(mc)
        except:
            mc = Chem.Mol(mol.ToBinary())
    if not mc.GetNumConformers():
        Chem.rdDepictor.Compute2DCoords(mc)

    drawer = rdMolDraw2D.MolDraw2DSVG(*molSize)
    drawer.DrawMolecules([mc])
    drawer.FinishDrawing()
    drawer.drawOptions().addStereoAnnotation=True
    svg = drawer.GetDrawingText().replace('svg:', '')
    svg = re.sub(pattern, '', svg)
    return svg
f=sys.argv[1]
fl = open(f,"r")
for l in fl:
    l=l[0:len(l)-1] #Remove the '\n' in the end
    sp=str.split(l.replace("\t"," "),' ')
    ss=sp[0]; name=sp[1]
    print(ss)
    #mc = Chem.MolFromPDBFile(ss,removeHs=False,proximityBonding=False,sanitize=True) #Change to MolFromSmiles when needed
    #mc = Chem.MolFromMol2File(ss,removeHs=False) #Change to MolFromSmiles when needed
    #mc = Chem.MolFromSmiles(Chem.MolToSmiles(mc)) #Change to MolFromSmiles when needed
    mc = Chem.MolFromSmiles(ss)
    svg=drawMol(mc,(600,600),False)
    spl=str.split(name,'.')
    imgf='.'.join(spl[0:max(len(spl),2)-1])
    fn=imgf+".svg"
    print(fn)
    outf=open(fn,"w")
    outf.write(svg+"\n")
    outf.close()
