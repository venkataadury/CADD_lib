import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj

# Perform automatic vacuum minimization of small molecules using Chimera (Open directly in chimera)

def fileen(fnam):
    f=open(fnam,'r')
    for x in f:
        ens=str.split(str.split(x,' ')[1],'\t')
        en=float(ens[0])
        hen=float(ens[1])
        break
    return en,hen
def filekey(fnam):
    e1,e2=fileen(fnam)
    return e1
#os.chdir("/home/venkata/CADD/betacatenin/results/size19")
file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb") and fn.startswith("result") and "emin" not in fn] #Select all PDB files in this folder
file_names=sorted(file_names,key=filekey,reverse=False)
#rc("open protein.gro")
i=0
ind=0
while ind<len(file_names):
        fn=file_names[ind]
        print("Loading:",fn)
        fnnew=fn.replace(".pdb","_autoemin.pdb")
        replyobj.status("Processing " + fn)  #show what file we're working on
        rc("open " + fn)
	#rc("align DNV ~DNV") # put ligand in front of remainder of molecule
        if not i:
            rc("focus: DNV") # center/zoom ligand
            i=1
        rc("addcharge nonstd :DNV 0 method gas")
        rc("minimize spec :DNV nsteps 20000 cgsteps 15000 cgstepsize 0.01 interval 5000 nogui true prep false")
        print("Will save as "+fnnew)
        rc("write format pdb 0 "+fnnew)
        #rc("preset apply publication 1") # make everything look nice
        #rc("pause")
        rc("close: DNV")
        ind+=1
	# save image to a file that ends in .png rather than .pdb
	#png_name = fn[:-3] + "png"
	#rc("copy file " + png_name + " supersample 3")
	#rc("close all")
