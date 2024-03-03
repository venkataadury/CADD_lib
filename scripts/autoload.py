import os
from chimera import runCommand as rc # use 'rc' as shorthand for runCommand
from chimera import replyobj

# OLD TEMPLATE - USE: autoplay.py #
os.chdir("/home/venkata/CADD/betacatenin/results/size19")
file_names = [fn for fn in os.listdir(".") if fn.endswith(".pdb")]
rc("open /home/venkata/CADD/betacatenin/betacat_reduced.pdb")
i=0
for fn in file_names:
        print("Loading:",fn)
	replyobj.status("Processing " + fn) # show what file we're working on
	rc("open " + fn)
	#rc("align DNV ~DNV") # put ligand in front of remainder of molecule
        if not i:
            rc("focus: DNV") # center/zoom ligand
            i=1
	rc("surf") # surface receptor
	#rc("preset apply publication 1") # make everything look nice
	rc("surftransp 15") # make the surface a little bit see-through
        rc("pause")
        rc("close: DNV")
	# save image to a file that ends in .png rather than .pdb
	#png_name = fn[:-3] + "png"
	#rc("copy file " + png_name + " supersample 3")
	#rc("close all")
