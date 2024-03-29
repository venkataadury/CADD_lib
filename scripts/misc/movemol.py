#!/usr/bin/python3
#REQUIRES molgeom.py

import os
import sys
import argparse
import molgeom


parser=argparse.ArgumentParser(prog="movemol.py")
parser.add_argument(dest="infile",action="store",help="The input pdb file")
parser.add_argument('-x',dest='xc',metavar='--move_x',action="store",help="Translation along x direction (Angstroms)",default=0.0)
parser.add_argument('-y',dest='yc',metavar='--move_y',action="store",help="Translation along y direction (Angstroms)",default=0.0)
parser.add_argument('-z',dest='zc',metavar='--move_z',action="store",help="Translation along z direction (Angstroms)",default=0.0)
parser.add_argument('-o',dest='outfile',metavar='--output',action="store",help="Output file name",default=None)

args=parser.parse_args().__dict__
molgeom.parmove(args)
