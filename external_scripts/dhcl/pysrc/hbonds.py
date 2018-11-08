from __future__ import with_statement
import math, sys, os, os.path
import Bio.PDB
#import Numeric
from dhcl.utils import *
from dhcl.pdb_data import prepareHbonds
from dhcl.hbonds import *
from collections import defaultdict
# Output DEBUG msgs
DEBUG = 0
import csv

if __name__=='__main__':
    pdb_file = sys.argv[1]
    my_bonds = hbonds(pdb_file, True)
    #hb_file = sys.argv[2]
    hb_file = sys.stdout
    writeHbonds(hb_file, hbonds(pdb_file))    
