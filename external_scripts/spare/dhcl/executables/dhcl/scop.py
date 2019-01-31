# Handling SCOP format data
from dhcl.utils import *
from collections import defaultdict
import sys

class ScopId(object):
    def __init__(self, klass, fold, sfam, fam):
        self.klass= klass
        self.sfam = sfam
        self.fam = fam
        self.fold = fold
    def upToSfam(self):
        return "%s.%s.%s" % (self.klass, self.fold, self.sfam)
    def __str__(self):
        return "%s.%s.%s.%s" % (self.klass, self.fold, self.sfam, self.fam)
# Read a CLA file and get dictionary of PDB identifiers with each chain mapped to SCOP domains
@handlify(mode='r')
def parseClaDict(fh):
    result = defaultdict(innerDefaultdict(dict))
    for line in fh:
        if line and not line.startswith('#'):
            fields = line.strip().split('\t')
            pdb_id, positions, scop = fields[1].strip(), fields[2].strip(), fields[3].strip()
            for position in positions.split(','):
                try:
                    chain, position = position.split(':')
                except:
                    #print pdb_id, chain, position
                    continue
                if not position:
                    position = '*'
                if not scop[0] in 'ijkl':
                    if scop.startswith('unassigned'):
                        continue
                    result[pdb_id][chain][position] = ScopId(*scop.split('.'))
                           
    return result

def findSfams(x):
    r = defaultdict(list)
    for k in x:
        scops = list()
        x_ = x[k]
        for k1 in x_:
            x__= x_[k1]
            for k2 in x__:
                scops.append( x__[k2].upToSfam() )
        scops.sort()
        r[ frozenset(scops) ].append(k)
    return r

if __name__ == '__main__':    
    x = parseClaDict(sys.stdin)
    print len(x)
    # Clustering (finding each superfamily combinations)
    r = defaultdict(list)
    for k in x:
        scops = list()
        x_ = x[k]
        for k1 in x_:
            x__= x_[k1]
            for k2 in x__:
                scops.append( x__[k2].upToSfam() )
        scops.sort()
        r[ frozenset(scops) ].append(k)
    y = [k for k in sorted(r) if len(r[k])>50]
    print len(y)
    for k in y:
        print k, len(r[k])
    # Read up on 'em in the PDB stats file
    

                
