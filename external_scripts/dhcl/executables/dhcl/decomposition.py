#####################################################################################################################################################################
#
# Grzegorz M Koczyk (2007-)
# 
# Decomposition of protein structures into loops
#
#####################################################################################################################################################################
from __future__ import with_statement
import sys, os, os.path
from dhcl.utils import *
from dhcl.pdb_data import *
from collections import defaultdict
DEBUG = 0
######################################################################################################################################################################
# FUNCTIONS - parsing data
@handlify(is_method=False, mode='r')
def readAtomDistsForResidues(fh, triple2seqno):
    ''' Utility function to read in residue contact distances data outputted by contacts program '''
    result = defaultdict( innerDefaultdict( innerDefaultdict( dict ) ) )
    for line in fh:
        if line:
            atomno1, resno1, icode1, resname1, chain1, atomname1, atomno2, resno2, icode2, resname2, chain2, atomname2, r = line.strip('\n').split('\t')
            resno1, resno2 = int(resno1), int(resno2)
            seqno1, seqno2 = triple2seqno[chain1][(resno1, icode1)], triple2seqno[chain2][(resno2, icode2)]
            r = float(r)
            r1 = result[chain1][seqno1][chain2]
            r2 = result[chain2][seqno2][chain1]
            if resno2 in r1:
                r1[resno2].no_contacts +=1
                r_old = r1[resno2].r
                r1[resno2].r = min( [r, r_old] )
            else:
                result[chain2][seqno2][chain1][seqno1] = result[chain1][seqno1][chain2][seqno2] = StructObject( resname1=resname1, resname2=resname2, r=r, no_contacts=1 )
    return result

MIN_LOOP_LEN = 15
MAX_LOOP_LEN = 45
MIN_INIT_DISTANCE = 2.5
MAX_INIT_DISTANCE = 5.0
DISTANCE_STEP = 1.0
MAX_FINAL_DISTANCE = 12.0
# Max gap length for stopping criterion
MAX_GAP_LEN = 10

class Loop(object):
    ''' Object representing a loop '''
    def __init__(self,resno1, resno2, distance, chain=' ', bounds=None):
        self.resno1, self.resno2 = (resno1,resno2) if resno2>resno1 else (resno2,resno1)
        if bounds is not None:
            self.bc = bounds[chain]
        else:
            self.bc = None
        self.distance = distance
        self.chain = chain
        self.use = True
    # Alternative accessors for first and last residue numbers
    def getBegin(self):
        return self.resno1
    begin = property(getBegin)
    def getEnd(self):
        return self.resno2
    end = property(getEnd)
    def encompasses(self, other):
        ''' Whether a loop encompasses another loop '''
        if self.chain!=other.chain:
            return False
        else:
            return (self.resno1<=other.resno1 and self.resno2>=other.resno2)
    def crosses(self, other):
        ''' Whether a loop overlaps with another loop '''
        if self.chain!=other.chain:
            return 0
        else:
            d = min([self.resno2, other.resno2]) - max([self.resno1, other.resno1])+1
            if d>0:
                return d
            else:
                return 0
    def __str__(self):
        if self.bc:
            beginr, begini = self.bc[self.begin]
            begin = "%s%s" % (beginr, begini.strip())
            endr, endi = self.bc[self.end]
            end = "%s%s" % (endr, endi.strip())
        # Otherwise use our numbering
        else:
            begin = self.begin
            end = self.end
        return "(%s:%s>%s:%s=%s)" % (begin, self.chain, end, self.chain, self.distance)
    def __len__(self):
        return self.resno2-self.resno1+1
    def getSize(self):
        return self.resno2-self.resno1+1
    size = property(getSize)
    
    @classmethod
    def parseString(klass, s, seqno2triple=None, triple2seqno=None):
        ''' Parse a string encoding this loop (see __str__ method for the input pattern)'''
        s = s.replace(' ', '')
        if len(s):
            loop_strs = s.split('/')
            for lstr in loop_strs:
                lstr = lstr[1:-1]
                bstr, tstr = lstr.split('=')
                distance = float(tstr)
                b1str, b2str = bstr.split('>')
                begin, chain = b1str.split(':')
                end, chain = b2str.split(':')
                chain = chain if chain else ' ' 
                if triple2seqno:
                    begin, end = triple2seqno[chain][parseResidueId(begin)], triple2seqno[chain][parseResidueId(end)]
                else:
                    pass
                yield klass(begin, end, distance, chain, seqno2triple)
    __repr__ = __str__
    
def decomposeStep(rdists, seqno2triple, r_min, r_max, no_len_bounds=False):
    ''' Decompose into loops at given (r_min, r_max) radius thresholds. '''
    for chain in seqno2triple:
        bchain = sorted(seqno2triple[chain])
        for i, resno1 in enumerate(bchain):
            r_ = rdists[chain][resno1][chain]
            for resno2 in bchain[i+1:]:
                if resno2 in r_ and (no_len_bounds or MIN_LOOP_LEN <= (resno2-resno1+1) and (resno2-resno1+1) <= MAX_LOOP_LEN):
                    r_dist = r_[resno2].r
                    if resno2 in r_ and (r_min <= r_dist) and ( r_dist <= r_max):
                        yield Loop(resno1, resno2, r_dist, chain, seqno2triple)

def curateSingle(loops):
    ''' Curate a single set of loops, so that only the non-overlapping tightest ones remain. '''
    loops = sorted(loops, key=lambda l: l.distance)
    keep = [True for l in loops]
    for i,loop in enumerate(loops):
        if not keep[i]:
            continue
        for j,loop_ in enumerate(loops[i+1:]):
            if loop.crosses(loop_)>5:
                keep[i+j+1] = False
    return [ l for i,l in enumerate(loops) if keep[i] ]

def curateByOld(loops, old_loops):
    ''' Curate a single set of loops, by other set of loops (removing overlapping segments). '''
    for loop in old_loops:
        yield loop
    for loop in loops:
        keep = True
        for loop_ in old_loops:
            if loop.crosses(loop_)>5:
                keep=False
                break
        if keep:
            yield loop

class Block(list):
    pass

def loopDecompose(rdists, seqno2triple):
    # Start with MIN_INIT_DISTANCE, MAX_FINAL_DISTANCE decomposition
    # Find the loops and reconstruct initial coverage
    whole_loops = list( decomposeStep(rdists, seqno2triple, MIN_INIT_DISTANCE, MAX_FINAL_DISTANCE) )
    if DEBUG:
        print "Whole", whole_loops
    loops = list( curateSingle(whole_loops) )
    if DEBUG:
        print loops
    # Partition loops by chain
    lbcs = defaultdict(list)
    for chain in seqno2triple:
        lbcs[chain] = list()
    for l in loops:
        lbcs[l.chain].append(l)
    # Also create smaller loops and partition them by chain
    smcs = defaultdict(list)
    for l in decomposeStep(rdists, seqno2triple, 4., 5.):
        if len(l)<30:
            smcs[l.chain].append(l)
    # For each chain
    for chain in lbcs:
        cloops = lbcs[chain] = sorted(lbcs[chain], key=lambda l: l.resno1)
        smaller_loops = sorted(smcs[chain], key=lambda l: l.distance)
        #####################################################################
        # Find blocks of large loops (>40AA)
        blocks = []
        block = Block()
        block.begin = block.end = min(seqno2triple[chain])
        new_cloops = []
        for l in cloops:
            if len(l)>=39:
                if block.begin > l.resno1:
                    block.begin = l.resno1
                if block.end < l.resno2:
                    block.end = l.resno2
                block.append(l)
            else:
                if block:
                    if block.end < l.resno1:
                        block.end = l.resno1
                    blocks.append(block)                 
                block = Block()
                block.begin = l.resno2
                block.end = l.resno2
                # This loop is saved for future (it is not part of a block to decompose)
                new_cloops.append(l)
        if block:
            block.end = max(seqno2triple[chain])
            blocks.append(block)
        if DEBUG:
            print blocks
            print smaller_loops
        ########################################################################
        # For each block of large loops - try to decompose it into sets of smaller ones
        for block in blocks:
            nsloops = list( curateSingle( l for l in smaller_loops if l.resno1>=block.begin and l.resno2<=block.end ) )
            if DEBUG:
                print "Decompose"
                print block, block.begin, block.end
                print nsloops
            old_but_used = []
            if len(nsloops)>0:
                for l in block:
                    append = True
                    for nl in nsloops:
                        if nl.crosses(l)>5:
                            append=False
                            break
                # Keep old loops if they are the only solution (no significant crosses between new ones and old one) :)
                if append:
                    old_but_used.append(l)
                for l in nsloops:
                    l.use = False
                new_cloops.extend(nsloops+old_but_used)
            else:
                new_cloops.extend(block)
        if DEBUG:
            print new_cloops
        ########################################################################
        # Curate ends
        #lbcs[chain] = new_cloops
        cloops = sorted(new_cloops, key=lambda l: l.resno1)
        if len(cloops)>1:
            # Curate N end - try to extend the original loop keeping to within 115% of tightness
            loop = cloops[0]
            if loop.use and len(loop)<20:
                begin = min(seqno2triple[chain])
                if len(cloops)>1:
                    if cloops[1].resno1 > loop.resno2:
                        end = cloops[1].resno1
                    else:
                        end = loop.resno2
                else:
                    end = max(seqno2triple[chain])
                cands = sorted( [ l for l in whole_loops if l.resno1>=begin and l.resno2<=end and l.distance<=loop.distance*1.15 and l.chain==chain],
                            key=lambda l: len(l), reverse=True )
                cloops[0] = cands[0]
            # Curate C end - try to extend the original loop keeping to within 115% of tightness
            loop = cloops[-1]
            if loop.use and len(loop)<20:
                end = max(seqno2triple[chain])
                if len(cloops)>1:
                    if cloops[-2].resno2 < loop.resno1:
                        begin = cloops[-2].resno2
                    else:
                        begin = loop.resno1
                else:
                    begin = min(seqno2triple[chain])
                cands = sorted( [ l for l in whole_loops if l.resno1>=begin and l.resno2<=end and l.distance<=loop.distance*1.15 and l.chain==chain],
                            key=lambda l: len(l), reverse=True )
                cloops[-1] = cands[0]
            lbcs[chain] = cloops
    ############################################################################
    # Results are returned by chain
    return lbcs

def prepareLoops(calpha_cs_fname, bounds_fname):
    ''' Prepares domains for chains on basis of the files. WARNING: returns Python data structure, not filename '''
    seqno2triple, triple2seqno = readChainBounds(bounds_fname)
    return loopDecompose( readAtomDistsForResidues(calpha_cs_fname, triple2seqno), seqno2triple )

# INPUT / OUTPUT OF LOOPS
import csv
@handlify(mode='w')
def writeLoopsTsv(fh, loops, uid=None):
    ''' Write loops into TSV format. Optionally give the unique PDB identifier (first column in the TSV file).'''
    writer = csv.writer(fh, dialect='excel-tab')
    for chain in sorted(loops):
        writer.writerow([uid, chain] +
                        ['/'.join( [ str(l) for l in loops[chain] ] ) ]
                        )
        
@handlify(mode='r')
def readLoopsTsv(fh, seqno2triple=None, triple2seqno=None):
    ''' Read loops in tab-separated format into a dictionary indexed by PDB of loop dictionaries indexed by chains.'''
    reader = csv.reader(fh, dialect='excel-tab')
    r = defaultdict( dict )
    for row in reader:
        uid, chain = row[0], row[1]
        r[uid][chain] = list(Loop.parseString(row[2], seqno2triple, triple2seqno))
    return r
