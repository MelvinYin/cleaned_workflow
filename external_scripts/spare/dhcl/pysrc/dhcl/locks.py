#######################################################################################################################################################################
#
# Grzegorz M Koczyk (2007-)
# Calculating van der Waals locks between regions of a protein (for _unrelated_ resource locking see: locking.py)
#
# NOTE: JMB_2001 = Berezovsky IN, Trifonov EN
#                  "Van der Waals locks: loop-n-lock structure of globular proteins."
#                  J Mol Biol. 2001; 307(5):1419-26. 
#
#######################################################################################################################################################################
from collections import defaultdict
import sys
from dhcl.utils import *
from dhcl.hierarchy import *
# Threshold of number of outside contacts to consider a block/cluster for locking
BLOCK_THRESHOLD = 200
# Threshold of number of contacts per reported lock
LOCK_THRESHOLD = 300
# Threshold of reported contacts per residue
RESIDUE_THRESHOLD = 100
# Minimum length of below threshold residue stretch to consider for separating a block
BLOCK_SPACING = 1
BLOCK_CONTACTS = 100
# Default sizes of locks to compute for
DEFAULT_LOCK_SIZES = range(3,6)
MIN_LOCK_SIZE = min(DEFAULT_LOCK_SIZES)
MAX_LOCK_SIZE = max(DEFAULT_LOCK_SIZES)

class Lock(object):
    ''' Object representing a van der Waals lock (two short stretches of contacting protein chain
    acting to stabilize the structure). '''
    def __init__(self, rset1, rset2, chain, no_contacts, seq2triple=None):
        ''' Constructs a lock object based on rset1,rset2 - sets of residue numbers, chain - chain identifiers, no_contacts - number
        of contacts made by lock residues, seq2triple - optional map of sequence numbers to (chain, residue number, insertion code) triplets. '''
        rset1, rset2 = ((rset1, rset2) if rset1[0]<=rset2[0] else (rset2,rset1))
        self.rset1 = tuple( sorted(rset1) )
        self.rset1_as_set = frozenset(self.rset1)
        self.rset2 = tuple( sorted(rset2) )
        self.rset2_as_set = frozenset(self.rset2)
        self.chain = chain
        self.no_contacts = no_contacts
        if seq2triple:
            self.bc = seq2triple[chain]
        else:
            self.bc = None
    def overlaps(self, other):
        ''' Whether a lock overlaps another lock'''
        if self.chain!=other.chain:
            return False
        else:
            return min( [len( self.rset1_as_set & other.rset1_as_set ), len(self.rset2_as_set & other.rset2_as_set)] )
    def __eq__(self, other):
        return self.rset1==other.rset1 and other.rset2==self.rset2
    def __hash__(self):
        return hash(self.rset1) ^ hash(self.rset2)
    def __cmp__(self,other):
        return cmp(self.rset1,other.rset1) or cmp(self.rset2,other.rset2)
    # TODO: Test this method
    @classmethod
    def parseString(klass, s, seqno2triple, triple2seqno):
        ''' Parse the string representation iteratively yielding locks encoded within it. '''
        s = s.replace(' ', '').replace('(', '').replace(')','')
        for lstr in s.split('/'):
            bstr, tstr = lstr.split('=')
            bstr_1, bstr_2 = bstr.split('>')
            try:
                rstr1, chain1 = bstr_1.split(':')
            except:
                rstr1 = bstr_1
                chain1 = ' '
            try:
                rstr2, chain2 = bstr_2.split(':')
            except:
                rstr2 = bstr_2
                chain2 = ' '
            chain = chain1
            rs1 = [ triple2seqno[chain][ parseResidueId(r) ] for r in rstr1.split(',') ]
            rs2 = [ triple2seqno[chain][ parseResidueId(r) ] for r in rstr2.split(',') ]
            no_contacts = int(tstr)
            chain = chain if chain else ' '
            yield klass(rs1, rs2, chain, no_contacts, seqno2triple)
    def __str__(self):
        ''' Encode into the string representation '''
        res1s = []
        res2s = []
        for res1 in self.rset1:
            # Hack to renumber back into the original structure scheme while writing the result out (if original numbering is available)
            if self.bc:
                resno, icode = self.bc[res1]
                res1s.append("%s%s" % (resno, icode.strip()) )
            # Otherwise use our numbering
            else:
                res1s.append("%s" % res1)
        for res2 in self.rset2:
            # Hack to renumber back into the original structure scheme while writing the result out (if original numbering is available)
            if self.bc:
                resno, icode = self.bc[res2]
                res2s.append("%s%s" % (resno, icode.strip()) )
            # Otherwise use our numbering
            else:
                res2s.append("%s" % res2)
        res1str = ','.join(res1s)
        res2str = ','.join(res2s)
        if self.chain!=' ':
            return "( %s:%s>%s:%s=%s )" % ( res1str, self.chain, res2str, self.chain, self.no_contacts)
        else:
            return "( %s>%s=%s )" % (res1str, res2str, self.no_contacts)
        
def continuous(rset):
    ''' Is a stretch of residue numbers continuous. '''
    last = None
    for r in rset:
        if last is not None and r!=last+1:
            return False
        last = r
    return True

def countContactsRes(rcs, resno1, chain1, resno2, chain2):
    ''' Count contacts between two residues, based on residue vdW data. '''
    try:
        return rcs[chain1][resno1][chain2][resno2].no_contacts
    except:
        return 0

def countContactsSets(first_set, second_set, chain, rcs):
    ''' Count contacts between two sets of residues, based on residue vdW data. '''
    r = 0
    for i in first_set:
        for j in second_set:
            if j<i-5 or j>i+5:
                try:
                    r += rcs[chain][i][chain][j].no_contacts
                except:                
                    pass    
    return r

# Block of residues, grouped for whatever reasons
class Block(list):
    def round(self, bcs=[]):
        return [ x for x in range(min(self), max(self)+1) if x in bcs ]

def blockContacts(block1, block2, chain, rcs):
    return countContactsSets(block1, block2, chain, rcs)

class Candidate(object):
    def __init__(self, rs1, rs2, conts):
        self.rset1 = rs1
        self.rset2 = rs2
        self.total = conts
        self.avg = conts/float( len(rs1) )
def findMaxLock(block1, block2, chain, rcs, seq2triple):
    ''' Find the best lock between two clusters of residues along the same chain '''
    conts1 = [ countContactsSets([x],block2, chain, rcs) for x in block1]
    conts2 = [ countContactsSets(block1,[y], chain, rcs) for y in block2]
    # On the basis of these - find best stretches of 3-5 residues
    cand_lock = None
    for lock_size in DEFAULT_LOCK_SIZES:
        for i1 in range(0, len(block1)-lock_size+1):
            rs1 = [ block1[x] for x in range(i1,i1+lock_size) ]
            for i2 in range(0, len(block2)-lock_size+1):
                rs2 = [  block2[x] for x in range(i2,i2+lock_size) ]
                total = countContactsSets( rs1, rs2, chain, rcs)
                # Criterion for choosing the best lock is avg number of residue contacts in the lock
                if cand_lock is None or (total/float(lock_size))>cand_lock.avg:                    
                    cand_lock = Candidate(rs1, rs2, total)
    # Rank candidates according to number of contacts
    if cand_lock:
        return Lock(cand_lock.rset1, cand_lock.rset2, chain, cand_lock.total, seq2triple)
    else:
        return None

def computeLocks(rcs, seq2triple):
    ''' Compute all locks in the protein. '''
    locksd = {}
    for chain in sorted(seq2triple):
        bcs = seq2triple[chain]
        # Compute total contacts for every residue
        totals = defaultdict(int)
        for i in sorted(bcs):
            # Compute total number of contacts between this residue and rest of the sequence
            # ignoring residues closer than 5 residues apart (turn of alpha-helix+1)
            for j in bcs:
                if j<i-5 or j>i+5:
                    totals[i]+= countContactsRes(rcs, i, chain, j, chain)
        # Create blocks of residues scoring above RESIDUE_THRESHOLD in contact numbers
        # Blocks are separate if there is more than a BLOCK_SPACING-residue gap between them
        idiff = 0
        last_i = None
        blocks = []
        for i in sorted(totals):
            # Increment in case of any discontinuities in the chain
            if (last_i is not None and last_i+1<i):
                idiff += i - last_i - 1
            if totals[i]<RESIDUE_THRESHOLD:
                idiff += 1
            else:
                # Skip to next block, if we are in an extended area below the threshold
                if len(blocks)==0 or idiff>=BLOCK_SPACING:
                    blocks.append( Block() )
                blocks[-1].append(i)
                idiff = 0
            last_i = i
        # Round out blocks with inlying residues scoring below threshold
        blocks = [ b.round(bcs) for b in blocks ]
        # Compute contacts among blocks (clusters)
        contacts = [ [j for j in range(0,len(blocks))] for i in range(0,len(blocks)) ]
        for i,block1 in enumerate(blocks):
            for j in range(i, len(blocks)):
                block2 = blocks[j]
                contacts[i][j] = contacts[j][i] = countContactsSets(block1, block2, chain, rcs)
        # For each possible pair of blocks scoring contacts above BLOCK_THRESHOLD
        result = []
        for i,block1 in enumerate(blocks):
            j = sorted( range(0,len(blocks)), key=lambda x: contacts[i][x])[-1]
            #for j in range(i+1, len(blocks)):
            if i != j:
                block2 = blocks[j]
                if contacts[i][j] > BLOCK_THRESHOLD:
                    #print "LOCK:", blocks[i], blocks[j]
                    lock = findMaxLock(block1, block2, chain, rcs, seq2triple)
                    if lock and lock.no_contacts>LOCK_THRESHOLD:
                        result.append(lock)
        locksd[chain] = sorted( frozenset(result) )
    return locksd

def prepareLocks(residue_cs_fname, bounds_fname):
    ''' Prepares domains for chains on basis of the files. WARNING: returns Python data structure, not filename '''
    seqno2triple, triple2seqno = readChainBounds(bounds_fname)
    return computeLocks( readResidueContacts(residue_cs_fname, triple2seqno), seqno2triple )

##########################################################################################################################################
# INPUT/OUTPUT METHODS
@handlify(mode='w')
def writeLocksTsv(fh, r, uid=None):
    for chain in sorted(r):
        print >>fh, "\t".join( [ uid, chain, " / ".join( str(lock) for lock in r[chain] ) ] )

@handlify(mode='r')
def readLocksTsv(fh, seq2triple, triple2seq):
    r = defaultdict(dict)
    for line in fh:
        uid, chain, locks_txt = line.split('\t')
        r[uid][chain] = locks = list( Lock.parseString(locks_txt,seq2triple, triple2seq) )
    return r
