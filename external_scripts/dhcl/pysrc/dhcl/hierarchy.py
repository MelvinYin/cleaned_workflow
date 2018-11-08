#######################################################################################################################################################################
#
# Grzegorz M Koczyk (2007-)
# Domain hierarchy module for:
#           - establishing split points based on van der Waals' interaction energies at a given potential barrier level.
#           - converting split points maxima/minima to a set of segments
#           - joining the segments to form domains
#
# NOTE: JBSD_1999 = "Hierarchy of the Interaction Energy Distribution in the Spatial Structure of Globular Proteins and the Problem of Domain Definition"
#                    Berezovsky IN, Namiot VA, Tumanyan VG, Esipova NG (1999) 
#                    Journal of Biomolecular Structure&Dynamics 17(1):133-155
#
#######################################################################################################################################################################
from __future__ import with_statement
import sys, os, os.path, time
from dhcl.utils import *
from dhcl.pdb_data import *
from collections import defaultdict
from operator import itemgetter, attrgetter
from copy import deepcopy, copy

######################################################################################################################################################################
# Module constants
# Turn on debugging comments
DEBUG = 0
# Benchmark timing 
TIMEBENCH=0
# Constants for prediction
# Critical percentage for merging domains based on majority of interactions
CRITICAL_MERGE_PERCENT = 0.7
# External/internal energy ratio for merging segments
CRITICAL_MERGE_ERATIO = 2.5
# Low energy segment external/internal energy ratio
LOWE_SEGMENT_ERATIO = 1.5
# Default deltaE values to compute segments and domain
DEFAULT_DELTAS = [ 0.3, 0.25, 0.2, 0.15, 0.1, 0.05,  ]

######################################################################################################################################################################
# Classes representing segments, domains etc.
class Segment(object):
    ''' A segment is a triplet containing begin, end and chain information for
    a continuous strip of residues. '''
    ######################################################################################################
    # Magic methods
    def __init__(self, begin, end, chain, bounds=None):
        # Bounds is specified optionally (allows to skip breaks in the protein chain while iterating over a segment)
        self.begin = begin
        self.end = end
        self.bc_parent = bounds
        if bounds is not None:
            self.bc = bounds[chain]
        else:
            self.bc = None
        self.chain = chain
    def __iter__(self):
        bc = self.bc
        for i in range(self.begin, self.end+1):
            if bc is None or i in bc:
                yield (i, self.chain)
    def __contains__(self, v):
        seqno, chain = v
        if self.chain!=chain:
            return False
        if (self.begin<=seqno) and (self.end>=seqno):
            return True
        return False
    def __cmp__(self, other):
        return cmp(self.begin,other.begin) and cmp(self.end, other.end) and cmp(self.chain,other.chain)
    def __eq__(self, other):
        return self.begin==other.begin and self.chain==other.chain and self.end==other.end
    def __hash__(self):
        return hash(self.begin) ^ hash(self.end) ^ hash(self.chain)
    def __str__(self):
        # Hack to renumber back into the original structure scheme while writing the result out (if original numbering is available)
        if self.bc:
            beginr, begini = self.bc[self.begin]
            begin = "%s%s" % (beginr, begini.strip())
            endr, endi = self.bc[self.end]
            end = "%s%s" % (endr, endi.strip())
        # Otherwise use our numbering
        else:
            begin = self.begin
            end = self.end
        if self.chain!=' ':
            return "%s:%s>%s:%s" % (begin, self.chain, end, self.chain)
        else:
            return "%s>%s" % (begin, end)
    __repr__ = __str__
    ######################################################################################################
    # Property attributes    
    def getSize(self):
        return len( list(x for x in self) )
        return (self.end - self.begin + 1)
    size = property(getSize, doc='Length of this segment (will take into account chain discontinuities).')
 
class Domain(set):
    ''' A domain is a set of segments, which can be iterated to yield all their respective segments. '''
    ######################################################################################################
    # Magic methods
    def __init__(self, *segments):
        self.segments = set(segments)
    def __iter__(self):
        ''' Iterate over residues in the domain '''
        for segment in sorted(self.segments, key=attrgetter('begin')):
            for x in segment:
                yield x
    def __contains__(self, v):
        for segment in self.segments:
            if v in segment:
                return True
        return False
    def __str__(self):
        return  "( %s )" % (', '.join(str(x) for x in self.sorted_segments),)
    __repr__=__str__
    def __hash__(self):
        return hash( frozenset(self.segments) )
    def __cmp__(self, other):
        return cmp( frozenset(self.segments), frozenset(other.segments) )
    # TODO: rewrite this (renumbering!)
    @classmethod
    def parseString(klass, s,seqno2triple=None, triple2seqno=None):
        ''' Parse a string encoding this domain (see __str__ method for the input pattern)'''
        dom_strs = s.replace(' ', '').replace('(','').replace(')','').split('/')
        for dom_str in dom_strs:
            #dom_str = dom_str[1:-1]
            segments = []
            for seg_str in dom_str.split(','):
                bstr, estr = seg_str.split('>')
                try:
                    begin, chain = bstr.split(':')
                except:
                    begin, chain = bstr, ' '
                try:
                    end, chain = estr.split(':')
                except:
                    end, chain = estr, ' '
                if triple2seqno:
                    begin, end = triple2seqno[chain][parseResidueId(begin)], triple2seqno[chain][parseResidueId(end)]
                else:
                    begin, end = int(begin), int(end)
                segments.append(Segment(begin, end, chain, seqno2triple))
            yield klass(*segments)
    ######################################################################################################
    # Property attributes
    def getSorted(self):
        ''' Return segments sorted by "begin"'''
        for segment in sorted(self.segments, key=attrgetter('begin')):
            yield segment
    sorted_segments = property(getSorted, doc='Iterator over segments sorted by beginning index')
    def getMerged(self):
        segs = []
        for s in self.sorted_segments:
            if len(segs) and s.begin==segs[-1].end+1:
                segs[-1] = Segment(segs[-1].begin, s.end, s.chain)
            else:
                segs.append(deepcopy(s))
        return Domain(*segs)
    merged = property(getMerged, doc='Same domain with adjacent segments merged')
    def getSize(self):        
        return sum(segment.size for segment in self.segments)
    size = property(getSize, doc='''Total (sum) length of this domain''')
    ######################################################################################################
    # Domain manipulation (addition, substraction)    
    def addSegment(self, segment):
        ''' Add a new segment '''
        self.append(segment)
    def mergeDomain(self, other):
        ''' Merge with segments from another domain '''
        self.segments |=  other.segments
    ######################################################################################################
    # Calculation of energies
    def internalEnergy(self, rcs):
        ''' Calculate internal energy of this domain, according to residue contact information in rcs.'''
        result = 0.
        selfrs = frozenset(x for x in self)
        seen = set()
        for r1 in selfrs:
            seen.add(r1)
            res1, chain1 = r1
            for r2 in selfrs:
                if r2 in seen:
                    continue
                res2, chain2 = r2
                try:
                    x = rcs[chain1][res1][chain2][res2].ev_lj
                    if result is None:
                        result = 0.
                    result+=x
                except KeyError:
                    pass #raise
        return result
    def externalEnergy(self, other, rcs):
        ''' Calculate the external interaction energy of this domain - with other domain,
        according to residue contact information in rcs.'''
        result = 0.
        rs1 = frozenset(x for x in self)
        rs2 = frozenset(x for x in other)
        for r1 in rs1:
            res1, chain1 = r1
            for r2 in rs2:
                res2, chain2 = r2
                try:
                    x = rcs[chain1][res1][chain2][res2].ev_lj
                    if result is None:
                        result = 0.
                    result+=x
                except KeyError:
                    pass
        return result
    def externalEnergies(self, others, rcs):
        ''' Calculate total interaction energy of this domain with all the others -
        according to residue contact information in rcs '''
        return sum(self.externalEnergy(other,rcs) for other in others)

######################################################################################################################################################################
# FUNCTIONS - parsing data
@handlify(is_method=False, mode='r')
def readResidueContacts(fh, triple2seqno):
    ''' Utility function to read in residue contact data outputted by vdw program '''
    result = defaultdict( innerDefaultdict( innerDefaultdict( dict ) ) )
    for line in fh:
        if line:
            chain1, resno1, icode1, resname1, chain2, resno2, icode2, resname2, no_contacts, ev_lj = line.strip('\n').split('\t')
            resno1, resno2, no_contacts = int(resno1), int(resno2), int(no_contacts)
            ev_lj = float(ev_lj)
            seqno1 = triple2seqno[chain1][(resno1, icode1)]
            seqno2 = triple2seqno[chain2][(resno2, icode2)]
            result[chain1][seqno1][chain2][seqno2] = StructObject(resname1=resname1, resname2=resname2, no_contacts=no_contacts, ev_lj=ev_lj)
            result[chain2][seqno2][chain1][seqno1] = StructObject(resname1=resname2, resname2=resname1, no_contacts=no_contacts, ev_lj=ev_lj)
    #raw_input('PAUSE')
    return result

######################################################################################################################################################################
# FUNCTIONS - calculating splitting points, energies, merging domains etc.
def calculateSplitEnergies(rcontacts, seqno2triple):
    ''' Calculating energy values for each possible split point in each of the chains. '''
    result = defaultdict( innerDefaultdict(float) )
    for chain_ in sorted(seqno2triple):
        result_for_chain = result[chain_]
        rc_chain = rcontacts[chain_]
        bounds_chain = seqno2triple[chain_]
        #print "Bounds", min(bounds_chain), max(bounds_chain)
        rnos = sorted(rc_chain)
        for divider_resno in sorted(bounds_chain):
            ev_until = 0.
            ev_after = 0.
            ev_between = 0.
            seen = set()
            for resno1 in sorted(rc_chain):
                rc_chain_resno1 = rc_chain[resno1][chain_]
                seen.add(resno1)
                for resno2 in sorted(rc_chain_resno1):
                    if resno2 in seen:
                        continue
                    if resno1<=divider_resno and resno2<=divider_resno:
                        ev_until += rc_chain_resno1[resno2].ev_lj                        
                    elif resno1>divider_resno and resno2>divider_resno:
                        ev_after += rc_chain_resno1[resno2].ev_lj
                    else:
                        ev_between += rc_chain_resno1[resno2].ev_lj
            result_for_chain[divider_resno]= StructObject(ev_until=ev_until, ev_after=ev_after, ev_between=ev_between, resno=divider_resno)
    if DEBUG:
        for chain in sorted(result):
            print "ENERGIES", chain
            for dno in sorted(result[chain]):
                print "%s\t%s" % (dno, result[chain][dno].ev_between)
    return result

MAX_FLOAT =  130000.
MIN_FLOAT = -130000.

def searchMin(i, splitvs, deltaE, pmin=MAX_FLOAT):
    ''' Search for next minimum, beginning at index i, in split points vector splitvs (minimum value is starting from MAX_FLOAT).
    Minimum must satisfy the deltaE energy threshold from the last maximum. '''
    min_cand = None
    min_ev = None
    for j in range(i, len(splitvs)):
        v = splitvs[j]
        ev = v.ev_between
        # Found the next minimum candidate
        if ev<=pmin:
            pmin = ev
            min_cand = j
        # Found next max. value satisfying deltaE
        elif (ev-pmin)>deltaE:
            return min_cand, pmin, j
    return min_cand, pmin, j
    
def searchMax(i, splitvs, deltaE, pmax=MIN_FLOAT):
    ''' Search for next maximum, beginning at index i, in split points vector splitvs (maximum value is starting from MIN_FLOAT).
    Maximum must satisfy the deltaE energy threshold from the last minimum. '''
    max_cand = None
    max_ev = None
    for j in range(i, len(splitvs)):
        v = splitvs[j]
        ev = v.ev_between
        # Found the next maximum candidate
        if ev>=pmax:
            pmax = ev
            max_cand = j
        # Found next min. value satisfying deltaE
        elif (pmax-ev)>deltaE:
            return max_cand, pmax, j
    return max_cand, pmax, j

def calculateSplitPoints(splitvs, chain, delta=0.15):
    ''' Converts the list of split points into a list of extrema satisfying the JBSD_1999 criteria for extremum selection. '''
    splitvs = sorted( splitvs[chain].values(), key=attrgetter('resno') )
    # Find the lowest minimum
    E =  -min(v.ev_between for v in splitvs)
    dE = (delta * E)
    if DEBUG:
        print dE
    # Find minima and maxima (at deltaE barrier level), using the approximate approach
    search_for_max = False
    segment_points = []
    if splitvs:
        splitvs[0].what='maximum'
        segment_points = [ splitvs[0], ]
    i = 1
    while i<len(splitvs)-1:
        if search_for_max:
            max_i, max_ev, i = searchMax(i, splitvs, dE)
            if DEBUG:
                print "FOUND: Max", max_i, max_ev, i, splitvs[max_i]
            splitvs[max_i].what='maximum'
            segment_points.append(splitvs[max_i])
            search_for_max = False
        else:
            min_i, min_ev, i = searchMin(i, splitvs, dE)
            if DEBUG:
                print "FOUND: Min", min_i, min_ev, i, splitvs[min_i]
            splitvs[min_i].what='minimum'
            segment_points.append(splitvs[min_i])
            search_for_max=True
    # If the last extremum is a minimum - add maximum in the end
    if segment_points and segment_points[-1].what=='minimum':
        x = copy(splitvs[-1])
        x.what='maximum'
        segment_points.append(x)
    # If the last extremum is a maximum - change maximum's position to 
    elif segment_points and segment_points[-1].what=='maximum' and segment_points[-1].resno!=splitvs[-1].resno:
        x = copy(splitvs[-1])
        x.what='maximum'
        segment_points[-1]=x
    if DEBUG:
        for sp in segment_points:
            print sp.__dict__
    return segment_points

from bisect import bisect
def splitPoints2Segments(split_points, residue_contacts, chain, seqno2triple):
    ''' Convert a list of split point extrema (for a given chain) into segments. Uses chain boundary data from residue contacts. '''
    r = []
    # No split points - no other checks are needed
    if len(split_points)==0:
        return r
    # Pairing maxima one by one
    prev_max = None
    seqsorted = sorted(seqno2triple[chain])
    for i,p in enumerate(split_points):
        if p.what=='maximum':
            if prev_max is not None:
                if not r:
                    coords = (prev_max.resno, p.resno)
                else:
                    coords = ( seqsorted[bisect(seqsorted, prev_max.resno)], p.resno)
                if DEBUG:
                    print coords
                r.append( Segment(coords[0], coords[1], chain, seqno2triple) )
            else:
                pass
            prev_max = p
    # If last segment was not created to span the entire rest of chain - let it do so
    # If nothing created - everything is one segment
    last_rno = seqsorted[-1]
    if not len(r):
        begin = seqsorted[0]
        r.append( Segment(begin, last_rno, chain, seqno2triple) )
    # Else check if one last segment needs to be created
    #elif r[-1].end != last_rno:
    #    # Find the next residue in the chain, after the end of the current segment (accounting for discontinuities)
    #    begin = seqsorted[ bisect(seqsorted, r[-1].end) ]
    #    r.append( Segment(begin, last_rno, chain, seqno2triple) )
    return r

# Defining an analogue of numpy's fromfunction working on Python lists of 2D
def fromfunction2d(func, lens):
    ''' Calls a function func over all the pairs of indices from lens=(l1,l2) grid.'''
    # Generate indices
    len1,len2 = lens
    r = []
    for i in range(0,len1):
        r_ = []
        r.append(r_)
        for j in range(0, len2):
            r_.append( func(i, j) )
    return r
# Defining an analogue of numpy's fromfunction working on Python lists of 2D (symmetric meaning f(i,j)==f(j,i)
def fromfunction2dSymm(func, lens):
    # Generate indices
    len1,len2 = lens
    r = [ range(0,len2) for i in range(0,len1) ]
    for i in range(0,len1):
        for j in range(i, len2):
            r[i][j] = r[j][i] = func(i,j)
    return r

def mergeSegments(segments, residue_contacts):
    if TIMEBENCH:
        t1 = time.time()
    ''' Merge segments according to rules from JBSD_1999 paper.'''
    # 0) Each segment becomes a domain
    domains = [Domain(segment) for segment in segments]
    if TIMEBENCH:
        print len(domains)
    # 1) For each segment - its internal energy, its interaction energies with the other segments,
    # and its integral energy of interactions with others. This yields a K*K matrix.    
    # Function calculating energies for every pair of segments:
    def energy(i,j):
        i,j = int(i), int(j)
        if i==j:
            x = -domains[i].internalEnergy(residue_contacts)
        else:
            x = -domains[i].externalEnergy(domains[j], residue_contacts)
        return x
    seg_e = fromfunction2dSymm( energy, [len(domains), len(domains)] )
    if TIMEBENCH:
        t1a = time.time()
    width_seg_e = len(seg_e)
    # 2) Decision which segments to merge. Forms an adjacency matrix for the merge graph.
    # Function computing the JBSD_1999 criteria for merging a pair of segments:
    def should_merge(i,j):
        i,j = int(i), int(j)
        if i==j:
            return False
        else:
            int_e = seg_e[i][i]
            int_e_j = seg_e[j][j]
            ext_e_i_wo_j = sum(seg_e[i][k] for k in range(0, len(domains)) if k!=i and k!=j)
            ext_e_i_wth_j = ext_e_i_wo_j + seg_e[i][j]
            ext_e_j_wo_i = sum(seg_e[j][k] for k in range(0, len(domains)) if k!=i and k!=j)
            ext_e_j_wth_i = ext_e_j_wo_i + seg_e[i][j]
            # If the number of domains is more than 2 !
            # Energy isolated domains are those with e_ii >= 3*sum(j!=i:e_ij). Potential domains i,k are merged if :
            #       a) e_ik >= sum(j!=i,k:e_ij) and e_ki >= sum(j!=i,k:e_kj) or,
            #       b) e_ik >= 0.7 * sum(j!=i:e_ij) or e_ki >= 0.7*sum(j!=k:e_kj)- more than 70% of external interactions of one domain pertains to the other domain
            if width_seg_e>2 and int_e < CRITICAL_MERGE_ERATIO * ext_e_i_wth_j and int_e_j < CRITICAL_MERGE_ERATIO * ext_e_j_wth_i:
                # First criterion - if both interactions exceed the other
                if seg_e[i][j] > ext_e_i_wo_j and seg_e[j][i] > ext_e_j_wo_i:
                    return True
                # Second criterion - if >= 70% of interactions are with each other
                elif seg_e[i][j] > CRITICAL_MERGE_PERCENT*ext_e_i_wth_j or seg_e[j][i] > CRITICAL_MERGE_PERCENT * ext_e_j_wth_i:
                    return True
            # Third criterion - low-energy segment is always joined to what it has best interactions with
            # Isolated segment with e_ii <= 1.5*sum(j!=i:e_ij) will be joined with the segment or potential domain with
            # which it has maximal interaction energy.
            if ext_e_i_wth_j>=0:
                try:
                    max_e = max(seg_e[i][k] for k in range(0,len(domains)) if k!=i and k!=j)
                    # If it had as good a choice before
                    #if len([seg_e[i][k_]==max_e for k_ in range(0,k-1)]):
                    #    return False
                except:
                    max_e = MIN_FLOAT
                if DEBUG:
                    print domains[i], "size:", domains[i].size
                # Low energy segment due to low internal energy
                if (int_e <= LOWE_SEGMENT_ERATIO * ext_e_i_wth_j or domains[i].size<=15) and seg_e[i][j]>=max_e:
                    return True
            else:
                try:
                    max_e = max(seg_e[i][k] for k in range(0,len(domains)) if k!=i and k!=j)
                except:
                    max_e = MIN_FLOAT
                # Low energy segment due to prevalence of repulsive term
                if seg_e[i][j]>=max_e and seg_e[i][j]>=0.:
                    return True
            return False
    to_merge = fromfunction2d(should_merge, [len(domains), len(domains)])
    if DEBUG:
        print "[",
        for row in to_merge:
            print row
        print "]"
        print "[",
        for row in seg_e:
            print row
        print "]"
        
    
    # 3) Merge the domains.
    # This requires finding all separated (unconnected with each other) components of the
    # merge graph (i.e. answers which segments will belong to same domains after merge).
    dom_inds = [ frozenset([i]) for i in range(0,len(domains)) ]
    # Make the graph undirected to_merge(i,j)=to_merge(j,i)
    for i in range(0, len(domains)):
        for j in range(0, len(domains)):
            if to_merge[i][j] or to_merge[j][i]:
                to_merge[i][j] = to_merge[j][i] = True

    # A simple depth-first search visitor keeping track of all reached vertices
    def visit(i, to_merge, visited, result=None):
        if result is None:
            result = set()
        result.add(i)
        for j, b_v in enumerate(to_merge[i]):
            if j not in visited and b_v:
                visited.add(j)
                visit(j, to_merge, visited, result)
        return frozenset(result)

    # DFS the graph noting index for each component in domains array
    comps = set()
    visited = set()
    for i in range(0, len(domains)):
        if i not in visited:
            comps.add( visit(i, to_merge, visited) )

    # Transform found components into full domains
    result = []
    for inds in comps:
        segments = []
        for i in inds:
            segments.extend(domains[i].segments)
        result.append( Domain(*segments) )
        
    # Return results sorted by domain beginning position
    if TIMEBENCH:
        t2 = time.time()
        print >> sys.stderr, "Time elapsed before merging decision: %s" % (t1a-t1,)
        print >> sys.stderr, "Total time elapsed in mergeSegments: %s" % (t2-t1,)
    return sorted( result, key=lambda x: list(x.sorted_segments)[0].begin)

def createDomainsForChains(rcs, seqno2triple, triple2seqno, deltas=DEFAULT_DELTAS):
    ''' Create a dictionary of segments and final domains for each chain and for each deltaE level '''
    r = defaultdict(dict)
    split_evs = calculateSplitEnergies(rcs, seqno2triple)
    for delta in deltas:
        for chain in sorted(split_evs):
            segment_points = calculateSplitPoints( split_evs, chain, delta )
            segments = splitPoints2Segments(segment_points, split_evs, chain, seqno2triple)
            final_domains = mergeSegments(segments, rcs)
            r[chain][delta] = final_domains
    return r

def createMDomainsForChains(rcs, seqno2triple, triple2seqno, deltas=DEFAULT_DELTAS):
    ''' Create a dictionary of segments and final domains for each chain and for each deltaE level '''
    r = defaultdict(dict)
    split_evs = calculateSplitEnergies(rcs, seqno2triple)
    for delta in deltas:
        fsegments = []
        for chain in sorted(split_evs):
            segment_points = calculateSplitPoints( split_evs, chain, delta )
            segments = splitPoints2Segments(segment_points, split_evs, chain, seqno2triple)
            fsegments.extend(segments)
            final_domains = mergeSegments(segments, rcs)
            r[chain][delta] = final_domains
        r['*'][delta] = mergeSegments(fsegments, rcs)
    return r

def prepareDomainsForChains(residue_cs_fname, bounds_fname, deltas=DEFAULT_DELTAS):
    ''' Prepares domains for chains on basis of the files. WARNING: returns Python data structure, not filename '''
    seqno2triple, triple2seqno = readChainBounds(bounds_fname)
    return createMDomainsForChains( readResidueContacts(residue_cs_fname, triple2seqno), seqno2triple, triple2seqno, deltas )

def doms2segments(doms):
    rdoms = defaultdict(dict)
    for chain in doms:
        cdoms = doms[chain]
        for delta in cdoms:
            dlist = cdoms[delta]
            segdoms = []
            for d in dlist:
                for s in d.sorted_segments:
                    segdoms.append(Domain(s))
            rdoms[chain][delta] = segdoms 
    return rdoms
#####################################################################################################################################################################
# Functions - output
def allInDefault(domain):
    ''' Checks if all the domains segments are in the null chain '''
    for segment in domain.segments:
        if segment.chain.strip():
            return False
    return True

from copy import deepcopy, copy
def reduceDomain(d):
    segs = []
    for s in d.sorted_segments:
        if len(segs) and s.begin==segs[-1].end+1:
            #print s.bc_parent
            segs[-1] = Segment(segs[-1].begin, s.end, s.chain, s.bc_parent)
        else:
            segs.append(copy(s))
    return Domain(*segs)

def dom2str(domain, merge=False):
    ''' Converts a domain to string representation '''
    if merge:
        domain = reduceDomain(domain)
    if allInDefault(domain):
        return "( %s )" % (', '.join( str(segment) for segment in domain.sorted_segments),)
    else:
        return "( %s )" % (', '.join( str(segment) for segment in domain.sorted_segments),)

def doms2str(domains,merge=False):
    ''' Converts a domain to string representation ''' 
    return '/'.join(dom2str(d, merge) for d in domains)

def seg2str(segment):
    ''' Converts a segment to string representation '''
    return dom2str( Domain(segment) )
def segs2str(segments):
    ''' Converts multiple segments to string representating the whole iterable. '''
    return '/'.join(seg2str(s) for s in segments)

import csv
@handlify(mode='w')
def writeDomainsTsv(fh, domains, uid):
    ''' Write domains into TSV format. Optionally give the unique PDB identifier (first column in the TSV file).'''
    writer = csv.writer(fh, dialect='excel-tab')
    for chain in sorted(domains):
        chd = domains[chain]
        for delta in sorted(chd, reverse=True):
            dlist = chd[delta]
            writer.writerow( [uid, chain, delta, '/'.join([str(domain) for domain in dlist])] )
@handlify(mode='r')
def readDomainsTsv(fh, seqno2triple=None, triple2seqno=None):
    ''' Read domains in tab-separated format into a dictionary of dicts of dicts indexed by PDB uid, chain and delta hierarchy level.'''
    reader = csv.reader(fh, dialect='excel-tab')
    r = defaultdict( innerDefaultdict(dict) )
    for row in reader:
        #print row
        uid, chain, delta = row[0], row[1], row[2]
        delta = float(delta)
        r[uid][chain][delta] = list(Domain.parseString(row[3], seqno2triple, triple2seqno))
    return r

@handlify(mode='w')
def writeEnergiesTsv(fh, domains, rcs, uid=None):
    ''' Write domain energies into TSV format. Optionally give the unique PDB identifier (first column in the TSV file).'''
    writer = csv.writer(fh, dialect='excel-tab')
    for chain in sorted(domains):        
        chd = domains[chain]
        for delta in sorted(chd, reverse=True):            
            dlist = chd[delta]
            # Re-sort by beginning of domain
            dlist = sorted( dlist, key=lambda x: list(x.sorted_segments)[0].begin)
            r = defaultdict(dict)
            for i, dom1 in enumerate(dlist):
                sdom1 = str(dom1)
                r[sdom1][sdom1] = dom1.internalEnergy(rcs)
                for dom2 in dlist[i+1:]:
                    sdom2 = str(dom2)                    
                    r[ sdom1 ][ sdom2 ] = r[ sdom2 ][ sdom1 ] = dom1.externalEnergy(dom2, rcs)
            sdlist = [str(x) for x in dlist]
            writer.writerow(["#STARTMATRIX"])
            writer.writerow( ["#HEADER", uid, chain, delta] )
            writer.writerow( ["#NAMES","-"]+sdlist )
            for k1 in sdlist:
                rk = r[k1]                
                writer.writerow( ["#MATRIX", k1] + [ str(rk[k2]) for k2 in sdlist ] )
            writer.writerow(["#ENDMATRIX"])

# TODO:
# Write readEnergiesTsv as well
