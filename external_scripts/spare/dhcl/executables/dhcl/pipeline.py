#####################################################################################################################################################################
#
# Grzegorz M Koczyk (2007-)
# 
# Pipeline building module for automatically generating domain hierarchy level information
#
#####################################################################################################################################################################
from __future__ import with_statement
from dhcl.pdb_data import *
from dhcl.utils import *
from dhcl.hierarchy import *
from dhcl.hbonds import *
import math, gzip, logging, time, traceback
#####################################################################################################################################################################
# Constants

#####################################################################################################################################################################
# Functions - computation
def doUntilVdws(fname):
    ''' Calculate residue contacts, chain boundaries and split point energies for a given PDB file.
    Prepare the file with recalculated hydrogen atoms, first. The results are returned as a tuple of
    the above-mentioned data objects.'''
    with tempDir() as tmp_dir:
        final_segments = []
        hfname = prepareWithHydrogens(fname, "query.pdb")
        atom_cs_fname, residue_cs_fname, bounds_fname = prepareVdws(hfname)
        seqno2triple, triple2seqno = readChainBounds(bounds_fname)
        rcs = readResidueContacts(residue_cs_fname, triple2seqno)
        return rcs, seqno2triple, triple2seqno
def doUntilSplitEnergies(fname):
    ''' Calculate residue contacts, chain boundaries and split point energies for a given PDB file.
    Prepare the file with recalculated hydrogen atoms, first. The results are returned as a tuple of
    the above-mentioned data objects.'''
    with tempDir() as tmp_dir:
        final_segments = []
        hfname = prepareWithHydrogens(fname, "query.pdb")
        atom_cs_fname, residue_cs_fname, bounds_fname = prepareVdws(hfname)
        seqno2triple, triple2seqno = readChainBounds(bounds_fname)
        rcs = readResidueContacts(residue_cs_fname, triple2seqno)
        return rcs, seqno2triple, triple2seqno, calculateSplitEnergies(rcs, seqno2triple)

def predictDomains(fname, deltas=DEFAULT_DELTAS, merge_chains=True):
    ''' Iteratively yield segments and final domains for each chain and for each deltaE level '''
    rcs, seqno2triple, triple2seqno, split_evs = doUntilSplitEnergies(fname)
    for delta in deltas:
        if merge_chains:
            fsegments = []
        for chain in sorted(split_evs):
            segment_points = calculateSplitPoints( split_evs, chain, delta )
            segments = splitPoints2Segments(segment_points, split_evs, chain, seqno2triple)
            if merge_chains:
                fsegments.extend(segments)
            final_domains = mergeSegments(segments, rcs)
            # Domain-ize segments ;)
            segments = [Domain(_) for _ in segments]
            yield chain, delta, segments, final_domains
        if merge_chains:
            chain = 'SUM'
            final_domains = mergeSegments(fsegments, rcs)
            # Domain-ize segments ;)
            segments = [Domain(_) for _ in fsegments]
            yield chain, delta, segments, final_domains

def predictDomainsForChainAndLevel(fname, chain, delta):
    ''' Predict domain just for the given chain and level '''
    rcs, seqno2triple, triple2seqno, split_evs = doUntilSplitEnergies(fname)
    segment_points = calculateSplitPoints( split_evs, chain, delta )
    segments = splitPoints2Segments(segment_points, split_evs, chain, seqno2triple)
    final_domains = mergeSegments(segments, rcs)
    segments = [Domain(_) for _ in segments]
    return segments, final_domains
#####################################################################################################################################################################
# Functions - computation of loops
try:
    from decomposition import *
    def doUntilContacts(fname):
        ''' Calculate Calpha-Calpha contacts and chain boundaries.
        Prepare the file, first. The results are returned as a tuple of
        the above-mentioned data objects.'''
        with tempDir() as tmp_dir:
            final_segments = []
            hfname = preparePdb(fname, "query.pdb")
            atom_cs_fname, bounds_fname = prepareContacts(hfname)
            seqno2triple, triple2seqno = readChainBounds(bounds_fname)
            rcs = readAtomDistsForResidues(atom_cs_fname, triple2seqno)                
            return rcs, seqno2triple, triple2seqno
    def predictLoops(fname):
        ''' Predict 15-45 size loops (returns of the trajectory) comprising this protein  '''
        rcs, seqno2triple, triple2seqno = doUntilContacts(fname)
        return loopDecompose(rcs, seqno2triple)
except ImportError:
    pass

try:
    from locks import *
    def predictLocks(fname):
        ''' Predict locks for a given PDB file  '''
        rcs, seqno2triple, triple2seqno = doUntilVdws(fname)
        return computeLocks(rcs, seqno2triple)
except ImportError:
    pass

def runAllInDir(fname, out_dir='.', pdb_id='query'):
    ''' Runs all the analyses in current directory out_dir using a given pdb input file - fname. '''
    fname, out_dir = os.path.abspath(fname), os.path.abspath(out_dir)
    os.chdir(out_dir)
    try:
        # Make the hydrogens
        hydr_fname = pdb_id + '.hhh.pdb'
        hydr2_fname = hydr_fname + '.gz'
        calpha_cs_fname = pdb_id + '.calpha.cs'
        bounds_fname = pdb_id + '.bounds'
        atom_cs_fname = pdb_id + '.atom.cs'
        residue_cs_fname = pdb_id + '.residue.cs'
        doms_fname = pdb_id + '.domains'
        seg_energies_fname = pdb_id + '.segment.energies'
        dom_energies_fname = pdb_id + '.domain.energies'
        cdoms_fname = pdb_id + '.complex.domains'
        cdom_energies_fname = pdb_id + '.complex.domain.energies'
        locks_fname = pdb_id + '.locks'
        loops_fname = pdb_id + '.loops'
        hb_fname = pdb_id + '.hbonds'
        # Prepare packed .pdb with hydrogens
        if  not hasContent(hydr_fname):
            prepareWithHydrogens(fname, hydr_fname)
        # Make the vdW contacts
        if   not hasContent(residue_cs_fname):
            prepareVdws(hydr_fname, atom_cs_fname, residue_cs_fname, bounds_fname)
            os.unlink(atom_cs_fname)
        # Make the Calpha-Calpha contacts
        if  not hasContent(calpha_cs_fname) or not hasContent(bounds_fname):
            prepareContacts(hydr_fname, calpha_cs_fname, bounds_fname)
        # Make the domains for individual chains
        doms = chains = None
        if  not hasContent(doms_fname):
            doms = prepareDomainsForChains(residue_cs_fname, bounds_fname)
            writeDomainsTsv( doms_fname, prepareDomainsForChains(residue_cs_fname, bounds_fname), pdb_id )
        # Make the energies file for individual segments
        if  not hasContent(seg_energies_fname):
            seqno2triple, triple2seqno = readChainBounds(bounds_fname)
            rcs = readResidueContacts(residue_cs_fname, triple2seqno)
            if doms is None:
                doms = prepareDomainsForChains(residue_cs_fname, bounds_fname)
            segdoms = doms2segments(doms)
            writeEnergiesTsv( seg_energies_fname, segdoms, rcs, pdb_id )
        # Make the energies file for whole domain interactions
        if  not hasContent(dom_energies_fname):
            seqno2triple, triple2seqno = readChainBounds(bounds_fname)
            rcs = readResidueContacts(residue_cs_fname, triple2seqno)
            if doms is None:                
                doms = prepareDomainsForChains(residue_cs_fname, bounds_fname)
            writeEnergiesTsv(dom_energies_fname, doms, rcs, pdb_id)
        # Make the loops
        if  not hasContent(loops_fname):
            writeLoopsTsv( loops_fname, prepareLoops(calpha_cs_fname, bounds_fname), pdb_id )
        # Make the locks
        if not hasContent(locks_fname):
            writeLocksTsv(locks_fname, prepareLocks(residue_cs_fname, bounds_fname), pdb_id )
        # Make the hydrogen bonds
        # pass
        # Compress the hydrogens .pdb file if it does not already exist
        if  not hasContent(hydr2_fname):
            subprocess.call('gzip %s' % hydr_fname, shell=True)
        else:
            os.unlink(hydr_fname)
        markAsDone(pdb_id)
    except Exception, e:
        raise

def makeResultDir(fname, out_dir, force=False, ignore_done=False, ino=None):
    ''' Calculates result for a single PDB id - if needed (result files do not already exist or force has been specified). Puts result files into pdb_id subdirectory of out_dir.'''
    pdb_id = fname2pdb(fname)
    pdb_fname = os.path.abspath(fname)
    if not os.path.exists(out_dir):
        os.mkdir(out_dir)
    t1 = t2 = time.time()
    with condMakeDir( os.path.join(out_dir,pdb_id) ) as work_dir:
        if existsDone(pdb_id) and not force and not ignore_done:
             t2 = time.time()
             logging.info('%d\tALREADY EXISTS: %s, TIME: %.1f' % (ino, pdb_id, t2-t1) )
        else:
            try:
                # Make the hydrogens
                hydr_fname = pdb_id + '.hhh.pdb'
                hydr2_fname = hydr_fname + '.gz'
                calpha_cs_fname = pdb_id + '.calpha.cs'
                bounds_fname = pdb_id + '.bounds'
                atom_cs_fname = pdb_id + '.atom.cs'
                residue_cs_fname = pdb_id + '.residue.cs'
                doms_fname = pdb_id + '.domains'
                seg_energies_fname = pdb_id + '.segment.energies'
                dom_energies_fname = pdb_id + '.domain.energies'
                cdoms_fname = pdb_id + '.complex.domains'
                cdom_energies_fname = pdb_id + '.complex.domain.energies'
                # Complex domains
                domst_fname = cdoms_fname
                domst_energies_fname = cdom_energies_fname
                segst_energies_fname = pdb_id + '.complex.segment.energies'
                
                locks_fname = pdb_id + '.locks'
                loops_fname = pdb_id + '.loops'
                hb_fname = pdb_id + '.hbonds'
                # Prepare packed .pdb with hydrogens
                if force or not hasContent(hydr_fname):
                    #print pdb_fname, hydr_fname
                    prepareWithHydrogens(pdb_fname, hydr_fname)
                # Make the vdW contacts
                if force or  not hasContent(residue_cs_fname):
                    prepareVdws(hydr_fname, atom_cs_fname, residue_cs_fname, bounds_fname)
                    os.unlink(atom_cs_fname)
                # Make the Calpha-Calpha contacts
                if force or not hasContent(calpha_cs_fname) or not hasContent(bounds_fname):
                    prepareContacts(hydr_fname, calpha_cs_fname, bounds_fname)
                # Make the domains for individual chains
                #doms = chains = None
                #if force or not hasContent(doms_fname):
                #    doms = prepareDomainsForChains(residue_cs_fname, bounds_fname)
                #    writeDomainsTsv( doms_fname, doms, pdb_id )
                # Make the energies file for individual segments
                #if force or not hasContent(seg_energies_fname):
                #    seqno2triple, triple2seqno = readChainBounds(bounds_fname)
                #    rcs = readResidueContacts(residue_cs_fname, triple2seqno)
                #    if doms is None:
                #        doms = prepareDomainsForChains(residue_cs_fname, bounds_fname)
                #    segdoms = doms2segments(doms)
                #    writeEnergiesTsv( seg_energies_fname, segdoms, rcs, pdb_id )
                # Make the energies file for whole domain interactions
                #if force or not hasContent(dom_energies_fname):
                #    seqno2triple, triple2seqno = readChainBounds(bounds_fname)
                #    rcs = readResidueContacts(residue_cs_fname, triple2seqno)
                #    if doms is None:                
                #        doms = prepareDomainsForChains(residue_cs_fname, bounds_fname)
                #    writeEnergiesTsv(dom_energies_fname, doms, rcs, pdb_id)
                # Extended: make Full Domains - this will likely supersed earlier results
                # Make the domains for total
                domst = chainst = None
                if force or not hasContent(domst_fname):
                    domst = prepareDomainsForChains(residue_cs_fname, bounds_fname)
                    writeDomainsTsv( domst_fname, prepareDomainsForChains(residue_cs_fname, bounds_fname), pdb_id )
                # Make the energies file for individual segments
                if force or not hasContent(segst_energies_fname):
                    seqno2triple, triple2seqno = readChainBounds(bounds_fname)
                    rcs = readResidueContacts(residue_cs_fname, triple2seqno)
                    if domst is None:
                        domst = prepareDomainsForChains(residue_cs_fname, bounds_fname)
                    segdomst = doms2segments(domst)
                    writeEnergiesTsv( segst_energies_fname, segdomst, rcs, pdb_id )
                # Make the energies file for whole domain interactions
                if force or not hasContent(domst_energies_fname):
                    seqno2triple, triple2seqno = readChainBounds(bounds_fname)
                    rcs = readResidueContacts(residue_cs_fname, triple2seqno)
                    if domst is None:                
                        domst = prepareDomainsForChains(residue_cs_fname, bounds_fname)
                    writeEnergiesTsv(domst_energies_fname, domst, rcs, pdb_id)

                # Make the loops
                if force or not hasContent(loops_fname):
                    writeLoopsTsv( loops_fname, prepareLoops(calpha_cs_fname, bounds_fname), pdb_id )
                # Make the hydrogen bonds
                if not hasContent(hb_fname):
                    prepareHbonds(hydr_fname, hb_fname)
                # Make the locks
                if not hasContent(locks_fname):
                    writeLocksTsv(locks_fname, prepareLocks(residue_cs_fname, bounds_fname), pdb_id )               
                # Compress the hydrogens .pdb file if it does not already exist
                if force or not hasContent(hydr2_fname):
                    subprocess.call('gzip %s' % hydr_fname, shell=True)
                else:
                    if os.path.exists(hydr_fname):
                        os.unlink(hydr_fname)
                markAsDone(pdb_id)
                t2 = time.time()
                logging.info('%d\tDONE: %s, TIME: %.1f' % (ino, pdb_id,t2-t1) )
                return True
            #except:
            #    pass
            except Exception, e:
                logging.error('%d\tERROR: occured with %s\n:%s' % (ino, pdb_id, ''.join(traceback.format_exception(sys.exc_type, sys.exc_value, sys.exc_traceback)) ) )
                return False
            #except:
            #    pass
            #    #raise

