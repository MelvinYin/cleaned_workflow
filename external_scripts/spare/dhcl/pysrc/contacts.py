#! /usr/bin/python2.5
# Calculate contact densities and individual pairs
# for a number of structures
from __future__ import with_statement
from __future__ import division
import sys, os, os.path
from optparse import OptionParser
from dhcl.pipeline import *
from dhcl.decomposition import *
from dhcl.locks import *
from dhcl.utils import *
import csv
import logging

@handlify(is_method=False, mode='w')
def writeResidueContacts(fh, rcs, seqno2triple):
    ''' Utility function to read in residue contact data outputted by vdw program '''
    result = defaultdict( innerDefaultdict( innerDefaultdict( dict ) ) )
    seen = set()
    for chain1 in sorted(rcs):
        rcs1 = rcs[chain1]
        seen2 = set()
        for seqno1 in sorted(rcs1):
            rcs12 = rcs1[seqno1]
            resno1, icode1 = seqno2triple[chain1][seqno1]
            seen2.add(seqno1)
            for chain2 in sorted(rcs12):
                if chain2 in seen:
                    continue
                rcs123 = rcs12[chain2]
                for seqno2 in sorted(rcs123):
                    if seqno2 in seen2 and chain1==chain2:
                        continue
                    v = rcs123[seqno2]
                    resno2, icode2 = seqno2triple[chain2][seqno2]
                    print >> fh, "\t".join([ str(s) for s in [chain1, resno1, icode1, v.resname1, chain2, resno2, icode2, v.resname2, v.no_contacts, v.ev_lj] ] )
        seen.add(chain1)

if __name__=='__main__':
    # Define a parser for command-line options and use it
    parser = OptionParser(usage="usage: %prog [options] <list of PDB format files>")
    parser.add_option('--with-ids', dest='with_ids', default=True, action='store_true',
                      help='attempt to coerce file names to PDB/sequence identifiers (in the output) [%default]')
    parser.add_option('-e', '--error', dest='error', default=None, action='store',
                      help='redirect error/warnings to a logfile [%default]')
    parser.add_option('--errors-append', dest='errors_append', default=False, action='store_true',
                      help='append errors to existing logfile [%default]')
    parser.add_option('-o', '--output', dest='odirname', default="output", action='store',
                      help='redirect output [%default]')
    parser.add_option('-l', '--listfile', dest='listfname', default=None, action='store',
                      help='read in a list of filenames to analyse (will then ignore the commandline list) [%default]')
    parser.add_option('-d', '--dir', dest='dir', default=None, action='store',
                      help='compute on a directory (will then ignore commandline or listfile) [%default]')
    parser.add_option('--suffix', dest='suffix', default=None, action='store',
                      help='use only on files with extension (used if dir option is specified) [%default]')
    parser.add_option('-v', '--verbose', dest='verbose', default=False, action='store_true',
                      help='give verbose output (includes errors/warnings from structure parsing etc.) [%default]')
    opts, args = parser.parse_args()
    odirname = os.path.abspath(opts.odirname)
    # Set up logging
    if opts.verbose:
        log_level=logging.DEBUG
    else:
        log_level=logging.WARNING
    log_format = '%(levelname)-8s %(message)s'
    log_mode = 'a' if opts.errors_append else 'w'
    if opts.error:
        logging.basicConfig(level=log_level,
                            format=log_format,
                            )
    else:
        logging.basicConfig(level=log_level,
                            format=log_format,
                            filename= opts.error,
                            filemode=log_mode,
                            )

    # Read in filenames either from commandline arguments or a list file
    if opts.dir:
        if not os.path.isdir(opts.dir):
            raise IOError, ''' Isn\'t a directory: %s ''' % opts.dir
        orig_fnames = [ os.path.join(opts.dir, fname) for fname in os.listdir(opts.dir) if opts.suffix is None or fname.endswith(opts.suffix) ]
    elif opts.listfname is None:
        orig_fnames = [ a for a in args]
    else:
        orig_fnames = parseSimpleList(opts.listfname)
    pdb_fnames = [ os.path.abspath(a) for a in orig_fnames ]
    # Coerce to identifiers if needed
    if opts.with_ids:
        orig_fnames = [ os.path.basename(a).split('.')[0] for a in orig_fnames]

    vals = []

    # Make output directory if it does not exists
    if os.path.exists(odirname):
        if not os.path.isdir(odirname):
            raise IOError('Not a directory %s' % odirname)
    else:
        os.mkdir(odirname)
    os.chdir(odirname)
    lengths_fn = "lengths.txt"
    summary_fn = "summary.txt"
    with open(lengths_fn, 'w') as lengths_fh:
        len_writer = csv.writer(lengths_fh, dialect='excel-tab')
        with open(summary_fn, 'w') as summary_fh:
            sum_writer =  csv.writer(summary_fh, dialect='excel-tab')
            for uid, fname in zip(orig_fnames, pdb_fnames):
                rcs, seq2triple, triple2seq = doUntilVdws(fname)
                sdikt = defaultdict(innerDefaultdict(int))
                # Write residue contacts
                idrcs_fn = uid + '.residue.cs'
                writeResidueContacts(idrcs_fn, rcs, seq2triple)
                # Write lengths                
                for chain in seq2triple:
                    len_writer.writerow([uid, chain, len(seq2triple[chain])])
                # Compute contact densities
                seen = set()
                totcon_all = 0
                totlen_all = 0
                for chain1 in sorted(rcs):
                    rcs1 = rcs[chain1]
                    length1 = len(seq2triple[chain1])
                    totlen_all += length1
                    totcon = 0
                    seqnoseen = set()
                    for seqno1 in sorted(rcs1):
                        rcs12 = rcs1[seqno1]
                        for chain2 in sorted(rcs12):
                            length2 = len(seq2triple[chain2])
                            if chain2 in seen:
                                 continue
                            rcs123 = rcs12[chain2]
                            for seqno2 in sorted(rcs123):
                                if seqno2 in seqnoseen and chain1==chain2:
                                    continue
                                v = rcs123[seqno2]
                                sdikt[chain1][chain2]+=v.no_contacts                                
                                totcon_all+=v.no_contacts
                        seqnoseen.add(seqno1)
                    seen.add(chain1)

                for chain1 in sorted(sdikt):
                    sdikt1 = sdikt[chain1]
                    length1 = len(seq2triple[chain1])
                    for chain2 in sorted(sdikt1):
                        length2 = len(seq2triple[chain2])
                        totcon = sdikt1[chain2]
                        if chain1==chain2:
                            sum_writer.writerow([ uid, chain1, chain2, totcon, (totcon/(length1)) ] )
                        else:
                            sum_writer.writerow( [uid, chain1, chain2, totcon, (totcon/(length1+length2))] )

                # Computing total contact density
                sum_writer.writerow([ uid, 'ALL', 'ALL', totcon_all, float(totcon_all)/totlen_all])
