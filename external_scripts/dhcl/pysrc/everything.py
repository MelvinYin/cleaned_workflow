#! /usr/bin/python2.5
# For every file do all calculations (segments&loops), saving them into a single file
# (with line prefixes)
# Make a prediction of domain structure for all chains of a single protein structure
import sys, os, os.path
from optparse import OptionParser
#from collections import defaultdict
from dhcl.pipeline import *
from dhcl.decomposition import *
from dhcl.utils import *
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)
import logging

if __name__=='__main__':
    # Define a parser for command-line options and use it
    parser = OptionParser(usage="usage: %prog [options] <list of PDB format files>")
    parser.add_option('--with-ids', dest='with_ids', default=False, action='store_true',
                      help='attempt to coerce file names to PDB/sequence identifiers (in the output) [%default]')
    parser.add_option('-e', '--error', dest='error', default=None, action='store',
                      help='redirect error/warnings to a logfile [%default]')
    parser.add_option('--errors-append', dest='errors_append', default=False, action='store_true',
                      help='append errors to existing logfile [%default]')
    parser.add_option('-o', '--output', dest='output', default=None, action='store',
                      help='redirect output [%default]')
    parser.add_option('-l', '--listfile', dest='listfname', default=None, action='store',
                      help='read in a list of filenames to analyse (will then ignore the commandline list) [%default]')
    parser.add_option('-s', '--segments', dest='segments', default=False, action='store_true',
                      help='output domains in a segmented form [%default]')
    parser.add_option('-d', '--dir', dest='dir', default=None, action='store',
                      help='compute on a directory (will then ignore commandline or listfile) [%default]')
    parser.add_option('--suffix', dest='suffix', default=None, action='store',
                      help='use only on files with extension (used if dir option is specified) [%default]')
    parser.add_option('--outdir', dest='outdir', default=None, action='store',
                      help='output in separate pieces to a new output directory (will override -o option) [%default]')
    parser.add_option('--out-suffix', dest='out_suffix', default='.dhcl.txt', action='store',
                      help='extension for output files (used only with --outdir option) [%default]')
    parser.add_option('-v', '--verbose', dest='verbose', default=False, action='store_true',
                      help='give verbose output (includes errors/warnings from structure parsing etc.) [%default]')
    opts, args = parser.parse_args()
    if opts.outdir:
        opts.outdir = os.path.abspath(opts.outdir)
        if not os.path.exists(opts.outdir):
            os.mkdir(opts.outdir)
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
    # Redirect output to file if needed
    if opts.output is None or opts.outdir is not None:
        fh = sys.stdout
    else:
        fh = open(opts.output, 'w').close()
    # Write out predictions
    r = defaultdict( innerDefaultdict(dict) )
    # Keep the ordering information from input
    ordering = []
    # Individual errors do not crash the whole program
    for i, pdb_fname in enumerate(pdb_fnames):
        try:
            logging.info("%s computing..." % orig_fnames[i])
            if opts.outdir is not None:
                root_fn = os.path.join( opts.outdir, os.path.basename(pdb_fname).split('.')[0]) 
                fname =  root_fn + opts.out_suffix
                hydr_pdb_fname = root_fn + '.pdb'
                # Continue to the next one if output data already exists
                if os.path.exists(fname) and endsWithDone(fname):
                    continue
                # If cleaned up PDB exists - do not clean up again
                if not os.path.exists(hydr_pdb_fname) or os.path.getsize(hydr_pdb_fname)==0:
                    prepareWithHydrogens(pdb_fname, hydr_pdb_fname)
                # Open the file, otherwise
                fh = open(fname, 'w')
            try:
                ordering.append(pdb_fname)
                
                r = predictLoops(pdb_fname)
                for chain in sorted(r):
                    print >> fh, '\t'.join([ "LOOPS", orig_fnames[i], chain, '/'.join(str(l) for l in r[chain]) ])
                for chain, delta, segments, domains in predictDomains(pdb_fname):
                    if opts.segments:
                        print >> fh, '\t'.join( stringize([ "DOMAINS", orig_fnames[i], chain, delta, doms2str(domains)]) )
                    else:
                        print >> fh, '\t'.join(stringize([  "DOMAINS", orig_fnames[i], chain, delta, doms2str(domains, True)]) )
                if opts.outdir is not None:
                    print >> fh, "#DONE"
            finally:
                if opts.outdir is not None:
                    fh.close()
            logging.info("%s computed OK" % orig_fnames[i])
        except:
            logging.error("Problem computing data for %s" % orig_fnames[i])
            raise
        
