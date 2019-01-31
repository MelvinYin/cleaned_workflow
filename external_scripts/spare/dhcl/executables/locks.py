#! /usr/bin/python2.5
# Compute locks on a subset of structures
import sys, os, os.path
from optparse import OptionParser
from dhcl.pipeline import *
from dhcl.decomposition import *
from dhcl.locks import *
from dhcl.utils import *
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
    parser.add_option('-d', '--dir', dest='dir', default=None, action='store',
                      help='compute on a directory (will then ignore commandline or listfile) [%default]')
    parser.add_option('--suffix', dest='suffix', default=None, action='store',
                      help='use only on files with extension (used if dir option is specified) [%default]')
    parser.add_option('-v', '--verbose', dest='verbose', default=False, action='store_true',
                      help='give verbose output (includes errors/warnings from structure parsing etc.) [%default]')
    opts, args = parser.parse_args()
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
    if opts.output is None:
        fh = sys.stdout
    else:
        fh = open(opts.output, 'w').close()

    # Keep the ordering information from input
    ordering = []
    # Individual errors do not crash the whole program
    for i, pdb_fname in enumerate(pdb_fnames):
        try:
            logging.info("%s computing..." % orig_fnames[i])
            r = predictLocks(pdb_fname)
            writeLocksTsv(fh, r, orig_fnames[i])
        except:
            logging.error("Problem computing data for %s" % orig_fnames[i])
            raise
        
