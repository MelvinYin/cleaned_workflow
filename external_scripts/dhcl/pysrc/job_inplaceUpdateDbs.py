#! /home/fimm/cbu/gkoczyk/bin/python2.5
from __future__ import with_statement, division
# For every file do all calculations (segments&loops), saving them into a single file
# (with line prefixes)
# Make a prediction of domain structure for all chains of a single protein structure
import sys, os, os.path, subprocess
if 'DHCL_EXEC_ROOT' in os.environ:
    sys.path.append(os.environ['DHCL_EXEC_ROOT'])
#this_dir = "/home/fimm/cbu/gkoczyk/dhcl_new/executables" #os.path.dirname(os.path.abspath(sys.argv[0]))
#if this_dir not in sys.path:
from optparse import OptionParser
from dhcl.pipeline import *
from dhcl.utils import * 
import logging, traceback
import time

if __name__=='__main__':
    # Define a parser for command-line options and use it
    parser = OptionParser(usage="usage: %prog [options] <input directory> <output directory>")
    parser.add_option('-d', '--pdb-dir', dest='pdb_dir', default=os.environ.get('DHCL_PDB_DIR',None), action='store',
                      help='directory of PDB files [%default]')
    parser.add_option('-s', '--suffix', dest='suffix', default='.gz', action='store',
                      help='suffix for PDB files [%default]')
    parser.add_option('-r', '--dhcl-root', dest='dhcl_root', default=os.environ.get('DHCL_ROOT', None), action='store',
                      help='root directory for DHcL database [%default]')
    parser.add_option('-l', '--listfile', dest='list_file', default=os.environ.get('DHCL_LIST_FILE', None), action='store',
                      help='file containing a list of PDB identifiers for PDB files [%default]')
    #parser.add_option('-o', '--out-dir', dest='out_dir', default=os.environ.get('DHCL_OUT_DIR', None), action='store',
    #                  help='output master directory for DHCL results [%default]')
    parser.add_option('-v', '--verbose', dest='verbose', default=False, action='store_true',
                      help='give verbose output (includes errors/warnings from structure parsing etc.) [%default]')
    opts, args = parser.parse_args()
    try:
        DHCL_ERRORS = args[0]
    except:
        DHCL_ERRORS = os.environ['DHCL_ERRORS']
    DHCL_ROOT = opts.dhcl_root
    # 'Absolutize' and set up directories
    #in_dir, out_dir = args
    opts.pdb_dir = os.path.abspath(opts.pdb_dir)
    # Set up logging
    if opts.verbose:
        log_level=logging.INFO
    else:
        log_level=logging.INFO
    log_format = '%(levelname)-8s %(message)s'
    logging.basicConfig(level=log_level,
                        format=log_format,
                        )
    if opts.list_file is not None:
        idset = parseSimpleSet(opts.list_file)
    else:
        idset = None

    # TODO: Code this properly - not hardcode dirs
    #DHCL_ROOT = r"C:\Mine\Work\Databases\dhcl_db"
    #DHCL_ERRORS = r"C:\Mine\Work\Databases\dhcl_db\status_dict.txt"
    FULLWORK = True
    if not os.path.exists(DHCL_ERRORS):
        open(DHCL_ERRORS, 'w').close()
    #    os.mkdir(DHCL_ERRORS)
        
    DHCL_DIRS = sorted( os.path.join(DHCL_ROOT, dirn) for dirn in os.listdir(DHCL_ROOT) if  os.path.isdir(os.path.join(DHCL_ROOT, dirn))  )
    for out_dir in DHCL_DIRS:
        if not os.path.exists(out_dir):
            os.mkdir(out_dir)
    mapdikt = dict( sum([  [ (fn,dirn) for fn in os.listdir(dirn)] for dirn in DHCL_DIRS ] , []) )

    dikt_fn = DHCL_ERRORS
    idikt = parseSimpleDict(dikt_fn)
    ERROR_CODE = 'ERROR'
    DONE_CODE = 'DONE'
    def markError(pdb_id):
        #error_fname = os.path.join(DHCL_ERRORS, pdb_id + '.error')
        idikt[pdb_id] = ERROR_CODE
        open(dikt_fn, 'a').write( "%s\t%s\n" % (pdb_id, ERROR_CODE) )

    def markDone(pdb_id):
        #error_fname = os.path.join(DHCL_ERRORS, pdb_id + '.error')
        idikt[pdb_id] = DONE_CODE
        open(dikt_fn, 'a').write( "%s\t%s\n" % (pdb_id, DONE_CODE) )

    def errorMarked(pdb_id):
        try:
            return idikt[pdb_id]==ERROR_CODE
        except:
            return False
    def doneMarked(pdb_id):
        try:
            return idikt[pdb_id]==DONE_CODE
        except:
            return False
        
    # For every PDB file in the directory (if it's in the listfile or there is no listfile) - compute (PDB with hydrogens, Ca-Ca contacts, vdW contacts etc. )
    for ino,fname in enumerate(sorted(iterateDir(opts.pdb_dir, suffix=opts.suffix))):
        #print fname
        #break
        pdb_id = fname2pdb(fname)
        if idset is not None and pdb_id not in idset:
            continue
        #if not os.path.isdir(fname):
        #    print fname
        #    print ino+1, pdb_id
        #    break
        #    continue
        if errorMarked(pdb_id) or doneMarked(pdb_id):
            continue
        query_fname = os.path.join(mapdikt.get(pdb_id, DHCL_DIRS[-1]), pdb_id, pdb_id + '.hhh.pdb.gz')
        if FULLWORK and not os.path.exists(query_fname):
            query_fname = fname
        if not os.path.exists(query_fname):
            continue
        #print query_fname
        if not makeResultDir(query_fname, mapdikt.get(pdb_id, DHCL_DIRS[-1]), force=False, ignore_done=True, ino=ino+1):
            markError(pdb_id)
        else:
            markDone(pdb_id)
    
