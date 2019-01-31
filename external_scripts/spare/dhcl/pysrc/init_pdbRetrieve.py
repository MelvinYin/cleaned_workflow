#! /usr/bin/python2.5
##########################################################################################################################
# DESCRIPTION: retrieves PDB files in a systematic way (can keep them compressed or uncompressed) - compression is done
# through bzip or gzip
##########################################################################################################################
# Default behavior for __future__ Python options
from __future__ import with_statement
from __future__ import division
##########################################################################################################################
# BASIC CONFIG
# Loading SiteProf data sources configuration
# see: scripts/config/siteprof_config.py for details
import os, os.path, sys
#sys.path.append(os.path.join(os.environ['SITEPROF'], 'config'))
#from siteprof_config import *
##########################################################################################################################
# SCRIPT
from dhcl.init_subs import defaultInitScript
from dhcl.sources import *

def makeSource(opts, log):
    return FtpPdbSource()    
def makeStorage(opts, log):
    return LocalPdbSource(opts.pdb_source)
def savePdb(uid, source, storage, opts, log,):
    return storage.store( uid, source.retrieve(uid) )

def addOptions(parser):
    parser.add_option('--pdb-source', dest='pdb_source', default=None,
                      help='Indicate a PDB source to put data into [%default]')
# COMPUTATION LOOP
defaultInitScript("dhcl.ini", source_f=makeSource, storage_f=makeStorage, compute_f=savePdb, error_f=None, logname=sys.stderr, adder_f=addOptions)
