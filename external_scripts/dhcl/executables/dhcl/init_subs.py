#####################################################################################################################
##
## The module is part of SiteProf2 framework. It serves to provide basic framework for 'database-filling' scripts.
## It may not be copied/changed without explicit permission of the author.
## Copyright: Grzegorz M. Koczyk (2005-2007)
##
#####################################################################################################################    
# Subroutines used to build init scripts
from __future__ import division
from __future__ import with_statement

import sys, os, os.path
from operator import itemgetter, attrgetter
from dhcl.utils import parseSimpleSet
from dhcl.utils import IniOptionParser
from dhcl.source_types import *
import logging, shutil

KEEP_LOGS=5

def useLog(f, rotate=True, format="%(levelname)-10s %(asctime)s %(message)s", level=logging.INFO, datefmt='%d/%m/%Y %H:%M:%S', name='siteprof'):
    ''' Tells the logging module to use a given logfile from now on - will rotate logs, if told to do so. Alternatively
    a stream can be given instead of filename (no rotation is performed then) '''
    logger = logging.getLogger(name)
    ft = logging.Formatter(format, datefmt)
    if isinstance(f, basestring):
        # Rotate old logs
        if rotate:
            dirn = os.path.dirname(f)
            base = os.path.basename(f)
            to_shift=[]
            for fname in os.listdir(dirn):
                if fname.startswith(base):
                    full_fname = os.path.join(dirn, fname)
                    try:
                        logn, which = fname.rsplit('.',1)
                        which = int(which)
                        if which>=KEEP_LOGS:
                            os.unlink(full_fname)
                        else:
                            to_shift.append([full_fname, os.path.join(dirn, "%s.%d" % (logn, which+1,)), which])
                    except:
                        pass
            to_shift.sort( key=itemgetter(2), reverse=True )
            for old, new, which in to_shift:
                shutil.move(old, new)
            if KEEP_LOGS>0 and os.path.exists(f):
                shutil.move(f, f+'.1')
        h = logging.FileHandler(f, 'w')                
    else:
        h = logging.StreamHandler(f)        
    h.setFormatter(ft)
    logger.setLevel(level)
    logger.addHandler(h)
    return logger

def parseInitOptions(ini_file, logname=sys.stderr, adder_f=None):
    parser = IniOptionParser(ini_file)
    # DRY-RUN - LOG ONLY - DO NOT COMPUTE
    parser.add_option('-d', '--dry-run', dest='dry_run', default=False, action='store_true',
                      help='Output log messages only - do not compute anything [%default]')
    # OPTIONS AFFECTING CHOICE OF UIDS TO COMPUTE
    parser.add_option('-s', '--single-id', dest='single_id', default=False,
                      help='Indicate a single id for which to compute [%default]')
    parser.add_option('-b', '--blacklist', dest='blacklist', default=False,
                      help='Indicate a blacklist file with identifiers to ignore [%default]')
    parser.add_option('-w', '--whitelist', dest='whitelist', default=False,
                      help='Indicate a whitelist file with identifiers to be used [%default]')
    parser.add_option('-f', '--force', dest='force', default=False, action='store_true',
                      help='Force recomputation, if result already exists [%default]')
    parser.add_option('--force-errors', dest='force_errors', default=False, action='store_true',
                      help='Force recomputation of old errors/empty results [%default]')
    # OPTIONS AFFECTING LOGGING
    parser.add_option('-l', '--logfile', dest='logfile', default=logname,
                      help='Direct output to SiteProf logfile [%default]')
    parser.add_option('-e', '--stderror', dest='stderr', default=False, action='store_true',
                      help='Tee output to standard error [%default]')
    parser.add_option('-q', '--quiet', dest='quiet', default=False, action='store_true',
                      help='Quiet mode - log only computed entries [%default]')
    if adder_f is not None:
        adder_f(parser)
    opts, args = parser.parse_args()
    return opts, args

def prepareLogging(opts):
    if opts.stderr:
        log = useLog( sys.stderr )
    else:
        log = useLog( opts.logfile )
    return log

def computeLoop(compute_f, source, storage, log, opts, error_f=None):
    # Preparing to do list
    if opts.blacklist:
        blacklist = parseSimpleSet(opts.blacklist)
    else:
        blacklist = set()                               
    if opts.whitelist:
        whitelist = parseSimpleSet(opts.whitelist)
    else:
        whitelist = set()
    source_all = frozenset(source.all())
    try:
        storage_all = frozenset(storage.all() )
    except:
        storage_all = frozenset([])
    try:
        storage_errors = frozenset(storage.errors())
    except:
        storage_errors = frozenset([])

    if opts.single_id:
        to_do = set([opts.single_id])
    else:    
        to_do = source_all.copy()
    if not opts.force:
        to_do -= storage_all
    if not opts.force and not opts.force_errors:
        to_do -= storage_errors
    if whitelist:
        to_do &= whitelist
    if blacklist:
        to_do -= blacklist
    # Actual compute loop
    for i, uid in enumerate( sorted( source_all ) ):    
        if uid not in to_do:
            if not opts.quiet:
                if uid in storage_all:
                    log.info("%s\t%s\t%s" % ( i+1, uid, "ALREADY EXISTS - SKIPPING") )
                elif uid in storage_errors:
                    if storage.query(uid)==STATUS.EMPTY:
                        log.info("%s\t%s\tMARKED AS EMPTY - SKIPPING" % (i+1,uid) )
                    else:
                        log.warning("%s\t%s\tMARKED AS ERROR - SKIPPING" % (i+1, uid))
        else:
            try:
                try:
                    if opts.dry_run:
                        log.info('"%s\t%s\t%s" % ( i+1, uid, "QUEUED (DRY RUN)")')
                    else:
                        compute_f( uid, source, storage, opts, log )                
                except KeyboardInterrupt:
                    raise
                except Exception, err:
                    if error_f is not None:
                        error = error_f(uid, err, storage, log)
                    else:
                        if isinstance(err, EmptyError):
                            raise
                        else:
                            error = GoneWrongError('Job could not be completed - unknown reason')
                            try:
                                if storage is not None:
                                    storage.store(uid, error)
                            except:
                                pass
                            raise
            except KeyboardInterrupt:
                sys.exit(1)
            except EmptyError:
                log.info("%s\t%s\tEMPTY" % (i+1,uid) )
            except:
                log.exception("%s\t%s\tERROR" % (i+1, uid))
            else:            
                log.info("%s\t%s\tDONE" % (i+1,uid) )


def defaultInitScript(ini_file, source_f, storage_f, compute_f, error_f=None, logname='defaultlog', adder_f=None):
    opts, args = parseInitOptions(ini_file, logname, adder_f)
    log = prepareLogging(opts)
    source = source_f(opts, log)
    storage = storage_f(opts, log)
    computeLoop(compute_f, source, storage, log, opts, error_f=error_f)
    return log, opts
