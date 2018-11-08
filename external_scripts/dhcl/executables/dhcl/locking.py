#####################################################################################################################
##
## The module is part of SAPLing framework. It serves to provide a resource (file) locking procedures and classes.
## It may not be copied/changed without explicit permission of the author.
## Copyright: Grzegorz M. Koczyk (2005-2007)
##
#####################################################################################################################
## FILE LOCKING
import time, sys, os, os.path, random

class BaseLockingIdiom(object):
    def __enter__(self):
        return self.lock()
    def __exit__(self, type, value, traceback):
        self.unlock()
        
# WINDOWS VERSION
if sys.platform.startswith('win'):
    import win32con, win32file, pywintypes
    import os
    LOCK_EX = win32con.LOCKFILE_EXCLUSIVE_LOCK
    LOCK_SH = 0 # the default
    LOCK_NB = win32con.LOCKFILE_FAIL_IMMEDIATELY
    # is there any reason not to reuse the following structure?
    _overlapped = pywintypes.OVERLAPPED()

    class FilenameLock(BaseLockingIdiom):
        ''' A lock object to restrict access to a filename. Works by creating a lockfile.
        Last lock to exist is responsible for destroying the lockfile. Should work over older NFS. '''
        def __init__(self, fname, wait_period=3):
            # Find the exact parameters of the file we are trying to protect
            self.fname = os.path.abspath(fname)
            self.wait_period = 3
            self.lock_fname = self.fname + '.lock'
            #print self.fname, self.lock_fname
        def __del__(self):
            self.unlock()
        def _getLock(self):
            self.fh = open(self.lock_fname, 'w')
            flags = LOCK_EX | LOCK_NB
            try:
                hfile = win32file._get_osfhandle(self.fh.fileno())
                win32file.LockFileEx(hfile, flags, 0, -0x10000, _overlapped)
                #print >> sys.stderr, "Acquired"
                return True
            except:
                #print >> sys.stderr, sys.exc_info()
                #sys.stderr.flush()
                self.fh.close()
                del self.fh
            return False
        def lock(self, nonblock=False):
            if getattr(self, 'has_lock', False):
                return
            try:            
                if nonblock:
                    if not self._getLock():
                        raise IOError, ''' The lock on %s is already acquired by different object ''' % self.fname
                else:
                    while not self._getLock():
                        time.sleep( self.wait_period )
                self.has_lock = True
            finally:
                try:
                    os.unlink(self.tmp_lock_fname)
                except:
                    pass
        def unlock(self):
            if getattr(self, 'has_lock', False):
                try:
                    hfile = win32file._get_osfhandle(self.fh.fileno())
                    win32file.UnlockFileEx(hfile, 0, -0x10000, _overlapped)
                    self.fh.close()
                    del self.fh
                    os.unlink(self.lock_fname)
                except:
                    pass
            self.has_lock = False
else:
    class FilenameLock(BaseLockingIdiom):
        ''' A lock object to restrict access to a filename. Works by creating a lockfile.
        Last lock to exist is responsible for destroying the lockfile. Should work over older NFS. '''
        def __init__(self, fname, wait_period=3):
            # Find the exact parameters of the file we are trying to protect
            self.fname = os.path.abspath(fname)
            self.wait_period = 3
            self.lock_fname = self.fname + '.lock'
            # My own unique temporary lockfile for purpose of multiple locks
            self.tmp_lock_fname = self.lock_fname + "." + os.uname()[1] + "." + str(os.getpid()) + "." + str(int(100000*random.random()))
        def __del__(self):
            self.unlock()
        def _getLock(self):
            try:
                n = os.stat(self.tmp_lock_fname)[3]
            except:
                return False
            try:
                os.link(self.tmp_lock_fname, self.lock_fname)
            except:
                return False
            try:
                m = os.stat(self.tmp_lock_fname)[3]
            except:
                return False
            if n+1 == m:
                return True
            return False
        def lock(self, nonblock=False):
            if getattr(self, 'has_lock', False):
                return
            try:
                open(self.tmp_lock_fname, 'w').close()       
                if nonblock:
                    if not self._getLock():
                        raise IOError, ''' The lock on %s is already acquired by different object ''' % self.fname
                else:
                    while not self._getLock():
                        time.sleep( self.wait_period )
                self.has_lock = True
            finally:
                try:
                    os.unlink(self.tmp_lock_fname)
                except:
                    pass
        def unlock(self):
            if getattr(self, 'has_lock', False):
                try:
                    os.unlink(self.lock_fname)
                except:
                    pass
            self.has_lock = False

# import fcntl
# class FcntlFilenameLock(object):
#     ''' A lock object to restrict access to a filename. Works by creating a lockfile.
#     Last lock to exist is responsible for destroying the lockfile. Will not work over older NFS. '''
#     def __init__(self, fname, wait_period=3):
#        # Find the exact parameters of the file we are trying to protect
#        self.fname = os.path.abspath(fname)
#        self.wait_period = 3
#        self.lock_fname = self.fname + '.lock'
#     def __del__(self):
#         self.unlock()
#     def _getLock(self):
#         try:
#             self.fh = open(self.lock_fname,'w')
#             fcntl.lockf(self.fh, fcntl.LOCK_EX|fcntl.LOCK_NB)
#             return True
#         except:
#             if hasattr(self, 'fh'):
#                 fcntl.lockf(self.fh, fcntl.LOCK_UN)
#                 del self.fh
#             return False       
#     def lock(self, nonblock=False):
#         if getattr(self, 'has_lock', False):
#             return
#         try:
#             if nonblock:
#                 if not self._getLock():
#                     raise IOError, ''' The lock on %s is already acquired by different object ''' % self.fname
#             else:
#                 while not self._getLock():
#                     time.sleep( self.wait_period )
#             self.has_lock = True
#         finally:
#             pass
#     def unlock(self):
#         if hasattr(self, 'fh'):
#             fcntl.lockf(self.fh, fcntl.LOCK_UN)
#             del self.fh
#         self.has_lock=False

 
