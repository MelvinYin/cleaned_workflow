#####################################################################################################################
##
## The module is part of SiteProf2 framework. It serves to provide data source handling methods and idioms.
## It may not be copied/changed without explicit permission of the author.
## Copyright: Grzegorz M. Koczyk (2005-2007)
##
#####################################################################################################################    
from __future__ import with_statement
# Definition of DataSources
# Specific sources defined
from dhcl.source_types import *
# PDB sources
from dhcl.pdb_gk import PdbFile, TmpPdbFile
from cStringIO import StringIO
import subprocess
import os, os.path
from contextlib import contextmanager, closing

# FTP data sources
from ftplib import FTP
@contextmanager
def ftp_connection(url, user='', passwd='', acct=''):
    ''' Context manager for self-cleaning FTP connection '''
    #print url
    ftp = FTP(url)
    #sys.stdout.flush()
    try:
        ftp.login(user, passwd, acct)
        yield ftp
    finally:
        ftp.quit()

from functools import partial
class FtpDataSource(PackagedInput, PackagedOutput):
    ''' FTP data source. Transforms requests to (or from) files on the server. '''
    def __init__(self, server, path, prefix='', ext='', user='', passwd='', acct=''):
        self.server = server
        self.path = path
        self.prefix = prefix
        self.ext = ext
        self.user, self.passwd, self.acct = user, passwd, acct
    def map_request(self, request):
        return self.prefix + request + self.ext
    def all(self):
        with ftp_connection(self.server) as ftp:
            ftp.cwd(self.path)
            #print self.path
            #print ftp.nlst()
            return sorted( x[len(self.prefix):-len(self.ext)] for x in ftp.nlst() if x.endswith(self.ext) and x.startswith(self.prefix) )
    def get_content(self, tocall):
        raise NotImplementedError(''' Abstract method ''')
    @contextmanager
    def retrieving(self):
        def _f():
            with ftp_connection(self.server) as ftp:
                ftp.cwd(self.path)
                orig_request = ( yield )
                #print orig_request
                request = self.map_request(orig_request)
                while 1:
                    r = self.get_content(partial(ftp.retrbinary, 'RETR %s' % request ) )
                    orig_request = ( yield r )
                    request = self.map_request(orig_request)
        try:
            r = _f()
            r.next()
            yield r
        finally:
            r.close()

class FtpPdbSource(FtpDataSource):
    def __init__(self, server='ftp.wwpdb.org', path='/pub/pdb/data/structures/all/pdb', prefix='pdb', ext='.ent.gz',
                 user='', passwd='', acct=''):
        FtpDataSource.__init__(self, server, path, prefix, ext, user, passwd, acct)
        self.uncompress_cmd = "gzip -d"
    def get_content(self, tocall):
        self.sfh = StringIO()
        def _handle(data):
            self.sfh.write(data)
        tocall(_handle)
        p = subprocess.Popen("%s" % (self.uncompress_cmd), stdin=subprocess.PIPE, stdout=subprocess.PIPE, bufsize=0, shell=True)
        data = p.communicate( self.sfh.getvalue() )[0]
        p.wait()
        result = TmpPdbFile.fromText( data )
        return result
  
class LocalPdbSource(DirectoryDataSource):
    def __init__(self, dirname, prefix='', suffix='.gz', pack=False):
        DirectoryDataSource.__init__(self, dirname=dirname, prefix=prefix, suffix=suffix, pack=pack)
    def readData(self, fh):
        return TmpPdbFile.fromFile(fh)
    def writeData(self, fh, value):
        if not isinstance(value, PdbFile):
            raise TypeError, ''' Value must be a PDBFile object. '''
        return value.write(fh)

