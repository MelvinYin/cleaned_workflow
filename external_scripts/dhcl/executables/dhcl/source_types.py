from __future__ import with_statement
# TODO: Fix the filename lock so it works on FAT partitions under Linux
from dhcl.locking import FilenameLock

# Definition of DataSources
class STATUS:
    ''' Status of a resource - either MISSING, READY, WAIT '''
    MISSING=0
    READY=1
    OK=1
    WAITING=2
    ABORTED=4
    ERROR=8
    EMPTY=16

import sys
from contextlib import contextmanager, closing
class DataSourceError(Exception):  pass
class RetrieveError(DataSourceError): pass
class StoreError(DataSourceError): pass
class QueryError(DataSourceError): pass
class WaitingError(DataSourceError): pass
class AbortedError(DataSourceError): pass
class GoneWrongError(DataSourceError): pass
class EmptyError(DataSourceError): pass

class DataSource(object):
    def __init__(self):
        raise NotImplementedError('Abstract class can not be instantiated')
    def all(self):
        ''' Sequentially yields all possible requests objects that it is able to service. By default, no errors are iterated. '''
        raise NotImplementedError('Object can not formulate all possible requests')
    def retrieve(self, request):
        ''' Retrieves a single value from itself.
        @param request - a request (e.g. a PDB identifier) '''
        raise NotImplementedError('Object does not provide read access')
    def delete(self, request):
        ''' Deletes a single value from itself.
        @param request - a request (e.g. a PDB identifier) '''
        raise NotImplementedError('Object does not provide delete operation')
    def store(self, request, value):
        ''' Stores a single value in itself.
        @param request - a request (e.g. a PDB identifier)
        @param value - value to store under this request
        '''        
        raise NotImplementedError('Object does not provide store operation')
    def query(self, request):
        ''' Return a status code (see STATUS class for more info)
        @param request - a request (e.g. a PDB identifier)
        '''
        raise NotImplementedError('Object is not able to answer status queries')
    def retrieving(self):
        ''' A consumer/producer serving requests sequentially and yielding results from retrieve operation (if any) '''
        raise NotImplementedError('Object is not able to retrieve continuously.')
    def storing(self):
        ''' A consumer/producer storing values for requests sequentially and yielding results of storage operation (if any) '''
        raise NotImplementedError('Object is not able to store continuously.')
    def deleting(self):
        ''' A consumer/producer deleting values for corresponding requests sequentially and yielding results of delete operation (if any) '''
        raise NotImplementedError('Object is not able to delete continuously.')

class SingleInput(DataSource):
    ''' A data source defining each retrieve in separate transaction. Classes seeking to define retrieving in context of a repeated, single
    retrieve operation should use this. '''
    @contextmanager
    def retrieving(self):
        def _f():
            request = (yield)
            while 1:
                request = ( yield self.retrieve(request) )        
        try:
            r = _f()
            r.next()
            yield r
        finally:
            r.close()

class PackagedInput(DataSource):
    ''' A data source packaging retrieve calls into a single transaction block.
    Retrieving is a primitive and storing is defined in its terms. '''
    def retrieve(self, request):
        with self.retrieving() as retrieve:
            return retrieve.send(request)

class SingleOutput(DataSource):
    ''' A data source defining each store in separate transaction. Classes seeking to define storing in context of a repeated, single
    store operation should use this. '''
    @contextmanager
    def storing(self):
        def _f():
            request, value = (yield)
            while 1:
                request, value = ( yield self.store(request, value) )        
        try:
            r = _f()
            r.next()
            yield r
        finally:
            r.close()

class PackagedOutput(DataSource):
    ''' A data source packaging store calls into a single transaction block.
    Storing is a primitive and store is defined in its terms.
    '''
    def store(self, request, value):
        with closing_consumer(self.storing()) as store:
            return store.send(request, value)

class SingleDelete(DataSource):
    ''' A data source defining each retrieve in separate transaction. Classes seeking to define retrieving in context of a repeated, single
    retrieve operation should use this. '''
    @contextmanager
    def deleting(self):
        def _f():
            request = (yield)
            while 1:
                request = ( yield self.delete(request) )        
        try:
            r = _f()
            r.next()
            yield r
        finally:
            r.close()

class PackagedDelete(DataSource):
    ''' A data source packaging retrieve calls into a single transaction block.
    Retrieving is a primitive and storing is defined in its terms. '''
    def delete(self, request):
        with closing_consumer(self.deleting()) as delete:
            return delete.send(request)


# TODO: This earlier definition is needed only so that variables are correctly seen in the latter definition. Fix it.
class MappingSource:
    ANY = 0
    DELETE = 1
    READ = 2
    WRITE = 4
    QUERY = 8

class MappingSource(DataSource):
    ANY = 0
    DELETE = 1
    READ = 2
    WRITE = 4
    QUERY = 8
    ''' An abstract parent class for defining access to a dictionary(like) data source '''
    def writeData(self, fh, value):
        ''' Writes data correspnding to a mapped request filehandle'''
        raise NotImplementedError, '''Abstract method'''
    def readData(self, fh):
        ''' Reads data corresponding to an already mapped request filehandle'''
        raise NotImplementedError, '''Abstract method'''
    def deleteData(self, mapped):
        ''' Deletes data corresponding to an already mapped request'''
        raise NotImplementedError, '''Abstract method'''
    def queryData(self, mapped):
        ''' Queries data corresponding to an already mapped request'''
        raise NotImplementedError, '''Abstract method'''
    @contextmanager
    def mapRequest(self, request, reason=MappingSource.ANY):
        ''' Maps a request to a - must be implemented as a context manager ! Method is responsible for any initializing/cleaning up (e.g. writeData expect a filehandle) '''
        yield request
    def retrieve(self, request):
        with self.mapRequest(request, MappingSource.READ) as mapped:
            return self.readData(mapped)
    def store(self, request, value):
        with self.mapRequest(request, MappingSource.WRITE) as mapped:
            return self.writeData(mapped, value)
    def delete(self, request):
        with self.mapRequest(request, MappingSource.DELETE) as mapped:
            return self.deleteData(mapped)
    def query(self, request):
        with self.mapRequest(request, MappingSource.QUERY) as mapped:
            return self.queryData(mapped)

class TypedMappingSource(MappingSource):
    def __init__(self, dtype):
        self.dtype = dtype
    def store(self, request, value):
        if not isinstance(value, self.dtype):
            raise TypeError, 'Values is not of type: %s' % (str(self.dtype),)
        MappingSource.store(request, value)
        
import os, os.path, gzip, shelve, bz2
class DirectoryDataSource(SingleInput, SingleOutput, MappingSource):
    ''' Default implementation of accessing directory data sources maps to filenames using dirname, prefix and suffix.
    The files can optionally be packed using either gzip or bz2 compressions. Directory is locked (waiting lock) during all operations.
    '''
    def __init__(self, dirname, prefix='', suffix='', pack=False, create=True, clevel=9):
        self.dirname = os.path.abspath(dirname)
        self.prefix = prefix
        self.suffix = suffix
        self.pack = pack
        self.clevel = clevel
        self.lock = FilenameLock(self.dirname)
        with self.lock:
            if create and not os.path.exists(self.dirname):
                os.makedirs(self.dirname)
            if suffix=='.dbm':
                raise RuntimeError("DBM files presently cannot be stored")
            else:
                self.errordb = os.path.join(dirname, '.errors.dbm')
                shelve.open(self.errordb, 'c', 2)                
    @contextmanager
    def mapRequest(self, request, reason=MappingSource.ANY):
        with self.lock:
            fname = os.path.join(self.dirname, self.prefix + request + self.suffix)
            if reason==MappingSource.READ:
                if fname.endswith('.gz') or self.pack=='gzip' or self.pack==True:
                    yield gzip.open(fname, 'r', self.clevel)
                elif fname.endswith('.bz2') or self.pack=='bz2':
                    yield bz2.BZ2File(fname, 'r', compresslevel=self.clevel)
                else:
                    yield open(fname, 'r')    
            elif reason==MappingSource.WRITE:
                if fname.endswith('.gz') or self.pack=='gzip' or self.pack==True:
                    yield gzip.open(fname, 'w', self.clevel)
                elif fname.endswith('.bz2') or self.pack=='bz2':
                    yield bz2.BZ2File(fname, 'w', compresslevel=self.clevel)
                else:
                    yield open(fname, 'w')    
            else:        
                yield fname
    def all(self):
        #print os.listdir(self.dirname)[:10]
        with self.lock:
            return sorted( x[len(self.prefix):-len(self.suffix)] for x in os.listdir(self.dirname)
                           if x.startswith(self.prefix) and x.endswith(self.suffix) )
    def errors(self):
        with self.lock:
            try:
                db = shelve.open(self.errordb, 'c', 2)
                return sorted(db.keys())
            finally:
                db.close()            
    def queryData(self, mapped):
        if os.path.exists( mapped ):
            return STATUS.READY
        else:
            return STATUS.MISSING

    def retrieve(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if str(request) in error_db:
                    raise error_db[str(request)]
        return super(DirectoryDataSource, self).retrieve(request)
    def store(self, request, value):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if str(request) in error_db:
                    del error_db[str(request)]
                if isinstance(value, DataSourceError):
                    error_db[str(request)]=value
                    return
        return super(DirectoryDataSource, self).store(request, value)
    def query(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if str(request) in error_db:
                    error = error_db[str(request)]
                    if isinstance(error, WaitingError):
                        return STATUS.WAITING
                    elif isinstance(error, AbortedError):
                        return STATUS.ABORTED
                    elif isinstance(error, EmptyError):
                        return STATUS.EMPTY
                    else:
                        return STATUS.ERROR
        return super(DirectoryDataSource, self).query(request)
    def delete(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if str(request) in error_db:
                    del error_db[str(request)]                
        return super(DirectoryDataSource, self).delete(request)
    def deleteData(self, mapped):
        with self.lock:
            return os.unlink(mapped)


class ShelvedDataSource(SingleInput, SingleOutput, MappingSource):
    ''' A shelved data source. Encapsulates requests with a locking mechanism, allowing for multiple access.  '''
    def __init__(self, fname):
        self.fname = os.path.abspath(fname)
        self.lock = FilenameLock(self.fname)
        with self.lock:
            db = shelve.open(self.fname, 'c', 2)
            self.errordb = self.fname + '.errors.dbm'
            db = shelve.open(self.errordb, 'c', 2)
    def all(self):
        with self.lock:
            try:
                db = shelve.open(self.fname, 'c', 2)
                return sorted(db.keys())
            finally:
                db.close()
    def errors(self):
        with self.lock:
            try:
                db = shelve.open(self.errordb, 'c', 2)
                return sorted(db.keys())
            finally:
                db.close()            
    def queryData(self, mapped):
        with self.lock:
            db, key = mapped
            if key in db:
                return STATUS.READY
            else:
                return STATUS.MISSING
    def readData(self, mapped):
        with self.lock:
            try:
                db,key = mapped
                return db[key]
            except:            
                raise RetrieveError('No data exists')
    def writeData(self, mapped, value):
        with self.lock:
            db,key=mapped
            db[key]=value
    def deleteData(self, mapped):
        with self.lock:            
            db,key=mapped
            del db[key]
    @contextmanager
    def mapRequest(self, request, reason=MappingSource.ANY):
        with self.lock:
            db = shelve.open(self.fname, 'c', 2)
            try:
                yield db, request
            finally:
                db.close()
    def retrieve(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if request in error_db:
                    raise error_db[request]
        return super(ShelvedDataSource, self).retrieve(request)
    def store(self, request, value):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if request in error_db:
                    del error_db[request]
                if isinstance(value, DataSourceError):
                    error_db[request]=value
                    return
        return super(ShelvedDataSource, self).store(request, value)
    def query(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if request in error_db:
                    error = error_db[request]
                    if isinstance(error, WaitingError):
                        return STATUS.WAITING
                    elif isinstance(error, AbortedError):
                        return STATUS.ABORTED
                    elif isinstance(error, EmptyError):
                        return STATUS.EMPTY
                    else:
                        return STATUS.ERROR
        return super(ShelvedDataSource, self).query(request)
    def delete(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if request in error_db:
                    del error_db[request]
                else:
                    return super(ShelvedDataSource, self).delete(request)
                
class TypedShelvedDataSource(ShelvedDataSource):
    def __init__(self, fname, dtype):
        ShelvedDataSource.__init__(self, fname)
        self.dtype = dtype        
    def writeData(self, mapped, value):
        if not isinstance(value, self.dtype):
            raise TypeError, 'Values is not DataSourceError or of type: %s' % (str(self.dtype),)
        return ShelvedDataSource.writeData(self, mapped, value)

try:
    from pysqlite2 import dbapi2 as sqlite3
except:
    import sqlite3
@contextmanager
def transaction(conn):
    try:
        yield conn.cursor()
    except:
        conn.rollback()
        raise
    else:
        conn.commit()
@contextmanager
def connection(*args, **kwargs):
    try:
        conn = sqlite3.connect(*args, **kwargs)
        yield conn 
    finally:
        try:
            conn.close()
        except:
            pass

@contextmanager
def opening_transaction(*args, **kwargs):
    with connection(*args, **kwargs) as conn:
        with transaction(conn) as cursor:
            yield cursor

from cPickle import loads, dumps
from zlib import compress, decompress
class SqliteDataSource(SingleInput, SingleOutput, MappingSource):
    ''' An SQLite table data source. Encapsulates requests with a locking mechanism, allowing for multiple access.  '''
    def __init__(self, fname, table):
        self.fname = os.path.abspath(fname)
        self.lock = FilenameLock(self.fname)
        self.table = table
        with self.lock:
            with opening_transaction(self.fname) as cursor:
                cursor.execute('create table if not exists %s (uid TEXT primary key, value BLOB)' % (self.table,))
            self.errordb = self.fname + '.errors.dbm'
            db = shelve.open(self.errordb, 'c', 2)
    def all(self):
        with self.lock:
            with opening_transaction(self.fname) as cursor:
                cursor.execute('select uid from %s' % (self.table,) )
                keys = [ r[0] for r in cursor ]
                return sorted(keys)
    def errors(self):
        with self.lock:
            try:
                db = shelve.open(self.errordb, 'c', 2)
                return sorted(db.keys())
            finally:
                db.close()            
    def queryData(self, mapped):
        cursor, key = mapped
        cursor.execute('select uid from %s where uid=?' % self.table, (key,) )
        if key in [ r[0] for r in cursor ]:
            return STATUS.READY
        else:
            return STATUS.MISSING
    def readData(self, mapped):        
        try:
            cursor,key=mapped
            cursor.execute('select uid, value from %s where uid=?' % self.table, (key,) )
            for row in cursor:
                return loads(decompress(row[1]))
            raise RetrieveError('No data exists')
        except:
            raise
    def writeData(self, mapped, value):
        cursor,key=mapped
        # compress(dumps(value,2))
        cursor.execute('insert or replace into %s values (?,?)' % self.table , (key, sqlite3.Binary(compress(dumps(value,2)))) )
    def deleteData(self, mapped):
        cursor,key=mapped
        # compress(dumps(value,2))
        cursor.execute('delete from %s where uid=?' % self.table , (key, ) )        
    @contextmanager
    def mapRequest(self, request, reason=MappingSource.ANY):
        request = unicode(request)
        with self.lock:
            with opening_transaction(self.fname) as cursor:
                yield cursor, request                    
    def retrieve(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if str(request) in error_db:
                    raise error_db[request]
        return super(SqliteDataSource, self).retrieve(request)
    def store(self, request, value):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if str(request) in error_db:
                    del error_db[str(request)]
                if isinstance(value, DataSourceError):
                    error_db[request]=value
                    return
        return super(SqliteDataSource, self).store(request, value)
    def query(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if str(request) in error_db:
                    error = error_db[str(request)]
                    if isinstance(error, WaitingError):
                        return STATUS.WAITING
                    elif isinstance(error, AbortedError):
                        return STATUS.ABORTED
                    elif isinstance(error, EmptyError):
                        return STATUS.EMPTY
                    else:
                        return STATUS.ERROR
        return super(SqliteDataSource, self).query(request)
    def delete(self, request):
        with self.lock:
            with closing( shelve.open(self.errordb, 'c', 2) ) as error_db:
                if request in error_db:
                    del error_db[request]
        return super(SqliteDataSource, self).delete(request)

class TypedSqliteDataSource(SqliteDataSource):
    def __init__(self, fname, table, dtype):
        SqliteDataSource.__init__(self, fname, table)
        self.dtype = dtype        
    def writeData(self, mapped, value):
        if not isinstance(value, self.dtype):
            raise TypeError, 'Values is not DataSourceError or of type: %s' % (str(self.dtype),)
        return SqliteDataSource.writeData(self, mapped, value)
