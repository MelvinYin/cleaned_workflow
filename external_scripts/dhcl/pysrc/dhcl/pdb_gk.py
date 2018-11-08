#####################################################################################################################
##
## The module is part of SiteProf2 framework, for processing structures& sequences of biological macromolecules.
## It may not be copied/changed without explicit permission of the author.
## Copyright: Grzegorz M. Koczyk (2005-2007)
##
#####################################################################################################################
# A module to handle retrieving of information in form of PDB files
from __future__ import with_statement
from dhcl.utils import provideTempdir, handlify, provideTempfile
import Bio.PDB
import os, os.path, sys, string, tempfile, subprocess, shutil
####################################################################################################################
# SELECTING PARTS OF A MODEL
class HetatmSelector(Bio.PDB.Select):
    def accept_residue(self, r):
        if r.get_id()[0].startswith('H_'):
            return True
        else:
            return False

class LigandSelector(Bio.PDB.Select):
    def __init__(self, *ligands):
        self.ligands = set(ligands)
    def accept_residue(self, r):
        if r.get_resname() in self.ligands:
            return True
        else:
            return False

class ChainSelector(Bio.PDB.Select):
    ''' Selects only residues belonging to a named chain. '''
    def __init__(self, *chains):
        self._chains = set(chains)
    def accept_residue(self, r):
        if r.get_parent().get_id() in self._chains:
            return True
        else:
            return False
class ResnameSelector(Bio.PDB.Select):
    ''' Selects only named residues. '''
    def __init__(self, *names):
        self._names = set(names)
    def accept_residue(self, r):
        if r.get_resname() in self._names:
            return True
        else:
            return False

class NotResnameSelector(Bio.PDB.Select):
    ''' Selects only named residues. '''
    def __init__(self, *names):
        self._names = set(names)
    def accept_residue(self, r):
        if r.get_resname() in self._names:
            return False
        else:
            return True

class OrSelector(Bio.PDB.Select):
    ''' Selects based on OR operation on component Selectors.'''
    def __init__(self, *args):
        self._subs = args
    def accept_residue(self, r):
        for s in self._subs:
            if s.accept_residue(r):
               return True
        return False
    def accept_chain(self, c):
        for s in self._subs:
            if s.accept_chain(c):
               return True
        return False
    def accept_model(self, m):
        for s in self._subs:
            if s.accept_model(m):
               return True
        return False
class AndSelector(Bio.PDB.Select):
    ''' Selects based on AND operation on component Selectors. '''
    def __init__(self, *args):
        self._subs = args
    def accept_residue(self, r):
        for s in self._subs:
            if not s.accept_residue(r):
               return False
        return True
    def accept_chain(self, c):
        for s in self._subs:
            if not s.accept_chain(c):
               return False
        return True           
    def accept_model(self, m):
        for s in self._subs:
            if not s.accept_model(m):
               return False
        return True

class NotSelector(Bio.PDB.Select):
    ''' Selects based on NOT operation on component Selectors. '''
    def __init__(self, *args):
        self._subs = args
    def accept_residue(self, r):
        for s in self._subs:
            if s.accept_residue(r):
               return False
        return True
    def accept_chain(self, c):
        for s in self._subs:
            if s.accept_chain(c):
               return False
        return True           
    def accept_model(self, m):
        for s in self._subs:
            if s.accept_model(m):
               return False
        return True

# Helper functions
from Bio.PDB import to_one_letter_code
def isHeteroResidue(residue):
    ''' Is a residue an heteroatom '''
    return residue.id[0]!=' ' #and not residue.resname in to_one_letter_code
def isAminoResidue(residue):
    ''' Is a residue an aminoacid '''
    return residue.id[0]==' ' and residue.resname in to_one_letter_code
def isWaterResidue(residue):
    ''' Is a residue water '''
    return residue.id[0]=='W'

def toOneLetter(residue):
    return to_one_letter_code.get(residue.resname, 'X')
class AminoSelector(Bio.PDB.Select):
    def accept_residue(self, r):
        return isAminoResidue(r)

####################################################################################################################
# PDB data structures usable to PyMol/MODELLER programs, as well as through BioPython
# Each such structure has an explicitly/implicitly associated file.
# An object representing a PDB file, suitable for both MODELLER and BioPython use
# Its structure can be modified, however
class PdbFile(object):
    def __init__(self, fname, chain=None):
        self.fname = os.path.abspath(fname)
        #print self.fname
        if not os.path.exists(self.fname):
            raise IOError(''' File %s does not exist ''' % fname)
    
    @classmethod
    def makeEmpty(klass, fname):
        self = klass.__new__(klass)
        self.fname = os.path.abspath(fname)
        return self
    @classmethod
    def fromStruct(klass, fname, struct, selector=None):
        if selector:
            io=Bio.PDB.PDBIO()
            io.set_structure(struct)
            io.save(fname, selector)        
            self = klass(fname)
        else:
            self = klass.makeEmpty(fname)
            self._struct = struct
            self.save()
        return self

    @classmethod
    def fromText(klass, fname, text):
        self = klass.makeEmpty(fname)
        wfh = open(self.fname, 'w')
        wfh.write(text)
        wfh.close()
        return self

    @classmethod
    def fromFile(klass, fname, source_fname):
        self = klass.makeEmpty(fname)
        if isinstance(source_fname, basestring):
            shutil.copy(source_fname, self.fname)
        else:
            with open(fname, 'w') as wfh:
                wfh.write(source_fname.read())                
        return self
    
    def __str__(self):
        return self.fname
    def getStructure(self):
        if not hasattr(self, '_struct'):
            p = Bio.PDB.PDBParser()
            self._struct = p.get_structure('protein', self.fname)
        return self._struct            
    struct = property(getStructure)
    def getFh(self):
        return open(self.fname, 'r')
    fh = property(getFh)

    def save(self):
        io=Bio.PDB.PDBIO()
        io.set_structure(self.struct)
        io.save(self.fname)        
        
    @handlify(is_method=True, mode='w')
    def write(self, fh):
        rfh = self.fh
        #print self, fh
        for line in rfh:
            fh.write(line)
        rfh.close()

import new
import time, os
# An object representing a PDB file, suitable for both MODELLER and BioPython use
# Optionally a selector is applied
class TmpPdbFile(PdbFile):
    def __init__(self, pdb_id, selector=None):
        raise RuntimeError, 'Abstract method'
        fh, self.fname = tempfile.mkstemp(suffix='.pdb')
        os.close(fh)
        self.pdb_id = pdb_id
        if selector is None:
            getPdb(pdb_id, self.fname)
        else:
            getSelectedPdb( pdb_id, self.fname, selector )

    @classmethod
    def makeEmpty(klass):
        self = klass.__new__(klass)
        fh, self.fname = tempfile.mkstemp('.pdb', 'tmp', os.getcwd())
        os.close(fh)
        return self
    @classmethod
    def fromStruct(klass, struct, selector=None):
        self = klass.makeEmpty()
        if selector:       
            io=Bio.PDB.PDBIO()            
            io.set_structure(struct)
            io.save(self.fname, selector)                
        else:
            self._struct = struct
            self.save()
        return self

    @classmethod
    def fromFile(klass, fname):
        self = klass.__new__(klass)
        fh, self.fname = tempfile.mkstemp('.pdb', 'tmp', os.getcwd())
        os.close(fh)
        # Filename
        if isinstance(fname, basestring):
            shutil.copy(fname, self.fname)
        # Filehandle
        else:
            with open(self.fname, 'w') as wfh:
                wfh.write(fname.read())                
        return self

    @classmethod
    def fromText(klass, text):
        self = klass.makeEmpty()
        wfh = open(self.fname, 'w')
        wfh.write(text)
        wfh.close()
        return self
    
    def finalize(self):
        if hasattr(self, 'fname'):
            try:
                #time.sleep(1)
                os.unlink(self.fname)
            except:
                pass        
    def __del__(self):
        self.finalize()

import csv
from dhcl.utils import handlify
class ResIndex(dict):
    def __init__(self, seq):
        if seq:
            for i, (v1, v2) in enumerate(seq):
                self[i+1] = (v1, v2)
    @classmethod
    @handlify(is_method=True)
    def readData(klass, fh):
        self = klass.__new__()
        reader = csv.reader(fh, dialect='excel-tab')
        for row in reader:
            self[ int(row[0]) ] = (row[1], row[2])
        return self
    @handlify(is_method=True, mode='r')
    def writeData(self, fh):
        writer = csv.writer(fh, dialect='excel-tab')
        for k, (v1, v2) in sorted( self.iteritems(), key=itemgetter(0) ):
            writer.writerow([k, v1, v2])
            
from operator import isNumberType, itemgetter
class PdbSeqRecord(dict):
    def getFlatHeader(self):
        r = {}
        for k,v in sorted(self.header.iteritems(), key=itemgetter(0)):
            if isinstance(v, dict):
                for k1,v1 in sorted(v.iteritems(), key=itemgetter(0)):
                    if isinstance(v1, dict):
                        for k2,v2 in sorted(v1.iteritems(), key=itemgetter(0)):
                            r[k2.upper()] = v2 if isNumberType(v2) else str(v2).upper().strip()
                    else:
                        r[k1.upper()] = v1 if isNumberType(v1) else str(v1).upper().strip()
            else:
                r[k.upper()] = v if isNumberType(v) else str(v).upper().strip()
        return r
    flat_header = property(getFlatHeader)

