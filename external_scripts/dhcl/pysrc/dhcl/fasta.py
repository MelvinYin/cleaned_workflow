#####################################################################################################################
##
## The module is part of SAPLing framework for processing structures& sequences of biological macromolecules.
## It may not be copied without explicit permission of the author.
## Copyright: Grzegorz M. Koczyk (2005-2007)
##
#####################################################################################################################
import re
###############################################################################################
## Extraction and manipulation of sequence identifiers
###############################################################################################
# Regex to match SwissProt identifier
__sprot_id = re.compile(r'^sp\|([A-Za-z0-9]+)\|([A-Za-z_0-9]+)')
# Regex to match identifier in customized PDB database
__gkpdb_id = re.compile(r'^([A-Za-z0-9]{4}_([A-Z\ ]))[\|]')
# Regex to match entries from NCBI sequence database
__gi_id = re.compile(r'^gi\|([0-9]+)')
# Regex to match 'db/id' type sequence identifiers
__slash_id = re.compile(r'^([^\s\|\/]+)[\/\|]([^\s\|]+)')
# Regex to match entries from lcl database
__lcl_id = re.compile(r'^([^\s\|]+)[\s\|]')

__space_regex = re.compile(r'\s+')
# Extraction of unique sequence identifier
def parseId(id):
    #print id
    match = __sprot_id.match(id)
    if match:
        return match.group(2)
    match = __gi_id.match(id)
    # Regex to match identifier in customized PDB database
    #print "ID", id
    match = __gkpdb_id.match(id) 
    if match:
        #print "MATCHED"
        return match.group(1)
    if match:
        return match.group(1)
    match = __slash_id.match(id)
    if match:    
        return match.group(2)
    match = __lcl_id.match(id)
    if match:
        return match.group(1)
    return id

class FastaSequence(object):
    def __init__(self, description='', text=''):
        self.description = description
        self._descriptions = description.split(chr(1)) # was '>'
        self.id = parseId(self.description)
        self.text = text
    def __len__(self):
        return len(self.text)
    def __iter__(self):
        ''' Iterate over the sequence text '''
        for s in self.text:
            yield s
    def at(self,position):
        ''' Deprecated method to be removed in 1.1'''
        return self.text[position-1]
    def __getitem__(self,position):
        return self.text[position]
    def __str__(self):
        return ">" + self.description + "\n" + self.text + "\n"
    def singleDescriptions(self):
        ''' Acess to the split up description for the sequence
        (useful if there are many copies with different annotations)'''
        return self._descriptions
    def getSeqId(self):
        return self.id

    seq_id = property(getSeqId)
    def getSeqDesc(self):
        return self.description
    def setSeqDesc(self, description):
        self.description = description
        self._descriptions = description.split(chr(1)) # was '>'
        self.id = parseId(self.description)
    seq_description = property(getSeqDesc, setSeqDesc)

    @classmethod
    @handlify(is_method=True)
    def readData(klass, datafile):
        ''' Can take either a filename, or a filelike object. Assumes a single FASTA sequence in file '''
        description, text, firstline = '', '', True
        for line in datafile:
            line = line.strip()
            if len(line):
                if line[0]=='>':
                    if firstline:
                        description = line[1:]
                    else:
                        raise RuntimeError, ''' Multiple/embedded definition lines found in FASTA file '''
                else:
                    text+=re.sub(r'\s+','',line)
            if line:
                firstline=False
        if not text:
            raise RuntimeError, ''' No data found - empty file '''
        return klass(description, text)

import subprocess
@handlify(is_method=False)
def readMultiFasta(fh):
    ''' Reads multiple FASTA format (multiple records, separated by line breaks)'''
    result = []
    for line in fh:
        line = line.strip()
        if len(line):
            if line[0]=='>':
                try:
                    result.append(FastaSequence(description, text))
                except NameError:
                    pass
                description = line[1:]
                text = ''
            else:
                text+=re.sub(r'\s+','',line)            
    return result

@handlify(is_method=False)
def xreadMultiFasta(fh):
    ''' Sequentially reads, yielding results. Requires multiple FASTA format (multiple records, separated by line breaks)'''
    result = []
    text = ''
    description = ''
    for line in fh:
        line = line.strip()
        if len(line):
            if line[0]=='>':            
                if len(text) or len(description):
                    yield FastaSequence(description, text)
                description = line[1:]
                text = ''
            else:
                text+= __space_regex.sub('',line)            
    yield FastaSequence(description, text)

