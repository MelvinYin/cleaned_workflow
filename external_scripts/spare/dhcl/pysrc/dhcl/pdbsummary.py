from __future__ import with_statement, division
import sys, os
from dhcl.utils import StructObject, handlify, innerDefaultdict
from collections import defaultdict

# Splitter functions for parsing my extracts from PDB headers
def splitter(s, delim=',', quoter='"'):
    v = ""
    in_quotes = False
    for l in s:
        #print l
        if l==quoter:
            in_quotes = False if in_quotes else True
        if l==delim and not in_quotes:
            yield v
            v = ""
        else:
            v+=l
    if v.strip():
        yield v

def isplitter(s):
    return s[0], s[2:]

# Parse the PDB annotations file
@handlify()
def parsePdbDataFile(fh):
    r = dict()
    for line in fh:
        if 'ERROR PARSING' in line:
            continue
        line = line.split('\t')
        #print line
        while len(line)>9:
            line = line[:3] + [ line[3]+line[4] ] + line[5:]
        pdb_id, moltype, method, description, resolution, sum_len, lengths, organisms, moldescs = line
        sum_len = int(sum_len)
        if resolution:
            resolution = float(resolution) if resolution and resolution!='None' else None
        else:
            resolution = None
        lengths =  dict( (k, int(v)) for k,v in ( isplitter(x) for x in lengths.split(',') ) )
        organisms = organisms
        #print moldescs        
        #for x in splitter( moldescs ):
        #    print "X1", x
        moldescs =  dict( (k, v.replace('"','')) for k,v in ( x.split(':') for x in splitter( moldescs ) ) )
        v = StructObject(pdb_id=pdb_id, moltype=moltype, method=method, description=description, sum_len=sum_len,
                         lengths=lengths, organisms=organisms, moldescs=moldescs, resolution=resolution )
        r[v.pdb_id] = v
    return r


# Parse the OGT-species associations from Alex's database
@handlify()
def parseOrganismOgts(fh):
    r = list()
    for line in fh:
        if line.startswith('organism'):
            continue
        #print line
        line = [ x.strip() for x in line.split('\t') ]
        organism, superkingdom, division, lineage, ogt, gb_accession, accession, loaded = line
        ogt = int(ogt)
        superkingdom = "Archea" if superkingdom=='A' else 'Bacteria'
        #print organism
        ognames = organism.split(' ')
        genus = ognames[0]
        species = ognames[1]
        species = "%s %s" % (genus, species)
        # Find the temperature class, for divisions
        if ogt>80:
            tempklass = "hyperthermophile"
        elif ogt>50:
            tempklass = "thermophile"
        elif ogt>23:
            tempklass = "mesophile"
        else:
            ogt = "psychrophile"
        v = StructObject(tempklass=tempklass, organism=organism, superkingdom=superkingdom, genus = genus, species = species,
                         lineage=lineage, ogt=ogt, gb_accession=gb_accession, accession=accession )
        r.append(v)
    return r

@handlify()
def parsePdbClusters(fh):
    r_clusts = r1 = defaultdict( innerDefaultdict(set) )
    r_pdb_ids = r2 = defaultdict(dict)
    for line in fh:
        clust_no, seq_no, pdb_id = line.strip().split()
        clust_no = int(clust_no)
        pdb_id,chain = pdb_id.split(':')
        pdb_id = pdb_id.lower()
        r1[clust_no][pdb_id].add(chain)
        r2[pdb_id][chain]=clust_no
    return r1, r2
