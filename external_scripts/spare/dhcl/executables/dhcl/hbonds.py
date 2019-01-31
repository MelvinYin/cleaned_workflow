from __future__ import with_statement
import math, sys, os, os.path
import Bio.PDB
#import Numeric
from dhcl.utils import *
from dhcl.pdb_data import prepareHbonds
from collections import defaultdict
# Output DEBUG msgs
DEBUG = 0
#DEBUG = 0

# def handlify(mode='r', close=False, is_method=False, opener_f = file):
#     ''' Decorator on a parse/read function.method, so that if a string is supplied it will attempt to open a file
#     and use it instead '''
#     # Intermediate function is needed so that decorator creation can be customised, by options above
#     def intermediate(func):
#         if is_method:
#             def result(self, fh, *args, **kwargs):
#                 close_ = close
#                 if isinstance(fh, basestring):
#                     fh = opener_f(fh, mode)
#                     close_ = True
#                 try:                    
#                     r = func(self, fh, *args, **kwargs)
#                 finally:
#                     if close_:
#                         fh.close()
#                 return r
#             return result
#         else:                
#             def result(fh, *args, **kwargs):
#                 close_ = close
#                 if isinstance(fh, basestring):
#                     fh = opener_f(fh, mode)
#                     close_ = True
#                 try:
#                     r = func(fh, *args, **kwargs)
#                 finally:
#                     if close_:
#                         fh.close()
#                 return r
#         return result
#     return intermediate

# # Dummy object corresponding to C struct construct
# class StructObject(object):
#     def __init__(self, **dictargs):    
#         for k, v in dictargs.iteritems():
#             object.__setattr__(self,k, v)
            
#     def __setattr__(self, k, v):
#         object.__setattr__(self, k, v)
#     def __str__(self):
#         return "( " + ", ".join(["%s: %s" % (k, self.__dict__[k]) for k in sorted(self.__dict__)]) + " )";
#     def __repr__(self):
#         return self.__str__()


# def innerDefaultdict(t=None):
#     def __inner():
#         return defaultdict(t)
#     return __inner
    
# # CLASSES CONTAINING DATA ABOUT ATOM TYPES FOR DETERMINATION OF HYDROGEN BONDS
# class AtomType(object):
#     """ Atom type defines following properties:
#     - atom_name - name of the atom (acceptor or donor depending on subclass)
#     - hybridization - hybridization of the atom (Nsp2, Nsp3, Osp2, Osp3, Ssp3)
#     - max_bonds - maximum number of bonds formed by this atom in its given role
#     - excludes - whether the atom excludes any other atom from forming a bond (default: None)
#     """
#     def __init__(self, atom_name, hybridization, radius, max_bonds = 1, excludes=None):
#         self.atom_name = atom_name
#         self.thybridization = hybridization
#         self.hybridization = hybridization[1:]
#         self.radius = radius
#         self.max_bonds = max_bonds
#         self.excludes = excludes

# class AtomDonor(AtomType):
#     """ Donor type defined additional properties:
#     - d - donor atom name - same as atom_name
#     - dd - antecedent atom
#     - dd1 - postcedent atom (for defining planes - sometimes it is a different atom, in case of terminal groups)
#     """
#     def __init__(self, atom_name, dd, dd1, hybridization, radius, max_bonds = 1, excludes=None, hydrogen=None):
#         super(AtomDonor, self).__init__(atom_name, hybridization, radius, max_bonds, excludes)
#         self.dd = dd
#         self.dd1 = dd1
#         self.hydrogen = hydrogen
#     def getD(self):
#         return self.atom_name
#     d = property(getD)

# class AtomAcceptor(AtomType):
#     """ Donor type defined additional properties:
#     - a - acceptor atom name - same as atom_name
#     - aa - antecedent atom
#     - aa2 - postcedent atom (for defining planes - sometimes it is a different atom, in case of terminal groups)
#     """
#     def __init__(self, atom_name, aa, aa2, hybridization, radius, max_bonds = 1, excludes=None):
#         super(AtomAcceptor, self).__init__(atom_name, hybridization, radius, max_bonds, excludes)
#         self.aa = aa
#         self.aa2 = aa2
#     def getA(self):
#         return self.atom_name
#     a = property(getA)


# peptide_donor = AtomDonor( atom_name="N", dd="CA", dd1="C_prev", hybridization="Nsp2", radius=1.9)
# peptide_acceptor = AtomAcceptor( atom_name="O", aa="C", aa2="CA", hybridization="Osp2", radius=1.6)

# def atomdikt(*it):
#     return dict( (x.atom_name, x) for x in it )
# # Donors according to (Stickle et al. 1992)
# donors = {
#     'ALA': atomdikt( peptide_donor ),
#     'ARG': atomdikt( peptide_donor,
#                      AtomDonor(atom_name='NE', dd='CD', dd1='CZ', hybridization='Nsp2', radius=1.9, max_bonds=1),
#                      AtomDonor(atom_name='NH1', dd='CZ', dd1='NE', hybridization='Nsp2', radius=1.9, max_bonds=2),
#                      AtomDonor(atom_name='NH2', dd='CZ', dd1='NE', hybridization='Nsp2', radius=1.9, max_bonds=2)
#                      ),
#     'ASP': atomdikt( peptide_donor ),
#     'ASN': atomdikt( peptide_donor,
#                      AtomDonor(atom_name='ND2', dd='CG', dd1='CB', hybridization='Nsp2', radius=1.9, max_bonds=2),
#                      ),
#     'CYS': atomdikt( peptide_donor ),
#     'GLU': atomdikt( peptide_donor ),
#     'GLN': atomdikt( peptide_donor,
#                      AtomDonor(atom_name='NE2', dd='CD', dd1='CG', hybridization='Nsp2', radius=1.9, max_bonds=2),
#                      ),
#     'GLY': atomdikt( peptide_donor ),
#     'HIS': atomdikt( peptide_donor,
#                      AtomDonor(atom_name='ND1', dd='CG', dd1='CE1', hybridization='Nsp2', radius=1.9, max_bonds=1, excludes='NE2'),
#                      AtomDonor(atom_name='NE2', dd='CE1', dd1='CD2', hybridization='Nsp2', radius=1.9, max_bonds=1, excludes='ND1'),
#                      ),
#     'ILE': atomdikt( peptide_donor ),
#     'LEU': atomdikt( peptide_donor ),
#     'LYS': atomdikt( peptide_donor,
#                      AtomDonor(atom_name='NZ', dd='CE', dd1='CD', hybridization='Nsp3', radius=2.1, max_bonds=3),
#                      ),
#     'MET': atomdikt( peptide_donor ),
#     'PHE': atomdikt( peptide_donor ),
#     'PRO': atomdikt( ),
#     'SER': atomdikt( peptide_donor,
#                      AtomDonor(atom_name='OG', dd='CB', dd1='CA', hybridization='Osp3', radius=1.7, max_bonds=1),
#                      ),
#     'THR': atomdikt( peptide_donor,
#                      # TODO: Check radius=1.9 with Igor - whether it is on purpose (article says r=1.7)
#                      AtomDonor(atom_name='OG1', dd='CB', dd1='CA', hybridization='Osp3', radius=1.7, max_bonds=1),
#                      ),
#     'TRP': atomdikt( peptide_donor,
#                      AtomDonor(atom_name='NE1', dd='CD1', dd1='CE2', hybridization='Nsp2', radius=1.9, max_bonds=1),
#                     ) ,
#     'TYR': atomdikt( peptide_donor,
#                      AtomDonor(atom_name='OH', dd='CZ', dd1='CE1', hybridization='Osp3', radius=1.7, max_bonds=1),
#                      ),
#     'VAL': atomdikt( peptide_donor ),
#     }

# acceptors = {
#     'ALA': atomdikt( peptide_acceptor ),
#     'ARG': atomdikt( peptide_acceptor ),
#     'ASP': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='OD1', aa='CG', aa2='CB', hybridization='Osp2', radius=1.6, max_bonds=1),
#                      AtomAcceptor(atom_name='OD2', aa='CG', aa2='CB', hybridization='Osp2', radius=1.6, max_bonds=1),
#                      ),
#     'ASN': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='OD1', aa='CG', aa2='CB', hybridization='Osp2', radius=1.6, max_bonds=2),
#                      ),
#     'CYS': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='SG', aa='CB', aa2='CA', hybridization='Ssp3', radius=2.1, max_bonds=2),
#                      ),
#     'GLU': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='OE1', aa='CD', aa2='CG', hybridization='Osp2', radius=1.6, max_bonds=2),
#                      AtomAcceptor(atom_name='OE2', aa='CD', aa2='CG', hybridization='Osp2', radius=1.6, max_bonds=2),
#                      ),
#     'GLN': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='OE1', aa='CD', aa2='CG', hybridization='Nsp2', radius=1.6, max_bonds=2),
#                      ),
#     'GLY': atomdikt( peptide_acceptor ),
#     'HIS': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='ND1', aa='CG', aa2='CE1', hybridization='Nsp2', radius=1.6, max_bonds=1, excludes='NE2'),
#                      AtomAcceptor(atom_name='NE2', aa='CE1', aa2='CD2', hybridization='Nsp2', radius=1.6, max_bonds=1, excludes='ND1'),
#                      ),
#     'ILE': atomdikt( peptide_acceptor ),
#     'LEU': atomdikt( peptide_acceptor ),
#     'LYS': atomdikt( peptide_acceptor ),
#     'MET': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='SD', aa='CG', aa2='CB', hybridization='Ssp3', radius=1.95, max_bonds=2),
#                      ),
#     'PHE': atomdikt( peptide_acceptor ),
#     'PRO': atomdikt( peptide_acceptor ),
#     'SER': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='OG', aa='CB', aa2='CA', hybridization='Osp3', radius=1.7, max_bonds=1),
#                      ),
#     'THR': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='OG1', aa='CB', aa2='CA', hybridization='Osp3', radius=1.7, max_bonds=1),
#                      ),
#     'TRP': atomdikt( peptide_acceptor ),
#     'TYR': atomdikt( peptide_acceptor,
#                      AtomAcceptor(atom_name='OH', aa='CZ', aa2='CE1', hybridization='Osp3', radius=1.7, max_bonds=1),
#                      ),
#     'VAL': atomdikt( peptide_acceptor ),
#     }

# def norm(v):
#     ''' Norm(length) of a vector '''
#     return math.sqrt( sum(v_*v_ for v_ in v) )

# def cross(v1, v2):
#     ''' Cross product of two vectors '''
#     x1, x2, x3 = v1
#     y1, y2, y3 = v2
#     return Numeric.array([
#         x2*y3 - x3*y2,
#         x3*y1 - x1*y3,
#         x1*y2 - x2*y1,
#         ])

# def dot(v1,v2):
#     return sum( v1_*v2_ for v1_,v2_ in zip(v1,v2) )

# def sgn(x):
#     if x<0:
#         return -1
#     else:
#         return 1
    
# def angle(v1, v2, abs=True):
#     m = sum(v1*v2)/(norm(v1)*norm(v2))
#     if m>1.:
#         m=1.0
#     if m<-1:
#         m=-1.0    
#     va = math.acos(m)
#     #if s<0:
#     #    va+=math.pi
#     #else:
#     #    va-=math.pi
#     if va<0:
#         va=math.pi - va
#     #if not abs:
#     #    if -va:
#     #        
#     #    s = sgn(dot(v1,v2))
#     #    #if va>math.pi/2:
#     #    #    va= -(math.pi-va)
#     #    va*=s        
#     return va

# def linkResidues(s):
#     for m in s:
#         for c in m:
#             last_r = None
#             for r in c:
#                 if DEBUG:
#                     print "## %s %s" % (c.id, r.id)
#                 r.previous = None
#                 if last_r is not None:
#                     if r.id[1]>last_r.id[1]+1:
#                         pass
#                     else:
#                         last_r.next = r
#                         r.previous = last_r
#                 r.next = None
#                 last_r = r
                
# def todeg(angle):
#     try:
#         return float("%.2f" % (angle*180 / math.pi,))
#     except:
#         return None

# def checkAtomPair(donor, donor_atom, acceptor, acceptor_atom):
#     distance = (donor_atom - acceptor_atom)[0]
#     # Distance condition must be satisfied
#     if distance >= donor.radius + acceptor.radius:
#         if DEBUG:
#             print "REJECTED %s:%s %s:%s r=%s" % (donor_atom.get_parent().id, donor_atom.id, acceptor_atom.get_parent().id, acceptor_atom.id, distance)
#         return None
#     donor_residue = donor_atom.get_parent()
#     acceptor_residue = acceptor_atom.get_parent()
#     if DEBUG:
#         print "%s:%s %s:%s r=%s" % (donor_atom.get_parent().id, donor_atom.id, acceptor_atom.get_parent().id, acceptor_atom.id, distance)
#         print "%s <-> %s" % (donor.hybridization, acceptor.hybridization)
#     if DEBUG:
#         print "Distance accepted"
#     # Defining atoms for angles, if any can not be found - abort
#     # DD
#     try:
#         dd = donor_residue[donor.dd]
#     except:
#         return None
#     if DEBUG:
#         print "DD: %s" % donor.dd
#     if donor.dd1=='C_prev':
#         try:
#             if DEBUG:
#                 print "DD - C_prev %s" % donor_residue.previous
#             dd1 = donor_residue.previous['C']
#         except:
#             dd1 = None
#             #return None
#     else:
#         try:
#             dd1 = donor_residue[donor.dd1]
#         except:
#             return None
#     if DEBUG:
#         print "DD1: %s" % donor.dd1
#     # AA
#     try:
#         aa = acceptor_residue[acceptor.aa]
#     except:
#         return None
#     # AA2
#     if acceptor.aa2=='N_next':
#         try:
#             aa2 = donor_residue.next['N']
#             # Hax to fix angles 
#             #aa= acceptor_residue['C']
#         except:
#             return None
#     else:
#         try:
#             aa2 = acceptor_residue[acceptor.aa2]
#         except:
#             return None
#     a = acceptor_atom
#     d = donor_atom
#     #########################################################
#     # Checking scalar angles
#     #
#     # 1) D-A-AA angle in (60-180 for sp3, 90-180 for sp2)
#     v1 = aa.coord - a.coord
#     v2 = d.coord - a.coord
#     d_a_aa_angle = angle(v1, v2)
#     if acceptor.hybridization.endswith('sp3'):
#         #if donor_residue.resname == 'LYS':
#         #    if d_a_aa_angle<=0:
#         #        if DEBUG:
#         #            print "Rejected due to D-A-AA %s>60" % todeg(d_a_aa_angle)
#         #        return None
#         if d_a_aa_angle <= math.pi/3:
#             if DEBUG:
#                 print "Rejected due to D-A-AA %s>60" % todeg(d_a_aa_angle)
#             return None
#     else:
#         if d_a_aa_angle <= math.pi/2:
#             if DEBUG:
#                 print "Rejected due to D-A-AA %s>90 (sp2)" % ( todeg(d_a_aa_angle), )
#             return None
#     if DEBUG:
#         print "Accepted due to D-A-AA %s>60" % todeg(d_a_aa_angle)
#     #
#     # 2) A-D-DD angle in (90-180 for both)
#     v1 = a.coord - d.coord
#     v2 = dd.coord - d.coord
#     if a is dd:
#         return None
#     a_d_dd_angle =  angle(v1, v2)
#     if (donor_residue.resname=='LYS' and not donor is peptide_donor and a_d_dd_angle<=math.pi/3) or (a_d_dd_angle <= math.pi/2):
#         if DEBUG:
#                 print "Rejected due to A-D-DD %s>90 (sp2)" % todeg(a_d_dd_angle)       
#         return None
#     if DEBUG:
#         print "CHECKING PLANARS..."
#     #########################################################
#     # Checking planar angles, if at least one partner is sp2
#     planar1=planar2=a_d_dd=None
#     if True or donor.hybridization.endswith('sp2') or acceptor.hybridization.endswith('sp2'):
#         # 1) Checking acceptor in donor plane (D-DD-DD1 vs. A-D-DD)
#         if dd1 is None:
#             d_dd_dd1 = None
#         else:
#             d_dd_dd1 = cross( d.coord - dd.coord, dd1.coord - dd.coord )
#         a_d_dd = cross( a.coord-d.coord, dd.coord-d.coord )
#         if d_dd_dd1 is None:
#             planar1 = 0.
#         else:
#             planar1 =  angle( d_dd_dd1, -a_d_dd, abs=False)
        
#         # Special case of LYS - 90 degrees of latitude        
#         if donor.hybridization.endswith('sp2'):
#             if donor_residue.resname == 'LYS' and not donor is peptide_donor:
#                 if planar1>math.pi/2 and -planar1<=-math.pi/2:
#                     if DEBUG:
#                         print "Rejected due to (planar1; LYS) A-D< %s" % todeg(planar1)
#                     return None
#             elif planar1>=math.pi/3 and -planar1>=-2*math.pi/3:
#                 if DEBUG:
#                     print "Rejected due to (planar1) A-D< %s" % todeg(planar1)
#                 return None
#             #if donor_residue.resname == 'LYS':
#             #    if planar1 == math.pi/2 or planar1 == -math.pi/2:
#             #        if DEBUG:
#             #            print "Rejected due to (planar1) A-D< %s" % todeg(planar1)
#             #        return None
#             # Else 60
#             #elif (planar1<=math.pi and planar1 >= math.pi/3) or (planar1>=-math.pi and planar1 <= -math.pi/3) or (planar1>=math.pi and planar1 <= math.pi+math.pi/3) or (planar1<=-math.pi and planar1 >= -math.pi-math.pi/3):
#              #   if DEBUG:
#              #       print "Rejected due to (planar1) A-D< %s" % todeg(planar1)
#              #   return None
#         # 2) Checking donor in acceptor plane (A-AA-AA2 vs. D-A-AA)
#         a_aa_aa2 = cross( a.coord - aa.coord, aa2.coord - aa.coord)
#         d_a_aa = cross(d.coord - a.coord, aa.coord - a.coord)
#         planar2 = angle(a_aa_aa2, d_a_aa, abs=False) 
#         # This is more lax - 90 degs always
#         if acceptor.hybridization.endswith('sp2'):
#             if planar2 >= math.pi/2 and -planar2 >= -math.pi/2:
#                 if DEBUG:
#                     print "Rejected due to (planar2) D-A< %s" % todeg(planar2)
#                 return None
#     if DEBUG:
#         print "ACCEPTED"
#     # All tests passed - we have a valid donor-acceptor pair !!!
#     return StructObject(donor=donor, donor_atom=donor_atom, acceptor=acceptor, acceptor_atom=acceptor_atom, distance=distance, planar1=todeg(planar1), planar2=todeg(planar2), a_d_dd=todeg(a_d_dd_angle), a_d_dd2=todeg(a_d_dd), d_a_aa=todeg(d_a_aa_angle))

# class Clist(list):
#     pass

# def checkResiduePair(residue1, residue2):
#     ''' Outputs a list of possible donor-acceptor pairs for this pair of residues. '''
#     if DEBUG:
#         print "CONSIDERING###", residue1.id[1], residue1.resname, residue2.id[1], residue2.resname
#     result = Clist()
#     acceptors1 = acceptors.get(residue1.resname, {})
#     acceptors2 = acceptors.get(residue2.resname, {})
#     donors1 = donors.get(residue1.resname, {})
#     donors2 = donors.get(residue2.resname, {})
#     # Donor-acceptor - first pairing
#     for donor in donors1.values():
#         for acceptor in acceptors2.values():
#             try:
#                 donor_atom = residue1[ donor.atom_name ]
#                 acceptor_atom = residue2[ acceptor.atom_name ]
#             except:
#                 continue
#             if donor_atom is acceptor_atom:
#                 continue
#             x = checkAtomPair(donor, donor_atom, acceptor, acceptor_atom)
#             if x is not None:
#                 x.type = ("N" if donor is peptide_donor else "R") + "_" + ("O" if acceptor is peptide_acceptor else "R")
#                 x.donor_residue = residue1
#                 x.acceptor_residue = residue2
#                 result.append(x)
#                 if DEBUG:
#                     print x
#     # Donor-acceptor - second pairing
#     for donor in donors2.values():
#         for acceptor in acceptors1.values():
#             try:
#                 donor_atom = residue2[ donor.atom_name ]
#                 acceptor_atom = residue1[ acceptor.atom_name ]
#             except:
#                 continue
#             if donor_atom is acceptor_atom:
#                 continue
#             x = checkAtomPair(donor, donor_atom, acceptor, acceptor_atom)
#             if x is not None:
#                 x.type = ("N" if donor is peptide_donor else "R") + "_" + ("O" if acceptor is peptide_acceptor else "R")
#                 x.donor_residue = residue2
#                 x.acceptor_residue = residue1
#                 result.append(x)
#                 if DEBUG:
#                     print x
#     return result

# def ordHbonds(hbs):
#     r = defaultdict( innerDefaultdict( innerDefaultdict( innerDefaultdict( innerDefaultdict( dict ) ) ) ) )
#     for hbond in hbs:
#         acc_atom = hbond.acceptor_atom
#         don_atom = hbond.donor_atom
#         dist = hbond.distance
#         res1 = hbond.donor_residue
#         res2 = hbond.acceptor_residue
#         icode1 = res1.id[2]
#         icode2 = res2.id[2]
#         chain1 = res1.get_parent().id
#         chain2 = res2.get_parent().id
#         x = StructObject(resname1 = res1.resname, chain1=chain1, chain2=chain2, resno1 = res1.id[1], atom1 = don_atom.id, icode1 = icode1, icode2 = icode2,
#                          resname2 = res2.resname, resno2 = res2.id[1], atom2 = acc_atom.id, distance = dist,
#                          d_a_aa = hbond.d_a_aa, a_d_dd = hbond.a_d_dd,
#                          planar1 = hbond.planar1, planar2 = hbond.planar2, type=hbond.type)
#         r [ chain1 ][ (x.resno1, x.icode1) ][ don_atom.id ][ chain2 ][ (x.resno2, x.icode2) ][ acc_atom.id ] = x    
#     return r

# def old_hbonds(pdb_fname, chain=None):
#     result = []
#     p = Bio.PDB.PDBParser()
#     s = p.get_structure('query', pdb_fname)
#     linkResidues(s)
#     m = s[0]
#     if chain:
#         #print m.__dict__
#         chain = sorted(m.child_dict.keys())[0]
#     seen = set()
#     for chain1 in m:
#         if chain is not None and chain1!=chain:
#             continue
#         for residue1 in chain1:            
#             if DEBUG:
#                 print "#" *100
#                 print "Considering %s" % ( residue1.id, )
#             for chain2 in m:
#                 if chain is not None and chain2!=chain:
#                     continue
#                 for residue2 in chain2:
#                     if residue2 not in seen:
#                         rp =  checkResiduePair(residue1, residue2)
#                         if rp:
#                             result.extend( rp )
#             seen.add(residue1)
#     return ordHbonds(result)

    
# def tosdeg(a):
#     if a is None:
#         return "None"
#     else:
#         return "%.2f" % (a*180/math.pi)

import csv
@handlify(mode='w')
def writeHbonds(fh, hbs, pdb_id='query'):
    writer = csv.writer(fh, dialect='excel-tab')
    for chain1 in sorted(hbs):
        hbs_ch1 = hbs[chain1]
        for resno1, icode1 in sorted(hbs_ch1):
            hbs_ch1_r1 = hbs_ch1[ (resno1, icode1) ]
            for donor_atom in sorted(hbs_ch1_r1):
                hbs_ch1_r1_a1 = hbs_ch1_r1[donor_atom]
                for chain2 in sorted(hbs_ch1_r1_a1):
                    hbs_ch1_r1_a1_ch2 = hbs_ch1_r1_a1[chain2]
                    for resno2, icode2 in hbs_ch1_r1_a1_ch2:
                        hbs_ch1_r1_a1_ch2_r2 = hbs_ch1_r1_a1_ch2[(resno2,icode2)]
                        for acceptor_atom in sorted(hbs_ch1_r1_a1_ch2_r2):
                            hb = hbs_ch1_r1_a1_ch2_r2[ acceptor_atom ]                            
                            writer.writerow([pdb_id, chain1, hb.resname1, resno1, icode1, donor_atom, chain2, hb.resname2, resno2, icode2, acceptor_atom, "%.2f" % hb.distance,
                                             "%.2f" % hb.a_d_dd, "%.2f" % hb.d_a_aa, "%.2f" % hb.planar1, "%.2f" % hb.planar2, hb.type])

def tofloat(x):
    try:
        return float(x)
    except:
        return None
    
@handlify(mode='r')
def readHbonds(fh):
    r = defaultdict( innerDefaultdict( innerDefaultdict( innerDefaultdict( innerDefaultdict( innerDefaultdict( dict ) ) ) ) ) )
    reader = csv.reader(fh, dialect='excel-tab')
    for row in reader:
        pdb_id, chain1, resname1, resno1, icode1, donor_atom, chain2,  resname2, resno2, icode2, acceptor_atom, distance, a_d_dd, d_a_aa, planar1, planar2, htype = row
        x = StructObject(chain1 = chain1, resname1 = resname1, resno1 = int(resno1), icode1 = icode1, donor_atom = donor_atom,
                         chain2 = chain2, resno2 = int(resno2), resname2 = resname2, icode2 = icode2, acceptor_atom = acceptor_atom,
                         distance = tofloat(distance), a_d_dd = tofloat(a_d_dd), d_a_aa = tofloat(d_a_aa), planar1 = tofloat(planar1), planar2 = tofloat(planar2), type=htype )
        r[pdb_id][ chain1 ][ (x.resno1, x.icode1) ][ x.donor_atom ][ chain2 ][ (x.resno2, x.icode2) ][ x.acceptor_atom ] = x 
    return r

def hbonds(pdb_fname):
    pdb_fname = os.path.abspath(pdb_fname)
    with tempDir() as tmp_dir:
        hb_fname = prepareHbonds(pdb_fname)
        return readHbonds(hb_fname)

################################################################
# ROSE FORMAT PARSING

@handlify()
def parseRoseData(fh):
    lc = 0
    r = defaultdict( innerDefaultdict( innerDefaultdict( dict ) ) )
    for line in fh:
        lc+=1
        if lc==1:
            continue
        resname1, resno1, atom1, resname2, resno2, atom2, distance, d_a_aa, a_d_dd, planar1, planar2 = ( line[1:4], int(line[4:8]), line[13:16],
                                                                                                         line[14:22], int(line[22:26]), line[31:34],
                                                                                                         float(line[34:41]),
                                                                                                         float(line[43:49]), float(line[49:55]),
                                                                                                         float(line[55:61]), float(line[61:67]) )
        atom1 = atom1.strip()
        atom2 = atom2.strip()
        x = StructObject(resname1 = resname1, resno1 = resno1, atom1 = atom1, resname2 = resname2, resno2 = resno2, atom2 = atom2, distance = distance, d_a_aa = d_a_aa, a_d_dd = a_d_dd,
                         planar1 = planar1, planar2 = planar2)
        r[resno1][resno2][atom1][atom2] = x
    return r

@handlify()
def parseIgorData(fh):
    r = defaultdict( innerDefaultdict( innerDefaultdict( dict ) ) )        
    for line in fh:
        atom1, resno1, atom2, resno2, d_a_aa, a_d_dd, planar1, planar2, distance = line.strip().split()
        resno1, resno2 = int(resno1), int(resno2)
        d_a_aa, a_d_dd, planar1, planar2, distance = float(d_a_aa), float(a_d_dd), float(planar1), float(planar2), float(distance)
        x = StructObject(resno1 = resno1, atom1=atom1, resno2=resno2, atom2=atom2, distance=distance, d_a_aa = d_a_aa, a_d_dd = a_d_dd, planar1=planar1, planar2=planar2)
        r[resno1][resno2][atom1][atom2] = x
    return r

if __name__=='__main__':
    pdb_file = sys.argv[1]
    my_bonds = hbonds(pdb_file, True)
    #hb_file = sys.argv[2]
    hb_file = sys.stdout
    writeHbonds(hb_file, hbonds(pdb_file))
    sys.exit(0)
    
    #x = readHbonds(hb_file)
    #igor_fn = sys.argv[2]
    #igor_bonds = parseIgorData(igor_fn)
    for chain1 in sorted(hbs):
        if chain1!=the_chain:
            continue
        hbs_ch1 = hbs[chain1]
        for resno1, icode1 in sorted(hbs_ch1):
            hbs_ch1_r1 = hbs_ch1[ (resno1, icode1) ]
            for donor_atom in sorted(hbs_ch1_r1):
                hbs_ch1_r1_a1 = hbs_ch1_r1[donor_atom]
                for chain2 in sorted(hbs_ch1_r1_a1):
                    if chain2!=chain1:
                        continue
                    hbs_ch1_r1_a1_ch2 = hbs_ch1_r1_a1[chain2]
                    for resno2, icode2 in hbs_ch1_r1_a1_ch2:
                        hbs_ch1_r1_a1_ch2_r2 = hbs_ch1_r1_a1_ch2[(resno2,icode2)]
                        for acceptor_atom in sorted(hbs_ch1_r1_a1_ch2_r2):
                            hb = hbs_ch1_r1_a1_ch2_r2[ acceptor_atom ]
                            
