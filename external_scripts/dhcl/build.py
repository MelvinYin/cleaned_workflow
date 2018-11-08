#! /usr/bin/python2.5
from __future__ import with_statement
import sys, os, os.path, shutil, subprocess
# Change this if other compilers should be used
C_COMPILER = "gcc"
CXX_COMPILER = "g++"
THIS_DIR = os.path.dirname( os.path.abspath( sys.argv[0] ) )

from contextlib import contextmanager
@contextmanager
def changedDir(new_dir):
    old_dir = os.getcwd()
    os.chdir(new_dir)
    yield new_dir
    os.chdir(old_dir)

def makeExeName(rootn):
    if 'win' in sys.platform:
        return rootn + '.exe'
    else:
        return rootn
#####################################################################################
# REBUILD&INSTALL
# Rebuild everything from source and install it
EXEC_ROOT_DIR = os.path.join(THIS_DIR, 'executables')
SRC_ROOT_DIR = os.path.join(THIS_DIR, 'src')
PYSRC_ROOT_DIR = os.path.join(THIS_DIR, 'pysrc')

if not os.path.exists(EXEC_ROOT_DIR):
    os.mkdir(EXEC_ROOT_DIR)
# vdw executable
vdw_src_dir = os.path.join(SRC_ROOT_DIR, 'vdw')
vdw_exec = makeExeName("vdw")
contacts_exec = makeExeName("contacts")
hbonds_exec = makeExeName("hbonds")
with changedDir(vdw_src_dir):
    new_vdw_exec = os.path.join(EXEC_ROOT_DIR, vdw_exec)
    if os.path.exists(new_vdw_exec) and max( os.path.getmtime(fname) for fname in os.listdir('.') if fname.endswith('.cpp') or fname.endswith('.hpp') ) <= os.path.getmtime(new_vdw_exec):
        print >> sys.stderr, "'vdw' already up to date. Skipping."
    else:
        cmd = "%s -O3 main.cpp -o %s" % (CXX_COMPILER, vdw_exec)
        print >> sys.stderr, cmd
        if subprocess.call(cmd, shell=True)==0:
            shutil.move(vdw_exec, EXEC_ROOT_DIR)
            print >> sys.stderr, "'vdw' built and moved to %s." % EXEC_ROOT_DIR
        else:
            print >> sys.stderr, "Compilation of 'vdw' failed."
    new_contacts_exec = os.path.join(EXEC_ROOT_DIR, contacts_exec)
    if os.path.exists(new_contacts_exec) and max( os.path.getmtime(fname) for fname in os.listdir('.') if fname.endswith('.cpp') or fname.endswith('.hpp') ) <= os.path.getmtime(new_contacts_exec):
        print >> sys.stderr, "'contacts' already up to date. Skipping."
    else:
        cmd = "%s -O3 ca_main.cpp -o %s" % (CXX_COMPILER, contacts_exec)
        print >> sys.stderr, cmd
        if subprocess.call(cmd, shell=True)==0:
            shutil.move(contacts_exec, EXEC_ROOT_DIR)
            print >> sys.stderr, "'contacts' built and moved to %s." % EXEC_ROOT_DIR
        else:
            print >> sys.stderr, "Compilation of 'contacts' failed."
    new_hbonds_exec = os.path.join(EXEC_ROOT_DIR, hbonds_exec)
    if os.path.exists(new_hbonds_exec) and max( os.path.getmtime(fname) for fname in os.listdir('.') if fname.endswith('.cpp') or fname.endswith('.hpp') ) <= os.path.getmtime(new_hbonds_exec):
        print >> sys.stderr, "'hbonds' already up to date. Skipping."
    else:
        cmd = "%s -O3 hbs_main.cpp -o %s" % (CXX_COMPILER, hbonds_exec)
        print >> sys.stderr, cmd
        if subprocess.call(cmd, shell=True)==0:
            shutil.move(hbonds_exec, EXEC_ROOT_DIR)
            print >> sys.stderr, "'hbonds' built and moved to %s." % EXEC_ROOT_DIR
        else:
            print >> sys.stderr, "Compilation of 'hbonds' failed."

# # prep23 executable
# prep_src_dir = os.path.join(SRC_ROOT_DIR, 'prep23')
# prep_exec = makeExeName("prep23")
# with changedDir(prep_src_dir):
#     new_prep_exec = os.path.join(EXEC_ROOT_DIR, prep_exec)
#     if os.path.exists(new_prep_exec) and max( os.path.getmtime(fname) for fname in os.listdir('.') if fname.endswith('.c') or fname.endswith('.h') ) <= os.path.getmtime(new_prep_exec):
#         print >> sys.stderr, "'prep23' already up to date. Skipping."
#     else:
#         cmd = "%s -O3 -c -g Prep23.c" % (C_COMPILER,)
#         print >> sys.stderr, cmd
#         if subprocess.call(cmd, shell=True)==0:
#             pass
#         else:
#             print >> sys.stderr, "Compilation of 'prep23' failed."
#         cmd = "%s -O3 -c -g Proteins.c" % (C_COMPILER,)
#         print >> sys.stderr, cmd
#         if subprocess.call(cmd, shell=True)==0:
#             pass
#         else:
#             print >> sys.stderr, "Compilation of 'prep23' failed."
#         cmd = "%s -O3 -lm -o %s Prep23.o Proteins.o" % (C_COMPILER, prep_exec)
#         print >> sys.stderr, cmd
#         if subprocess.call(cmd, shell=True)==0:
#             shutil.move(prep_exec, EXEC_ROOT_DIR)
#             shutil.copy('Hi.dat', EXEC_ROOT_DIR)
#             print >> sys.stderr, "'prep23' built and moved to %s." % EXEC_ROOT_DIR
#         else:
#             print >> sys.stderr, "Compilation of 'prep23' failed."

# Python programs
# Simply move these to executable directory
py_src_dir = PYSRC_ROOT_DIR
with changedDir(py_src_dir):
    for fname in os.listdir('.'):
        #if fname.endswith('.py'):
        if os.path.isdir(fname):
            if os.path.exists( os.path.join(EXEC_ROOT_DIR, fname) ):
                shutil.rmtree( os.path.join(EXEC_ROOT_DIR, fname) )
            shutil.copytree(fname, os.path.join(EXEC_ROOT_DIR, fname) )
        else:
            shutil.copy(fname, os.path.join(EXEC_ROOT_DIR, fname) )
    print "Python sources installed."
