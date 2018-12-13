# Description
Loops are first obtained from the structure of a few input samples, using dhcl. The sequences of these loops are used as seed sequences, and PSSMs are created based on a set of sequence data from a superfamily, using either meme or converge. For those derived from meme, a filter for profiles with low e-value or entropy is carried out. Subsequently for both meme and converge, a screen for correlated profiles (as determined using mast) is performed. The profiles and superfamily sequence data are then run through mast to determine ideal combinations of profiles that best fit each of the sequences. Combinations with fewer than a certain number of sequences (5) are removed as noise, while the rest are clustered using their levenshtein distance and the agglomerative clustering method in sklearn. Clustering continues until we have a certain number of clusters (20). For each cluster, we then obtain the combination that is most representative, based on minimising the levenshtein distance between other combinations in the cluster, weighted by the number of sequences for each combination. These combinations are named cluster centroids, and profiles that make up these centroids are assembled, while those that are not present in any centroid are dropped. The assembled profiles are then run again through mast with the input superfamily sequence data, and the combinations clustered as in previous. These clusters form our output.

To test its performance, we obtain sequence data from each family in a particular superfamily instead. Sample structures from each family is also obtained. The family sequences are labelled and merged together, and the label is used at the end, when clusters have been determined, to compare how well the clusters correspond to actual family designations. 

# Build
The current implementation requires use of three external tools, meme_suite, dhcl, and converge from pipeline. Modification of these tools is kept to a minimum to provide for better compatibility. The build process is therefore somewhat convoluted as each has to be built separately.

dhcl requires p2.7 to build and run, so that will need to be installed as a separate env and the exec path provided. 

# Skip Start
If using anaconda, a p2.7 env can be added by:
conda create --name dhcl_p python=2.7
conda activate dhcl_p
conda install biopython
conda deactivate
The new conda env will be called dhcl_p. If using your own env, please change the exec path in the build process below, and in config.py in src. Remember to install biopython as well, tried on 1.70. 

Once the new env is installed, find the path to the python exec, it is usually in the form /home/<YOUR_USR>/anaconda3/envs/dhcl_p/bin/python, depending on where your anaconda is installed. 
# Skip End

Tested on Ubuntu 18.04

First, cd to the project folder/external_scripts, there should be three zip files.

> Build meme
tar xzf meme-5.0.1_1.tar.gz
mv meme-5.0.1 meme
cd meme
./configure --prefix=$PWD --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make install
mv src/ceqlogo bin/ceqlogo
cd ..

> Build converge
unzip pipeline.zip
cd pipeline
make converge
cd ..

> Build dhcl
unzip dhcl.zip
cd dhcl
/home/<YOUR_USR>/anaconda3/envs/dhcl_p/bin/python build.py
cd ..

Edits to dhcl:
1. everything.py
delete <from dhcl.pdb import *>
add
<
import warnings
from Bio.PDB.PDBExceptions import PDBConstructionWarning
warnings.simplefilter("ignore", PDBConstructionWarning)
>
2. executables/hprep.py
change:
from Bio.PDB import to_one_letter_code
to
from Bio.Data.SCOPData import protein_letters_3to1 as to_one_letter_code

<TODO: REMOVE IMPORT BIO.PDB FROM HPREP.PY

The build process is now complete. The remainder of the program are done in python, and does not require compilation.
# BUILD END

# Set Parameters
Next, if this is not done, go to config.py and edit p2_7_env such that it points to your dhcl_p exec. If you are not sure, simply paste the path into your bash and see if a python console is launched. 

While there, change num_processor used if you wish. This affects runtime for the motif-finding part (both meme and converge). Memory use also increases at least linearly for meme at least, so do adjust accordingly. 

When testing, setting seeds_divisor to a positive int decreases the number of profiles used from dhcl output by that factor. Setting seq_divisor decreases the number of sequences profiles are derived and compared against. Setting either or both to None turns off the dividing function. 

The program currently works only with sequence data labelled with their family initials. If this has been generated, set input_seqdir to None and input_seqs to the sequence path. Otherwise, set input_seqdir as the folder containing the individual family sequence data files. 

PSSMs can be obtained using either converge from pipeline or meme from meme_suite. For the former, go to project_main.py, under set_switches, and set True to all listed under converge, and False for those under meme. For the latter, do the reverse. Do not set both to be True. 
# Set End

# Add Data
Two types of files are required, the superfamily sequence data, and structural .pdb files for which profiles are generated. 

For now, labelled superfamily sequence data is used, to assess classification accuracy. 

For the former, go to sfld website, and download the sequence .fasta files for each family. Label the files with their family initials (Enolase => EL.fasta), then place them in files/sfld_datasets.

For the latter, go to sfld website, and copy a few pdb IDs from as many families as possible in the superfamily. Then, go to pdb website and download the .pdb and corresponding .fasta files. Move the .pdb files into input_pdb and fasta into input_fasta. 

When testing, one can reduce the number of .pdb files in the pdb folder, so dhcl doesn't have to process all of them. The fasta files can remain, only the relevant ones will be called. 
# Add End

# Run
The program requires an env of p3.6+. Activate the env, go to project folder (with project_main.py), and run python -B project_main.py. Output is in output folder. 

# Ignore subsequent if the above works for you, and if you're just using it as it is. 

# Notes
BUG: Logos currently don't correspond to the combination profiles when using MEME.

> Edits made to external_scripts are:
In converge_pssm.c, change delta=30 in all occurrences, so profiles are generated in blocks of 30, without overlap.

> Known bugs
converge_pssm.c has a known segfault, where random_index_30, defined in PSSM.h, is of size [1][30], but in line:
 S += Information[random_index_30[randomization][i]] * PSSM[random_index_30[randomization][i]][q-'A'];
in convergePSSM.c, is called with a first index larger than 1. Doesn't really matter since it's a random index anyway, but it can lead to compiler segfault-ing when any unused variable in the code is removed. 

There is also a bug regarding F/H/M, where there's a wrong definition/use of one, somewhere in the heavily nested parts of convergePSSM.c, but I cannot remember where it is now, and I recall previous tests not changing the output (significantly) when the bug is corrected. Should only lead to wrong output and not a compile or runtime error. 

Other found bugs and inconsequential lines did not significantly affect output and are left in. A proper refactor will probably necessitate a rewrite, and so is deferred until the final program is stable.

> Oddities in external_scripts
> pipeline
Order of #include matter in a few scripts, switching the header file (or header-of-header) import order can lead to segfaults due to variable not defined.

Reason why pipeline needs to be cd-ed in to run using bash, instead of running from outside by providing path of exec and files, is because BLOSUM filepath is hard-coded in and so needs to be present in the working directory. There may be other reasons. 

> dhcl
There are times when dhcl simply returns very few loops (2-4 instead of usual ~20) for a particular structure, no idea why and not that interested in finding out, just take note when running with very few .pdb files. 

Ediitng dhcl is necessary because dhcl is apparently written for a version of Biopython I cannot install. Don't bother attempting to get it to work as well, biopython apparently changes its api even at the 0.01 versioning level, making it troublesome to find out when the change actually takes place.

Be careful when installing older versions of biopython as well, I have not replicated this but doing a pip install biopython==1.60 on a p3 conda env leads to a successful install (not supposed to since biopython should work only on p2=>p3.5, if documentation is accurate), but this leads to a corruption of conda paths, and trying to call conda leads to ImportError on conda.cli. This necessitates a reinstallation of the entire conda and not only that particular env. 

Manual editing of dhcl is a temporary measure, once it is verified that dhcl can be pre-built, or that the changes can be implemented in src without needing to change the code post-build, the edits can be removed. 

> meme
mast output changes in mast.txt depending on how many input sequences are used. The additional content (at low seq numbers) are irrelevant to us, and are cropped out in clean_pssm.py.

Instructions as per meme_suite documentation in building, is to run make test in between make and make install. I'm getting failed tests (5 at last run), but there's no documentation on what to do with them or why they occur, and the program (meme/mast at least) appears to run fine, so it's ignored for now. 

> Generic
gcc required for meme_suite, pipeline, and possibly dhcl compilation, not yet tested on a fresh Ubuntu to see which minimal version works. pipeline requires gcc, c99, known to fail on Visual and/or Clang. 

> REQUIREMENT LIST
main: p3.7
dhcl_p: p2.7
conda biopython==1.70




								
