Build process: 

Original file sent should contain the following:
src
project_main.py
files, containing input_fasta, input_pdb, sfld_datasets
external_scripts, and within it the 3 zip files

Build script should do:
cd external_scripts
>build meme
tar xzf meme-5.0.1_1.tar.gz
cd meme-5.0.1
./configure --prefix=$PWD --with-url=http://meme-suite.org/ --enable-build-libxml2 --enable-build-libxslt
make
make install
cd ..

> build dhcl
> Install conda env first
unzip dhcl.zip
cd dhcl
<ask user to write here their python directory>
<something like> 
conda create --name dhcl_p python=2.7
/home/melvin/anaconda3/envs/dhcl_p/bin/python build.py
<or if they do conda activate, then python build.py>
cd ..

> build converge
unzip pipeline.zip
cd pipeline
make converge
cd ..

That settles the external script building. 




Notes:
Not using unpatched meme_5.0.2 because ./configure line does not appear to work with a changed prefix.

Todo:
check for gcc, advise on that. 
Change num_processor to a cmd arg? Or at least justify why not in readme
