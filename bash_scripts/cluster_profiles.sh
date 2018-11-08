cp -f ./files/meme_format3.txt ./external_scripts/meme
cd ./external_scripts/meme
mast -remcorr meme_format3.txt consolidated_cropped.fasta
cd ..
cd ..
rm -r ./files/mast_out
cp -f -r ./external_scripts/meme/mast_out ./files
cp -f ./files/mast_out/mast.txt ./files
python cluster_final.py
