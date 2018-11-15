import shutil
import os
import subprocess

COPY_INIT_SEQS=False
RUN_DHCL=False
EXTRACT_CONSENSUS=False
RUN_MEME=False
BUILD_STARTER=False
CLEAN_MEME=False
MERGE_MEME=False
SCREEN_EVALUE=False
SCREEN_ENTROPY=False
SCREEN_CORRELATED=False
SCREEN_NON_COMBI=False
MAST_COMBI=True
CLUSTER_COMBI=True
CREATE_CLUSTER_MOTIFS=True
CREATE_CLUSTER_LOGOS=True
DELETE_INIT_FASTA=False
DELETE_INTERMEDIATE=False

# For testing purposes
REDUCE_DHCL=False
SHRINK_INPUT=False

LOG="files/log.txt"
TRASH="files/_trash"
INPUT_SEQS="files/input_seqs_copy.fasta"
INPUT_PDB="files/input_pdb_test"

shutil.rmtree(TRASH, ignore_errors=True)
os.mkdir(TRASH)
open(LOG, 'w').close()

# exec > $LOG

# Use a reduced seq file for testing
if SHRINK_INPUT:
    print("SHRINK_INPUT:")
    from shrink_input_for_test import main
    main(dict(seqs=INPUT_SEQS, output=INPUT_SEQS, divideby=100))
    print("Success!\n")
# Input: ./files/input_seqs.fasta
# Output: ./files/input_seqs.fasta


# Copy initial seq files
if COPY_INIT_SEQS:
    print("COPY_INIT_SEQS:")
    shutil.copy(INPUT_SEQS, './external_scripts/pipeline')
    shutil.copy('./files/single_seq.fasta', './external_scripts/meme')
    shutil.copy(INPUT_SEQS, './external_scripts/meme')
    print("Success!\n")

# Input: ./files/input_seqs.fasta    ./files/single_seq.fasta
# Output: ./external_scripts/pipeline/input_seqs.fasta    ./external_scripts/meme/input_seqs.fasta
#         ./external_scripts/meme/single_seq.fasta


# Run dhcl
# WARNING: TAKES 0.5 HOUR
if RUN_DHCL:
    print("RUN_DHCL:")
    shutil.rmtree("files/from_dhcl", ignore_errors=True)
    os.mkdir("files/from_dhcl")
    command = 'source activate p2.7 ' \
              '&& python ./external_scripts/dhcl/executables/everything.py ' \
              '-d {} --outdir files/from_dhcl ' \
              '&& source deactivate'.format(INPUT_PDB)
    subprocess.run(command.split())
    print("Success!\n")
# Input: ./files/input_pdb
# Output: ./files/from_dhcl


# Match dhcl loop index to seqs and extract loop seq
if EXTRACT_CONSENSUS:
    print("EXTRACT_CONSENSUS:")
    python src/process_dhcl_output_meme.py \
    -dhcl_dir files/from_dhcl -fasta_dir files/input_fasta -output files/init_seed_seqs.fasta
    print("Success!\n")

# Input: ./files/from_dhcl    ./files/input_fasta
# Output: ./files/init_seed_seqs.fasta


# Reduce number of DHCL loops, for testing
if $REDUCE_DHCL
then
    print "REDUCE_DHCL:"
    python src/reduce_dhcl_test.py \
    -consensus files/consensus_seqs.txt -output files/consensus_seqs.txt
    print "Success!"
    print ""
fi
# Input: ./files/consensus_seqs.txt
# Output: ./files/consensus_seqs.txt


# Get meme profile per dhcl loop.
if $RUN_MEME
then
    print "RUN_MEME:"
    mv -f ./files/meme_full $TRASH
    mkdir ./files/meme_full
    python src/run_meme_on_dhcl_loops.py -consensus files/consensus_seqs.txt \
    -seqs $INPUT_SEQS -output_folder ./files/meme_full
    print "Success!"
    print ""
fi
# Input: ./files/consensus_seqs.txt    ./files/input_seqs.fasta
# Output: ./files/meme_full


# Build meme_starter
if $BUILD_STARTER
then
    print "BUILD_STARTER:"
    ./external_scripts/meme/bin/meme -text -protein -w 30 -p 7 -nmotifs 2 -nostatus \
    $INPUT_SEQS &>>./files/meme_starter.txt
    print "Success!"
    print ""
fi
# Input: ./files/input_seqs.txt
# Output: ./files/meme_starter.txt


# Clean meme files if len(seqs) is small
if $CLEAN_MEME
then
    print "CLEAN_MEME:"
    python src/meme_cleaner.py -input $INPUT_SEQS -output $INPUT_SEQS
    for file in ./files/meme_full/*
    do
        python src/meme_cleaner.py -input $file -output $file
    done
    print "Success!"
    print ""
fi



# Merge meme
if $MERGE_MEME
then
    print "MERGE_MEME:"
    python src/meme_merger.py -meme_folder files/meme_full \
    -output files/meme_merged.txt -meme_starter files/meme_starter.txt
    print "Success!"
    print ""
fi
# Input: ./files/meme_full    ./files/meme_starter.txt
# Output: ./files/meme_merged.txt


# Screen evalue
if $SCREEN_EVALUE
then
    print "SCREEN_EVALUE:"
    cp -f ./files/meme_merged.txt    ./files/meme.txt
    python src/remove_motifs_with_low_evalue.py \
    -meme files/meme.txt -output files/meme_evalue_screened.txt
    mv -f ./files/meme.txt $TRASH
    print "Success!"
fi
# Input: ./files/meme_merged.txt
# Output: ./files/meme_evalue_screened.txt


# Screen entropy
if $SCREEN_ENTROPY
then
    print "SCREEN_ENTROPY:"
    cp -f ./files/meme_evalue_screened.txt    ./files/meme.txt
    python src/screen_motif_entropy.py -meme files/meme.txt -output files/meme_format.txt
    mv -f ./files/meme.txt $TRASH
    print "Success!"
fi
# Input: ./files/meme_evalue_screened.txt
# Output: ./files/meme_format.txt


# Screen correlated profiles
if $SCREEN_CORRELATED
then
    print "SCREEN_CORRELATED:"
    cp -f ./files/meme_format.txt ./external_scripts/meme
    cd ./external_scripts/meme
    mast -remcorr meme_format.txt consolidated_single.fasta &>> ../.$LOG
    cd ..
    cd ..
    mv -f ./external_scripts/meme/mast_out ./files/mast_single
    cp -f ./files/mast_single/mast.txt ./files
    python src/mast_remove_profiles.py \
    -meme_in files/meme_format.txt -meme_out files/meme_format2.txt &>> $LOG
    mv -f ./external_scripts/meme/meme_format.txt $TRASH
    mv -f ./files/mast.txt $TRASH
    print "Success!"
fi
# Input: ./files/meme_format.txt    ./external_scripts/meme/consolidated_single.fasta
# Output: ./files/mast_single    ./files/meme_format2.txt


# Screen non-combi profiles
if $SCREEN_NON_COMBI
then
    print "SCREEN_NON_COMBI:"
    cp -f ./files/meme_format2.txt ./external_scripts/meme
    # Can add profile logo creation here, on meme_format2.txt
    cd ./external_scripts/meme
    mast -remcorr meme_format2.txt $INPUT_SEQS &>> ../.$LOG
    cd ..
    cd ..
    mv -f ./external_scripts/meme/mast_out ./files/mast_nocorr
    cp -f ./files/mast_nocorr/mast.txt ./files
    python src/cluster_get_profiles.py &>> $LOG
    python src/mast_remove_profiles_using_pkl.py &>> $LOG
    mv -f ./external_scripts/meme/meme_format2.txt $TRASH
    mv -f ./files/mast.txt $TRASH
    print "Success!"
fi
# Input: ./files/meme_format2.txt    ./external_scripts/meme/consolidated.fasta
# Output: ./files/mast_nocorr    ./files/meme_format3.txt


# MAST on final set of profiles
# WARNING: TAKES 1 HOURS
if $MAST_COMBI
then
    print "MAST_COMBI:"
    cp -f ./files/meme_format3.txt ./external_scripts/meme
    cd ./external_scripts/meme
    mast -remcorr meme_format3.txt $INPUT_SEQS &>> ../.$LOG
    cd ..
    cd ..
    mv ./external_scripts/meme/mast_out ./files/mast_onlycombi
    cp -r ./files/mast_onlycombi ./output/mast
    mv -f ./external_scripts/meme/meme_format3.txt $TRASH
    print "Success!"
fi
# Input: ./files/meme_format3.txt    ./external_scripts/meme/consolidated.fasta
# Output: ./files/mast_onlycombi    ./output/mast


# Cluster profiles
if $CLUSTER_COMBI
then
    print "CLUSTER_COMBI:"
    cp -f ./files/mast_onlycombi/mast.txt ./files
    python src/cluster_final.py
    mv -f ./files/mast.txt $TRASH
    print "Success!"
fi
# Input: ./files/mast_onlycombi/mast.txt
# Output: ./output/cluster_description.txt    ./files/cluster_centroids.pkl

# Create separate cluster motif files
if $CREATE_CLUSTER_MOTIFS
then
    print "CREATE_CLUSTER_MOTIFS:"
    mkdir -p ./files/motifs
    python src/split_motifs_into_cluster_motifs.py
    print "Success!"
fi
# Input: ./files/meme_format3.txt    ./files/cluster_centroids.pkl
# Output: multiple ./files/motifs/motifs_in_cluster_{}.txt

# Create logos for each cluster
if $CREATE_CLUSTER_LOGOS
then
    print "CREATE_CLUSTER_LOGOS:"
    mkdir -p ./output/logos
    for filename in ./files/motifs/*; do
        filename=${filename:15}
        cluster_i=${filename:18:-4}
        BLAH=1
        mkdir -p ./output/logos/cluster_$cluster_i
        while ./external_scripts/meme/ceqlogo -i$BLAH ./files/motifs/$filename -o ./output/logos/cluster_$cluster_i/logo_$BLAH.png -f PNG
        do
            BLAH=$((BLAH+1))
        done
    done
    python src/rename_cluster_logos.py
    print "Success!"
fi
# Input: multiple ./files/motifs/motifs_in_cluster_{}.txt
# Output: multiple ./output/logos/cluster_{}/logo_{}.png

# Delete init .fasta files
if $DELETE_INIT_FASTA
then
    print "DELETE_INIT_FASTA:"
    mv -f ./external_scripts/pipeline/$INPUT_SEQS $TRASH
    mv -f ./external_scripts/meme/consolidated_single.fasta $TRASH
    mv -f ./external_scripts/meme/$INPUT_SEQS $TRASH
    print "Success!"
fi

# Delete all intermediate files
if $DELETE_INTERMEDIATE
then
    [ -f ./files/init_seed_seqs.fasta ] && mv -f ./files/init_seed_seqs.fasta $TRASH
    [ -f ./files/meme_format.txt ] && mv -f ./files/meme_format.txt $TRASH
    [ -f ./files/mast_single ] && mv -rf ./files/mast_single $TRASH
    [ -f ./files/meme_format2.txt ] && mv -f ./files/meme_format2.txt $TRASH
    [ -f ./files/mast_nocorr ] && mv -rf ./files/mast_nocorr $TRASH
    [ -f ./files/meme_format3.txt ] && mv -f ./files/meme_format3.txt $TRASH
    [ -f ./files/mast_onlycombi ] && mv -rf ./files/mast_onlycombi $TRASH
    [ -f ./files/cluster_centroids.pkl ] && mv -f ./files/cluster_centroids.pkl $TRASH
    [ -f ./files/motifs ] && mv -rf ./files/motifs $TRASH
    [ -f $LOG ] && mv -f $LOG $TRASH
fi