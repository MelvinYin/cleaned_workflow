# chmod +rwx full_script.sh

MEME_BASH=false # FALSE
RUN_MEME2=false # FALSE
SCREEN_CORRELATED2=false # FALSE

COPY_INIT_SEQS=false
RUN_DHCL=true
EXTRACT_CONSENSUS=true
RUN_MEME=true
BUILD_STARTER=true
CLEAN_MEME=true
MERGE_MEME=false
SCREEN_EVALUE=false
SCREEN_ENTROPY=false
SCREEN_CORRELATED=false
SCREEN_NON_COMBI=false
MAST_COMBI=false
CLUSTER_COMBI=false
CREATE_CLUSTER_MOTIFS=false
CREATE_CLUSTER_LOGOS=false
DELETE_INIT_FASTA=false
DELETE_INTERMEDIATE=false

# For testing purposes
REDUCE_DHCL=false
SHRINK_INPUT=false

LOG=./files/log.txt
TRASH=./files/_trash
INPUT_SEQS=./files/input_seqs_copy.fasta
INPUT_PDB=./files/input_pdb_test

mkdir -p $TRASH
[ -f $LOG ] && mv $LOG $TRASH
touch $LOG

# exec > $LOG

# Use a reduced seq file for testing
if $SHRINK_INPUT
then
    echo "SHRINK_INPUT:"
    python src/shrink_input_for_test.py \
    -seqs $INPUT_SEQS -output $INPUT_SEQS -divideby 100
    echo "Success!"
    echo ""
fi
# Input: ./files/input_seqs.fasta
# Output: ./files/input_seqs.fasta


# Copy initial seq files
if $COPY_INIT_SEQS
then
    echo "COPY_INIT_SEQS:"
    cp -f $INPUT_SEQS ./external_scripts/pipeline
    cp -f ./files/single_seq.fasta ./external_scripts/meme
    cp -f $INPUT_SEQS ./external_scripts/meme
    echo "Success!"
    echo ""
fi
# Input: ./files/input_seqs.fasta    ./files/single_seq.fasta
# Output: ./external_scripts/pipeline/input_seqs.fasta    ./external_scripts/meme/input_seqs.fasta
#         ./external_scripts/meme/single_seq.fasta


# Run dhcl
# WARNING: TAKES 0.5 HOUR
if $RUN_DHCL
then
    echo "RUN_DHCL:"
    source activate p2.7
    mkdir -p ./files/from_dhcl
    python ./external_scripts/dhcl/executables/everything.py \
    -d $INPUT_PDB --outdir ./files/from_dhcl
    source deactivate
    echo "Success!"
    echo ""
fi
# Input: ./files/input_pdb
# Output: ./files/from_dhcl


# Match dhcl loop index to seqs and extract loop seq
if $EXTRACT_CONSENSUS
then
    echo "EXTRACT_CONSENSUS:"
    python src/process_dhcl_output_meme.py \
    -dhcl_dir files/from_dhcl -fasta_dir files/input_fasta -output files/init_seed_seqs.fasta
    echo "Success!"
    echo ""
fi
# Input: ./files/from_dhcl    ./files/input_fasta
# Output: ./files/init_seed_seqs.fasta


# Reduce number of DHCL loops, for testing
if $REDUCE_DHCL
then
    echo "REDUCE_DHCL:"
    python src/reduce_dhcl_test.py \
    -consensus files/consensus_seqs.txt -output files/consensus_seqs.txt
    echo "Success!"
    echo ""
fi
# Input: ./files/consensus_seqs.txt
# Output: ./files/consensus_seqs.txt


# Get meme profile per dhcl loop. 
if $RUN_MEME
then
    echo "RUN_MEME:"
    mv -f ./files/meme_full $TRASH
    mkdir ./files/meme_full
    python src/run_meme_on_dhcl_loops.py -consensus files/consensus_seqs.txt \
    -seqs $INPUT_SEQS -output_folder ./files/meme_full
    echo "Success!"
    echo ""
fi
# Input: ./files/consensus_seqs.txt    ./files/input_seqs.fasta
# Output: ./files/meme_full


# Build meme_starter
if $BUILD_STARTER
then
    echo "BUILD_STARTER:"
    ./external_scripts/meme/bin/meme -text -protein -w 30 -p 7 -nmotifs 2 -nostatus \
    $INPUT_SEQS &>>./files/meme_starter.txt
    echo "Success!"
    echo ""
fi
# Input: ./files/input_seqs.txt
# Output: ./files/meme_starter.txt


# Clean meme files if len(seqs) is small
if $CLEAN_MEME
then
    echo "CLEAN_MEME:"
    python src/meme_cleaner.py -input $INPUT_SEQS -output $INPUT_SEQS
    for file in ./files/meme_full/*
    do
        python src/meme_cleaner.py -input $file -output $file
    done
    echo "Success!"
    echo ""
fi
# Input: ./files/meme_full    ./files/meme_starter.txt
# Output: ./files/meme_full    ./files/meme_starter.txt


# Merge meme
if $MERGE_MEME
then
    echo "MERGE_MEME:"
    python src/meme_merger.py -meme_folder files/meme_full \
    -output files/meme_merged.txt -meme_starter files/meme_starter.txt
    echo "Success!"
    echo ""
fi
# Input: ./files/meme_full    ./files/meme_starter.txt
# Output: ./files/meme_merged.txt


# Screen evalue
if $SCREEN_EVALUE
then
    echo "SCREEN_EVALUE:"
    cp -f ./files/meme_merged.txt    ./files/meme.txt
    python src/remove_motifs_with_low_evalue.py \
    -meme files/meme.txt -output files/meme_evalue_screened.txt
    mv -f ./files/meme.txt $TRASH
    echo "Success!"
fi
# Input: ./files/meme_merged.txt
# Output: ./files/meme_evalue_screened.txt


# Screen entropy
if $SCREEN_ENTROPY
then
    echo "SCREEN_ENTROPY:"
    cp -f ./files/meme_evalue_screened.txt    ./files/meme.txt
    python src/screen_motif_entropy.py -meme files/meme.txt -output files/meme_format.txt
    mv -f ./files/meme.txt $TRASH
    echo "Success!"
fi
# Input: ./files/meme_evalue_screened.txt
# Output: ./files/meme_format.txt


# Screen correlated profiles
if $SCREEN_CORRELATED
then
    echo "SCREEN_CORRELATED:"
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
    echo "Success!"
fi
# Input: ./files/meme_format.txt    ./external_scripts/meme/consolidated_single.fasta
# Output: ./files/mast_single    ./files/meme_format2.txt


# Screen non-combi profiles
if $SCREEN_NON_COMBI
then
    echo "SCREEN_NON_COMBI:"
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
    echo "Success!"
fi
# Input: ./files/meme_format2.txt    ./external_scripts/meme/consolidated.fasta
# Output: ./files/mast_nocorr    ./files/meme_format3.txt


# MAST on final set of profiles
# WARNING: TAKES 1 HOURS
if $MAST_COMBI
then
    echo "MAST_COMBI:"
    cp -f ./files/meme_format3.txt ./external_scripts/meme
    cd ./external_scripts/meme
    mast -remcorr meme_format3.txt $INPUT_SEQS &>> ../.$LOG
    cd ..
    cd ..
    mv ./external_scripts/meme/mast_out ./files/mast_onlycombi
    cp -r ./files/mast_onlycombi ./output/mast
    mv -f ./external_scripts/meme/meme_format3.txt $TRASH
    echo "Success!"
fi
# Input: ./files/meme_format3.txt    ./external_scripts/meme/consolidated.fasta
# Output: ./files/mast_onlycombi    ./output/mast


# Cluster profiles
if $CLUSTER_COMBI
then
    echo "CLUSTER_COMBI:"
    cp -f ./files/mast_onlycombi/mast.txt ./files
    python src/cluster_final.py
    mv -f ./files/mast.txt $TRASH
    echo "Success!"
fi
# Input: ./files/mast_onlycombi/mast.txt
# Output: ./output/cluster_description.txt    ./files/cluster_centroids.pkl

# Create separate cluster motif files
if $CREATE_CLUSTER_MOTIFS
then
    echo "CREATE_CLUSTER_MOTIFS:"
    mkdir -p ./files/motifs
    python src/split_motifs_into_cluster_motifs.py
    echo "Success!"
fi
# Input: ./files/meme_format3.txt    ./files/cluster_centroids.pkl
# Output: multiple ./files/motifs/motifs_in_cluster_{}.txt

# Create logos for each cluster
if $CREATE_CLUSTER_LOGOS
then
    echo "CREATE_CLUSTER_LOGOS:"
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
    echo "Success!"
fi
# Input: multiple ./files/motifs/motifs_in_cluster_{}.txt
# Output: multiple ./output/logos/cluster_{}/logo_{}.png

# Delete init .fasta files
if $DELETE_INIT_FASTA
then
    echo "DELETE_INIT_FASTA:"
    mv -f ./external_scripts/pipeline/$INPUT_SEQS $TRASH
    mv -f ./external_scripts/meme/consolidated_single.fasta $TRASH
    mv -f ./external_scripts/meme/$INPUT_SEQS $TRASH
    echo "Success!"
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

##########################################################################################

# Screen correlated profiles again
if $SCREEN_CORRELATED2
then
    echo "SCREEN_CORRELATED2:"
    cp -f ./files/meme_format2.txt ./external_scripts/meme
    cd ./external_scripts/meme
    mast -remcorr meme_format2.txt consolidated_single.fasta &>> ../.$LOG
    cd ..
    cd ..
    mv -r ./files/mast_single $TRASH
    mv -f ./external_scripts/meme/mast_out ./files/mast_single
    cp -f ./files/mast_single/mast.txt ./files
    mv ./files/meme_format.txt $TRASH
    mv -f ./files/meme_format2.txt ./files/meme_format.txt
    python src/mast_remove_profiles.py &>> $LOG
    mv -f ./external_scripts/meme/meme_format2.txt $TRASH
    mv -f ./files/mast.txt $TRASH
    echo "Success!"
fi
# Input: ./files/meme_format.txt    ./external_scripts/meme/consolidated_single.fasta
# Output: ./files/mast_single    ./files/meme_format2.txt

# Create meme bash
if $MEME_BASH
then
    echo "MEME_BASH:"
    python src/create_bash_script_for_meme.py \
    -consensus files/consensus_seqs.txt -bash_output files/call_meme.sh \
    -seqs $INPUT_SEQS -output_folder ./meme_output
    echo "Success!"
fi
# Input: ./files/consensus_seqs.txt    ./files/consolidated.fasta
# Output: ./files/call_meme.sh

# Run meme script
# WARNING: TAKES 2 HOURS
if $RUN_MEME2
then
    echo "RUN_MEME2:"
    chmod +rwx ./files/call_meme.sh
    ./files/call_meme.sh
    mv -f ./external_scripts/meme/meme_output ./files/meme_full
    echo "Success!"
fi
# Input: ./files/call_meme.sh    ./files/consolidated.fasta
# Output: ./files/meme_full







