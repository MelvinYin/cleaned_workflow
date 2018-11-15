import subprocess
import os
import multiprocessing

with open("test_output.txt", "w") as file:
    test = subprocess.Popen(["../external_scripts/meme/bin/meme",
                             "../files/consolidated_test.fasta", "-text",
                             "-protein", "-cons",
                             "KVVPVAGHDSMLLNLSGAHGPLFTRNILIL", "-w", "30",
                             "-nmotifs", "1"],
                            stdout=file)
    test.wait()