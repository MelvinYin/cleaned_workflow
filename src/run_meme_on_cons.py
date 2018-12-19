import sys
from subprocess import Popen

def run(kwargs):
    consensus_filename = kwargs['consensus']
    seq_filename = kwargs['seqs']
    output_folder = kwargs['output_folder']
    num_p = kwargs['num_p']
    meme_dir = kwargs['meme_dir']
    with open(consensus_filename, 'r') as cons_file:
        consensus_seqs = []
        for line in cons_file:
            consensus_seqs.append(line.strip())
    i = 0
    while consensus_seqs:
        popens = []
        while len(popens) < num_p and consensus_seqs:
            cons = consensus_seqs.pop(0)
            i += 1
            filename = f"{output_folder}/meme_out{i}.txt"
            open(filename, 'w').close()
            command = f"{meme_dir}/meme {seq_filename} -text -protein " \
                f"-cons {cons} -w 30 -nmotifs 1 -nostatus &>> {filename}"
            popen = Popen(command, shell=True, executable='/bin/bash')
            popens.append(popen)
        while popens:
            popen = popens.pop(0)
            popen.wait()

def run2(kwargs):
    consensus_filename = kwargs['consensus']
    seq_filename = kwargs['seqs']
    output_folder = kwargs['output_folder']
    num_p = kwargs['num_p']
    meme_dir = kwargs['meme_dir']
    with open(consensus_filename, 'r') as cons_file:
        consensus_seqs = []
        for line in cons_file:
            consensus_seqs.append(line.strip())
    # create command
    filename = f"{output_folder}/meme_out1.txt"
    command = f"{meme_dir}/meme {seq_filename} -text -protein "
    for cons in consensus_seqs:
        command += f"-cons {cons} "
    command += f"-w 30 -nmotifs 1 -nostatus -p {num_p} &>> {filename}"
    popen = Popen(command, shell=True, executable='/bin/bash')
    popen.wait()
    return

def main(kwargs):
    # consensus seqs output_folder
    run(kwargs)
    return
