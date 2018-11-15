import sys
from utils import read_cmd_args
import subprocess

def run(kwargs):
    consensus_filename = kwargs['consensus']
    seq_filename = kwargs['seqs']
    output_folder = kwargs['output_folder']
    with open(consensus_filename, 'r') as cons_file:
        consensus_seqs = []
        for line in cons_file:
            consensus_seqs.append(line.strip())
    file_popens = []
    for i, cons in enumerate(consensus_seqs):
        file = open(output_folder + "/meme_out{}.txt".format(i), 'w')
        popen = subprocess.Popen(["./external_scripts/meme/bin/meme",
                             seq_filename, "-text", "-protein", "-cons",
                             cons, "-w", "30", "-nmotifs", "1"], stdout=file)
        file_popens.append((file, popen))
    while file_popens:
        file, popen = file_popens.pop(0)
        popen.wait()
        file.close()

def main(args):
    kwargs = read_cmd_args(args, 'consensus seqs output_folder')
    run(kwargs)
    return

if __name__ == '__main__':
    main(sys.argv)