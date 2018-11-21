import sys
from utils import read_cmd_args

def reduce_dhcl(kwargs):
    consensus_filename = kwargs['input']
    output = kwargs['output']
    denominator = kwargs['denominator']
    cons_seqs = []
    with open(consensus_filename, 'r') as rfile:
        for i, line in enumerate(rfile):
            if i % denominator == 0:
                cons_seqs.append(line)

    with open(output, 'w') as wfile:
        for line in cons_seqs:
            wfile.write(line)
    return True

def main(kwargs):
    reduce_dhcl(kwargs)
    return

if __name__ == '__main__':
    kwargs = read_cmd_args(sys.argv, 'consensus denominator output')
    main(kwargs)