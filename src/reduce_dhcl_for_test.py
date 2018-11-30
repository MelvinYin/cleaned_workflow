import sys

def reduce_dhcl(kwargs):
    consensus_filename = kwargs['input']
    output = kwargs['output']
    divisor = kwargs['divisor']
    cons_seqs = []
    with open(consensus_filename, 'r') as rfile:
        for i, line in enumerate(rfile):
            if i % divisor == 0:
                cons_seqs.append(line)

    with open(output, 'w') as wfile:
        for line in cons_seqs:
            wfile.write(line)
    return True

def main(kwargs):
    # consensus denominator output
    reduce_dhcl(kwargs)
    return
