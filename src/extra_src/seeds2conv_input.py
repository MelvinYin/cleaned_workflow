

def main(kwargs):
    seedseq_filename = kwargs['seed_seqs']
    output = kwargs['output']
    to_write = ""

    with open(seedseq_filename, 'r') as rfile:
        for line in rfile:
            to_write += line.strip()

    with open(output, 'w') as wfile:
        wfile.write(">RANDOM\n")
        for i in range(len(to_write) // 60):
            wfile.write(to_write[i*60:(i+1)*60] + "\n")
        if not (len(to_write) % 60 == 0):
            wfile.write(to_write[(len(to_write) // 60) * 60:])
