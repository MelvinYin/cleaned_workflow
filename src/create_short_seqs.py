def main(kwargs):
    input_seqs = kwargs['input']
    output = kwargs['output']
    length = kwargs['length']
    lines_to_copy = []
    with open(input_seqs, 'r') as rfile:
        counter = 0
        for line in rfile:
            if line.startswith(">"):
                counter += 1
            if counter > length:
                break
            lines_to_copy.append(line)

    with open(output, 'w') as wfile:
        for line in lines_to_copy:
            wfile.write(line)
    return