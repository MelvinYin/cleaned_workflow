import sys

def cropped_seqs(seqs_filename, divisor):
    to_keep = []
    with open(seqs_filename, 'r') as file:
        selected = True
        seq_count = 0
        for line in file:
            if line.startswith(">"):
                seq_count += 1
                if seq_count % divisor == 0:
                    selected = True
                else:
                    selected = False
            if selected:
                # Breaking up by list so it isn't one large str
                to_keep.append(line)
    return to_keep

def write_to_file(filelines, output_filename):
    with open(output_filename, 'w') as file:
        for line in filelines:
            file.write(line)
    return

def main(kwargs):
    filelines = cropped_seqs(kwargs['seqs'], int(kwargs['divisor']))
    write_to_file(filelines, kwargs['output'])
    return