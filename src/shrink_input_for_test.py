import sys

def select_to_keep(seqs_filename, divide_by):
    to_keep = []
    with open(seqs_filename, 'r') as file:
        selected = True
        seq_count = 0
        for line in file:
            if line.startswith(">"):
                seq_count += 1
                if seq_count % divide_by == 0:
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
    filelines = select_to_keep(kwargs['seqs'], int(kwargs['divideby']))
    write_to_file(filelines, kwargs['output'])
    return