import os
import shutil
import re

def move_replace(input_path, output_dir):
    if not (os.path.isdir(input_path) or os.path.isfile(input_path)):
        return
    if os.path.isdir(output_dir):
        filename = input_path.rsplit("/", maxsplit=1)[1]
        output_path = f"{output_dir}/{filename}"
        if os.path.isdir(output_path):
            shutil.rmtree(output_path)
        elif os.path.isfile(output_path):
            os.remove(output_path)
    shutil.move(input_path, output_dir)
    return

def meme_rewritter(profiles, fname, to_keep=True, output=None, sub=False):
    if output is None:
        output = fname
    motif_count = 0
    to_write = []
    with open(fname, 'r') as rfile:
        deleting = False
        for i, line in enumerate(rfile):
            if line.startswith("MOTIF") and deleting:
                deleting = False
            if deleting:
                continue
            if line.startswith("MOTIF"):
                motif_count += 1
                if to_keep and motif_count not in profiles:
                    deleting = True
                    continue
                elif not to_keep and motif_count in profiles:
                    deleting = True
                    continue
            if sub:
                line = re.sub("MEME-[0-9]+", "MEME-{}".format(motif_count), line)
            to_write.append(line)

    with open(output, 'w') as wfile:
        for line in to_write:
            wfile.write(line)
    return

def check_fasta_validity(filename):
    assert os.path.isfile(filename)
    with open(filename) as file:
        try:
            first_line = next(file)
        except StopIteration:
            print(f"Empty fasta file in {filename}")
            raise
        assert first_line.startswith(">")
        assert first_line.strip()[1:]
    return

def conv_meme_parser(fname):
    assert os.path.isfile(fname)
    start = []
    pssms = []
    pssm_start = False
    current_pssm = []
    with open(fname, 'r') as file:
        for line in file:
            if line.startswith("MOTIF"):
                pssm_start = True
            if not pssm_start:
                start.append(line)
                continue
            if line.startswith('MOTIF') and current_pssm:
                assert len(current_pssm) > 30, current_pssm
                pssms.append(current_pssm)
                current_pssm = [line]
                continue
            current_pssm.append(line)
        assert len(current_pssm) > 30  # width of motif
        pssms.append(current_pssm)
    assert start
    assert pssms
    return start, pssms

def meme_format_parser(filename):
    assert os.path.isfile(filename)
    start = []
    pssms = []
    end = []
    at_pssm = False
    at_end = False
    with open(filename, 'r') as file:
        current_pssm = []
        for line in file:
            # Check if we are in start
            if line.startswith("MOTIF") and not at_pssm:
                at_pssm = True
            if not at_pssm:
                start.append(line)

            # Split into pssms
            if line.startswith("MOTIF") and current_pssm:
                assert len(current_pssm) > 30
                pssms.append(current_pssm)
                current_pssm = [line]
            elif at_pssm:
                current_pssm.append(line)

            # Check if we are at end
            if line.startswith("Stopped"):
                at_pssm = False
                assert len(current_pssm) > 30
                pssms.append(current_pssm)
                at_end = True
            if at_end:
                end.append(line)
    assert start
    assert pssms
    assert end
    return start, pssms, end

def get_composition_meme(start_lines):
    assert start_lines
    composition = []
    at_composition = False
    for line in start_lines:
        # Check if we are in composition
        if line.startswith("Letter frequencies in dataset"):
            at_composition = True
        if line.startswith("Background letter frequencies"):
            break
        if at_composition:
            freq_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
            for (f1, f2) in freq_re:
                composition.append(float(f1) if f1 else float(f2))
    assert len(composition) == 20   # num alphabets
    return composition

def get_composition_conv(start_lines):
    assert start_lines
    composition = []
    at_composition = False
    for line in start_lines:
        if line.startswith("Background"):
            at_composition = True
            continue
        if at_composition:
            freq_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
            if freq_re: # skip any \n after frequencies
                for (f1, f2) in freq_re:
                    freq = f1 if f1 else f2
                    composition.append(float(freq))
    assert len(composition) == 20   # num alphabets
    return composition

