def obtain_lines(fname):
    composition = ""
    pssms = []
    in_composition = False
    current_pssm = []
    with open(fname, 'r') as file:
        for line in file:
            if not in_composition \
                and line.startswith("Letter frequencies"):
                in_composition = True
                continue
            if in_composition and line.startswith("Background letter"):
                in_composition = False
                continue
            if in_composition:
                composition += line
                continue
            if line.startswith("letter-probability matrix"):
                current_pssm.append(line)
                continue
            if current_pssm and line.startswith("------------"):
                pssms.append(current_pssm)
                current_pssm = []
                continue
            if current_pssm:
                current_pssm.append(line)
    return composition, pssms

def format_output_lines(composition, pssms):
    output = []
    output.append("MEME version 4\n\n")
    output.append("ALPHABET= ACDEFGHIKLMNPQRSTVWY\n\n")
    output.append("Background letter frequencies\n")
    output += composition
    output.append("\n")
    for i, pssm in enumerate(pssms):
        output.append(f"MOTIF MEME-{i}\n")
        output += pssm
        output.append("\n")
    return output

def main(kwargs):
    input_fname = kwargs['input']
    output = kwargs['output']
    composition, pssms = obtain_lines(input_fname)
    output_lines = format_output_lines(composition, pssms)
    with open(output, 'w') as file:
        for line in output_lines:
            file.write(line)
    return


