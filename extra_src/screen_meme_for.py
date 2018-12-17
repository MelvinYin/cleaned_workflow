from utils import meme_format_parser, meme_pssm_get_entropy_bits, \
    meme_pssm_get_evalue

def main(kwargs):
    entropy_threshold = None
    evalue_threshold = None
    if "entropy_bits_threshold" in kwargs:
        entropy_threshold = kwargs['entropy_bits_threshold']
    if "evalue_threshold" in kwargs:
        evalue_threshold = kwargs['evalue_threshold']

    start, pssms, end = meme_format_parser(kwargs['memefile'])
    to_del = []
    for i, pssm in enumerate(pssms):
        if entropy_threshold:
            entropy = meme_pssm_get_entropy_bits(pssm)
            if entropy < entropy_threshold:
                to_del.append(i)
        if evalue_threshold:
            evalue = meme_pssm_get_evalue(pssm)
            if evalue > evalue_threshold:
                to_del.append(i)
    for i in to_del[::-1]:
        del pssms[i]
    with open(kwargs['memefile'], 'w') as file:
        for line in start:
            file.write(line)
        for pssm in pssms:
            for line in pssm:
                file.write(line)
        for line in end:
            file.write(line)
    return