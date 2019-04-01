from utils import meme_format_parser, meme_pssm_get_entropy_bits, \
    meme_pssm_get_evalue

def main(kwargs):
    entropy_threshold = 40
    start, pssms, end = meme_format_parser(kwargs['memefile'])
    to_del = []
    for i, pssm in enumerate(pssms):
        entropy = meme_pssm_get_entropy_bits(pssm)
        if entropy < entropy_threshold:
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
