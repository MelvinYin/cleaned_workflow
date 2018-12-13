import os
import re


class PSSM_meme:
    def __init__(self, **kwargs):
        self.fname = kwargs['filename']
        self.start, self.pssms, self.end = self.meme_format_parser()

    def relabel_pssms(self):
        for i, pssm in enumerate(self.pssms):
            for j, line in enumerate(pssm):
                if re.search("MEME-[0-9]+", line):
                    line = re.sub("MEME-[0-9]+", f"MEME-{i}", line)
                pssm[j] = line
            self.pssms[i] = pssm
        return True

    def get_evalue(self):
        evalues = []
        for pssm in self.pssms:
            evalue = None
            for line in pssm:
                if line.startswith("MOTIF"):
                    evalue = re.search(
                        "E-value = ([0-9]+[\.]?[0-9]?e[+-][0-9]+)", line)
                    evalue = float(evalue.group(1))
                    break
            # evalue can be == 0 so avoid assert evalue
            assert evalue is not None
            evalues.append(evalue)
        return evalues

    def get_entropy_bits(self):
        bits = []
        for pssm in self.pssms:
            entropy = None
            for line in pssm:
                if re.match("\(([0-9]+\.[0-9]+) bits\)", line):
                    entropy = float(re.match("\(([0-9]+\.[0-9]+) bits\)", line)
                                    .group(1))
                    break
            assert entropy is not None
            bits.append(entropy)
        return bits

    def get_composition(self):
        assert self.start
        composition = []
        at_composition = False
        for line in self.start:
            # Check if we are in composition
            if line.startswith("Letter frequencies in dataset"):
                at_composition = True
            if line.startswith("Background letter frequencies"):
                break
            if at_composition:
                freq_re = re.findall("(1.[0]+)|(0.[0-9]+)", line)
                for (f1, f2) in freq_re:
                    composition.append(float(f1) if f1 else float(f2))
        assert len(composition) == 20  # num alphabets
        return composition

    def meme_format_parser(self):
        assert os.path.isfile(self.fname)
        start = []
        pssms = []
        end = []
        at_start = True
        at_pssm = False
        at_end = False
        with open(self.fname, 'r') as file:
            current_pssm = []
            for line in file:
                # Check if we are in start
                if line.startswith("MOTIF") and at_start:
                    at_start = False
                    at_pssm = True
                if at_start:
                    start.append(line)

                # Split into pssms
                if line.startswith("MOTIF") and current_pssm:
                    assert len(current_pssm) > 30
                    pssms.append(current_pssm)
                    current_pssm = [line]
                elif line.startswith("Stopped"):
                    at_end = True
                    at_pssm = False
                    assert len(current_pssm) > 30
                    pssms.append(current_pssm)
                elif at_pssm:
                    current_pssm.append(line)

                # Check if we are at end
                if at_end:
                    end.append(line)
        assert start
        assert pssms
        assert end
        return start, pssms, end

    def delete_pssm(self, to_delete):
        # Deletion is only by order and not label
        for i in sorted(to_delete, reverse=True):
            del self.pssms[i]
        return

    def keep_pssm(self, to_keep):
        for i in range(len(self.pssms))[::-1]:
            if i not in to_keep:
                del self.pssms[i]
        return

    def output(self, fname=None):
        if not fname:
            fname = self.fname
        with open(fname, "w") as file:
            for line in self.start:
                file.write(line)
            for pssm in self.pssms:
                for line in pssm:
                    file.write(line)
            for line in self.end:
                file.write(line)
        return

class PSSM_conv:
    def __init__(self, kwargs):
        pass
