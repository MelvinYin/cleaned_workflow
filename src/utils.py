import os
import shutil

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

def meme_rewritter(profiles, fname, to_keep=True, output=None):
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
            to_write.append(line)

    with open(output, 'w') as wfile:
        for line in to_write:
            wfile.write(line)
    return