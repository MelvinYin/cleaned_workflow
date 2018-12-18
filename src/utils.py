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

