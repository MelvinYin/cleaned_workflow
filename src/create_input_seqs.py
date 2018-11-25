import os

def main(kwargs):
    dataset_dir = kwargs['input_dir']
    output = kwargs['output']
    edited_lines = []
    for filename in os.listdir(dataset_dir):
        header = filename.split(".")[0]
        with open(f"{dataset_dir}/{filename}") as file:
            for line in file:
                if line.startswith(">"):
                    line = ">" + header + line[1:]
                edited_lines.append(line)

    with open(output, "w") as file:
        for line in edited_lines:
            file.write(line)
