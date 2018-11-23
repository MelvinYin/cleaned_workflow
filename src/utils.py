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

# def rename(path, new_name):
#     assert os.path.isfile(path) or os.path.isdir(path)
#     dir_path = path.rsplit("/", maxsplit=1)[0]
#     new_path = dir_path + "/" + new_name
#     os.rename(path, new_path)
#     return