import os
import shutil

def move_replace(input_path, output_dir, output_filename=None):
    if not (os.path.isdir(input_path) or os.path.isfile(input_path)):
        return
    filename = input_path.rsplit("/", maxsplit=1)[-1]
    if os.path.isdir(input_path):
        if output_filename:
            output_path = output_dir + output_filename
        else:
            output_path = output_dir + filename
        if os.path.isdir(output_path):
            shutil.rmtree(output_path)
        shutil.move(input_path, output_dir)

        if output_filename:
            os.rename(output_dir + filename, output_dir + output_filename)
    else:
        if output_filename:
            output_path = output_dir + output_filename
        else:
            output_path = output_dir + filename
        if os.path.isfile(output_path):
            os.remove(output_path)
        shutil.move(input_path, output_path)
    return
