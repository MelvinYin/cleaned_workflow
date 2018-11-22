import os
import shutil

def read_cmd_args(sys_args, commands_str):
    kwargs = dict()
    commands = commands_str.split(" ")
    for i, arg in enumerate(sys_args):
        for command in commands:
            if arg == '-{}'.format(command):
                kwargs[command] = sys_args[i + 1]
    assert len(kwargs) == len(commands), "Command-line args missing, " \
                                         "present=<{}>".format(list(kwargs.keys()))
    return kwargs

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

def to_trash(file, trash_dir):
    return move_replace(file, self.dir.trash)