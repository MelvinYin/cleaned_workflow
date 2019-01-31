import os
import shutil
import functools

# def move_replace(input_path, output_dir):
#     if not (os.path.isdir(input_path) or os.path.isfile(input_path)):
#         return
#     if os.path.isdir(output_dir):
#         filename = input_path.rsplit("/", maxsplit=1)[1]
#         output_path = f"{output_dir}/{filename}"
#         if os.path.isdir(output_path):
#             shutil.rmtree(output_path)
#         elif os.path.isfile(output_path):
#             os.remove(output_path)
#     shutil.move(input_path, output_dir)
#     return

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

def move_replace(input_path, output_path):
    assert (os.path.isdir(input_path) or os.path.isfile(input_path))
    if os.path.isfile(output_path):
        os.remove(output_path)
    shutil.move(input_path, output_path)
    return

def move_into(input_path, output_dir):
    assert (os.path.isdir(input_path) or os.path.isfile(input_path))
    assert os.path.isdir(output_dir)
    filename = input_path.rsplit("/", maxsplit=1)[1]
    output_path = f"{output_dir}/{filename}"
    if os.path.isdir(output_path):
        shutil.rmtree(output_path)
    elif os.path.isfile(output_path):
        os.remove(output_path)
    shutil.move(input_path, output_dir)
    return

def check_inout(input_file=(), input_dir=(), output=()):
    def decorator(func):
        for file in input_file:
            try:
                assert os.path.isfile(file)
            except AssertionError:
                print(f"In {func}, input file {file} does not exist.")
                raise
        for dir in input_dir:
            try:
                assert os.path.isdir(dir)
            except AssertionError:
                print(f"In {func}, input dir {dir} does not exist.")
                raise
        for item in output:
            if os.path.isfile(item):
                print(f"In {func}, output file {item} exists, will remove.")
                os.remove(item)
            elif os.path.isdir(item):
                print(f"In {func}, output dir {item} exists, will remove.")
                shutil.rmtree(item)
        @functools.wraps(func)
        def decorated(*args, **kwargs):
            ret_value = func(*args, **kwargs)
            for item in output:
                try:
                    assert (os.path.isfile(item) or os.path.isdir(item))
                except AssertionError:
                    print(f"In {func}, output {item} is not produced.")
                    raise
            return ret_value
        return decorated
    return decorator
