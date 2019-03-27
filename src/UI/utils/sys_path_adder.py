import sys
import inspect
src_path = inspect.currentframe().f_code.co_filename.rsplit("/", maxsplit=2)[0]
sys.path.append(src_path)

def folders_to_add(folders, suffix=None):
    sys.path.append(src_path)
    if suffix:
        _src_path = src_path + "/" + suffix
    else:
        _src_path = src_path + "/"
    for folder in folders:
        sys.path.append(_src_path + folder)
