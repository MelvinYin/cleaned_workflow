import inspect
import os
import sys

src_path = inspect.currentframe().f_code.co_filename.rsplit("/", maxsplit=1)[0]
sys.path.append(src_path)

def folders_to_add(folders, suffix=None):
    sys.path.append(src_path)
    if suffix:
        _src_path = src_path + "/" + suffix
    else:
        _src_path = src_path + "/"
    for folder in folders:
        sys.path.append(_src_path + folder)

def _convert_url_to_bokeh(url):
    _url = os.path.join(os.path.basename(os.path.dirname(__file__)),
                        'static', url)
    # _url = f"./static/{url}"
    return _url
