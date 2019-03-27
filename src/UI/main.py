
from utils.sys_path_adder import folders_to_add
from bokeh.plotting import show, curdoc

folders_to_add(['bokeh_ui', 'config', 'p_model', 'tests', 'gs_components',
                'utils'])
folders_to_add(['figures', 'widgets'], suffix='bokeh_ui')

from ui import UI

from ui_config import UISpecs
import random
random.seed(1)

def _callback():
    print("call")
    return


def main():
    ui = UI(_callback, UISpecs())
    # show(ui.ui_layout)
    curdoc().add_root(ui.layout)

main()
# from bokeh.models import Image
# from bokeh.plotting import figure
# import numpy as np
#
# img = np.full([300], fill_value=200)
#
# _image2 = figure(x_range=(0, 10), y_range=(0, 10))
# _image2.image_url(url=['output_80.png'], x=0, y=10, w=10, h=10)
#
# show(_image2)