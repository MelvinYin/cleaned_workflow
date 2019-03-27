from utils.sys_path_adder import folders_to_add
from bokeh.plotting import show, curdoc

folders_to_add(['bokeh_ui', 'utils'])
folders_to_add(['figures'], suffix='bokeh_ui')

from ui import UI
from ui_config import UISpecs

def _callback():
    print("call")
    return

def main():
    ui = UI(_callback, UISpecs())
    # show(ui.ui_layout)
    curdoc().add_root(ui.layout)

main()