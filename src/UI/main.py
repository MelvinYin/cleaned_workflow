from utils.sys_path_adder import folders_to_add
from bokeh.plotting import show, curdoc
from scorer import Scorer

folders_to_add(['bokeh_ui', 'utils'])
folders_to_add(['figures'], suffix='bokeh_ui')

from ui import UI
from ui_config import UISpecs

scorer = Scorer()

def _callback(args):
    input_seq = args.strip().replace("\n", "")
    output = scorer.score_seq(input_seq)
    to_return = dict()
    to_return['class_probs'] = output
    return to_return

def main():
    ui = UI(_callback, UISpecs())
    # show(ui.layout)
    curdoc().add_root(ui.layout)

main()


# from bokeh.models import ColumnDataSource, OpenURL, TapTool
# from bokeh.plotting import figure, output_file, show
#
# output_file("openurl.html")
#
# p = figure(plot_width=400, plot_height=400,
#            tools="tap", title="Click the Dots")
#
# source = ColumnDataSource(data=dict(
#     x=[1, 2, 3, 4, 5],
#     y=[2, 5, 8, 2, 7],
#     color=["navy", "orange", "olive", "firebrick", "gold"]
#     ))
#
# p.circle('x', 'y', color='color', size=20, source=source)
#
# url = "http://www.colors.commutercreative.com/@color/"
# taptool = p.select(type=TapTool)
# taptool.callback = OpenURL(url=url)
#
# show(p)
