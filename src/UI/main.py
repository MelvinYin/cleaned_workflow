from bokeh.plotting import show, curdoc

# This need to be above the rest, so sys_path is updated.
from ui_utils import folders_to_add

folders_to_add(['bokeh_ui', 'utils'])
folders_to_add(['figures'], suffix='bokeh_ui')

from scorer import Scorer
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

