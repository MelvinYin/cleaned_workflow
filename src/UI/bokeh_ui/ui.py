from bokeh.layouts import column, row, Spacer
from figures import ConsoleOutput, TextInputComponent, ButtonComponent, \
    TextBoxComponent, ConsoleTextConsoleRow, SingleImageComponent, \
    MultiImageComponent, ButtonURLComponent
from bokeh.models.callbacks import CustomJS
from bokeh.models import OpenURL
from ui_config import _convert_url_to_bokeh

class UI:
    def __init__(self, callback, specs):
        self.callback = callback
        self.specs = specs
        self._console = ConsoleOutput(specs.console)
        self._ti = TextInputComponent(specs.ti)

        self._app_title = TextBoxComponent(specs.app_title)
        self._input_header = TextBoxComponent(specs.input_header)
        self._prob_header = TextBoxComponent(specs.prob_header)
        self._placeholder1_header = TextBoxComponent(specs.tmp_output_header)
        self._alignment_header = TextBoxComponent(specs.alignment_header)

        self._family_prob_1 = ConsoleTextConsoleRow(specs.con_text_con_1)
        self._family_prob_2 = ConsoleTextConsoleRow(specs.con_text_con_2)
        self._family_prob_3 = ConsoleTextConsoleRow(specs.con_text_con_3)
        self._family_prob_4 = ConsoleTextConsoleRow(specs.con_text_con_4)
        self._family_prob_5 = ConsoleTextConsoleRow(specs.con_text_con_5)
        self._button = ButtonComponent(specs.button, self._button_callback)

        self._alignment_img = ButtonURLComponent(specs.mast_img,
                                              self._url_callback)
        self._profile_logo_header = TextBoxComponent(specs.logo_header)
        self._logos_fig = MultiImageComponent(specs.logos_fig)

        self._logo_descr_1 = TextBoxComponent(specs.logo_descr_1)
        self._logo_descr_2 = TextBoxComponent(specs.logo_descr_2)
        self._logo_descr_3 = TextBoxComponent(specs.logo_descr_3)
        self._logo_descr_4 = TextBoxComponent(specs.logo_descr_4)
        self._logo_descr_5 = TextBoxComponent(specs.logo_descr_5)
        self._logo_descr_6 = TextBoxComponent(specs.logo_descr_6)

        self.layout = self._plot()

    def _url_callback(self):
        url = "mast.html"
        url = _convert_url_to_bokeh(url)
        args = dict(url=url)
        obj = CustomJS(args=args, code='window.open(url);')
        return obj

    def _button_callback(self):
        # check for input values in self._ti
        print("button callback called")
        values_to_update = self.callback(self._ti.current_value)
        fam_probs = values_to_update['class_probs']
        self._family_prob_1.figure_update(fam_probs[0])
        self._family_prob_2.figure_update(fam_probs[1])
        self._family_prob_3.figure_update(fam_probs[2])
        self._family_prob_4.figure_update(fam_probs[3])
        self._family_prob_5.figure_update(fam_probs[4])
        return

    def _plot(self):
        header_row = row(Spacer(width=400), self._app_title.figure)
        title_row = row(Spacer(width=40), self._input_header.figure,
                        Spacer(width=110), self._prob_header.figure,
                        Spacer(width=250), self._placeholder1_header.figure)
        family_prob_table = column(self._family_prob_1.figure,
                                  self._family_prob_2.figure,
                                  self._family_prob_3.figure,
                                  self._family_prob_4.figure,
                                  self._family_prob_5.figure)
        left_text_input_col = column(self._ti.widget, self._button.widget)
        right_console_col = column(Spacer(height=5), family_prob_table)
        logo_header_row = row(Spacer(width=190),
                              self._profile_logo_header.figure)
        logo_descr_col = column(Spacer(height=20), self._logo_descr_1.figure,
                                Spacer(height=30),
                                self._logo_descr_2.figure, Spacer(height=30),
                                self._logo_descr_3.figure, Spacer(height=30),
                                self._logo_descr_4.figure, Spacer(height=30),
                                self._logo_descr_5.figure, Spacer(height=30),
                                self._logo_descr_6.figure)
        logos_set = row(logo_descr_col, Spacer(width=30),
                        self._logos_fig.figure)


        # alignment_col = column(Spacer(height=17), logo_header_row,
        #                        Spacer(height=10), logos_set, Spacer(height=17),
        #                        self._alignment_header.figure, Spacer(height=17),
        #                        self._alignment_img.widget)

        alignment_row = column(row(Spacer(width=300),
                                   self._alignment_header.figure),
                               Spacer(height=10),
                               row(Spacer(width=300),
                                   self._alignment_img.widget))

        logos_col = column(Spacer(height=17), logo_header_row,
                               Spacer(height=10), logos_set, Spacer(height=50),
                               alignment_row)
        fig_row = row(left_text_input_col, Spacer(width=210),
                      right_console_col, Spacer(width=50), logos_col)
        layout = column(header_row, Spacer(height=30), title_row, fig_row)
        return layout
