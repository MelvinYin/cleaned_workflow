from bokeh.layouts import column, row, Spacer
from figures import ConsoleOutput, TextInputComponent, ButtonComponent, \
    TextBoxComponent, ConsoleTextConsoleRow, SingleImageComponent, MultiImageComponent

class UI:
    def __init__(self, callback, specs):
        self.callback = callback
        self.specs = specs
        self._console = ConsoleOutput(specs.console)
        self._ti = TextInputComponent(specs.ti)
        self._button = ButtonComponent(specs.button, self._button_callback)
        self._app_title = TextBoxComponent(specs.app_title)
        self._input_header = TextBoxComponent(specs.input_header)
        self._prob_header = TextBoxComponent(specs.prob_header)
        self._alignment_header = TextBoxComponent(specs.alignment_header)

        self._family_prob_1 = ConsoleTextConsoleRow(specs.con_text_con_1)
        self._family_prob_2 = ConsoleTextConsoleRow(specs.con_text_con_2)
        self._family_prob_3 = ConsoleTextConsoleRow(specs.con_text_con_3)
        self._family_prob_4 = ConsoleTextConsoleRow(specs.con_text_con_4)
        self._family_prob_5 = ConsoleTextConsoleRow(specs.con_text_con_5)

        self._alignment_img = SingleImageComponent(specs.mast_img)
        self._profile_logo_header = TextBoxComponent(specs.logo_header)
        self._profile_logos = MultiImageComponent(specs.logos_specs)

        self.layout = self._plot()

    def _button_callback(self):
        # check for input values in self._ti
        self._console.figure_update(self._ti.current_value)
        return

    def _plot(self):
        header_row = row(Spacer(width=400), self._app_title.figure)
        title_row = row(Spacer(width=40), self._input_header.figure,
                        Spacer(width=110), self._prob_header.figure,
                        Spacer(width=250), self._alignment_header.figure)
        family_prob_table = column(self._family_prob_1.figure,
                                  self._family_prob_2.figure,
                                  self._family_prob_3.figure,
                                  self._family_prob_4.figure,
                                  self._family_prob_5.figure)
        left_text_input_col = column(self._ti.widget, self._button.widget)
        right_console_col = column(Spacer(height=5), family_prob_table)
        logo_header_row = row(Spacer(width=190),
                              self._profile_logo_header.figure)
        alignment_col = column(Spacer(height=17), self._alignment_img.figure,
                               Spacer(height=17), logo_header_row,
                               Spacer(height=10), self._profile_logos.figure)
        fig_row = row(left_text_input_col, Spacer(width=210),
                      right_console_col, Spacer(width=50), alignment_col)
        layout = column(header_row, Spacer(height=30), title_row, fig_row)
        return layout
