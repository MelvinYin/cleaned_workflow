from bokeh.layouts import column, row, Spacer
from figures import ConsoleOutput, TextInputComponent, ButtonComponent, \
    TextBoxComponent, ConsoleTextConsoleRow, SingleImageComponent, \
    RBGImageComponent

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
        self._logo_fig_1 = RBGImageComponent(specs.logo_fig_1)
        self._logo_fig_2 = RBGImageComponent(specs.logo_fig_2)
        self._logo_fig_3 = RBGImageComponent(specs.logo_fig_3)
        self._logo_fig_4 = RBGImageComponent(specs.logo_fig_4)
        self._logo_fig_5 = RBGImageComponent(specs.logo_fig_5)
        self._logo_fig_6 = RBGImageComponent(specs.logo_fig_6)

        self._logo_descr_1 = TextBoxComponent(specs.logo_descr_1)
        self._logo_descr_2 = TextBoxComponent(specs.logo_descr_2)
        self._logo_descr_3 = TextBoxComponent(specs.logo_descr_3)
        self._logo_descr_4 = TextBoxComponent(specs.logo_descr_4)
        self._logo_descr_5 = TextBoxComponent(specs.logo_descr_5)
        self._logo_descr_6 = TextBoxComponent(specs.logo_descr_6)

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
        logo_1_row = column(row(Spacer(width=100), self._logo_descr_1.figure),
                            self._logo_fig_1.figure)
        logo_2_row = column(row(Spacer(width=100), self._logo_descr_2.figure),
                            self._logo_fig_2.figure)
        logo_3_row = column(row(Spacer(width=100), self._logo_descr_3.figure),
                            self._logo_fig_3.figure)
        logo_4_row = column(row(Spacer(width=100), self._logo_descr_4.figure),
                            self._logo_fig_4.figure)
        logo_5_row = column(row(Spacer(width=100), self._logo_descr_5.figure),
                            self._logo_fig_5.figure)
        logo_6_row = column(row(Spacer(width=100), self._logo_descr_6.figure),
                            self._logo_fig_6.figure)
        alignment_col = column(Spacer(height=17), self._alignment_img.figure,
                               Spacer(height=17), logo_header_row,
                               Spacer(height=10), logo_1_row, logo_2_row,
                               logo_3_row, logo_4_row, logo_5_row, logo_6_row)
        fig_row = row(left_text_input_col, Spacer(width=210),
                      right_console_col, Spacer(width=50), alignment_col)
        layout = column(header_row, Spacer(height=30), title_row, fig_row)
        return layout
