from enum import Enum, auto

from ui_utils import _convert_url_to_bokeh

class FigureNames(Enum):
    console_output = auto()
    con_text_con_1 = auto()

class GenericSpecs:
    def __init__(self, width=1, height=1, text='label', style=None):
        self.width = width
        self.height = height
        self.text = text
        if style is not None:
            self.style = style
        else:
            # self.style = dict(border='2px solid rgb(200, 200, 200)')
            self.style = dict()

class SingleLineConsoleSpecs:
    def __init__(self, textbox_width=100, textbox_height=10,
                 html_height=None, html_width=None,
                 text='SingleLineConsole'):
        self.textbox_width = textbox_width
        self.textbox_height = textbox_height
        if html_height is None:
            self.html_height = textbox_height
        else:
            self.html_height = html_height
        if html_width is None:
            self.html_width = textbox_width
        else:
            self.html_width = html_width
        self.text = text
        _style = dict(height='{}px'.format(self.textbox_height),
                      width='{}px'.format(self.html_width))
        self.style = _style

class ConsoleTextConsoleSpecs:
    def __init__(self, console_1=SingleLineConsoleSpecs(),
                 text_inbetween=GenericSpecs(),
                 console_2=SingleLineConsoleSpecs()):
        self.name = FigureNames.con_text_con_1
        self.console_1 = console_1
        self.text_inbetween = text_inbetween
        self.console_2 = console_2

class UISpecs:
    def __init__(self):
        self.ti = GenericSpecs(width=100, height=200, text='Input Sequence')
        self.input_button = GenericSpecs(width=100, height=100, text="Submit")

        app_title_style = dict(color='rgb(0, 0, 0)')
        app_title_style['text-decoration'] = 'underline'
        app_title_style['font-size'] = '200%'

        subheader_style = dict(color='rgb(0, 0, 0)')
        subheader_style['text-decoration'] = 'underline'
        subheader_style['font-size'] = '120%'

        self.app_title = GenericSpecs(width=1000, height=10,
                                      text="Enolase Family Predictor",
                                      style=app_title_style)
        self.input_header = GenericSpecs(width=200, height=10,
                                         text="Input Sequence",
                                         style=subheader_style)
        self.prob_header = GenericSpecs(width=200, height=10,
                                        text="Family Prediction",
                                        style=subheader_style)
        self.alignment_header = GenericSpecs(width=200, height=10,
                                             text="Profile Alignment",
                                             style=subheader_style)

        t1 = GenericSpecs(width=50, height=10, text="=>")

        c11 = SingleLineConsoleSpecs(text='Subgroup 1')
        c12 = SingleLineConsoleSpecs(text='20.0%')

        c21 = SingleLineConsoleSpecs(text='Subgroup 2')
        c22 = SingleLineConsoleSpecs(text='20.0%')

        c31 = SingleLineConsoleSpecs(text='Subgroup 3')
        c32 = SingleLineConsoleSpecs(text='20.0%')

        c41 = SingleLineConsoleSpecs(text='Subgroup 4')
        c42 = SingleLineConsoleSpecs(text='20.0%')

        c51 = SingleLineConsoleSpecs(text='Subgroup 5')
        c52 = SingleLineConsoleSpecs(text='20.0%')

        self.con_text_con_1 = ConsoleTextConsoleSpecs(c11, t1, c12)
        self.con_text_con_2 = ConsoleTextConsoleSpecs(c21, t1, c22)
        self.con_text_con_3 = ConsoleTextConsoleSpecs(c31, t1, c32)
        self.con_text_con_4 = ConsoleTextConsoleSpecs(c41, t1, c42)
        self.con_text_con_5 = ConsoleTextConsoleSpecs(c51, t1, c52)
        self.show_align = GenericSpecs(width=100, height=100,
                                     text="Show Alignment")
        self.width = 100
        self.height = 100



