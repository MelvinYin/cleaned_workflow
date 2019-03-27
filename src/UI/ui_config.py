from enum import Enum, auto

def _convert_url_to_bokeh(url):
    _url = os.path.join(os.path.basename(os.path.dirname(__file__)),
                        'static', url)
    return _url

class FigureNames(Enum):
    console_output = auto()
    con_text_con_1 = auto()

class ConsoleOutputSpecs:
    def __init__(self):
        self.name = FigureNames.console_output
        self.title = "Console"
        self.text = 'Initial<p>'
        self.width = 50
        self.height = 20
        self.textbox_width = 400
        self.textbox_height = 500
        # Division by 2 so it fits well and within what bokeh uses.
        self.html_height = int(self.textbox_height / 2)
        _style = dict(border='2px solid rgb(200, 200, 200)',
                      height='{}px'.format(self.html_height))
        _style['overflow-y'] = 'auto'
        self.style = _style

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

class MastSpecs:
    def __init__(self, width=500, height=100,
                 image_fname='mast_alignment.png', style=None):
        image_fname = _convert_url_to_bokeh(image_fname)
        self.width = width
        self.height = height
        self.image_fname = image_fname
        self.x_range = (0, 10)
        self.y_range = (0, 5)
        self.img_xy = (0, 4)
        self.img_wh = (10, 4)
        if style is not None:
            self.style = style
        else:
            self.style = dict()

import os
class LogosSpecs:
    def __init__(self, width=500, height=400,
                 images=('logos_1.png', 'logos_2.png', 'logos_3.png',
                         'logos_4.png', 'logos_5.png', 'logos_6.png')):
        images = tuple([_convert_url_to_bokeh(url) for url in images])
        self.width = width
        self.height = height
        self.images = images
        self.x_range = (-1, 11)
        self.y_range = (-1, 26)
        template_img_xy = (0, 4)
        self.img_wh = (10, 4)
        self.img_xys = list([(template_img_xy[0],
                              self.img_wh[1]*i+template_img_xy[1])
                             for i in range(len(images))])
        # _style = dict(border='2px solid rgb(200, 200, 200)')
        # self.style = _style

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
        # _style = dict(border='2px solid rgb(200, 200, 200)',
        #               height='{}px'.format(self.textbox_height),
        #               width='{}px'.format(self.html_width))
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
        self.console = ConsoleOutputSpecs()
        self.ti = GenericSpecs(width=100, height=200, text='Input Sequence')
        self.button = GenericSpecs(width=100, height=100, text="Submit")

        app_title_style = dict(color='rgb(0, 0, 0)')
        app_title_style['text-decoration'] = 'underline'
        app_title_style['font-size'] = '200%'

        input_header_style = dict(color='rgb(0, 0, 0)')
        input_header_style['text-decoration'] = 'underline'
        input_header_style['font-size'] = '120%'

        prob_header_style = dict(color='rgb(0, 0, 0)')
        prob_header_style['text-decoration'] = 'underline'
        prob_header_style['font-size'] = '120%'

        alignment_header_style = dict(color='rgb(0, 0, 0)')
        alignment_header_style['text-decoration'] = 'underline'
        alignment_header_style['font-size'] = '120%'

        logos_header_style = dict(color='rgb(0, 0, 0)')
        logos_header_style['text-decoration'] = 'underline'
        logos_header_style['font-size'] = '120%'

        self.app_title = GenericSpecs(width=1000, height=10,
                                      text="Enolase Family Predictor",
                                      style=app_title_style)
        self.input_header = GenericSpecs(width=200, height=10,
                                         text="Input Sequence",
                                         style=input_header_style)
        self.prob_header = GenericSpecs(width=200, height=10,
                                        text="Family Prediction",
                                        style=prob_header_style)
        self.alignment_header = GenericSpecs(width=200, height=10,
                                             text="Profile Alignment",
                                             style=alignment_header_style)
        self.logo_header = GenericSpecs(width=200, height=10,
                                        text="Profile Logos",
                                        style=logos_header_style)

        t1 = GenericSpecs(width=50, height=10, text="=>")

        c11 = SingleLineConsoleSpecs(text='Enolase')
        c12 = SingleLineConsoleSpecs(text='98.2%')

        c21 = SingleLineConsoleSpecs(text='Galacterate Dehydratase')
        c22 = SingleLineConsoleSpecs(text='1.0%')

        c31 = SingleLineConsoleSpecs(text='Glucarate Dehydratase')
        c32 = SingleLineConsoleSpecs(text='0.6%')

        c41 = SingleLineConsoleSpecs(text='Mannonate Dehydratase')
        c42 = SingleLineConsoleSpecs(text='0.1%')

        c51 = SingleLineConsoleSpecs(text='Mandelate Racemase')
        c52 = SingleLineConsoleSpecs(text='0.1%')

        self.con_text_con_1 = ConsoleTextConsoleSpecs(c11, t1, c12)
        self.con_text_con_2 = ConsoleTextConsoleSpecs(c21, t1, c22)
        self.con_text_con_3 = ConsoleTextConsoleSpecs(c31, t1, c32)
        self.con_text_con_4 = ConsoleTextConsoleSpecs(c41, t1, c42)
        self.con_text_con_5 = ConsoleTextConsoleSpecs(c51, t1, c52)
        self.mast_img = MastSpecs()
        self.width = 100
        self.height = 100
        self.logos_specs = LogosSpecs()


