from bokeh.layouts import row, column
from bokeh.models.widgets import Paragraph, Div
from bokeh.models.layouts import WidgetBox, Spacer
from bokeh.models.widgets import Button
from bokeh.models.widgets.inputs import TextAreaInput
from bokeh.plotting import figure
from bokeh.models import WheelZoomTool

class SingleImageComponent:
    def __init__(self, specs):
        self.figure = self._set_figure(specs)

    def _set_figure(self, specs):
        fig = figure(x_range=specs.x_range, y_range=specs.y_range,
                     width=specs.width, height=specs.height)
        fig.image_url(url=[specs.image_fname], x=specs.img_xy[0],
                      y=specs.img_xy[1], w=specs.img_wh[0], h=specs.img_wh[1])
        return fig

class MultiImageComponent:
    def __init__(self, specs):
        self.figure = self._set_figure(specs)

    def _set_figure(self, specs):
        fig = figure(x_range=specs.x_range, y_range=specs.y_range,
                     width=specs.width, height=specs.height,
                     active_scroll='wheel_zoom')
        for i in range(len(specs.images)):
            fig.image_url(url=[specs.images[i]], x=specs.img_xys[i][0],
                          y=specs.img_xys[i][1], w=specs.img_wh[0],
                          h=specs.img_wh[1])
        return fig

class TextBoxComponent:
    def __init__(self, specs):
        self.figure = self._set_TB(specs)

    def _set_TB(self, specs):
        TB = Div(text=specs.text)
        TB.width = specs.width
        TB.height = specs.height
        TB.style = specs.style
        return TB

class ConsoleTextConsoleRow:
    def __init__(self, specs):
        self.name = specs.name
        self.specs = specs
        self._console_1 = SingleLineConsole(specs.console_1)
        self._text_inbetween = TextBoxComponent(specs.text_inbetween)
        self._console_2 = SingleLineConsole(specs.console_2)
        self.figure = self._build_layout()

    def _build_layout(self):
        layout = row(self._console_1.figure, Spacer(width=10),
                     self._text_inbetween.figure,
                     self._console_2.figure)

        return layout

    def figure_update(self, lines):
        self._console_1.figure_update(lines[0])
        self._console_2.figure_update(lines[1])
        return


class SingleLineConsole:
    def __init__(self, specs):
        self.specs = specs
        self.figure = self._build_paragraph()

    def _build_paragraph(self):
        paragraph = Div(width=self.specs.textbox_width,
                        height=self.specs.textbox_height,
                        text=self.specs.text, style=self.specs.style)
        return paragraph

    def figure_update(self, add_line):
        # Can't get bokeh div to scroll to end, it'll always reset to top
        # even if scrollHeight==scrollTop, etc. Keep it like this for now.
        self.figure.text = add_line
        return True

class TextInputComponent:
    def __init__(self, specs):
        # specs=namedtuple(typename=Header, width, height, placeholder)
        self.widget = self._set_TI(specs)
        self.current_value = 'Init'

    def _set_TI(self, specs):
        TI = TextAreaInput()
        TI.width = specs.width
        TI.height = specs.height
        TI.cols = 20

        TI.placeholder = specs.text
        TI.on_change("value", self._ti_callback)
        TI.title = None
        return TI

    def _ti_callback(self, attr, old, new):
        self.current_value = new
        return

class ButtonComponent:
    def __init__(self, specs, widget_callback):
        self.widget_callback = widget_callback
        self.widget = self._set_button(specs)

    def _set_button(self, specs):
        button = Button()
        button.label = specs.text
        button.width = specs.width
        button.height = specs.height
        button.on_click(self.widget_callback)
        return button

class ConsoleOutput:
    def __init__(self, specs):
        self.name = specs.name
        self.specs = specs
        self._paragraph = self._build_paragraph()
        self.figure = self._set_textbox(specs)
        self._rollover_count = 20

    def _build_paragraph(self):
        paragraph = Div(width=self.specs.textbox_width,
                        height=self.specs.textbox_height,
                        text=self.specs.text, style=self.specs.style)
        return paragraph

    def _set_textbox(self, specs):
        fig = column(row(Spacer(height=specs.height)),
                     row(Spacer(width=specs.width), self._paragraph))
        return fig

    def figure_update(self, add_line):
        # Can't get bokeh div to scroll to end, it'll always reset to top
        # even if scrollHeight==scrollTop, etc. Keep it like this for now.
        self._paragraph.text = add_line + "<br />" + self._paragraph.text
        return True