import json
import os
from Orange.widgets.widget import OWWidget, Input, Output
from Orange.widgets import gui
from Orange.data import Table
from PyQt5.QtWidgets import QWidget, QLabel, QHBoxLayout
from PyQt5.QtCore import Qt, pyqtSignal, QRect
from PyQt5.QtGui import QPainter, QPen

# Minimal implementation of an interval slider with two handles.
# This implementation is inspired by the approach described in the Orange Data Mining blog post :contentReference[oaicite:0]{index=0}.
class IntervalSlider(QWidget):
    intervalChanged = pyqtSignal(float, float)
    
    def __init__(self, low, high, minimum=0.0, maximum=100.0, parent=None):
        super().__init__(parent)
        self.minimum = minimum
        self.maximum = maximum
        self._low = low
        self._high = high
        self.setMinimumWidth(200)  # fixed width for alignment
        self.setFixedHeight(30)
        self.setMouseTracking(True)
        self.dragging = None

    @property
    def lowValue(self):
        return self._low

    @property
    def highValue(self):
        return self._high

    def paintEvent(self, event):
        painter = QPainter(self)
        rect = self.rect()
        # Define the groove rectangle with a margin on the sides.
        groove_rect = QRect(rect.left()+10, rect.center().y()-5, rect.width()-20, 10)
        # Draw the groove.
        painter.setPen(QPen(Qt.gray, 1))
        painter.setBrush(Qt.lightGray)
        painter.drawRect(groove_rect)
        # Calculate handle positions based on current low/high values.
        pos_low = groove_rect.left() + (self._low - self.minimum) / (self.maximum - self.minimum) * groove_rect.width()
        pos_high = groove_rect.left() + (self._high - self.minimum) / (self.maximum - self.minimum) * groove_rect.width()
        # Draw the selected interval between handles.
        selection_rect = QRect(int(pos_low), groove_rect.top(), int(pos_high - pos_low), groove_rect.height())
        painter.setBrush(Qt.blue)
        painter.drawRect(selection_rect)
        # Draw the two handles.
        handle_radius = 8
        painter.setBrush(Qt.white)
        painter.setPen(QPen(Qt.black, 1))
        painter.drawEllipse(int(pos_low - handle_radius), groove_rect.center().y()-handle_radius, handle_radius*2, handle_radius*2)
        painter.drawEllipse(int(pos_high - handle_radius), groove_rect.center().y()-handle_radius, handle_radius*2, handle_radius*2)

    def mousePressEvent(self, event):
        pos = event.pos()
        rect = self.rect()
        groove_rect = QRect(rect.left()+10, rect.center().y()-5, rect.width()-20, 10)
        pos_low = groove_rect.left() + (self._low - self.minimum) / (self.maximum - self.minimum) * groove_rect.width()
        pos_high = groove_rect.left() + (self._high - self.minimum) / (self.maximum - self.minimum) * groove_rect.width()
        # Determine if the click is near one of the handles.
        if abs(pos.x() - pos_low) < 10:
            self.dragging = 'low'
        elif abs(pos.x() - pos_high) < 10:
            self.dragging = 'high'

    def mouseMoveEvent(self, event):
        if self.dragging is None:
            return
        rect = self.rect()
        groove_rect = QRect(rect.left()+10, rect.center().y()-5, rect.width()-20, 10)
        pos = event.pos().x()
        # Convert the mouse position to a value.
        value = self.minimum + (pos - groove_rect.left()) / groove_rect.width() * (self.maximum - self.minimum)
        value = max(self.minimum, min(self.maximum, value))
        if self.dragging == 'low' and value < self._high:
            self._low = value
            self.intervalChanged.emit(self._low, self._high)
            self.update()
        elif self.dragging == 'high' and value > self._low:
            self._high = value
            self.intervalChanged.emit(self._low, self._high)
            self.update()

    def mouseReleaseEvent(self, event):
        self.dragging = None

class ADMETFilter(OWWidget):
    category = "Chemoinformatics"
    name = "ADMET Filter"
    description = "Filters compounds based on ADMET properties using configurable ranges."
    icon = "icons/admet_filter.png"
    priority = 82

    class Inputs:
        data = Input("Data", Table)

    class Outputs:
        filtered_data = Output("Filtered Data", Table)

    def __init__(self):
        super().__init__()
        # Use a different attribute name to avoid conflicting with the input signal.
        self.input_data = None  
        self.filter_config = []
        self.filter_widgets = []

        # Load filter configuration from JSON file
        json_path = os.path.join(os.path.dirname(__file__), "admet_filter.json")
        self.load_filter_config(json_path)

        # Create UI controls for each filter property using interval sliders
        self.create_filter_controls()

    def load_filter_config(self, filepath):
        """Load the filter configuration from a JSON file."""
        if not os.path.exists(filepath):
            self.error("Filter configuration file not found.")
            return
        try:
            with open(filepath, 'r') as f:
                self.filter_config = json.load(f)
        except Exception as e:
            self.error(f"Error loading filter config: {e}")

    def create_filter_controls(self):
        """Create UI controls (interval sliders with value labels) for each property defined in the filter config."""
        for prop in self.filter_config:
            prop_name = prop["Property"]
            display_name = prop.get("name", prop_name)
            ranges = prop.get("ranges", {}).get("range", {})
            min_val = ranges.get("min", 0.0)
            max_val = ranges.get("max", 100.0)

            # Use type checking: if both min_val and max_val are ints, format as integer;
            # otherwise, format as float.
            if isinstance(min_val, int) and isinstance(max_val, int):
                fmt = lambda x: f"{int(x)}"
            else:
                fmt = lambda x: f"{x:.2f}"
            
            # Create a horizontal box for each property.
            box = gui.hBox(self.controlArea, margin=0)
            gui.label(box, self, f"{display_name}:", tooltip=prop_name)
            
            # Create a container widget with a QHBoxLayout for labels and slider.
            container = QWidget()
            container_layout = QHBoxLayout(container)
            container_layout.setSpacing(5)
            # Left label showing current low value (fixed width for alignment)
            left_label = QLabel(fmt(min_val))
            left_label.setFixedWidth(40)
            left_label.setAlignment(Qt.AlignCenter)
            container_layout.addWidget(left_label)
            # The interval slider (fixed width)
            slider = IntervalSlider(low=min_val, high=max_val, minimum=min_val, maximum=max_val, parent=self)
            slider.setFixedWidth(200)
            container_layout.addWidget(slider)
            # Right label showing current high value (fixed width)
            right_label = QLabel(fmt(max_val))
            right_label.setFixedWidth(40)
            right_label.setAlignment(Qt.AlignCenter)
            container_layout.addWidget(right_label)
            
            # Connect the slider's signal to update the value labels and trigger filtering.
            slider.intervalChanged.connect(lambda low, high, left=left_label, right=right_label, fmt=fmt: (
                left.setText(fmt(low)),
                right.setText(fmt(high)),
                self.apply_filters()
            ))
            
            # Add the container to the main control area.
            box.layout().addWidget(container)
            
            # Store the slider widget for later use in filtering.
            self.filter_widgets.append({
                "property": prop_name,
                "slider": slider
            })

    @Inputs.data
    def set_data(self, data):
        """Receive input data and apply filters."""
        self.input_data = data
        self.apply_filters()

    def apply_filters(self):
        """Apply the current filters to the data and send the filtered output."""
        if self.input_data is None:
            self.Outputs.filtered_data.send(None)
            return

        # Ensure that we have an Orange Table. If self.input_data is a list,
        # convert it using the domain of the first instance.
        if not hasattr(self.input_data, 'domain'):
            if isinstance(self.input_data, list) and len(self.input_data) > 0 and hasattr(self.input_data[0], 'domain'):
                data_table = Table(self.input_data[0].domain, self.input_data)
            else:
                self.Outputs.filtered_data.send(None)
                return
        else:
            data_table = self.input_data

        filtered = []
        # Iterate over each instance in the data table.
        for instance in data_table:
            valid = True
            # Check each filter condition.
            for filt in self.filter_widgets:
                prop_name = filt["property"]
                slider = filt["slider"]
                # Find the variable (attribute or meta) with the matching name.
                var = next(
                    (v for v in data_table.domain.variables + data_table.domain.metas if v.name == prop_name),
                    None
                )
                if var is None:
                    valid = False
                    break
                try:
                    value = float(instance[var])
                    if not (slider.lowValue <= value <= slider.highValue):
                        valid = False
                        break
                except (ValueError, TypeError):
                    valid = False
                    break
            if valid:
                filtered.append(instance)

        # Construct a new table from the filtered instances.
        new_table = Table(data_table.domain, filtered) if filtered else None
        self.Outputs.filtered_data.send(new_table)

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(ADMETFilter).run()
