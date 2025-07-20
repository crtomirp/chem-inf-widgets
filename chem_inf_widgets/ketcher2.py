import os
import logging
from typing import List, Dict

from Orange.widgets import gui, widget
from Orange.widgets.settings import Setting
from PyQt5.QtWidgets import QLabel, QSizePolicy
from PyQt5.QtWebEngineWidgets import QWebEngineView, QWebEnginePage
from PyQt5.QtWebChannel import QWebChannel
from PyQt5.QtCore import QUrl, Qt, QObject, pyqtSlot, pyqtSignal
from Orange.data import Table, Domain, StringVariable

# Configure logging for debugging
logging.basicConfig(level=logging.INFO)
logger = logging.getLogger(__name__)

class Bridge(QObject):
    """Bridge exposed to JS via QWebChannel; emits SMILES strings."""
    smilesReceived = pyqtSignal(str)

    @pyqtSlot(str)
    def receiveSmiles(self, smiles: str) -> None:
        """Slot called from JS to deliver the SMILES string."""
        self.smilesReceived.emit(smiles)

class CustomWebEnginePage(QWebEnginePage):
    """Capture JavaScript console messages for debugging."""
    def javaScriptConsoleMessage(self, level, message, line, source):
        logger.debug(f"JS Console [{level}] {message} ({source}:{line})")

class OWKetcherMolecularSketcher(widget.OWWidget):
    name = "MolSKetcher II"
    description = "Ketcher-based molecular editor: export SMILES."
    icon = "icons/ketcher.png"
    outputs = [("Compounds", Table)]
    priority = 12
    want_main_area = True
    resizing_enabled = True

    # Stored SMILES data
    data: List[Dict[str, str]] = Setting([])

    def __init__(self):
        super().__init__()
        self.bridge = None
        self.channel = None

        self._setup_ketcher_view()
        self._setup_ui()
        # Set a fixed width for the control panel
        self.controlArea.setFixedWidth(250)
        self.data.clear()

    def _setup_ketcher_view(self) -> None:
        """Load standalone Ketcher HTML and set up QWebChannel bridge."""
        html_path = os.path.join(os.path.dirname(__file__), "ketcher", "standalone", "indexqt.html")
        if not os.path.isfile(html_path):
            raise FileNotFoundError(f"Ketcher HTML not found: {html_path}")

        self.web_view = QWebEngineView(self.mainArea)
        self.web_view.setPage(CustomWebEnginePage(self.web_view))

        # Initialize WebChannel
        self.channel = QWebChannel(self.web_view.page())
        self.bridge = Bridge()
        self.channel.registerObject('bridge', self.bridge)
        self.web_view.page().setWebChannel(self.channel)
        self.bridge.smilesReceived.connect(self._process_smiles)

        # Load the HTML
        self.web_view.load(QUrl.fromLocalFile(html_path))
        self.web_view.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        self.mainArea.layout().addWidget(self.web_view)
        self.web_view.loadFinished.connect(self._on_load_finished)

    def _setup_ui(self) -> None:
        """Build UI controls: count label, SMILES display, info label, and action buttons."""
        # Total count label
        self.count_label = QLabel("Total: 0")
        self.count_label.setWordWrap(True)
        self.controlArea.layout().addWidget(self.count_label)

        # SMILES label
        self.smiles_label = QLabel("Current SMILES: None")
        self.smiles_label.setWordWrap(True)
        self.controlArea.layout().addWidget(self.smiles_label)

        # Info label
        self.info_label = QLabel()
        self.info_label.setObjectName("status")
        self.info_label.setWordWrap(True)
        self.controlArea.layout().addWidget(self.info_label)

        # Actions box with fixed width
        actions = gui.widgetBox(self.controlArea, "Actions", orientation=Qt.Vertical)
        actions.setFixedWidth(200)  # limit actions panel width
        self.add_btn = gui.button(actions, self, "Add SMILES", callback=self._add_smiles)
        self.clear_btn = gui.button(actions, self, "Clear All", callback=self._clear)
        self.clear_btn.setEnabled(False)

    def _on_load_finished(self, ok: bool) -> None:
        """Called when the web view finishes loading."""
        if not ok:
            self.error("Failed to load Ketcher interface.")
        else:
            logger.info("Ketcher interface loaded.")

    def _add_smiles(self) -> None:
        """Trigger Ketcher to send SMILES via the QWebChannel bridge."""
        js = """
        (async ()=>{
          try {
            const s = await window.ketcher.getSmiles();
            window.bridge.receiveSmiles(String(s));
          } catch(e) {
            window.bridge.receiveSmiles('');
          }
        })();
        """
        self.web_view.page().runJavaScript(js)

    def _process_smiles(self, smiles: str) -> None:
        """Handle SMILES returned from Ketcher: update table and UI with count."""
        if not smiles or not smiles.strip():
            self.error("Draw a molecule first.")
            return

        # Store data and compute count
        self.data.append({"SMILES": smiles})
        count = len(self.data)

        # Update UI
        self.count_label.setText(f"Total: {count}")
        self.smiles_label.setText(f"Current SMILES: {smiles}")
        self.info_label.setText(f"✓ Added: {smiles} (Total: {count})")
        self.info_label.setStyleSheet("color: #4CAF50;")
        self.clear_btn.setEnabled(True)

        # Send updated table
        self._update_output()

    def _update_output(self) -> None:
        """Construct a simple Orange Table with only SMILES metadata."""
        if not self.data:
            self.send("Compounds", None)
            return

        var = StringVariable("SMILES")
        domain = Domain(attributes=[], metas=[var])
        rows = [[d["SMILES"]] for d in self.data]
        table = Table.from_list(domain, rows)
        table.name = "Compounds"
        self.send("Compounds", table)

    def _clear(self) -> None:
        """Clear all stored SMILES and reset the UI."""
        self.data.clear()
        self._update_output()
        self.clear_btn.setEnabled(False)
        self.count_label.setText("Total: 0")
        self.smiles_label.setText("Current SMILES: None")
        self.info_label.setText("✓ All cleared")
        self.info_label.setStyleSheet("color: #4CAF50;")

if __name__ == "__main__":
    from Orange.widgets.utils.widgetpreview import WidgetPreview
    WidgetPreview(OWKetcherMolecularSketcher).run()

