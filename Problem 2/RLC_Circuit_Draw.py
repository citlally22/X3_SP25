import sys
from PyQt5.QtWidgets import (
    QApplication, QMainWindow, QFileDialog,
    QGraphicsView, QGraphicsScene, QGraphicsItem, QGraphicsTextItem
)
from PyQt5.QtGui import QPainter, QPen, QBrush, QPainterPath
from PyQt5.QtCore import Qt, QPointF, QRectF
import math


# -------------------------
# Base class for components
# -------------------------
class WireComponentBase(QGraphicsItem):
    """Base class for all circuit components defined between two nodes."""

    def __init__(self, p1, p2, parent=None):
        """
        Initialize component endpoints.

        Args:
            p1: Tuple or QPointF, first node position.
            p2: Tuple or QPointF, second node position.
            parent: Optional parent QGraphicsItem.
        """
        super().__init__(parent)
        self.p1 = QPointF(*p1)
        self.p2 = QPointF(*p2)
        self.setFlag(QGraphicsItem.ItemIsSelectable, True)

    def boundingRect(self):
        """
        Defines the bounding rectangle around the component for redrawing.
        Returns:
            QRectF: Expanded bounding box around the item.
        """
        return QRectF(self.p1, self.p2).normalized().adjusted(-20, -20, 20, 20)


class StraightWire(WireComponentBase):
    """Draws a straight wire between two points."""

    def paint(self, painter, option, widget):
        painter.setPen(QPen(Qt.black, 2))
        painter.drawLine(self.p1, self.p2)


class ZigzagResistor(WireComponentBase):
    """Draws a resistor using a zigzag symbol."""

    def paint(self, painter, option, widget):
        painter.setRenderHint(QPainter.Antialiasing)
        dx, dy = self.p2.x() - self.p1.x(), self.p2.y() - self.p1.y()
        total_length = math.hypot(dx, dy)
        ux, uy = dx / total_length, dy / total_length
        perp = QPointF(-uy, ux)

        segments = 6  # Number of zigzag segments
        zigzag_length = total_length * 0.6
        seg_len = zigzag_length / segments
        amplitude = 5  # Vertical amplitude of zigzags

        start = self.p1
        zigzag_start = start + QPointF(ux * (total_length - zigzag_length) / 2, uy * (total_length - zigzag_length) / 2)

        path = QPainterPath(start)
        path.lineTo(zigzag_start)

        for i in range(1, segments + 1):
            base = zigzag_start + QPointF(ux * seg_len * i, uy * seg_len * i)
            sign = 1 if i % 2 else -1
            path.lineTo(base + perp * (amplitude * sign))

        path.lineTo(self.p2)
        painter.setPen(QPen(Qt.black, 2))
        painter.drawPath(path)


class PlateCapacitor(WireComponentBase):
    """Draws a parallel-plate capacitor."""

    def paint(self, painter, option, widget):
        painter.setRenderHint(QPainter.Antialiasing)
        dx, dy = self.p2.x() - self.p1.x(), self.p2.y() - self.p1.y()
        length = math.hypot(dx, dy)
        ux, uy = dx / length, dy / length
        perp = QPointF(-uy, ux)
        mid = (self.p1 + self.p2) * 0.5

        plate_width = 15
        spacing = 8  # Gap between capacitor plates

        # Draw wires to the plates
        painter.setPen(QPen(Qt.black, 2))
        painter.drawLine(self.p1, mid - QPointF(ux * spacing, uy * spacing))
        painter.drawLine(mid + QPointF(ux * spacing, uy * spacing), self.p2)

        # Draw the plates
        p1a = mid - QPointF(ux * spacing, uy * spacing) - perp * plate_width
        p1b = mid - QPointF(ux * spacing, uy * spacing) + perp * plate_width
        p2a = mid + QPointF(ux * spacing, uy * spacing) - perp * plate_width
        p2b = mid + QPointF(ux * spacing, uy * spacing) + perp * plate_width

        painter.drawLine(p1a, p1b)
        painter.drawLine(p2a, p2b)


class CoilInductor(WireComponentBase):
    """Draws a coil-style inductor using arcs."""

    def paint(self, painter, option, widget):
        painter.setRenderHint(QPainter.Antialiasing)
        dx, dy = self.p2.x() - self.p1.x(), self.p2.y() - self.p1.y()
        total_length = math.hypot(dx, dy)
        ux, uy = dx / total_length, dy / total_length
        perp = QPointF(-uy, ux)

        arc_count = 4
        coil_length = total_length * 0.7
        spacing = coil_length / arc_count
        radius = spacing / 2
        start = self.p1
        coil_start = start + QPointF(ux * (total_length - coil_length) / 2, uy * (total_length - coil_length) / 2)

        painter.setPen(QPen(Qt.black, 2))
        painter.drawLine(start, coil_start)

        # Draw arcs for the coil
        for i in range(arc_count):
            center = coil_start + QPointF(ux * spacing * (i + 0.5), uy * spacing * (i + 0.5))
            arc_rect = QRectF(center.x() - radius, center.y() - radius, 2 * radius, 2 * radius)
            angle = math.degrees(math.atan2(dy, dx))

            painter.save()
            painter.translate(center)
            painter.rotate(angle)
            painter.translate(-center)
            painter.drawArc(arc_rect, 0, 180 * 16)  # Qt uses 1/16 degree units
            painter.restore()

        coil_end = coil_start + QPointF(ux * coil_length, uy * coil_length)
        painter.drawLine(coil_end, self.p2)


class PowerSourceItem(WireComponentBase):
    """Draws a voltage source as a circle with + and - signs."""

    def paint(self, painter, option, widget):
        painter.setRenderHint(QPainter.Antialiasing)
        mid = (self.p1 + self.p2) * 0.5
        dx, dy = self.p2.x() - self.p1.x(), self.p2.y() - self.p1.y()
        length = math.hypot(dx, dy)
        ux, uy = dx / length, dy / length
        radius = 10

        # Draw wires and the circular voltage source
        painter.setPen(QPen(Qt.black, 2))
        painter.drawLine(self.p1, mid - QPointF(ux * radius, uy * radius))
        painter.drawEllipse(mid, radius, radius)
        painter.drawLine(mid + QPointF(ux * radius, uy * radius), self.p2)

        # Draw + and - signs for polarity
        font = painter.font()
        font.setBold(True)
        painter.setFont(font)
        painter.drawText(mid + QPointF(-6, -radius - 4), "+")
        painter.drawText(mid + QPointF(-6, radius + 12), "-")


# -------------------------
# Main GUI class
# -------------------------
class CircuitDiagramApp(QMainWindow):
    """Main application window for rendering the circuit diagram."""

    def __init__(self):
        super().__init__()
        self.setWindowTitle("RLC Circuit Visualizer")
        self.setGeometry(100, 100, 800, 600)

        self.view = QGraphicsView()
        self.scene = QGraphicsScene()
        self.view.setScene(self.scene)
        self.setCentralWidget(self.view)

        self._setup_menu()
        self.load_default_circuit()

    def _setup_menu(self):
        """Set up the File menu with load option."""
        menu = self.menuBar().addMenu("File")
        open_action = menu.addAction("Load Circuit File")
        open_action.triggered.connect(self.load_circuit_file)

    def load_default_circuit(self):
        """Load the default circuit from 'Circuit.txt'."""
        try:
            with open("Circuit.txt", "r") as f:
                self.parse_and_draw(f.read())
        except FileNotFoundError:
            print("Error: Circuit.txt not found.")

    def load_circuit_file(self):
        """Open a file dialog to choose and load a circuit file."""
        fname, _ = QFileDialog.getOpenFileName(self, "Open Circuit File", "", "Text files (*.txt);;All files (*)")
        if fname:
            with open(fname, "r") as f:
                self.parse_and_draw(f.read())

    def parse_and_draw(self, text):
        """
        Parse a circuit description and render it visually.

        Args:
            text (str): Text content from a circuit file.
        """
        self.scene.clear()
        nodes = {}

        # Pass 1: Parse and draw all node points
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if parts[0] == "node":
                _, label, x, y = parts[:4]
                nodes[label] = (float(x), float(y))
                self.scene.addEllipse(float(x) - 3, float(y) - 3, 6, 6, QPen(Qt.black), QBrush(Qt.black))

        # Pass 2: Parse and draw all components
        for line in text.splitlines():
            line = line.strip()
            if not line or line.startswith("#"):
                continue
            parts = line.split()
            if len(parts) < 4:
                continue
            element, n1, n2, val = parts[:4]
            if n1 not in nodes or n2 not in nodes:
                continue

            p1, p2 = nodes[n1], nodes[n2]

            # Instantiate appropriate component
            if element == "wire":
                item = StraightWire(p1, p2)
                unit = ""
            elif element == "resistor":
                item = ZigzagResistor(p1, p2)
                unit = " Ω"
            elif element == "capacitor":
                item = PlateCapacitor(p1, p2)
                unit = " F"
            elif element == "inductor":
                item = CoilInductor(p1, p2)
                unit = " H"
            elif element == "voltageSource":
                item = PowerSourceItem(p1, p2)
                unit = " V"
            else:
                continue

            self.scene.addItem(item)

            # Add value label next to component
            if element != "wire":
                mid = QPointF((p1[0] + p2[0]) / 2, (p1[1] + p2[1]) / 2)
                value_label = QGraphicsTextItem(f"{val}{unit}")
                value_label.setPos(mid + QPointF(10, -10))
                self.scene.addItem(value_label)

                # Add component symbol (R, C, L, V)
                symbol = "L" if element == "inductor" else element[0].upper()
                symbol_item = QGraphicsTextItem(symbol)
                symbol_item.setPos(mid + QPointF(-25, 10))
                font = symbol_item.font()
                font.setBold(True)
                symbol_item.setFont(font)
                self.scene.addItem(symbol_item)

        # Fit the scene to the view
        self.view.fitInView(self.scene.itemsBoundingRect(), Qt.KeepAspectRatio)


if __name__ == "__main__":
    app = QApplication(sys.argv)
    window = CircuitDiagramApp()
    window.show()
    sys.exit(app.exec_())


