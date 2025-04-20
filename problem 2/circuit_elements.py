# circuit_elements.py
from PyQt5.QtWidgets import QGraphicsItem
from PyQt5.QtGui import QPainter, QPen
from PyQt5.QtCore import Qt, QRectF

class Node:
    """A class to represent a node in the circuit with a position."""
    def __init__(self, id, x, y):
        self.id = id
        self.x = float(x)
        self.y = float(y)

class CircuitElement(QGraphicsItem):
    """Base class for circuit elements, determining orientation and bounding box."""
    def __init__(self, node1, node2):
        super().__init__()
        self.node1 = node1
        self.node2 = node2
        self.is_horizontal = abs(node1.x - node2.x) > abs(node1.y - node2.y)
        self.length = 50 if self.is_horizontal else 30
        self.width = 20

    def boundingRect(self):
        """Defines the bounding rectangle for the element."""
        if self.is_horizontal:
            x = min(self.node1.x, self.node2.x)
            y = min(self.node1.y, self.node2.y) - self.width / 2
            return QRectF(x, y, abs(self.node1.x - self.node2.x), self.width)
        else:
            x = min(self.node1.x, self.node2.x) - self.width / 2
            y = min(self.node1.y, self.node2.y)
            return QRectF(x, y, self.width, abs(self.node1.y - self.node2.y))

class VoltageSource(CircuitElement):
    """Represents a voltage source, drawn as a circle with +/- signs."""
    def paint(self, painter, option, widget):
        try:
            painter.setPen(QPen(Qt.black, 2))
            mid_x = int((self.node1.x + self.node2.x) / 2)
            mid_y = int((self.node1.y + self.node2.y) / 2)
            if self.is_horizontal:
                painter.drawEllipse(mid_x - 10, mid_y - 10, 20, 20)
                painter.drawLine(mid_x - 5, mid_y, mid_x + 5, mid_y)  # Plus
                painter.drawLine(mid_x, mid_y - 5, mid_x, mid_y + 5)
                painter.drawLine(mid_x + 15, mid_y, mid_x + 25, mid_y)  # Minus
                painter.drawLine(int(self.node1.x), int(self.node1.y), mid_x - 10, mid_y)
                painter.drawLine(mid_x + 25, mid_y, int(self.node2.x), int(self.node2.y))
            else:
                painter.drawEllipse(mid_x - 10, mid_y - 10, 20, 20)
                painter.drawLine(mid_x - 5, mid_y, mid_x + 5, mid_y)
                painter.drawLine(mid_x, mid_y - 5, mid_x, mid_y + 5)
                painter.drawLine(mid_x, mid_y + 15, mid_x, mid_y + 25)
                painter.drawLine(int(self.node1.x), int(self.node1.y), mid_x, mid_y - 10)
                painter.drawLine(mid_x, mid_y + 25, int(self.node2.x), int(self.node2.y))
        except Exception as e:
            print(f"Error drawing VoltageSource: {e}")
            raise

class Inductor(CircuitElement):
    """Represents an inductor, drawn as a series of arcs."""
    def paint(self, painter, option, widget):
        try:
            painter.setPen(QPen(Qt.black, 2))
            mid_x = int((self.node1.x + self.node2.x) / 2)
            mid_y = int((self.node1.y + self.node2.y) / 2)
            if self.is_horizontal:
                for i in range(3):
                    x = mid_x - 15 + i * 10
                    painter.drawArc(x, mid_y - 5, 10, 10, 0, 180 * 16)
                painter.drawLine(int(self.node1.x), int(self.node1.y), mid_x - 15, mid_y)
                painter.drawLine(mid_x + 15, mid_y, int(self.node2.x), int(self.node2.y))
            else:
                for i in range(3):
                    y = mid_y - 15 + i * 10
                    painter.drawArc(mid_x - 5, y, 10, 10, 90 * 16, 180 * 16)
                painter.drawLine(int(self.node1.x), int(self.node1.y), mid_x, mid_y - 15)
                painter.drawLine(mid_x, mid_y + 15, int(self.node2.x), int(self.node2.y))
        except Exception as e:
            print(f"Error drawing Inductor: {e}")
            raise

class Resistor(CircuitElement):
    """Represents a resistor, drawn as a rectangle."""
    def paint(self, painter, option, widget):
        try:
            painter.setPen(QPen(Qt.black, 2))
            mid_x = int((self.node1.x + self.node2.x) / 2)
            mid_y = int((self.node1.y + self.node2.y) / 2)
            if self.is_horizontal:
                painter.drawRect(mid_x - 15, mid_y - 5, 30, 10)
                painter.drawLine(int(self.node1.x), int(self.node1.y), mid_x - 15, mid_y)
                painter.drawLine(mid_x + 15, mid_y, int(self.node2.x), int(self.node2.y))
            else:
                painter.drawRect(mid_x - 5, mid_y - 15, 10, 30)
                painter.drawLine(int(self.node1.x), int(self.node1.y), mid_x, mid_y - 15)
                painter.drawLine(mid_x, mid_y + 15, int(self.node2.x), int(self.node2.y))
        except Exception as e:
            print(f"Error drawing Resistor: {e}")
            raise

class Capacitor(CircuitElement):
    """Represents a capacitor, drawn as parallel lines."""
    def paint(self, painter, option, widget):
        try:
            painter.setPen(QPen(Qt.black, 2))
            mid_x = int((self.node1.x + self.node2.x) / 2)
            mid_y = int((self.node1.y + self.node2.y) / 2)
            if self.is_horizontal:
                painter.drawLine(mid_x - 10, mid_y - 10, mid_x - 10, mid_y + 10)
                painter.drawLine(mid_x + 10, mid_y - 10, mid_x + 10, mid_y + 10)
                painter.drawLine(int(self.node1.x), int(self.node1.y), mid_x - 10, mid_y)
                painter.drawLine(mid_x + 10, mid_y, int(self.node2.x), int(self.node2.y))
            else:
                painter.drawLine(mid_x - 10, mid_y - 10, mid_x + 10, mid_y - 10)
                painter.drawLine(mid_x - 10, mid_y + 10, mid_x + 10, mid_y + 10)
                painter.drawLine(int(self.node1.x), int(self.node1.y), mid_x, mid_y - 10)
                painter.drawLine(mid_x, mid_y + 10, int(self.node2.x), int(self.node2.y))
        except Exception as e:
            print(f"Error drawing Capacitor: {e}")
            raise

class Wire(QGraphicsItem):
    """Represents a wire, drawn as a straight line between two nodes."""
    def __init__(self, node1, node2):
        super().__init__()
        self.node1 = node1
        self.node2 = node2

    def boundingRect(self):
        """Defines the bounding rectangle for the wire."""
        x = min(self.node1.x, self.node2.x)
        y = min(self.node1.y, self.node2.y)
        w = abs(self.node1.x - self.node2.x)
        h = abs(self.node1.y - self.node2.y)
        return QRectF(x, y, w, h)

    def paint(self, painter, option, widget):
        try:
            painter.setPen(QPen(Qt.black, 2))
            painter.drawLine(int(self.node1.x), int(self.node1.y), int(self.node2.x), int(self.node2.y))
        except Exception as e:
            print(f"Error drawing Wire: {e}")
            raise