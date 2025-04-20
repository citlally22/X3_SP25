# draw_circuit.py
import sys
import os

# macOS-specific settings to improve rendering
os.environ["QT_MAC_WANTS_LAYER"] = "1"
os.environ["QT_OPENGL"] = "software"

from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsView, QGraphicsScene
from PyQt5.QtGui import QPen
from PyQt5.QtCore import Qt
from circuit_parser import parse_circuit_file


class CircuitWindow(QMainWindow):
    """A window to display the circuit diagram using Graphics View Framework."""

    def __init__(self):
        super().__init__()
        print("Starting CircuitWindow initialization...")

        try:
            # Parse the circuit file
            print("Parsing circuit file...")
            self.title, self.nodes, self.elements, self.wires = parse_circuit_file("circuit.txt")
            print(
                f"Parsed: Title={self.title}, Nodes={len(self.nodes)}, Elements={len(self.elements)}, Wires={len(self.wires)}")

            # Setup the window
            self.setWindowTitle(self.title)
            self.setGeometry(100, 100, 800, 600)
            print("Window setup complete.")

            # Setup Graphics View and Scene
            self.scene = QGraphicsScene()
            self.view = QGraphicsView(self.scene, self)
            # Set scene rect to encompass the circuit (nodes span x: 50-350, y: 150)
            self.scene.setSceneRect(0, 100, 400, 100)  # Adjusted to focus on the circuit area
            self.view.fitInView(0, 100, 400, 100, Qt.KeepAspectRatio)  # Ensure the circuit is visible
            self.setCentralWidget(self.view)
            print("Graphics View and Scene setup complete.")

            # Draw the circuit
            self.draw_circuit()
            print("Circuit drawing complete.")

        except Exception as e:
            print(f"Error during initialization: {e}")
            sys.exit(1)

    def draw_circuit(self):
        """Draws the circuit elements, wires, and nodes on the scene."""
        try:
            print("Drawing elements...")
            for element in self.elements:
                print(f"Adding element: {type(element).__name__}")
                self.scene.addItem(element)

            print("Drawing wires...")
            for wire in self.wires:
                print(f"Adding wire between nodes {wire.node1.id} and {wire.node2.id}")
                self.scene.addItem(wire)

            print("Drawing nodes...")
            for node in self.nodes.values():
                print(f"Adding node {node.id} at ({node.x}, {node.y})")
                self.scene.addEllipse(int(node.x) - 3, int(node.y) - 3, 6, 6, QPen(Qt.black), Qt.black)

        except Exception as e:
            print(f"Error during drawing: {e}")
            sys.exit(1)


if __name__ == "__main__":
    print("Starting application...")
    app = QApplication(sys.argv)
    window = CircuitWindow()
    print("Showing window...")
    window.show()
    print("Entering event loop...")
    # Force a repaint to ensure the scene is rendered
    window.repaint()
    sys.exit(app.exec_())