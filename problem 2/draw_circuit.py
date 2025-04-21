# draw_circuit.py
import sys
import os

# macOS-specific settings to improve rendering
os.environ["QT_MAC_WANTS_LAYER"] = "1"
os.environ["QT_OPENGL"] = "software"

from PyQt5.QtWidgets import QApplication, QMainWindow, QGraphicsView, QGraphicsScene, QGraphicsTextItem
from PyQt5.QtGui import QPen, QFont
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
            self.setCentralWidget(self.view)
            print("Graphics View and Scene setup complete.")

            # Draw the circuit
            self.draw_circuit()

            # Adjust the scene rect and view to fit the entire circuit
            bounding_rect = self.scene.itemsBoundingRect()
            print(
                f"Items bounding rect: x={bounding_rect.x()}, y={bounding_rect.y()}, width={bounding_rect.width()}, height={bounding_rect.height()}")
            # Add padding around the circuit for better visibility
            padding = 30
            self.scene.setSceneRect(
                bounding_rect.x() - padding,
                bounding_rect.y() - padding,
                bounding_rect.width() + 2 * padding,
                bounding_rect.height() + 2 * padding
            )
            self.view.fitInView(self.scene.sceneRect(), Qt.KeepAspectRatio)
            # Increase the zoom level to make the circuit larger
            self.view.scale(12.0, 12.0)
            print("Scene rect and view adjusted.")

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
                # Add node label with adjusted position
                label = QGraphicsTextItem(f"N{node.id}")
                label.setFont(QFont("Arial", 10))  # Increased font size
                # Adjust label position based on node location to avoid overlap
                if node.id == "1":
                    label.setPos(int(node.x) - 25, int(node.y) - 15)  # Top-left
                elif node.id == "2":
                    label.setPos(int(node.x) + 10, int(node.y) - 15)  # Top-right
                elif node.id == "3":
                    label.setPos(int(node.x) + 10, int(node.y) + 5)  # Bottom-right
                elif node.id == "4":
                    label.setPos(int(node.x) - 25, int(node.y) + 5)  # Bottom-left
                self.scene.addItem(label)

            print("Circuit drawing complete.")
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
    window.repaint()
    sys.exit(app.exec_())