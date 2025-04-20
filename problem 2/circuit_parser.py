# circuit_parser.py
import sys
from circuit_elements import Node, VoltageSource, Inductor, Resistor, Capacitor, Wire

def parse_circuit_file(filename):
    """Parses a circuit description file and returns the title, nodes, elements, and wires."""
    nodes = {}
    elements = []
    wires = []
    title = "Circuit Diagram"

    try:
        with open(filename, 'r') as file:
            for line in file:
                line = line.strip()
                if not line or line.startswith('#'):
                    continue
                parts = [part.strip().strip('"') for part in line.split(',')]
                keyword = parts[0].upper()

                if keyword == "TITLE":
                    title = parts[1]
                    print(f"Title set to: {title}")
                elif keyword == "NODE":
                    if len(parts) != 4:
                        raise ValueError(f"Invalid NODE format: {line}")
                    node_id, x, y = parts[1], parts[2], parts[3]
                    nodes[node_id] = Node(node_id, x, y)
                    print(f"Added node {node_id} at ({x}, {y})")
                elif keyword == "ELEMENT":
                    if len(parts) != 4:
                        raise ValueError(f"Invalid ELEMENT format: {line}")
                    element_type, node1_id, node2_id = parts[1].upper(), parts[2], parts[3]
                    if node1_id not in nodes or node2_id not in nodes:
                        raise ValueError(f"Node not found in ELEMENT: {line}")
                    node1 = nodes[node1_id]
                    node2 = nodes[node2_id]
                    if element_type == "VOLTAGE_SOURCE":
                        elements.append(VoltageSource(node1, node2))
                        print(f"Added VoltageSource between nodes {node1_id} and {node2_id}")
                    elif element_type == "INDUCTOR":
                        elements.append(Inductor(node1, node2))
                        print(f"Added Inductor between nodes {node1_id} and {node2_id}")
                    elif element_type == "RESISTOR":
                        elements.append(Resistor(node1, node2))
                        print(f"Added Resistor between nodes {node1_id} and {node2_id}")
                    elif element_type == "CAPACITOR":
                        elements.append(Capacitor(node1, node2))
                        print(f"Added Capacitor between nodes {node1_id} and {node2_id}")
                elif keyword == "WIRE":
                    if len(parts) != 3:
                        raise ValueError(f"Invalid WIRE format: {line}")
                    node1_id, node2_id = parts[1], parts[2]
                    if node1_id not in nodes or node2_id not in nodes:
                        raise ValueError(f"Node not found in WIRE: {line}")
                    node1 = nodes[node1_id]
                    node2 = nodes[node2_id]
                    wires.append(Wire(node1, node2))
                    print(f"Added Wire between nodes {node1_id} and {node2_id}")
    except FileNotFoundError:
        print(f"Error: {filename} not found.")
        sys.exit(1)
    except Exception as e:
        print(f"Error parsing file: {e}")
        sys.exit(1)

    return title, nodes, elements, wires