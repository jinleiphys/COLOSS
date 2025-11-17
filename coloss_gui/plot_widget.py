"""
Plot widget for displaying COLOSS results using matplotlib
"""

from PySide6.QtWidgets import QWidget, QVBoxLayout, QHBoxLayout, QPushButton, QComboBox, QLabel
from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT as NavigationToolbar
from matplotlib.figure import Figure
import matplotlib.pyplot as plt
import numpy as np
import os


class PlotWidget(QWidget):
    """Widget for plotting COLOSS results"""

    def __init__(self):
        super().__init__()
        self.current_data = {}
        self.working_directory = None  # Store working directory for refresh
        self.init_ui()

    def init_ui(self):
        """Initialize the plot widget"""
        layout = QVBoxLayout(self)

        # Control panel
        control_layout = QHBoxLayout()

        self.plot_type = QComboBox()
        self.plot_type.addItems([
            "Angular Distribution (fort.67)",
            "S-Matrix Elements (fort.60)",
            "Scattering Amplitude (fort.61)",
            "All Available Data"
        ])
        self.plot_type.currentIndexChanged.connect(self.update_plot)
        control_layout.addWidget(QLabel("Plot Type:"))
        control_layout.addWidget(self.plot_type)

        control_layout.addStretch()

        refresh_btn = QPushButton("Refresh")
        refresh_btn.clicked.connect(self.load_results)
        control_layout.addWidget(refresh_btn)

        clear_btn = QPushButton("Clear")
        clear_btn.clicked.connect(self.clear_plot)
        control_layout.addWidget(clear_btn)

        layout.addLayout(control_layout)

        # Matplotlib figure
        self.figure = Figure(figsize=(8, 6), dpi=100)
        self.canvas = FigureCanvas(self.figure)
        self.toolbar = NavigationToolbar(self.canvas, self)

        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        # Initial plot
        self.clear_plot()

    def load_results(self, working_dir=None):
        """Load results from COLOSS output files"""
        try:
            # Store working directory for future refreshes
            if working_dir is not None:
                self.working_directory = working_dir

            # Use stored working directory if available
            if self.working_directory is None:
                self.working_directory = os.getcwd()

            # Clear previous data
            self.current_data.clear()

            # Load specific COLOSS output files (fort.* files)
            fort_files = {
                'fort.60': 'S-matrix',
                'fort.61': 'Scattering amplitude',
                'fort.67': 'Angular distribution',
                'fort.10': 'General output',
            }

            for fort_file, description in fort_files.items():
                filepath = os.path.join(self.working_directory, fort_file)
                if os.path.exists(filepath):
                    try:
                        data = self.read_fort_file(filepath)
                        if data is not None and len(data) > 0:
                            self.current_data[fort_file] = {
                                'data': data,
                                'description': description
                            }
                    except Exception as e:
                        print(f"Error reading {fort_file}: {e}")

            # Also look for .dat files
            for filename in os.listdir(self.working_directory):
                if filename.endswith('.dat') and not filename.startswith('.'):
                    filepath = os.path.join(self.working_directory, filename)
                    try:
                        # Check if file is text (not binary)
                        with open(filepath, 'r') as f:
                            f.read(100)  # Try to read some text
                        data = np.loadtxt(filepath)
                        if len(data) > 0:
                            self.current_data[filename] = {
                                'data': data,
                                'description': filename
                            }
                    except:
                        pass  # Skip binary or unreadable files

            self.update_plot()

        except Exception as e:
            print(f"Error loading results: {e}")

    def read_fort_file(self, filepath):
        """Read a COLOSS fort.* file, skipping header lines"""
        try:
            # Read file and skip lines starting with '&'
            data_lines = []
            with open(filepath, 'r') as f:
                for line in f:
                    line = line.strip()
                    if line and not line.startswith('&'):
                        # Try to extract numerical data
                        parts = line.split()
                        # Find numerical parts (skip text like "(L S J):")
                        nums = []
                        for part in parts:
                            try:
                                nums.append(float(part))
                            except ValueError:
                                continue
                        if nums:
                            data_lines.append(nums)

            if data_lines:
                # Convert to numpy array, padding with NaN if rows have different lengths
                max_len = max(len(row) for row in data_lines)
                padded_data = []
                for row in data_lines:
                    padded_row = row + [np.nan] * (max_len - len(row))
                    padded_data.append(padded_row)
                return np.array(padded_data)
            else:
                return None
        except Exception as e:
            print(f"Error reading {filepath}: {e}")
            return None

    def update_plot(self):
        """Update the plot based on selected type"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)

        plot_type_index = self.plot_type.currentIndex()

        if not self.current_data:
            ax.text(0.5, 0.5, 'No data loaded\nRun a calculation to see results',
                   ha='center', va='center', fontsize=12, color='gray')
            self.canvas.draw()
            return

        try:
            if plot_type_index == 0:  # Angular Distribution (fort.67)
                if 'fort.67' in self.current_data:
                    data = self.current_data['fort.67']['data']
                    ax.plot(data[:, 0], data[:, 1], 'b-', linewidth=2)
                    ax.set_xlabel('Angle (degrees)', fontsize=12)
                    ax.set_ylabel('Cross Section (mb)', fontsize=12)
                    ax.set_title('Angular Distribution', fontsize=14, fontweight='bold')
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, 'fort.67 not found\nRun a calculation first',
                           ha='center', va='center', fontsize=12, color='orange')

            elif plot_type_index == 1:  # S-Matrix Elements (fort.60)
                if 'fort.60' in self.current_data:
                    data = self.current_data['fort.60']['data']
                    # Data format: Real(S), Imag(S), L, S, J
                    L_values = data[:, 2]
                    real_S = data[:, 0]
                    imag_S = data[:, 1]

                    ax.plot(L_values, real_S, 'bo-', label='Real(S)', markersize=4)
                    ax.plot(L_values, imag_S, 'ro-', label='Imag(S)', markersize=4)
                    ax.set_xlabel('Angular Momentum L', fontsize=12)
                    ax.set_ylabel('S-Matrix Element', fontsize=12)
                    ax.set_title('S-Matrix Elements', fontsize=14, fontweight='bold')
                    ax.legend()
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, 'fort.60 not found\nRun a calculation first',
                           ha='center', va='center', fontsize=12, color='orange')

            elif plot_type_index == 2:  # Scattering Amplitude (fort.61)
                if 'fort.61' in self.current_data:
                    data = self.current_data['fort.61']['data']
                    # Data format: Real(f), Imag(f), L, S, J
                    L_values = data[:, 2]
                    real_f = data[:, 0]
                    imag_f = data[:, 1]

                    ax.plot(L_values, real_f, 'go-', label='Real(f)', markersize=4)
                    ax.plot(L_values, imag_f, 'mo-', label='Imag(f)', markersize=4)
                    ax.set_xlabel('Angular Momentum L', fontsize=12)
                    ax.set_ylabel('Scattering Amplitude', fontsize=12)
                    ax.set_title('Scattering Amplitude', fontsize=14, fontweight='bold')
                    ax.legend()
                    ax.grid(True, alpha=0.3)
                else:
                    ax.text(0.5, 0.5, 'fort.61 not found\nRun a calculation first',
                           ha='center', va='center', fontsize=12, color='orange')

            elif plot_type_index == 3:  # All Available Data
                if not self.current_data:
                    ax.text(0.5, 0.5, 'No data available',
                           ha='center', va='center', fontsize=12, color='gray')
                else:
                    # List all available data files
                    info_text = "Available data files:\n\n"
                    for filename, file_info in self.current_data.items():
                        desc = file_info.get('description', filename)
                        data = file_info['data']
                        info_text += f"• {filename}: {desc}\n"
                        info_text += f"  Shape: {data.shape}\n\n"

                    ax.text(0.5, 0.5, info_text,
                           ha='center', va='center', fontsize=10,
                           family='monospace', verticalalignment='center')
                    ax.set_title('Available Data Files', fontsize=14, fontweight='bold')
                    ax.axis('off')

        except Exception as e:
            ax.text(0.5, 0.5, f'Error plotting data:\n{str(e)}',
                   ha='center', va='center', fontsize=10, color='red')

        self.canvas.draw()

    def clear_plot(self):
        """Clear the current plot"""
        self.figure.clear()
        ax = self.figure.add_subplot(111)
        ax.text(0.5, 0.5, 'No plot data',
               ha='center', va='center', fontsize=14, color='gray')
        ax.set_xlim(0, 1)
        ax.set_ylim(0, 1)
        ax.axis('off')
        self.canvas.draw()
        self.current_data.clear()
