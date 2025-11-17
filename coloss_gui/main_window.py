"""
Main window for COLOSS GUI application
"""

from PySide6.QtWidgets import (
    QMainWindow, QWidget, QVBoxLayout, QHBoxLayout, QSplitter,
    QTabWidget, QToolBar, QStatusBar, QFileDialog, QMessageBox
)
from PySide6.QtCore import Qt, QTimer
from PySide6.QtGui import QAction, QKeySequence
import os

from input_panel import InputPanel
from plot_widget import PlotWidget
from log_widget import LogWidget
from runner import ColossRunner
from styles import apply_modern_style
from path_utils import get_repo_root, find_executable, get_default_test_directory, get_executable_info


class MainWindow(QMainWindow):
    """Main application window with modern layout"""

    def __init__(self):
        super().__init__()
        self.coloss_runner = ColossRunner()
        self.current_file = None
        self.working_directory = None  # Store current working directory
        self.is_running = False  # Track whether COLOSS is running

        # Detect repository root and default paths
        self.repo_root = get_repo_root()
        self.default_test_dir = get_default_test_directory(self.repo_root)

        self.init_ui()
        self.setup_connections()

    def init_ui(self):
        """Initialize the user interface"""
        self.setWindowTitle("COLOSS - Complex-scaled Optical and couLOmb Scattering Solver")
        self.setGeometry(100, 100, 1600, 900)

        # Apply modern styling (default: light theme, change to "dark" for dark theme)
        apply_modern_style(self)

        # Create central widget with splitter layout
        central_widget = QWidget()
        self.setCentralWidget(central_widget)

        main_layout = QHBoxLayout(central_widget)
        main_layout.setContentsMargins(0, 0, 0, 0)

        # Create main splitter (horizontal)
        splitter = QSplitter(Qt.Horizontal)

        # Left panel: Input forms
        self.input_panel = InputPanel()
        splitter.addWidget(self.input_panel)

        # Right panel: Tabs for plotting and output
        right_widget = QWidget()
        right_layout = QVBoxLayout(right_widget)
        right_layout.setContentsMargins(0, 0, 0, 0)

        self.right_tabs = QTabWidget()
        self.right_tabs.setDocumentMode(True)

        # Plot tab
        self.plot_widget = PlotWidget()
        self.right_tabs.addTab(self.plot_widget, "Plot")

        # Log tab
        self.log_widget = LogWidget()
        self.right_tabs.addTab(self.log_widget, "Output Log")

        right_layout.addWidget(self.right_tabs)
        splitter.addWidget(right_widget)

        # Set initial splitter sizes (40% input, 60% output/plot)
        splitter.setSizes([640, 960])

        main_layout.addWidget(splitter)

        # Create status bar
        self.status_bar = QStatusBar()
        self.setStatusBar(self.status_bar)
        self.status_bar.showMessage("Ready")

        # Create menu bar and toolbar
        self.create_menus()
        self.create_toolbar()

    def create_menus(self):
        """Create menu bar with all menu items"""
        menubar = self.menuBar()

        # File menu
        file_menu = menubar.addMenu("&File")

        new_action = QAction("&New", self)
        new_action.setShortcut(QKeySequence.New)
        new_action.triggered.connect(self.new_file)
        file_menu.addAction(new_action)

        open_action = QAction("&Open...", self)
        open_action.setShortcut(QKeySequence.Open)
        open_action.triggered.connect(self.open_file)
        file_menu.addAction(open_action)

        save_action = QAction("&Save", self)
        save_action.setShortcut(QKeySequence.Save)
        save_action.triggered.connect(self.save_file)
        file_menu.addAction(save_action)

        save_as_action = QAction("Save &As...", self)
        save_as_action.setShortcut(QKeySequence.SaveAs)
        save_as_action.triggered.connect(self.save_file_as)
        file_menu.addAction(save_as_action)

        file_menu.addSeparator()

        quit_action = QAction("&Quit", self)
        quit_action.setShortcut(QKeySequence.Quit)
        quit_action.triggered.connect(self.close)
        file_menu.addAction(quit_action)

        # Run menu
        run_menu = menubar.addMenu("&Run")

        run_action = QAction("&Run COLOSS", self)
        run_action.setShortcut("Ctrl+R")
        run_action.triggered.connect(self.run_coloss)
        run_menu.addAction(run_action)

        stop_action = QAction("&Stop", self)
        stop_action.setShortcut("Ctrl+C")
        stop_action.triggered.connect(self.stop_coloss)
        run_menu.addAction(stop_action)

        # View menu
        view_menu = menubar.addMenu("&View")

        light_theme_action = QAction("Light Theme", self)
        light_theme_action.triggered.connect(lambda: apply_modern_style(self, "light"))
        view_menu.addAction(light_theme_action)

        dark_theme_action = QAction("Dark Theme", self)
        dark_theme_action.triggered.connect(lambda: apply_modern_style(self, "dark"))
        view_menu.addAction(dark_theme_action)

        # Help menu
        help_menu = menubar.addMenu("&Help")

        about_action = QAction("&About", self)
        about_action.triggered.connect(self.show_about)
        help_menu.addAction(about_action)

    def create_toolbar(self):
        """Create toolbar with common actions"""
        toolbar = QToolBar("Main Toolbar")
        toolbar.setMovable(False)
        self.addToolBar(toolbar)

        # Run button
        run_action = QAction("Run", self)
        run_action.triggered.connect(self.run_coloss)
        toolbar.addAction(run_action)

        # Stop button
        stop_action = QAction("Stop", self)
        stop_action.triggered.connect(self.stop_coloss)
        toolbar.addAction(stop_action)

        toolbar.addSeparator()

        # Open button
        open_action = QAction("Open", self)
        open_action.triggered.connect(self.open_file)
        toolbar.addAction(open_action)

        # Save button
        save_action = QAction("Save", self)
        save_action.triggered.connect(self.save_file)
        toolbar.addAction(save_action)

    def setup_connections(self):
        """Setup signal-slot connections"""
        # Connect runner signals
        self.coloss_runner.output_ready.connect(self.handle_output)
        self.coloss_runner.error_ready.connect(self.handle_error)
        self.coloss_runner.finished.connect(self.handle_finished)
        self.coloss_runner.started.connect(self.handle_started)

    def new_file(self):
        """Create a new input file"""
        self.current_file = None
        self.status_bar.showMessage("New file created")

    def open_file(self):
        """Open an existing input file"""
        file_path, _ = QFileDialog.getOpenFileName(
            self,
            "Open COLOSS Input File",
            self.default_test_dir,
            "Input Files (*.in);;All Files (*)"
        )

        if file_path:
            try:
                with open(file_path, 'r') as f:
                    content = f.read()

                # Parse and load into input panel (simplified - would need full parser)
                self.current_file = file_path
                self.working_directory = os.path.dirname(file_path)
                self.status_bar.showMessage(f"Opened: {file_path}")

                QMessageBox.information(
                    self,
                    "File Opened",
                    f"Loaded input file:\n{os.path.basename(file_path)}\n\n"
                    f"Note: Use 'Generate Input File' to create new input from current parameters."
                )

            except Exception as e:
                QMessageBox.critical(self, "Error", f"Failed to open file:\n{str(e)}")

    def save_file(self):
        """Save the current input file"""
        if self.current_file:
            self.save_to_file(self.current_file)
        else:
            self.save_file_as()

    def save_file_as(self):
        """Save the input file with a new name"""
        file_path, _ = QFileDialog.getSaveFileName(
            self,
            "Save COLOSS Input File",
            self.default_test_dir,
            "Input Files (*.in);;All Files (*)"
        )

        if file_path:
            self.save_to_file(file_path)
            self.current_file = file_path

    def save_to_file(self, file_path):
        """Save input to specified file"""
        try:
            input_text = self.input_panel.get_input_text()
            with open(file_path, 'w') as f:
                f.write(input_text)

            self.status_bar.showMessage(f"Saved: {file_path}")
            QMessageBox.information(self, "Success", f"Input file saved:\n{file_path}")

        except Exception as e:
            QMessageBox.critical(self, "Error", f"Failed to save file:\n{str(e)}")

    def run_coloss(self):
        """Run COLOSS calculation"""
        if self.is_running:
            QMessageBox.warning(self, "Warning", "A calculation is already running!")
            return

        # Find COLOSS executable
        exe_path, message, found = get_executable_info("COLOSS", self.repo_root)

        if not found:
            QMessageBox.critical(self, "COLOSS Not Found", message)
            return

        # Save input to temporary file if not already saved
        if not self.current_file:
            temp_dir = self.default_test_dir
            temp_file = os.path.join(temp_dir, "temp_input.in")
            self.save_to_file(temp_file)
            self.current_file = temp_file

        # Set working directory
        if self.working_directory is None:
            self.working_directory = os.path.dirname(self.current_file)

        # Get potential file path if external potential is being used
        pot_file_path = self.input_panel.get_pot_file_path()

        # Clear log
        self.log_widget.clear()
        self.log_widget.append_info(f"Running COLOSS: {exe_path}")
        self.log_widget.append_info(f"Input file: {self.current_file}")
        self.log_widget.append_info(f"Working directory: {self.working_directory}")
        if pot_file_path:
            self.log_widget.append_info(f"External potential: {pot_file_path}")
        self.log_widget.append_info("-" * 50)

        # Switch to log tab
        self.right_tabs.setCurrentWidget(self.log_widget)

        # Run COLOSS
        self.coloss_runner.run(exe_path, self.current_file, self.working_directory, pot_file_path)

    def stop_coloss(self):
        """Stop the running calculation"""
        if self.is_running:
            self.coloss_runner.stop()
            self.log_widget.append_warning("Calculation stopped by user")
            self.is_running = False
            self.status_bar.showMessage("Stopped")

    def handle_started(self):
        """Handle calculation start"""
        self.is_running = True
        self.status_bar.showMessage("Running COLOSS...")

    def handle_output(self, text):
        """Handle output from COLOSS"""
        if text:
            self.log_widget.append_output(text)

    def handle_error(self, text):
        """Handle error output from COLOSS"""
        if text:
            self.log_widget.append_error(text)

    def handle_finished(self, exit_code):
        """Handle calculation completion"""
        self.is_running = False

        if exit_code == 0:
            self.log_widget.append_success("Calculation completed successfully!")
            self.status_bar.showMessage("Calculation completed")

            # Load results into plot widget
            if self.working_directory:
                self.plot_widget.load_results(self.working_directory)
                self.right_tabs.setCurrentWidget(self.plot_widget)
        else:
            self.log_widget.append_error(f"Calculation failed with exit code {exit_code}")
            self.status_bar.showMessage(f"Failed with exit code {exit_code}")

    def show_about(self):
        """Show about dialog"""
        QMessageBox.about(
            self,
            "About COLOSS GUI",
            "<h3>COLOSS GUI</h3>"
            "<p>A modern graphical interface for COLOSS coupled-channel calculations.</p>"
            "<p>Version 1.0.0</p>"
            "<p>COLOSS: Coupled-channel calculations for nuclear reactions</p>"
        )
