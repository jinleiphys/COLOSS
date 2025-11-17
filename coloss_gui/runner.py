"""
COLOSS runner for executing calculations
"""

from PySide6.QtCore import QObject, QProcess, Signal
import os
import shutil


class ColossRunner(QObject):
    """Class to run COLOSS calculations"""

    output_ready = Signal(str)
    error_ready = Signal(str)
    finished = Signal(int)
    started = Signal()

    def __init__(self):
        super().__init__()
        self.process = None

    def run(self, executable, input_file, working_dir=None, pot_file_path=None):
        """Run COLOSS with the given input file in specified working directory

        Args:
            executable: Path to COLOSS executable
            input_file: Path to input file
            working_dir: Working directory for execution
            pot_file_path: Path to external potential file (optional)
        """
        if self.process is not None and self.process.state() == QProcess.Running:
            self.error_ready.emit("A calculation is already running")
            return

        self.process = QProcess()
        self.process.setProcessChannelMode(QProcess.MergedChannels)

        # Connect signals
        self.process.readyReadStandardOutput.connect(self.handle_output)
        self.process.readyReadStandardError.connect(self.handle_error)
        self.process.finished.connect(self.handle_finished)
        self.process.started.connect(self.handle_started)

        # Set working directory (use specified dir or default to input file's dir)
        if working_dir is None:
            working_dir = os.path.dirname(input_file)
        self.process.setWorkingDirectory(working_dir)

        # Copy external potential file to working directory if provided
        if pot_file_path and os.path.exists(pot_file_path):
            try:
                dest_path = os.path.join(working_dir, os.path.basename(pot_file_path))
                shutil.copy2(pot_file_path, dest_path)
                self.output_ready.emit(f"Copied potential file: {os.path.basename(pot_file_path)}")
            except Exception as e:
                self.error_ready.emit(f"Error copying potential file: {str(e)}")
                return

        # Start the process with input redirection
        self.process.start(executable, [])

        # Write input file content to stdin
        if os.path.exists(input_file):
            with open(input_file, 'r') as f:
                input_content = f.read()
            self.process.write(input_content.encode())
            self.process.closeWriteChannel()

    def stop(self):
        """Stop the running calculation"""
        if self.process is not None and self.process.state() == QProcess.Running:
            self.process.kill()

    def handle_output(self):
        """Handle standard output from the process"""
        if self.process:
            data = self.process.readAllStandardOutput()
            text = bytes(data).decode('utf-8', errors='replace')
            self.output_ready.emit(text.strip())

    def handle_error(self):
        """Handle standard error from the process"""
        if self.process:
            data = self.process.readAllStandardError()
            text = bytes(data).decode('utf-8', errors='replace')
            self.error_ready.emit(text.strip())

    def handle_finished(self, exit_code, exit_status):
        """Handle process completion"""
        self.finished.emit(exit_code)

    def handle_started(self):
        """Handle process start"""
        self.started.emit()
