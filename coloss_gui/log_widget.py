"""
Log widget for displaying COLOSS output with syntax highlighting
"""

from PySide6.QtWidgets import QTextEdit, QWidget, QVBoxLayout
from PySide6.QtGui import QTextCharFormat, QColor, QFont, QTextCursor
from PySide6.QtCore import Qt, Signal


class LogWidget(QWidget):
    """Widget for displaying colored log output"""

    def __init__(self):
        super().__init__()
        self.init_ui()

    def init_ui(self):
        """Initialize the log widget"""
        layout = QVBoxLayout(self)
        layout.setContentsMargins(0, 0, 0, 0)

        self.text_edit = QTextEdit()
        self.text_edit.setReadOnly(True)
        self.text_edit.setLineWrapMode(QTextEdit.NoWrap)

        # Use monospace font
        font = QFont("Menlo", 10)
        if not font.exactMatch():
            font = QFont("Monaco", 10)
        if not font.exactMatch():
            font = QFont("Courier New", 10)
        self.text_edit.setFont(font)

        layout.addWidget(self.text_edit)

        # Define text formats for different log levels
        self.format_normal = QTextCharFormat()
        self.format_normal.setForeground(QColor("#333333"))

        self.format_info = QTextCharFormat()
        self.format_info.setForeground(QColor("#0066CC"))
        self.format_info.setFontWeight(QFont.Bold)

        self.format_warning = QTextCharFormat()
        self.format_warning.setForeground(QColor("#FF8800"))
        self.format_warning.setFontWeight(QFont.Bold)

        self.format_error = QTextCharFormat()
        self.format_error.setForeground(QColor("#CC0000"))
        self.format_error.setFontWeight(QFont.Bold)

        self.format_success = QTextCharFormat()
        self.format_success.setForeground(QColor("#009900"))
        self.format_success.setFontWeight(QFont.Bold)

    def append_output(self, text):
        """Append normal output text"""
        self._append_text(text, self.format_normal)

    def append_info(self, text):
        """Append info message"""
        self._append_text(f"[INFO] {text}", self.format_info)

    def append_warning(self, text):
        """Append warning message"""
        self._append_text(f"[WARNING] {text}", self.format_warning)

    def append_error(self, text):
        """Append error message"""
        self._append_text(f"[ERROR] {text}", self.format_error)

    def append_success(self, text):
        """Append success message"""
        self._append_text(f"[SUCCESS] {text}", self.format_success)

    def _append_text(self, text, text_format):
        """Append text with specific format"""
        cursor = self.text_edit.textCursor()
        cursor.movePosition(QTextCursor.End)
        cursor.insertText(text + '\n', text_format)
        self.text_edit.setTextCursor(cursor)
        self.text_edit.ensureCursorVisible()

    def clear(self):
        """Clear the log"""
        self.text_edit.clear()

    def get_text(self):
        """Get all log text"""
        return self.text_edit.toPlainText()
