#!/usr/bin/env python3
"""
COLOSS GUI - Modern interface for COLOSS coupled-channel calculations
"""

import sys
from PySide6.QtWidgets import QApplication
from PySide6.QtCore import Qt
from main_window import MainWindow


def main():
    """Main entry point for the COLOSS GUI application"""
    # Enable high DPI scaling
    QApplication.setHighDpiScaleFactorRoundingPolicy(
        Qt.HighDpiScaleFactorRoundingPolicy.PassThrough
    )

    app = QApplication(sys.argv)
    app.setApplicationName("COLOSS")
    app.setOrganizationName("COLOSS Team")

    # Set application-wide font
    from PySide6.QtGui import QFont
    font = QFont("SF Pro Display", 10)
    app.setFont(font)

    window = MainWindow()
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()
