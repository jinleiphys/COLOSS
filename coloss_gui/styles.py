"""
Modern styling for the COLOSS GUI with beautiful gradients and modern colors
"""

LIGHT_THEME = """
/* Main Window - Soft gradient background */
QMainWindow {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #f0f4ff, stop:0.5 #e8f0fe, stop:1 #f0f8ff);
}

/* Base Widget Styling */
QWidget {
    background-color: #ffffff;
    color: #1a1a2e;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
    font-size: 13px;
}

/* Tab Widget - Modern glass-morphism style */
QTabWidget::pane {
    border: none;
    background-color: rgba(255, 255, 255, 0.95);
    border-radius: 16px;
    padding: 4px;
    margin-top: -1px;
}

QTabWidget::tab-bar {
    alignment: left;
}

/* Tab Bar - Beautiful gradient tabs */
QTabBar::tab {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                stop:0 #f8fafc, stop:1 #e2e8f0);
    color: #64748b;
    border: none;
    border-radius: 10px;
    padding: 12px 24px;
    margin-right: 6px;
    margin-top: 6px;
    font-weight: 600;
    min-width: 90px;
}

QTabBar::tab:selected {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #667eea, stop:0.5 #764ba2, stop:1 #f093fb);
    color: white;
    font-weight: 700;
}

QTabBar::tab:hover:!selected {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                stop:0 #e0e7ff, stop:1 #c7d2fe);
    color: #667eea;
}

/* Buttons - Vibrant gradient style */
QPushButton {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #667eea, stop:1 #764ba2);
    color: white;
    border: none;
    border-radius: 10px;
    padding: 12px 24px;
    font-weight: 700;
    font-size: 13px;
    min-height: 22px;
}

QPushButton:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #7c3aed, stop:1 #8b5cf6);
}

QPushButton:pressed {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #5b21b6, stop:1 #6d28d9);
    padding: 13px 23px 11px 25px;
}

QPushButton:disabled {
    background: #e2e8f0;
    color: #94a3b8;
}

/* Input Fields - Clean with subtle shadow */
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
    border: 2px solid #e2e8f0;
    border-radius: 10px;
    padding: 10px 14px;
    background-color: #f8fafc;
    selection-background-color: #667eea;
    selection-color: white;
    min-height: 26px;
    color: #1a1a2e;
}

QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
    border: 2px solid #667eea;
    background-color: #ffffff;
}

QLineEdit:hover, QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover {
    border: 2px solid #cbd5e1;
    background-color: #ffffff;
}

/* ComboBox - Modern dropdown */
QComboBox::drop-down {
    border: none;
    padding-right: 15px;
    width: 24px;
}

QComboBox::down-arrow {
    image: none;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 7px solid #64748b;
    margin-right: 6px;
}

QComboBox:hover::down-arrow {
    border-top: 7px solid #667eea;
}

QComboBox QAbstractItemView {
    border: 2px solid #e2e8f0;
    border-radius: 10px;
    background-color: #ffffff;
    selection-background-color: #667eea;
    selection-color: white;
    padding: 6px;
}

/* Text Edit - Code-style background */
QTextEdit {
    border: 2px solid #e2e8f0;
    border-radius: 10px;
    background-color: #f8fafc;
    padding: 12px;
    color: #1a1a2e;
    font-family: 'SF Mono', Monaco, 'Cascadia Code', 'Roboto Mono', Consolas, 'Courier New', monospace;
}

/* Group Box - Modern card with gradient title */
QGroupBox {
    border: 2px solid #e2e8f0;
    border-radius: 12px;
    margin-top: 20px;
    padding-top: 20px;
    font-weight: 600;
    background-color: #ffffff;
}

QGroupBox::title {
    subcontrol-origin: margin;
    left: 20px;
    padding: 6px 16px;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #667eea, stop:1 #764ba2);
    color: white;
    border-radius: 8px;
    font-size: 13px;
    font-weight: 800;
}

/* Scroll Area */
QScrollArea {
    border: none;
    background-color: transparent;
}

/* Scrollbar - Modern thin design */
QScrollBar:vertical {
    border: none;
    background-color: transparent;
    width: 12px;
    margin: 0px;
}

QScrollBar::handle:vertical {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #a5b4fc, stop:1 #c4b5fd);
    border-radius: 6px;
    min-height: 40px;
}

QScrollBar::handle:vertical:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #818cf8, stop:1 #a78bfa);
}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0px;
}

/* Status Bar - Gradient background */
QStatusBar {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #f8fafc, stop:1 #f1f5f9);
    color: #475569;
    border-top: 2px solid #e2e8f0;
    padding: 6px;
    font-weight: 500;
}

/* Menu Bar */
QMenuBar {
    background-color: #ffffff;
    color: #1a1a2e;
    border-bottom: 2px solid #e2e8f0;
    padding: 6px;
}

QMenuBar::item {
    padding: 10px 16px;
    background-color: transparent;
    border-radius: 8px;
    font-weight: 600;
}

QMenuBar::item:selected {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                stop:0 #e0e7ff, stop:1 #ddd6fe);
    color: #5b21b6;
}

/* Menu */
QMenu {
    background-color: #ffffff;
    border: 2px solid #e2e8f0;
    border-radius: 12px;
    padding: 8px;
}

QMenu::item {
    padding: 10px 28px;
    border-radius: 8px;
    font-weight: 500;
}

QMenu::item:selected {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #667eea, stop:1 #764ba2);
    color: white;
}

/* Toolbar */
QToolBar {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                stop:0 #ffffff, stop:1 #f8fafc);
    border: none;
    spacing: 10px;
    padding: 10px;
}

/* Checkbox - Modern style */
QCheckBox {
    spacing: 8px;
    color: #1a1a2e;
    font-weight: 500;
}

QCheckBox::indicator {
    width: 20px;
    height: 20px;
    border-radius: 6px;
    border: 2px solid #cbd5e1;
    background-color: #f8fafc;
}

QCheckBox::indicator:hover {
    border: 2px solid #667eea;
}

QCheckBox::indicator:checked {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #667eea, stop:1 #764ba2);
    border: 2px solid #667eea;
    image: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTIiIGhlaWdodD0iOSIgdmlld0JveD0iMCAwIDEyIDkiIGZpbGw9Im5vbmUiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+CjxwYXRoIGQ9Ik0xIDRMNC41IDdMMTEgMSIgc3Ryb2tlPSJ3aGl0ZSIgc3Ryb2tlLXdpZHRoPSIyIiBzdHJva2UtbGluZWNhcD0icm91bmQiLz4KPC9zdmc+Cg==);
}

/* Label - Enhanced readability */
QLabel {
    color: #1a1a2e;
    font-weight: 500;
}

/* Splitter */
QSplitter::handle {
    background-color: #e2e8f0;
    width: 2px;
}

QSplitter::handle:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #a5b4fc, stop:1 #c4b5fd);
}
"""

DARK_THEME = """
/* Main Window - Deep space gradient */
QMainWindow {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:1,
                                stop:0 #0f172a, stop:0.5 #1e1b4b, stop:1 #1e293b);
}

/* Base Widget Styling - Dark */
QWidget {
    background-color: #1e293b;
    color: #e2e8f0;
    font-family: -apple-system, BlinkMacSystemFont, "Segoe UI", "Helvetica Neue", Arial, sans-serif;
    font-size: 13px;
}

/* Tab Widget - Dark glass */
QTabWidget::pane {
    border: none;
    background-color: rgba(30, 41, 59, 0.95);
    border-radius: 16px;
    padding: 4px;
    margin-top: -1px;
}

QTabWidget::tab-bar {
    alignment: left;
}

/* Tab Bar - Neon gradient tabs */
QTabBar::tab {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                stop:0 #334155, stop:1 #1e293b);
    color: #94a3b8;
    border: none;
    border-radius: 10px;
    padding: 12px 24px;
    margin-right: 6px;
    margin-top: 6px;
    font-weight: 600;
    min-width: 90px;
}

QTabBar::tab:selected {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #8b5cf6, stop:0.5 #a78bfa, stop:1 #c4b5fd);
    color: white;
    font-weight: 700;
}

QTabBar::tab:hover:!selected {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                stop:0 #4c1d95, stop:1 #5b21b6);
    color: #c4b5fd;
}

/* Buttons - Vibrant neon gradient */
QPushButton {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #8b5cf6, stop:1 #a78bfa);
    color: white;
    border: none;
    border-radius: 10px;
    padding: 12px 24px;
    font-weight: 700;
    font-size: 13px;
    min-height: 22px;
}

QPushButton:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #a78bfa, stop:1 #c4b5fd);
}

QPushButton:pressed {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #6d28d9, stop:1 #7c3aed);
    padding: 13px 23px 11px 25px;
}

QPushButton:disabled {
    background: #334155;
    color: #64748b;
}

/* Input Fields - Dark with glow */
QLineEdit, QSpinBox, QDoubleSpinBox, QComboBox {
    border: 2px solid #334155;
    border-radius: 10px;
    padding: 10px 14px;
    background-color: #0f172a;
    selection-background-color: #8b5cf6;
    selection-color: white;
    min-height: 26px;
    color: #e2e8f0;
}

QLineEdit:focus, QSpinBox:focus, QDoubleSpinBox:focus, QComboBox:focus {
    border: 2px solid #8b5cf6;
    background-color: #1e293b;
}

QLineEdit:hover, QSpinBox:hover, QDoubleSpinBox:hover, QComboBox:hover {
    border: 2px solid #475569;
    background-color: #1e293b;
}

/* ComboBox - Dark dropdown */
QComboBox::drop-down {
    border: none;
    padding-right: 15px;
    width: 24px;
}

QComboBox::down-arrow {
    image: none;
    border-left: 5px solid transparent;
    border-right: 5px solid transparent;
    border-top: 7px solid #94a3b8;
    margin-right: 6px;
}

QComboBox:hover::down-arrow {
    border-top: 7px solid #a78bfa;
}

QComboBox QAbstractItemView {
    border: 2px solid #334155;
    border-radius: 10px;
    background-color: #1e293b;
    selection-background-color: #8b5cf6;
    selection-color: white;
    padding: 6px;
}

/* Text Edit - Dark terminal */
QTextEdit {
    border: 2px solid #334155;
    border-radius: 10px;
    background-color: #0f172a;
    padding: 12px;
    color: #e2e8f0;
    font-family: 'SF Mono', Monaco, 'Cascadia Code', 'Roboto Mono', Consolas, 'Courier New', monospace;
}

/* Group Box - Dark card with neon title */
QGroupBox {
    border: 2px solid #334155;
    border-radius: 12px;
    margin-top: 20px;
    padding-top: 20px;
    font-weight: 600;
    background-color: #1e293b;
}

QGroupBox::title {
    subcontrol-origin: margin;
    left: 20px;
    padding: 6px 16px;
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #8b5cf6, stop:1 #a78bfa);
    color: white;
    border-radius: 8px;
    font-size: 13px;
    font-weight: 800;
}

/* Scroll Area - Dark */
QScrollArea {
    border: none;
    background-color: transparent;
}

/* Scrollbar - Neon design */
QScrollBar:vertical {
    border: none;
    background-color: transparent;
    width: 12px;
    margin: 0px;
}

QScrollBar::handle:vertical {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #8b5cf6, stop:1 #a78bfa);
    border-radius: 6px;
    min-height: 40px;
}

QScrollBar::handle:vertical:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #a78bfa, stop:1 #c4b5fd);
}

QScrollBar::add-line:vertical, QScrollBar::sub-line:vertical {
    height: 0px;
}

/* Status Bar - Dark gradient */
QStatusBar {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #0f172a, stop:1 #1e293b);
    color: #94a3b8;
    border-top: 2px solid #334155;
    padding: 6px;
    font-weight: 500;
}

/* Menu Bar - Dark */
QMenuBar {
    background-color: #1e293b;
    color: #e2e8f0;
    border-bottom: 2px solid #334155;
    padding: 6px;
}

QMenuBar::item {
    padding: 10px 16px;
    background-color: transparent;
    border-radius: 8px;
    font-weight: 600;
}

QMenuBar::item:selected {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                stop:0 #4c1d95, stop:1 #5b21b6);
    color: #c4b5fd;
}

/* Menu - Dark */
QMenu {
    background-color: #1e293b;
    border: 2px solid #334155;
    border-radius: 12px;
    padding: 8px;
}

QMenu::item {
    padding: 10px 28px;
    border-radius: 8px;
    font-weight: 500;
}

QMenu::item:selected {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #8b5cf6, stop:1 #a78bfa);
    color: white;
}

/* Toolbar - Dark */
QToolBar {
    background: qlineargradient(x1:0, y1:0, x2:0, y2:1,
                                stop:0 #1e293b, stop:1 #0f172a);
    border: none;
    spacing: 10px;
    padding: 10px;
}

/* Checkbox - Dark neon */
QCheckBox {
    spacing: 8px;
    color: #e2e8f0;
    font-weight: 500;
}

QCheckBox::indicator {
    width: 20px;
    height: 20px;
    border-radius: 6px;
    border: 2px solid #475569;
    background-color: #0f172a;
}

QCheckBox::indicator:hover {
    border: 2px solid #8b5cf6;
}

QCheckBox::indicator:checked {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #8b5cf6, stop:1 #a78bfa);
    border: 2px solid #8b5cf6;
    image: url(data:image/svg+xml;base64,PHN2ZyB3aWR0aD0iMTIiIGhlaWdodD0iOSIgdmlld0JveD0iMCAwIDEyIDkiIGZpbGw9Im5vbmUiIHhtbG5zPSJodHRwOi8vd3d3LnczLm9yZy8yMDAwL3N2ZyI+CjxwYXRoIGQ9Ik0xIDRMNC41IDdMMTEgMSIgc3Ryb2tlPSJ3aGl0ZSIgc3Ryb2tlLXdpZHRoPSIyIiBzdHJva2UtbGluZWNhcD0icm91bmQiLz4KPC9zdmc+Cg==);
}

/* Label - Dark */
QLabel {
    color: #e2e8f0;
    font-weight: 500;
}

/* Splitter - Dark */
QSplitter::handle {
    background-color: #334155;
    width: 2px;
}

QSplitter::handle:hover {
    background: qlineargradient(x1:0, y1:0, x2:1, y2:0,
                                stop:0 #8b5cf6, stop:1 #a78bfa);
}
"""


def apply_modern_style(widget, theme="light"):
    """Apply beautiful modern styling to a widget"""
    if theme == "dark":
        widget.setStyleSheet(DARK_THEME)
    else:
        widget.setStyleSheet(LIGHT_THEME)
