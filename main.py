"""
main.py
-------
Entry point for GRACE.
Run with: python main.py
"""

import sys
import os
import multiprocessing

from PyQt6.QtWidgets import QApplication

from app_state import AppState
from ui.main_window import MainWindow
from ui.style import STYLESHEET


def main():
    # Required for PyInstaller frozen multiprocessing on Windows
    multiprocessing.freeze_support()

    app = QApplication(sys.argv)
    app.setApplicationName("GRACE")
    app.setApplicationVersion("1.0.0")

    # Load IBM Plex fonts if bundled in assets/fonts/
    from PyQt6.QtGui import QFontDatabase
    _fonts_dir = os.path.join(os.path.dirname(os.path.abspath(__file__)), "assets", "fonts")
    if os.path.isdir(_fonts_dir):
        for _fname in os.listdir(_fonts_dir):
            if _fname.lower().endswith((".ttf", ".otf")):
                QFontDatabase.addApplicationFont(os.path.join(_fonts_dir, _fname))

    app.setStyleSheet(STYLESHEET)

    # High DPI — handled automatically in PyQt6

    state  = AppState()
    window = MainWindow(state)
    window.show()

    sys.exit(app.exec())


if __name__ == "__main__":
    main()