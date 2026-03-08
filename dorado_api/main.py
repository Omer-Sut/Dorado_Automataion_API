import sys
from PySide6.QtWidgets import (QApplication)
from app_window import AppWindow

if __name__ == "__main__":
    app = QApplication(sys.argv)
    w = AppWindow()
    w.resize(600, 500)
    w.show()
    sys.exit(app.exec())