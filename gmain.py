import sys
import matplotlib.pyplot as plt

sys.path.append("ui")

from PySide6.QtWidgets import QApplication, QMainWindow
from PySide6.QtGui import QResizeEvent
from PySide6.QtGui import QRegularExpressionValidator
from PySide6.QtCore import QRegularExpression

from ui.main_window_ui import Ui_MainWindow
from motion_constraint import MotionConstraint
from s_curve_partial import SCurvePartial
from s_curve_full import SCurveFull

class Window(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setupUi(self)
        self.connectSignalsSlots()

        self.inputInitialVelocity.setValidator(QRegularExpressionValidator(QRegularExpression(r'\d+\.\d+|\d+(?![\d.])')))
        self.inputMaxAcceleration.setValidator(QRegularExpressionValidator(QRegularExpression(r'\d+\.\d+|\d+(?![\d.])')))
        self.inputMaxVelocity.setValidator(QRegularExpressionValidator(QRegularExpression(r'\d+\.\d+|\d+(?![\d.])')))
        self.inputMaxJerk.setValidator(QRegularExpressionValidator(QRegularExpression(r'\d+\.\d+|\d+(?![\d.])')))
        self.inputDisplacement.setValidator(QRegularExpressionValidator(QRegularExpression(r'\d+\.\d+|\d+(?![\d.])')))
        self.inputTimeConstraint.setValidator(QRegularExpressionValidator(QRegularExpression(r'\d+\.\d+|\d+(?![\d.])')))

    def connectSignalsSlots(self):
        self.buttonSimulate.clicked.connect(self.simulate)
        self.cbTimeConstraint.stateChanged.connect(self.timeConstraintChanged)

    def timeConstraintChanged(self):
        self.inputTimeConstraint.setEnabled(self.cbTimeConstraint.isChecked())
        self.inputTimeConstraint.setText("0")
        self.inputTimeConstraint.setFocus()
    
    def simulate(self):
        v0 = float(self.inputInitialVelocity.text())
        a = float(self.inputMaxAcceleration.text())
        v = float(self.inputMaxVelocity.text())
        j = float(self.inputMaxJerk.text())
        s = float(self.inputDisplacement.text())

        if (self.cbTimeConstraint.isChecked()):
            t = float(self.inputTimeConstraint.text())
        else:
            t = 0

        plot1 = self.comboPlot1.currentIndex()
        plot2 = self.comboPlot2.currentIndex()

        curve = SCurveFull(MotionConstraint(v0, v, a, j, s, t), 400)

        if curve.solve():
            self.plotWidget.canvas.ax[0].clear()
            self.plotWidget.canvas.ax[1].clear()

            if (plot1 > 0):
                if (plot1 == 1):
                    curve.plotS(self.plotWidget.canvas.ax[0], "Displacement")
                elif (plot1 == 2):
                    curve.plotV(self.plotWidget.canvas.ax[0], "Velocity")
                elif (plot1 == 3):
                    curve.plotA(self.plotWidget.canvas.ax[0], "Acceleration")

            if (plot2 > 0):

                if (plot2 == 1):
                    curve.plotS(self.plotWidget.canvas.ax[1], "Displacement")
                elif (plot2 == 2):
                    curve.plotV(self.plotWidget.canvas.ax[1], "Velocity")
                elif (plot2 == 3):
                    curve.plotA(self.plotWidget.canvas.ax[1], "Acceleration")

            plt.subplots_adjust(hspace = 0.4)
            self.plotWidget.canvas.draw()

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Window()
    win.show()
    sys.exit(app.exec())