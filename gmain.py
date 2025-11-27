import sys
import matplotlib.pyplot as plt

sys.path.append("ui")

from PySide6.QtWidgets import QApplication, QMainWindow
from PySide6 import QtGui
from PySide6.QtGui import QResizeEvent
from PySide6.QtGui import QRegularExpressionValidator
from PySide6.QtCore import QRegularExpression

from ui.main_window_ui import Ui_MainWindow
from motion_constraint import MotionConstraint
from motion import Motion

class Window(QMainWindow, Ui_MainWindow):
    def __init__(self, parent=None):
        super().__init__(parent)

        self.setupUi(self)
        self.connectSignalsSlots()
        self.setWindowIcon(QtGui.QIcon('ui/resources/icon.png'))

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

        samplingPoints = self.cbPlotSamplingPoints.isChecked()

        motion = Motion(MotionConstraint(v0, v, a, j, s, t), 400)

        if motion.simulate():
            self.plotWidget.canvas.ax[0].clear()
            self.plotWidget.canvas.ax[1].clear()

            if (plot1 > 0):
                if (plot1 == 1):
                    motion.getCurve().plotS(self.plotWidget.canvas.ax[0], "Displacement", samplingPoints)
                elif (plot1 == 2):
                    motion.getCurve().plotV(self.plotWidget.canvas.ax[0], "Velocity", samplingPoints)
                elif (plot1 == 3):
                    motion.getCurve().plotA(self.plotWidget.canvas.ax[0], "Acceleration", samplingPoints)

            if (plot2 > 0):

                if (plot2 == 1):
                    motion.getCurve().plotS(self.plotWidget.canvas.ax[1], "Displacement", samplingPoints)
                elif (plot2 == 2):
                    motion.getCurve().plotV(self.plotWidget.canvas.ax[1], "Velocity", samplingPoints)
                elif (plot2 == 3):
                    motion.getCurve().plotA(self.plotWidget.canvas.ax[1], "Acceleration", samplingPoints)

            plt.subplots_adjust(hspace = 0.4)
            self.plotWidget.canvas.draw()

            total_time = motion.getCurve().getTotalTime()
            self.label_result_total_time.setText("%0.4f s, %0.4f ms" % (total_time, total_time * 1000))
            self.label_result_max_velocity.setText("%0.4f units/s" % (motion.getCurve().getMaxVelocity()))
            self.label_result_max_acceleration.setText("%0.4f uints/s^2" % (motion.getCurve().getMaxAcceleration()))
            self.label_result_motion_profile.setText(motion.getCurve().getProfile())
        else:
            self.label_result_total_time.setText("N/A")
            self.label_result_max_velocity.setText("N/A")
            self.label_result_max_acceleration.setText("N/A")
            self.label_result_motion_profile.setText("N/A")

if __name__ == "__main__":
    app = QApplication(sys.argv)
    win = Window()
    win.showMaximized()
    sys.exit(app.exec())