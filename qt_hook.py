from PyQt4 import QtGui, QtCore
#from windowUi import Ui_MainWindow
 
class FingerTabBarWidget(QtGui.QTabBar):
    def __init__(self, parent=None, *args, **kwargs):
        self.tabSize = QtCore.QSize(kwargs.pop('width',100), kwargs.pop('height',25))
        QtGui.QTabBar.__init__(self, parent, *args, **kwargs)
                 
    def paintEvent(self, event):
        painter = QtGui.QStylePainter(self)
        option = QtGui.QStyleOptionTab()
 
        for index in range(self.count()):
            self.initStyleOption(option, index)
            tabRect = self.tabRect(index)
            tabRect.moveLeft(10)
            painter.drawControl(QtGui.QStyle.CE_TabBarTabShape, option)
            painter.drawText(tabRect, QtCore.Qt.AlignVCenter |\
                             QtCore.Qt.TextDontClip, \
                             self.tabText(index));
        painter.end()
    def tabSizeHint(self,index):
        return self.tabSize
 
# Shamelessly stolen from this thread:
#   http://www.riverbankcomputing.com/pipermail/pyqt/2005-December/011724.html
class FingerTabWidget(QtGui.QTabWidget):
    """A QTabWidget equivalent which uses our FingerTabBarWidget"""
    def __init__(self, parent=None):
        QtGui.QTabWidget.__init__(self,parent=parent)
        self.setTabBar(FingerTabBarWidget(self))

class MainWidget(QtGui.QMainWindow):

    # Option #3 - Custom signal
    def __init__(self):
        QtGui.QMainWindow.__init__(self)
        # This is always the same
        #self.ui=Ui_MainWindow()
        #self.ui.setupUi(self)

if __name__ == '__main__':
    import sys
    app = QtGui.QApplication(sys.argv)
    root = MainWidget()
    FingerTabWidget = FingerTabWidget(parent=root)
    root.show()
    sys.exit(app.exec_())