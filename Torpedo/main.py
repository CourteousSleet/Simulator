from classes import *


def __main__():
    app = QApplication([])
    window = OceanWidget()
    window.show()
    app.exec_()


if __name__ == '__main__':
    __main__()
