QT += widgets
requires(qtConfig(filedialog))
qtHaveModule(printsupport): QT += printsupport

HEADERS       = imageviewer.h \
    tempwindow.h
SOURCES       = imageviewer.cpp \
                main.cpp \
                tempwindow.cpp

# install
target.path = $$[QT_INSTALL_EXAMPLES]/widgets/widgets/imageviewer
INSTALLS += target
