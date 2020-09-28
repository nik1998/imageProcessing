#include "imageviewer.h"

#include <QGuiApplication>


TempWindow::TempWindow(QWidget *parent)
   : QMainWindow(parent), tempLabel(new QLabel),
     scrollArea(new QScrollArea)
{
    tempLabel->setBackgroundRole(QPalette::Base);
    tempLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    tempLabel->setScaledContents(true);

    scrollArea->setBackgroundRole(QPalette::Dark);
    scrollArea->setWidget(tempLabel);
    scrollArea->setVisible(true);
    setCentralWidget(scrollArea);
    resize(800,600);
}

TempWindow::~TempWindow()
{

}
