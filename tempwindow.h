#ifndef TEMPWINDOW_H
#define TEMPWINDOW_H

#include <QLabel>
#include <QMainWindow>
#include <QScrollArea>

class TempWindow : public QMainWindow
{
    Q_OBJECT
    public:
    TempWindow(QWidget *parent = nullptr);
    ~TempWindow();
    public:
        QLabel * tempLabel;
        QScrollArea *scrollArea;
};

#endif // TEMPWINDOW_H
