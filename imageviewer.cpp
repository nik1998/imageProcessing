
#include "imageviewer.h"
#include "cmath"
#include <QtWidgets>
#include <QRgb>
#include <QtMath>
#include <cstdio>
#include <iostream>

#if defined(QT_PRINTSUPPORT_LIB)
#include <QtPrintSupport/qtprintsupportglobal.h>
#if QT_CONFIG(printdialog)
#include <QPrintDialog>
#endif
#endif

ImageViewer::ImageViewer(QWidget *parent)
   : QMainWindow(parent), imageLabel(new QLabel),
     scrollArea(new QScrollArea), scaleFactor(1)
{
    imageLabel->setBackgroundRole(QPalette::Base);
    imageLabel->setSizePolicy(QSizePolicy::Ignored, QSizePolicy::Ignored);
    imageLabel->setScaledContents(true);

    scrollArea->setBackgroundRole(QPalette::Dark);
    scrollArea->setWidget(imageLabel);
    scrollArea->setVisible(false);
    setCentralWidget(scrollArea);

    createActions();

    resize(QGuiApplication::primaryScreen()->availableSize() * 3 / 5);
}

//! [0]
//! [2]

bool ImageViewer::loadFile(const QString &fileName)
{
    QImageReader reader(fileName);
    reader.setAutoTransform(true);
    const QImage newImage = reader.read();
    if (newImage.isNull()) {
        QMessageBox::information(this, QGuiApplication::applicationDisplayName(),
                                 tr("Cannot load %1: %2")
                                 .arg(QDir::toNativeSeparators(fileName), reader.errorString()));
        return false;
    }
//! [2]

    setImage(newImage);

    setWindowFilePath(fileName);

    const QString message = tr("Opened \"%1\", %2x%3, Depth: %4")
        .arg(QDir::toNativeSeparators(fileName)).arg(image.width()).arg(image.height()).arg(image.depth());
    statusBar()->showMessage(message);
    return true;
}
void ImageViewer::createCosImage()
{
    int height = 600;
    int width = 1000;
    QImage im(width,height,QImage::Format_RGB32);
    double w = 0;
    for(int i=0;i<width;i++)
    {
        for(int j=0;j<height;j++)
        {
            int v =abs(static_cast<int>(255*cos(w*i/width)));
            im.setPixelColor(i,j,qRgb(v,v,v));
        }
        w+=1/2.0;
    }
    setImage(im);
}
QVector<QVector<int>> getGrayMatrix(QImage &image)
{
    int n =image.width();
    int m = image.height();
    QVector<QVector<int>>a(m);
    for (int ii = 0; ii < m; ii++) {
        a[ii].resize(n);
    }
    for (int ii = 0; ii < m; ii++) {
        for (int jj = 0; jj < n; jj++)
        {
            a[ii][jj]=0;
        }
    }
    for (int ii = 0; ii < image.height(); ii++) {
        uchar* scan = image.scanLine(ii);
        int depth =4;
        for (int jj = 0; jj < image.width(); jj++) {
            QRgb* rgbpixel = reinterpret_cast<QRgb*>(scan + jj*depth);
            int gray = qGray(*rgbpixel);
            a[ii][jj]=gray;
            *rgbpixel = QColor(gray, gray, gray).rgba();
        }
    }
    return a;
}
void ImageViewer::correlation()
{
    QImage image = imageLabel->pixmap()->toImage();
    int n =image.width();
    int m = image.height();
    QVector<QVector<int>>a = getGrayMatrix(image);
    QImage filter = imageLabel->pixmap()->copy(70,60,70,60).toImage();
    QVector<QVector<int>>f=getGrayMatrix(filter);
    QImage res(image.width()-filter.width()+1,image.height()-filter.height()+1,QImage::Format_RGB32);
    QVector<QVector<int>>r = getGrayMatrix(res);
    double maxx=0;
    double su=0;
    for(int i=0;i<m;i++)
    {
        for(int j=0;j<n;j++)
        {
            su+=a[i][j];
        }
    }
    su=su/n/m;
    for (int ii = filter.height()/2; ii < m-filter.height()/2; ii++)
    {
        for (int jj =filter.width()/2; jj < n-filter.width()/2; jj++)
        {
            double pr  = 0;
            r[ii-filter.height()/2][jj-filter.width()/2]=0;
            for(int i=0;i<filter.height();i++)
            {
                for(int j=0;j<filter.width();j++)
                {
                   pr += (f[i][j] - su)*(a[ii + i-filter.height()/2][jj +j-filter.width()/2] -su);
                }
            }
            r[ii-filter.height()/2][jj-filter.width()/2]=pr/255/255;
            maxx = max(maxx,abs(pr/255/255));
            //res.setPixelColor(jj-filter.width()/2,ii-filter.height()/2,qRgb(pr,pr,pr));
        }
    }
    for(int i=0;i<res.height();i++)
    {
        for(int j=0;j<res.width();j++)
        {
            r[i][j]=abs(r[i][j]*255/maxx);
            //cout<<r[i][j];
            res.setPixelColor(j,i,qRgb(r[i][j],r[i][j],r[i][j]));
        }
    }
    imageLabel->setPixmap(QPixmap::fromImage(image));
    tempWindow.tempLabel->setPixmap(QPixmap::fromImage(res));
    tempWindow.tempLabel->adjustSize();
    tempWindow.show();
    frWindow.tempLabel->setPixmap(QPixmap::fromImage(filter));
    frWindow.tempLabel->adjustSize();
    frWindow.show();
}
void ImageViewer::sobel()
{
    QImage image = imageLabel->pixmap()->toImage();
    int msobel[3][3]={{-1,0,1},{-2,0,2},{-1,0,1}};

    QVector<QVector<int>>a=getGrayMatrix(image);
    int n =image.width();
    int m = image.height();
    QImage res(image.width()-2,image.height()-2,QImage::Format_RGB32);
    for (int ii = 1; ii < m-1; ii++)
    {
        for (int jj = 1; jj < n-1; jj++)
        {
            int x=0;
            int y=0;
            for(int i=0;i<3;i++)
            {
                for(int j=0;j<3;j++)
                {
                    x+=a[ii+i-1][jj+j-1]*msobel[i][j];
                    y+=a[ii+i-1][jj+j-1]*msobel[j][i];
                }
            }
            int ans = static_cast<int>(sqrt(x*x+y+y));
            res.setPixelColor(jj-1,ii-1,qRgb(ans,ans,ans));
        }
    }
    imageLabel->setPixmap(QPixmap::fromImage(image));
    tempWindow.tempLabel->setPixmap(QPixmap::fromImage(res));
    tempWindow.tempLabel->adjustSize();
    tempWindow.show();
}
void ImageViewer::fourier()
{
   // open();
    QImage image = imageLabel->pixmap()->toImage();
    int n = 1;
    while (n < image.width())  n <<= 1;
    int m = 1;
    while (m < image.height())  m <<= 1;
    QVector<QVector<base>>a(m);
    for (int ii = 0; ii < m; ii++) {
        a[ii].resize(n);
    }
    for (int ii = 0; ii < m; ii++) {
        for (int jj = 0; jj < n; jj++)
        {
            a[ii][jj]=0;
        }
    }
    QImage res(image.width(),image.height(),QImage::Format_RGB32);
    for (int ii = 0; ii < image.height(); ii++) {
        uchar* scan = image.scanLine(ii);
        int depth =4;
        for (int jj = 0; jj < image.width(); jj++) {
            QRgb* rgbpixel = reinterpret_cast<QRgb*>(scan + jj*depth);
            int gray = qGray(*rgbpixel);
            a[ii][jj]=gray;
            *rgbpixel = QColor(gray, gray, gray).rgba();
            res.setPixelColor(jj,ii,qRgb(255,255,255));
        }
    }
    imageLabel->setPixmap(QPixmap::fromImage(image));
    QVector<pair<int,int>> ans =tenHighFrequencies(a, res);
    tempWindow.tempLabel->setPixmap(QPixmap::fromImage(res));
    tempWindow.tempLabel->adjustSize();
    tempWindow.show();
    QImage f(2500,2500,QImage::Format_RGB32);
    f.fill(QColor(Qt::white).rgb());
    for(int i=0;i<10;i++)
    {
        for(int j=0;j<10;j++)
        {
            for(int ci=25;ci<225;ci++)
            {
                for(int cj=25;cj<225;cj++)
                {
                    int in = i*10+j;
                    int v =abs(static_cast<int>(255*cos(ans[in].first*(ci/250.0)+ans[in].second*(cj/250.0))));
                    f.setPixelColor(250*j + cj,250*i + ci,qRgb(v,v,v));
                }
            }
        }
    }
    frWindow.tempLabel->setPixmap(QPixmap::fromImage(f));
    frWindow.tempLabel->adjustSize();
    frWindow.show();
}


void ImageViewer::fft (QVector<base> & a, bool invert) {
    int n = static_cast<int>(a.size());
    if (n == 1)  return;
    QVector<base> a0 (n/2),  a1 (n/2);
    for (int i=0, j=0; i<n; i+=2, ++j) {
        a0[j] = a[i];
        a1[j] = a[i+1];
    }
    fft (a0, invert);
    fft (a1, invert);

    double ang = 2*M_PI/n * (invert ? -1 : 1);
    base w (1),  wn (cos(ang), sin(ang));
    for (int i=0; i<n/2; ++i) {
        a[i] = a0[i] + w * a1[i];
        a[i+n/2] = a0[i] - w * a1[i];
        if (invert)
        {
            a[i] /= 2;
            a[i+n/2] /= 2;
        }
        w *= wn;
    }
}
template <typename T>
QVector<QVector<T>> ImageViewer::transpose(QVector<QVector<T>> v)
{
    QVector<QVector<T>> ans(v[0].size());
    for(int i=0;i<v[0].size();i++)
    {
        ans[i].resize(v.size());
        for(int j =0;j<v.size();j++)
        {
            ans[i][j]=v[j][i];
        }
    }
    return ans;
}
QVector<pair<int,int>> ImageViewer::tenHighFrequencies (QVector<QVector<base>> f,QImage &res) {
    for(int i=0;i<f.size();i++)
    {
        fft(f[i],false);
    }
    f = transpose<base>(f);
    for(int i=0;i<f[0].size();i++)
    {
        fft(f[i],false);
    }
    QVector<pair<int,int>> ans;
    auto cmp = [](base a, base b) {
        double va = a.real()*a.real()+a.imag()*a.imag();
        double vb = b.real()*b.real()+b.imag()*b.imag();
        if(vb==va)
        {
            return a.real()<b.real();
        }
        return va<vb;
    };
    std::set<base,decltype(cmp)>s(cmp);
    for(int i=0;i<f.size();i++)
    {
        for(int j=0;j<f[0].size();j++)
        {
            s.insert(f[i][j]);
            if(s.size()>100)
            {
                s.erase(s.begin());
            }
        }
    }

    for(int i=0;i<f.size();i++)
    {
        for(int j=0;j<f[0].size();j++)
        {
            if(s.find(f[i][j])==s.end())
            {
                f[i][j]=0;
            }
            else
            {
                ans.push_back(make_pair(i,j));
            }
        }
    }
    for(int i=0;i<f[0].size();i++)
    {
        fft(f[i],true);
    }
    f = transpose<base>(f);
    for(int i=0;i<f.size();i++)
    {
        fft(f[i],true);
    }
    for (int ii = 0; ii < res.height(); ii++) {
        for (int jj = 0; jj < res.width(); jj++) {
            int gray = static_cast<int>(f[ii][jj].real());
            //cout<< gray<<" ";
            res.setPixelColor(jj,ii,qRgb(gray, gray, gray));
        }
        //cout<< "\n";
    }
    return ans;
}

void ImageViewer::setImage(const QImage &newImage)
{
    image = newImage;
    imageLabel->setPixmap(QPixmap::fromImage(image));

    scrollArea->setVisible(true);
    printAct->setEnabled(true);
    fitToWindowAct->setEnabled(true);
    updateActions();

    if (!fitToWindowAct->isChecked())
        imageLabel->adjustSize();
}


bool ImageViewer::saveFile(const QString &fileName)
{
    QImageWriter writer(fileName);

    if (!writer.write(image)) {
        QMessageBox::information(this, QGuiApplication::applicationDisplayName(),
                                 tr("Cannot write %1: %2")
                                 .arg(QDir::toNativeSeparators(fileName)), writer.errorString());
        return false;
    }
    const QString message = tr("Wrote \"%1\"").arg(QDir::toNativeSeparators(fileName));
    statusBar()->showMessage(message);
    return true;
}

//! [1]

static void initializeImageFileDialog(QFileDialog &dialog, QFileDialog::AcceptMode acceptMode)
{
    static bool firstDialog = true;

    if (firstDialog) {
        firstDialog = false;
        const QStringList picturesLocations = QStandardPaths::standardLocations(QStandardPaths::PicturesLocation);
        dialog.setDirectory(picturesLocations.isEmpty() ? QDir::currentPath() : picturesLocations.last());
    }

    QStringList mimeTypeFilters;
    const QByteArrayList supportedMimeTypes = acceptMode == QFileDialog::AcceptOpen
        ? QImageReader::supportedMimeTypes() : QImageWriter::supportedMimeTypes();
    for (const QByteArray &mimeTypeName : supportedMimeTypes)
        mimeTypeFilters.append(mimeTypeName);
    mimeTypeFilters.sort();
    dialog.setMimeTypeFilters(mimeTypeFilters);
    dialog.selectMimeTypeFilter("image/jpeg");
    if (acceptMode == QFileDialog::AcceptSave)
        dialog.setDefaultSuffix("jpg");
}

void ImageViewer::open()
{
    QFileDialog dialog(this, tr("Open File"));
    initializeImageFileDialog(dialog, QFileDialog::AcceptOpen);

    while (dialog.exec() == QDialog::Accepted && !loadFile(dialog.selectedFiles().first())) {}
}


void ImageViewer::saveAs()
{
    QFileDialog dialog(this, tr("Save File As"));
    initializeImageFileDialog(dialog, QFileDialog::AcceptSave);

    while (dialog.exec() == QDialog::Accepted && !saveFile(dialog.selectedFiles().first())) {}
}

//! [5]
void ImageViewer::print()
//! [5] //! [6]
{
    Q_ASSERT(imageLabel->pixmap());
#if QT_CONFIG(printdialog)
//! [6] //! [7]
    QPrintDialog dialog(&printer, this);
//! [7] //! [8]
    if (dialog.exec()) {
        QPainter painter(&printer);
        QRect rect = painter.viewport();
        QSize size = imageLabel->pixmap()->size();
        size.scale(rect.size(), Qt::KeepAspectRatio);
        painter.setViewport(rect.x(), rect.y(), size.width(), size.height());
        painter.setWindow(imageLabel->pixmap()->rect());
        painter.drawPixmap(0, 0, *imageLabel->pixmap());
    }
#endif
}
//! [8]

void ImageViewer::copy()
{
#ifndef QT_NO_CLIPBOARD
    QGuiApplication::clipboard()->setImage(image);
#endif // !QT_NO_CLIPBOARD
}

#ifndef QT_NO_CLIPBOARD
static QImage clipboardImage()
{
    if (const QMimeData *mimeData = QGuiApplication::clipboard()->mimeData()) {
        if (mimeData->hasImage()) {
            const QImage image = qvariant_cast<QImage>(mimeData->imageData());
            if (!image.isNull())
                return image;
        }
    }
    return QImage();
}
#endif // !QT_NO_CLIPBOARD

void ImageViewer::paste()
{
#ifndef QT_NO_CLIPBOARD
    const QImage newImage = clipboardImage();
    if (newImage.isNull()) {
        statusBar()->showMessage(tr("No image in clipboard"));
    } else {
        setImage(newImage);
        setWindowFilePath(QString());
        const QString message = tr("Obtained image from clipboard, %1x%2, Depth: %3")
            .arg(newImage.width()).arg(newImage.height()).arg(newImage.depth());
        statusBar()->showMessage(message);
    }
#endif // !QT_NO_CLIPBOARD
}

//! [9]
void ImageViewer::zoomIn()
//! [9] //! [10]
{
    scaleImage(1.25);
}

void ImageViewer::zoomOut()
{
    scaleImage(0.8);
}

//! [10] //! [11]
void ImageViewer::normalSize()
//! [11] //! [12]
{
    imageLabel->adjustSize();
    scaleFactor = 1.0;
}
//! [12]

//! [13]
void ImageViewer::fitToWindow()
//! [13] //! [14]
{
    bool fitToWindow = fitToWindowAct->isChecked();
    scrollArea->setWidgetResizable(fitToWindow);
    if (!fitToWindow)
        normalSize();
    updateActions();
}
//! [14]


//! [15]
void ImageViewer::about()
//! [15] //! [16]
{
    QMessageBox::about(this, tr("About Image Viewer"),
            tr("<p>The <b>Image Viewer</b> example shows how to combine QLabel "
               "and QScrollArea to display an image. QLabel is typically used "
               "for displaying a text, but it can also display an image. "
               "QScrollArea provides a scrolling view around another widget. "
               "If the child widget exceeds the size of the frame, QScrollArea "
               "automatically provides scroll bars. </p><p>The example "
               "demonstrates how QLabel's ability to scale its contents "
               "(QLabel::scaledContents), and QScrollArea's ability to "
               "automatically resize its contents "
               "(QScrollArea::widgetResizable), can be used to implement "
               "zooming and scaling features. </p><p>In addition the example "
               "shows how to use QPainter to print an image.</p>"));
}
//! [16]

//! [17]
void ImageViewer::createActions()
//! [17] //! [18]
{
    QMenu *fileMenu = menuBar()->addMenu(tr("&File"));

    QAction *openAct = fileMenu->addAction(tr("&Open..."), this, &ImageViewer::open);
    openAct->setShortcut(QKeySequence::Open);

    saveAsAct = fileMenu->addAction(tr("&Save As..."), this, &ImageViewer::saveAs);
    saveAsAct->setEnabled(false);

    printAct = fileMenu->addAction(tr("&Print..."), this, &ImageViewer::print);
    printAct->setShortcut(QKeySequence::Print);
    printAct->setEnabled(false);

    fileMenu->addSeparator();

    QAction *exitAct = fileMenu->addAction(tr("E&xit"), this, &QWidget::close);
    exitAct->setShortcut(tr("Ctrl+Q"));

    QMenu *editMenu = menuBar()->addMenu(tr("&Edit"));

    copyAct = editMenu->addAction(tr("&Copy"), this, &ImageViewer::copy);
    copyAct->setShortcut(QKeySequence::Copy);
    copyAct->setEnabled(false);

    QAction *pasteAct = editMenu->addAction(tr("&Paste"), this, &ImageViewer::paste);
    pasteAct->setShortcut(QKeySequence::Paste);

    frequentlyAct = editMenu->addAction(tr("&Frequently strips"), this, &ImageViewer::createCosImage);
    frequentlyAct->setEnabled(true);
    fourierAct = editMenu->addAction(tr("&Fourier transformation"), this, &ImageViewer::fourier);
    fourierAct->setEnabled(true);
    sobelAct= editMenu->addAction(tr("&Sobel"), this, &ImageViewer::sobel);
    sobelAct->setEnabled(true);
    correlationFunctionAct = editMenu->addAction(tr("&Correlation function"), this, &ImageViewer::correlation);
    correlationFunctionAct->setEnabled(true);

    QMenu *viewMenu = menuBar()->addMenu(tr("&View"));

    zoomInAct = viewMenu->addAction(tr("Zoom &In (25%)"), this, &ImageViewer::zoomIn);
    zoomInAct->setShortcut(QKeySequence::ZoomIn);
    zoomInAct->setEnabled(false);

    zoomOutAct = viewMenu->addAction(tr("Zoom &Out (25%)"), this, &ImageViewer::zoomOut);
    zoomOutAct->setShortcut(QKeySequence::ZoomOut);
    zoomOutAct->setEnabled(false);

    normalSizeAct = viewMenu->addAction(tr("&Normal Size"), this, &ImageViewer::normalSize);
    normalSizeAct->setShortcut(tr("Ctrl+S"));
    normalSizeAct->setEnabled(false);

    viewMenu->addSeparator();

    fitToWindowAct = viewMenu->addAction(tr("&Fit to Window"), this, &ImageViewer::fitToWindow);
    fitToWindowAct->setEnabled(false);
    fitToWindowAct->setCheckable(true);
    fitToWindowAct->setShortcut(tr("Ctrl+F"));

    QMenu *helpMenu = menuBar()->addMenu(tr("&Help"));

    helpMenu->addAction(tr("&About"), this, &ImageViewer::about);
    helpMenu->addAction(tr("About &Qt"), &QApplication::aboutQt);
}
//! [18]

//! [21]
void ImageViewer::updateActions()
//! [21] //! [22]
{
    saveAsAct->setEnabled(!image.isNull());
    copyAct->setEnabled(!image.isNull());
    zoomInAct->setEnabled(!fitToWindowAct->isChecked());
    zoomOutAct->setEnabled(!fitToWindowAct->isChecked());
    normalSizeAct->setEnabled(!fitToWindowAct->isChecked());
}
//! [22]

//! [23]
void ImageViewer::scaleImage(double factor)
//! [23] //! [24]
{
    Q_ASSERT(imageLabel->pixmap());
    scaleFactor *= factor;
    imageLabel->resize(scaleFactor * imageLabel->pixmap()->size());

    adjustScrollBar(scrollArea->horizontalScrollBar(), factor);
    adjustScrollBar(scrollArea->verticalScrollBar(), factor);

    zoomInAct->setEnabled(scaleFactor < 3.0);
    zoomOutAct->setEnabled(scaleFactor > 0.333);
}
//! [24]

//! [25]
void ImageViewer::adjustScrollBar(QScrollBar *scrollBar, double factor)
//! [25] //! [26]
{
    scrollBar->setValue(int(factor * scrollBar->value()
                            + ((factor - 1) * scrollBar->pageStep()/2)));
}
//! [26]
