#include <QApplication>
#include "Segmenter.h"

int main(int argc, char** argv) {
    QApplication a(argc, argv);
    Segmenter w;
    w.show();

    return a.exec();
}
