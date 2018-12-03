#include <QtWidgets/QApplication>
#include "Mainwindow.h"
int main(int argc, char *argv[])
{
	QApplication a(argc, argv);
	Mainwindow w;
	return a.exec();
}