#include "tester.h"
#include <QApplication>


int main(int argc, char **argv) {
  QApplication app(argc, argv);
  Tester tester;
  tester.show();
  tester.testEnvironment();
  tester.testExecutables();
  tester.verdict();
  return app.exec();
}
