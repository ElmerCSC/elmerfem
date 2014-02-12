#ifndef TESTER_H
#define TESTER_H

#include <QtGui>
#include "ui_mainform.h"

class Tester : public QWidget
{
 Q_OBJECT

 public:
  Tester(QWidget *parent = 0);
  void testEnvironment();
  void testExecutables();
  void verdict();

 private:
  QString get(const QString &variable) const;
  bool testDir(const QString &variable, QLabel *label) const;
  bool testFile(const QString &value, QLabel *label) const;
  bool testPath(const QString &value, QLabel *label) const;
  bool testLdLibraryPath(const QString &value, QLabel *label) const;
  void testSolver();
  void testPost();
  void testGrid();
  void testTetgen();

  Ui::mainForm ui;
  QString elmerHome;
  QString elmerGuiHome;
  bool ok;
  QTextEdit *e;
};

#endif // TESTER_H
