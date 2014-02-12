#include "tester.h"

Tester::Tester(QWidget *parent)
  : QWidget(parent)
{
  ui.setupUi(this);

  connect(ui.closeButton, SIGNAL(clicked()), this, SLOT(close()));

  setWindowTitle("ElmerGUI installation tester");
  setWindowIcon(QIcon(":/img/ElmerGUItester.ico"));

  elmerHome = get("ELMER_HOME");
  elmerGuiHome = get("ELMERGUI_HOME");
  ok = true;
  e = ui.verdict;
}

QString Tester::get(const QString &variable) const
{
  QString value(getenv(qPrintable(variable)));

#ifdef Q_WS_WIN
  while(value.endsWith("\\"))
    value.chop(1);
#else
  while(value.endsWith("/"))
    value.chop(1);
#endif

  return value;
}

bool Tester::testDir(const QString &variable, QLabel *label) const
{
  QString value(get(variable));

  label->setText(value);
  label->setAutoFillBackground(true);

  if(!value.isEmpty() && QDir(value).exists()) {
    label->setPalette(QPalette(Qt::green));
    return true;
  }

  label->setPalette(QPalette(Qt::red));
  return false;
}

bool Tester::testFile(const QString &value, QLabel *label) const
{
  label->setText(value);
  label->setAutoFillBackground(true);

  if(!value.isEmpty() && QFile(value).exists()) {
    label->setPalette(QPalette(Qt::green));
    return true;
  }

  label->setPalette(QPalette(Qt::red));
  return false;
}

bool Tester::testPath(const QString &value, QLabel *label) const
{
  QString path(get("PATH"));

#ifdef Q_WS_WIN
  QStringList splitPath(path.toUpper().split(";"));
#else
  QStringList splitPath(path.split(":"));
#endif

  label->setText(value);
  label->setAutoFillBackground(true);

#ifdef Q_WS_WIN
  if(splitPath.contains(value.toUpper())) {
    label->setPalette(QPalette(Qt::green));
    return true;
  }
#else
  if(splitPath.contains(value)) {
    label->setPalette(QPalette(Qt::green));
    return true;
  }
#endif
  
  label->setPalette(QPalette(Qt::red));  
  return false;
}

bool Tester::testLdLibraryPath(const QString &value, QLabel *label) const
{
  QString ldLibraryPath(get("LD_LIBRARY_PATH"));

#ifdef Q_WS_WIN
  QStringList splitLdLibraryPath(ldLibraryPath.toUpper().split(";"));
#else
  QStringList splitLdLibraryPath(ldLibraryPath.split(":"));
#endif

  label->setText(value);
  label->setAutoFillBackground(true);

#ifdef Q_WS_WIN
  if(splitLdLibraryPath.contains(value.toUpper())) {
    label->setPalette(QPalette(Qt::green));
    return true;
  }
#else
  if(splitLdLibraryPath.contains(value)) {
    label->setPalette(QPalette(Qt::green));
    return true;
  }
#endif
  
  label->setPalette(QPalette(Qt::red));  
  return false;
}

void Tester::testEnvironment()
{
  ok &= testDir("ELMER_HOME", ui.elmerHomeResult);
  ok &= testDir("ELMERGUI_HOME", ui.elmerGuiHomeResult);
  ok &= testDir("ELMER_POST_HOME", ui.elmerPostHomeResult);

#ifdef Q_WS_WIN
  ui.ldLibraryPathLabel->setText("PATH");
  ok &= testPath(elmerHome + "\\bin", ui.pathResult);
  ok &= testPath(elmerHome + "\\lib", ui.ldLibraryPathResult);
#else
  ok &= testPath(elmerHome + "/bin", ui.pathResult);
  ok &= testLdLibraryPath(elmerHome + "/lib", ui.ldLibraryPathResult);
#endif
}

void Tester::testExecutables()
{
#ifdef Q_WS_WIN
  ok &= testFile(elmerHome + "\\bin\\ElmerSolver.exe", ui.elmerSolverResult);
  ok &= testFile(elmerGuiHome + "\\ElmerGUI.exe", ui.elmerGuiResult);
  ok &= testFile(elmerHome + "\\bin\\ElmerPost.exe", ui.elmerPostResult);
  ok &= testFile(elmerHome + "\\bin\\ElmerGrid.exe", ui.elmerGridResult);
#else
  ok &= testFile(elmerHome + "/bin/ElmerSolver", ui.elmerSolverResult);
  ok &= testFile(elmerGuiHome + "/ElmerGUI", ui.elmerGuiResult);
  ok &= testFile(elmerHome + "/bin/ElmerPost", ui.elmerPostResult);
  ok &= testFile(elmerHome + "/bin/ElmerGrid", ui.elmerGridResult);
#endif
}

void Tester::verdict()
{
  repaint();

  if(ok) {
    e->append("Elmer seems to be installed correctly on this system");
    e->append("");
    e->append("Performing some additional tests:");
    
    testSolver();
    testPost();
    testGrid();
    testTetgen();

    return;
  }

  e->append("Elmer seems to be installed incorrectly on this system");
  e->append("1) Make sure that ELMER_HOME has been set up properly");
  e->append("2) Set ELMERGUI_HOME to ELMER_HOME/bin");
  e->append("3) Set ELMER_POST_HOME to ELMER_HOME/share/elmerpost");
  e->append("4) Make sure that ELMER_HOME/bin is in PATH");
#ifdef Q_WS_WIN
  e->append("5) Make sure that ELMER_HOME/lib is in PATH");
#else
  e->append("5) Make sure that ELMER_HOME/lib is in LD_LIBRARY_PATH");
#endif
  e->append("6) Executables should be found from ELMER_HOME/bin");
}

void Tester::testSolver()
{
  e->append("");
  e->append("Checking whether ElmerSolver starts...");

  QString solverName("ElmerSolver");
  QStringList solverArgs;
  solverArgs << "-v";

  QProcess *solver = new QProcess(this);
  solver->start(solverName, solverArgs);

  if(!solver->waitForStarted()) {
    e->append("ERROR: ElmerSolver does not start");
    return;
  }

  if(!solver->waitForFinished()) {
    e->append("ERROR: ElmerSolver does not finish");
    return;
  }

  QString str(solver->readAllStandardOutput());
  str.replace("\r", "");
  QStringList list = str.split("\n");

  foreach(QString line, list) {
    if(line.contains("Library version"))
      e->append(line.replace("MAIN:", "").trimmed());
  }

  e->append("OK: ElmerSolver seems functional");

  repaint();
}

void Tester::testPost()
{
  e->append("");
  e->append("Checking whether ElmerPost starts...");

  QString postName("ElmerPost");
  QStringList postArgs;
  postArgs << "-v";

  QProcess *post = new QProcess(this);
  post->start(postName, postArgs);

  if(!post->waitForStarted()) {
    e->append("ERROR: ElmerPost does not start");
    return;
  }

  if(!post->waitForFinished()) {
    e->append("ERROR: ElmerPost does not finish");
    return;
  }

  QString str(post->readAllStandardOutput());
  str.replace("\r", "");
  str.replace("\n", "");
  str = str.trimmed();
  if(!str.isEmpty())
    e->append(str);

  e->append("OK: ElmerPost seems functional");

  repaint();
}

void Tester::testGrid()
{
  e->append("");
  e->append("Checking whether ElmerGrid starts...");

  QString gridName("ElmerGrid");

  QProcess *grid = new QProcess(this);
  grid->start(gridName);

  if(!grid->waitForStarted()) {
    e->append("ERROR: ElmerGrid does not start");
    return;
  }

  if(!grid->waitForFinished()) {
    e->append("ERROR: ElmerGrid does not finish");
    return;
  }

  QString str(grid->readAllStandardOutput());
  e->append("OK: ElmerGrid seems functional");

  repaint();
}

void Tester::testTetgen()
{
  e->append("");
  e->append("Checking whether Tetgen plugin loads...");

  QLibrary tetgen("tetplugin");

  if(!tetgen.load()) {
    e->append("INFO: Tetgen plugin is not available");
    return;
  }

  if(!tetgen.resolve("CreateObjectOfTetgenio")) {
    e->append("WARNING: libtet version mismatch");
    return;
  }
  
  e->append("OK: Tetgen plugin seems functional");
}
