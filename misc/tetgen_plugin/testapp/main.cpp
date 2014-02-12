#include <QtCore>
#include "tetgen.h"

typedef tetgenio* (*tetgenio_t)();

typedef void (*delegate_t)(int, char *tetgenbehavior, char*, tetgenio*,
			   tetgenio*, tetgenio*, tetgenio*);

int main(int argc, char **argv)
{
  QCoreApplication app(argc, argv);

  QLibrary plugin("tetplugin");

  if(!plugin.load())
    qFatal(qPrintable(plugin.errorString()));

  tetgenio_t ptetgenio = (tetgenio_t)plugin.resolve("CreateObjectOfTetgenio");

  if(!ptetgenio)
    qFatal(qPrintable(plugin.errorString()));

  delegate_t pdelegate = (delegate_t)plugin.resolve("delegate_tetrahedralize");

  if(!pdelegate)
    qFatal(qPrintable(plugin.errorString()));
  
  tetgenio *in = (ptetgenio)();
  tetgenio *out = (ptetgenio)();
  delegate_t delegate_tetrahedralize = pdelegate;
  
  in->initialize();
  in->load_poly((char *)"example");
  delegate_tetrahedralize(1, NULL, (char *)"JApq1.414", in, out, NULL, NULL);
}
