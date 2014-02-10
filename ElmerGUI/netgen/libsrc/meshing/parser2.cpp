#include <mystdlib.h>
#include "meshing.hpp"

#ifdef WIN32
#define COMMASIGN ':'
#else
#define COMMASIGN ','
#endif


namespace netgen
{


void LoadMatrixLine (istream & ist, DenseMatrix & m, int line)
{
  char ch;
  int pnum;
  float f;

  ist >> ch;
  while (ch != '}')
    {
      ist.putback (ch);
      ist >> f;
      ist >> ch;
      ist >> pnum;

      if (ch == 'x' || ch == 'X')
	m.Elem(line, 2 * pnum - 1) = f;
      if (ch == 'y' || ch == 'Y')
	m.Elem(line, 2 * pnum) = f;

      ist >> ch;
      if (ch == COMMASIGN)
	ist >> ch;
    }
}


void netrule :: LoadRule (istream & ist)
{
  char buf[256];
  char ch;
  Point2d p;
  INDEX_2 lin;
  int i, j;
  DenseMatrix tempoldutonewu(20, 20), tempoldutofreearea(20, 20),
    tempoldutofreearealimit(20, 20);

  tempoldutonewu = 0;
  tempoldutofreearea = 0;
  tempoldutofreearealimit = 0;

  noldp = 0;
  noldl = 0;

  ist.get (buf, sizeof(buf), '"');
  ist.get (ch);
  ist.get (buf, sizeof(buf), '"');
  ist.get (ch);

  // if(name != NULL) 
  delete [] name;
  name = new char[strlen (buf) + 1];
  strcpy (name, buf);
  //(*testout) << "name " << name << endl;
  //  (*mycout) << "Rule " << name << " found." << endl;

  do
    {
      ist >> buf;

      //(*testout) << "buf " << buf << endl;

      if (strcmp (buf, "quality") == 0)

	{
	  ist >> quality;
	}

      else if (strcmp (buf, "mappoints") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ')'

	      points.Append (p);
	      noldp++;

	      tolerances.SetSize (noldp);
	      tolerances.Elem(noldp).f1 = 1.0;
	      tolerances.Elem(noldp).f2 = 0;
	      tolerances.Elem(noldp).f3 = 1.0;

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      ist >> tolerances.Elem(noldp).f1;
		      ist >> ch;  // ','
		      ist >> tolerances.Elem(noldp).f2;
		      ist >> ch;  // ','
		      ist >> tolerances.Elem(noldp).f3;
		      ist >> ch;  // '}'
		    }
		  else if (ch == 'd')
		    {
		      //            delpoints.Append (noldp);
		      ist >> ch; // 'e'
		      ist >> ch; // 'l'
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}


      else if (strcmp (buf, "maplines") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> lin.I1();
	      ist >> ch;    // ','
	      ist >> lin.I2();
	      ist >> ch;    // ')'


	      //(*testout) << "read line " << lin.I1() << " " << lin.I2() << endl;
	      lines.Append (lin);
	      linevecs.Append (points.Get(lin.I2()) - points.Get(lin.I1()));
	      noldl++;
	      linetolerances.SetSize (noldl);
	      linetolerances.Elem(noldl).f1 = 0;
	      linetolerances.Elem(noldl).f2 = 0;
	      linetolerances.Elem(noldl).f3 = 0;

	      //(*testout) << "mapl1" << endl; 
	      ist >> ch;
	      while (ch != ';')
		{
		  //(*testout) << "working on character \""<<ch<<"\""<< endl;
		  if (ch == '{')
		    {
		      ist >> linetolerances.Elem(noldl).f1;
		      ist >> ch;  // ','
		      ist >> linetolerances.Elem(noldl).f2;
		      ist >> ch;  // ','
		      ist >> linetolerances.Elem(noldl).f3;
		      ist >> ch;  // '}'
		    }
		  else if (ch == 'd')
		    {
		      dellines.Append (noldl);
		      ist >> ch; // 'e'
		      ist >> ch; // 'l'
		      //(*testout) << "read del" << endl;
		    }

		  ist >> ch;
		  //(*testout) << "read character \""<<ch<<"\""<< endl;
		}

	      ist >> ch;
	      //(*testout) << "read next character \""<<ch<<"\""<< endl;
	    }
	  

	  ist.putback (ch);
	}

      else if (strcmp (buf, "newpoints") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ')'

	      points.Append (p);

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      LoadMatrixLine (ist, tempoldutonewu,
				      2 * (points.Size()-noldp) - 1);

		      ist >> ch; // '{'
		      LoadMatrixLine (ist, tempoldutonewu,
				      2 * (points.Size()-noldp));
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "newlines") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> lin.I1();
	      ist >> ch;    // ','
	      ist >> lin.I2();
	      ist >> ch;    // ')'

	      lines.Append (lin);
	      linevecs.Append (points.Get(lin.I2()) - points.Get(lin.I1()));

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "freearea") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ')'

	      freezone.Append (p);
	      freezonelimit.Append (p);

	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      LoadMatrixLine (ist, tempoldutofreearea,
				      2 * freezone.Size() - 1);

		      ist >> ch; // '{'
		      LoadMatrixLine (ist, tempoldutofreearea,
				      2 * freezone.Size());
		    }

		  ist >> ch;
		}

	      ist >> ch;
	    }

	  for (i = 1; i <= tempoldutofreearealimit.Height(); i++)
	    for (j = 1; j <= tempoldutofreearealimit.Width(); j++)
	      tempoldutofreearealimit.Elem(i,j) =
		tempoldutofreearea.Elem(i,j);


	  ist.putback (ch);
	}    
      else if (strcmp (buf, "freearea2") == 0)
	{
	  ist >> ch;
	  int freepi = 0;
	  tempoldutofreearealimit = 0;

	  while (ch == '(')
	    {
	      freepi++;

	      ist >> p.X();
	      ist >> ch;    // ','
	      ist >> p.Y();
	      ist >> ch;    // ')'

	      freezonelimit.Elem(freepi) = p;
	  
	      ist >> ch;
	      while (ch != ';')
		{
		  if (ch == '{')
		    {
		      LoadMatrixLine (ist, tempoldutofreearealimit,
				      2 * freepi - 1);

		      ist >> ch; // '{'
		      LoadMatrixLine (ist, tempoldutofreearealimit,
				      2 * freepi);
		    }

		  ist >> ch;
		}
	  
	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "elements") == 0)
	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      elements.Append (Element2d());

	      ist >> elements.Last().PNum(1);
	      ist >> ch;    // ','
	  
	      if (ch == COMMASIGN)
		{
		  ist >> elements.Last().PNum(2);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  ist >> elements.Last().PNum(3);
		  ist >> ch;    // ','
		}
	      if (ch == COMMASIGN)
		{
		  elements.Last().SetType (QUAD);
		  ist >> elements.Last().PNum(4);
		  ist >> ch;    // ','
		  
		  // const Element2d & el = elements.Last();
		  /*
		  orientations.Append (threeint(el.PNum(1), el.PNum(2), el.PNum(3)));
		  orientations.Append (threeint(el.PNum(2), el.PNum(3), el.PNum(4)));
		  orientations.Append (threeint(el.PNum(3), el.PNum(4), el.PNum(1)));
		  orientations.Append (threeint(el.PNum(4), el.PNum(1), el.PNum(2)));
		  */
		}

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "orientations") == 0)

	{
	  ist >> ch;

	  while (ch == '(')
	    {
	      //        threeint a = threeint();
	      orientations.Append (threeint());

	      ist >> orientations.Last().i1;
	      ist >> ch;    // ','
	      ist >> orientations.Last().i2;
	      ist >> ch;    // ','
	      ist >> orientations.Last().i3;
	      ist >> ch;    // ','

	      ist >> ch;
	      while (ch != ';')
		{
		  ist >> ch;
		}

	      ist >> ch;
	    }

	  ist.putback (ch);
	}

      else if (strcmp (buf, "endrule") != 0)
	{
	  PrintSysError ("Parser error, unknown token ", buf);
	}
    }
  while (!ist.eof() && strcmp (buf, "endrule") != 0);

  //(*testout) << "loadr1" << endl;

  oldutonewu.SetSize (2 * (points.Size() - noldp), 2 * noldp);
  oldutofreearea.SetSize (2 * freezone.Size(), 2 * noldp);
  oldutofreearealimit.SetSize (2 * freezone.Size(), 2 * noldp);

  for (i = 1; i <= oldutonewu.Height(); i++)
    for (j = 1; j <= oldutonewu.Width(); j++)
      oldutonewu.Elem(i, j) = tempoldutonewu.Elem(i, j);

  for (i = 1; i <= oldutofreearea.Height(); i++)
    for (j = 1; j <= oldutofreearea.Width(); j++)
      oldutofreearea.Elem(i, j) = tempoldutofreearea.Elem(i, j);

  for (i = 1; i <= oldutofreearea.Height(); i++)
    for (j = 1; j <= oldutofreearea.Width(); j++)
      oldutofreearealimit.Elem(i, j) = tempoldutofreearealimit.Elem(i, j);

  freesetinequ.SetSize (freezone.Size());


  //(*testout) << "loadr2" << endl;

  {
    char ok;
    int minn;
    ARRAY<int> pnearness (noldp);

    for (i = 1; i <= pnearness.Size(); i++)
      pnearness.Elem(i) = 1000;

    for (j = 1; j <= 2; j++)
      pnearness.Elem(GetPointNr (1, j)) = 0;

    //(*testout) << "loadr3" << endl;
    do
      {
	ok = 1;

	for (i = 1; i <= noldl; i++)
	  {
	    minn = 1000;
	    for (j = 1; j <= 2; j++)
	      minn = min2 (minn, pnearness.Get(GetPointNr (i, j)));

	    for (j = 1; j <= 2; j++)
	      if (pnearness.Get(GetPointNr (i, j)) > minn+1)
		{
		  ok = 0;
		  pnearness.Elem(GetPointNr (i, j)) = minn+1;
		}
	  }
      }
    while (!ok);
    //(*testout) << "loadr4" << endl;

    lnearness.SetSize (noldl);

    for (i = 1; i <= noldl; i++)
      {
	lnearness.Elem(i) = 0;
	for (j = 1; j <= 2; j++)
	  lnearness.Elem(i) += pnearness.Get(GetPointNr (i, j));
      }
  }
  //(*testout) << "loadr5" << endl;

  oldutofreearea_i.SetSize (10);
  for (i = 0; i < oldutofreearea_i.Size(); i++)
    {
      oldutofreearea_i[i] = new DenseMatrix (oldutofreearea.Height(), oldutofreearea.Width());
      DenseMatrix & mati = *oldutofreearea_i[i];
      for (j = 0; j < oldutofreearea.Height(); j++)
	for (int k = 0; k < oldutofreearea.Width(); k++)
	  mati(j,k) = 1.0 / (i+1) * oldutofreearea(j,k) + (1 - 1.0/(i+1)) * oldutofreearealimit(j,k);
    }

  //(*testout) << "loadr6" << endl;
}




extern const char * triarules[];
extern const char * quadrules[];

void Meshing2 :: LoadRules (const char * filename)
{
  char buf[256];
  istream * ist;
  //char *tr1 = NULL;
  string tr1;

  /*
  ifstream ist (filename);
  if (!ist.good())
    {
      cerr << "Rule description file " << filename << " not found" << endl;
      exit (1);
    }
  */


  if (filename)
    {
      //      (*mycout) << "rule-filename = " << filename << endl;
      ist = new ifstream (filename);
    }
  else 
    {
      /* connect tetrules to one string */
      const char ** hcp;

      if (!mparam.quad)
	{
	  hcp = triarules;
	  PrintMessage (3, "load internal triangle rules");
	}
      else
	{
	  hcp = quadrules;
	  PrintMessage (3, "load internal quad rules");
	  // LoadRules ("rules/quad.rls");
	}

      size_t len = 0;
      while (*hcp)
	{
	  //	  (*testout) << "POS2 *hcp " << *hcp << endl;
	  len += strlen (*hcp);
	  hcp++;
	}
      //tr1 = new char[len+1];
      //tr1[0] = 0;
      tr1.reserve(len+1);


      if (!mparam.quad)
	hcp = triarules;
      else
	hcp = quadrules;


      //char * tt1 = tr1;
      while (*hcp)
	{
	  //strcat (tt1, *hcp);
	  //tt1 += strlen (*hcp);
	  tr1.append(*hcp);
	  hcp++;
	}
      
#ifdef WIN32
      // VC++ 2005 workaround
	  for(string::size_type i=0; i<tr1.size(); i++)
	if(tr1[i] == ',')
	  tr1[i] = ':';
#endif

      ist = new istringstream (tr1);
    }


  if (!ist->good())
    {
      cerr << "Rule description file " << filename << " not found" << endl;
      delete ist;
      exit (1);
    }
    
  while (!ist->eof())
    {
      buf[0] = 0;
      (*ist) >> buf;

      if (strcmp (buf, "rule") == 0)
	{
	  //(*testout) << "found rule" << endl;
	  netrule * rule = new netrule;
	  //(*testout) << "fr1" << endl;
	  rule -> LoadRule(*ist);
	  //(*testout) << "fr2" << endl;
	  
	  rules.Append (rule);
	}
      //(*testout) << "loop" << endl;
    }
  //(*testout) << "POS3" << endl;

  delete ist;
  //delete [] tr1;
}

}
