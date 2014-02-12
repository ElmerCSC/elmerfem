%{
#include <mystdlib.h>
#include <myadt.hpp> 

#include <linalg.hpp> 
#include <csg.hpp>




// extern SYMBOLTABLE<Solid*> solids;
namespace netgen {
extern CSGeometry * parsegeom;
}
using namespace netgen;

#include "geometry.h"


ARRAY<char*> parsestrings;
int linenum;
%}

dig     [0-9]
id      [a-zA-Z][a-zA-Z0-9]*
num1    [-+]?{dig}+\.?([eE][-+]?{dig}+)?
num2    [-+]?{dig}*\.{dig}+([eE][-+]?{dig}+)?
number  {num1}|{num2}
%x      incl
%x      comment

%%
algebraic3d     { return TOK_RECO; }
solid           { return TOK_SOLID; }
tlo		{ return TOK_TLO; }
and             { return TOK_AND; }
or              { return TOK_OR; }
not             { return TOK_NOT; }
translate	{ return TOK_TRANSLATE; }
multitranslate	{ return TOK_MULTITRANSLATE; }
rotate		{ return TOK_ROTATE; }
multirotate	{ return TOK_MULTIROTATE; }
sphere          { return TOK_SPHERE; }
cylinder        { return TOK_CYLINDER; }
cone		{ return TOK_CONE; }
plain           { return TOK_PLAIN; }
plane		{ return TOK_PLAIN; }
tube	 	{ return TOK_TUBE; }
gencyl		{ return TOK_GENCYL; }
orthobrick	{ return TOK_ORTHOBRICK; }
polyhedron	{ return TOK_POLYHEDRON; }
revolution	{ return TOK_REVOLUTION; }

singular	{ return TOK_SINGULAR; }
edge		{ return TOK_EDGE; }
point		{ return TOK_POINT; }

identify	{ return TOK_IDENTIFY; }
closesurfaces 	{ return TOK_CLOSESURFACES; }
closeedges 	{ return TOK_CLOSEEDGES; }
periodic	{ return TOK_PERIODIC; }
boundarycondition { return TOK_BOUNDARYCONDITION; }
boundingbox	{ return TOK_BOUNDINGBOX; }

{number}        { yylval.val = atof (YYText()); return NUM; }
{id}+           {
                  yylval.chptr = new char [YYLeng()+1];
		  parsestrings.Append (yylval.chptr);
                  strcpy (yylval.chptr, YYText());
                  if (parsegeom->GetSolid (yylval.chptr))
                    return IDENTSOLID;
                  else
                    return IDENT;
                }
[ \t]             /* eat up ws */
.               { return int(*YYText()); }
\n		{ linenum++; }
"##".*\n        { linenum++; cout << (YYText()+2) ; }  /* line comment */
"#".*\n        { linenum++; }  /* line comment */


%%

extern FlexLexer * lexer;

int yylex ()
  {
  return lexer -> yylex();
  }

extern "C" int yywrap ()
  {
  return 1;
  }
