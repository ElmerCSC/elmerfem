#ifndef _CSGPARSER_HPP
#define _CSGPARSER_HPP




//namespace netgen
//{
  enum TOKEN_TYPE
    { 
      TOK_MINUS = '-', TOK_LP = '(', OK_RP = ')', TOK_LSP = '[', TOK_RSP = ']',
      TOK_EQU = '=', TOK_COMMA = ',', TOK_SEMICOLON = ';',
      TOK_NUM = 100, TOK_STRING, TOK_NAMED_SOLID, TOK_PRIMITIVE, 
      TOK_OR, TOK_AND, TOK_NOT, 
      TOK_SINGULAR, TOK_EDGE, TOK_POINT, TOK_FACE, TOK_IDENTIFY, TOK_CLOSESURFACES,
      TOK_CLOSEEDGES, TOK_PERIODIC,
      TOK_SOLID, TOK_RECO, TOK_TLO, TOK_CURVE2D, TOK_CURVE3D, TOK_BOUNDINGBOX,
      TOK_BOUNDARYCONDITION, TOK_BOUNDARYCONDITIONNAME,
      TOK_DEFINE, TOK_CONSTANT,
      TOK_END };

  struct kwstruct
  {
    TOKEN_TYPE kw; 
    const char * name;
  };

  enum PRIMITIVE_TYPE
    {
      TOK_SPHERE = 1, TOK_CYLINDER, TOK_PLANE, TOK_ELLIPTICCYLINDER, 
      TOK_ELLIPSOID, TOK_CONE, 
      TOK_ORTHOBRICK, TOK_POLYHEDRON, 
      TOK_TORUS,
      TOK_TUBE, TOK_GENCYL, TOK_EXTRUSION, TOK_REVOLUTION,

      TOK_TRANSLATE, TOK_MULTITRANSLATE, TOK_ROTATE, TOK_MULTIROTATE
    };

  struct primstruct
  {
    PRIMITIVE_TYPE kw; 
    const char * name;
  };


  class CSGScanner
  {
    TOKEN_TYPE token;
    PRIMITIVE_TYPE prim_token;
    double num_value;
    string string_value;
    
    int linenum;
    istream * scanin;

  public:

    CSGScanner (istream & ascanin);

    TOKEN_TYPE GetToken() const
    { return token; }

    double GetNumValue() const
    { return num_value; }

    const string & GetStringValue() const
    { return string_value; }

    char GetCharValue() const
    { return string_value[0]; }

    PRIMITIVE_TYPE GetPrimitiveToken() const
    { return prim_token; }
  
    void ReadNext();

    /*
    CSGScanner & Parse (char ch);
    CSGScanner & Parse (int & i);
    CSGScanner & Parse (double & d);
    CSGScanner & Parse (Point<3> & p);
    CSGScanner & Parse (Vec<3> & p);
    */
    void Error (const string & err);
  };


  
  CSGScanner & operator>> (CSGScanner & scan, char ch);
  CSGScanner & operator>> (CSGScanner & scan, double & d);
  CSGScanner & operator>> (CSGScanner & scan, int & i);
  CSGScanner & operator>> (CSGScanner & scan, Point<3> & p);
  CSGScanner & operator>> (CSGScanner & scan, Vec<3> & v);
  


//}









#endif // _CSGPARSER_HPP
