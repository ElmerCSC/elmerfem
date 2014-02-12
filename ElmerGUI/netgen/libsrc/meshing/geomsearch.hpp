#ifndef FILE_GEOMSEARCH
#define FILE_GEOMSEARCH

/**************************************************************************/
/* File:   geomsearch.hh                                                  */
/* Author: Johannes Gerstmayr                                             */
/* Date:   19. Nov. 97                                                    */
/**************************************************************************/

class FrontPoint3;
class FrontFace;
class MiniElement2d;

  /// class for quick access of 3D-elements; class cannot delete elements, but only append
class GeomSearch3d
{

public:
  ///
  GeomSearch3d();
  ///
  virtual ~GeomSearch3d();

  ///
  void Init (ARRAY <FrontPoint3,PointIndex::BASE> *pointsi, ARRAY <FrontFace> *facesi);

  ///get elements max extension
  void ElemMaxExt(Point3d& minp, Point3d& maxp, const MiniElement2d& elem);
  
  ///get minimum coordinates of two points ->p2
  void MinCoords(const Point3d& p1, Point3d& p2);

  ///get minimum coordinates of two points ->p2
  void MaxCoords(const Point3d& p1, Point3d& p2);

  ///create a hashtable from an existing array of triangles
  ///sizei = number of pieces in one direction
  void Create();

  ///add new element to Hashtable
  void AddElem(const MiniElement2d& elem, INDEX elemnum);

  ///GetLocal faces in sphere with radius xh and middlepoint p
  void GetLocals(ARRAY<MiniElement2d> & locfaces,  ARRAY<INDEX> & findex,
		 INDEX fstind, const Point3d& p0, double xh);

private:
  
  ARRAY <FrontFace> *faces; // Pointers to Arrays in Adfront
  ARRAY <FrontPoint3,PointIndex::BASE> *points;

  ARRAY <ARRAY <int>*> hashtable;

  Point3d minext; //extension of Hashdomain
  Point3d maxext;
  Point3d maxextreal;
  Vec3d elemsize;  //size of one Hash-Element

  threeint size; // size of Hashtable in each direction
  int reset;
  int hashcount;
};

#endif





















































