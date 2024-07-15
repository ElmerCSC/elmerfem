/* 3d.h */


extern void	glaux_error(char *);
extern void	diff3(GLdouble [3], GLdouble [3], GLdouble [3]);
extern void	add3(GLdouble [3], GLdouble [3], GLdouble [3]);
extern void	scalarmult(GLdouble, GLdouble [3], GLdouble [3]);
extern GLdouble	dot3(GLdouble [3], GLdouble [3]);
extern GLdouble	length3(GLdouble [3]);
extern GLdouble	dist3(GLdouble [3], GLdouble [3]);
extern void	copy3(GLdouble [3], GLdouble [3]);
extern void	crossprod(GLdouble [3], GLdouble [3], GLdouble [3]);
extern void	normalize(GLdouble [3]);
extern void	print3(GLdouble [3]);
extern void	printmat3(GLdouble [3][3]);
extern void	identifymat3(GLdouble [3][3]);
extern void	copymat3(GLdouble *, GLdouble *);
extern void	xformvec3(GLdouble [3], GLdouble [3][3], GLdouble [3]);

extern void m_resetmatrixstack(void);
extern void m_xformpt(GLdouble [3], GLdouble [3], GLdouble [3], GLdouble [3]);
extern void m_xformptonly(GLdouble [3], GLdouble [3]);
extern void m_pushmatrix(void);
extern void m_popmatrix(void);
extern void m_shear(GLdouble, GLdouble, GLdouble);
extern void m_translate(GLdouble, GLdouble, GLdouble);
extern void m_scale(GLdouble, GLdouble, GLdouble);
extern void m_rotate(GLdouble, char);


