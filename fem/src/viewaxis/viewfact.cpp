/* This subroutine computes the viewfactors for an axisymmetric
   geometry. The code is written by Juha Katajam�ki while
   working for CSC. */

#include <math.h>
#include <stdio.h>
#include <stdlib.h>

#include "viewfact.h"
#include "../../config.h"
 
static Real r1, r2, z1, z2;     /* Katseltavan pinnan koordinaatit */
static Real r3, r4, z3, z4;     /* Katselevan pinnan koordinaatit */
static Real r12, r34, z12, z34; /* Keskilinjojen koordinaatit */
static Real zd1, zd3, zd;
static Real rd1, rd3;

static Real g1, g3, d1, d3, t, rratio;

static const Real eps = 1e-7, eps2 = 1.0e-14; /* eps*eps; */
static const Real delta = 1e-6; /* Suurin kosinien ero, joka aiheuttaa */
/* integroinnin */
static int nsurf,nsurfShade,inode,jnode;
static Real * coord;
static int *surfEltop, *surfEltopShade, * shadeParent;

static int compact = 1, verify = 0, selfshading = 1;

extern "C" void STDCALLBULL viewfactorsaxis
  (int *n,int *surf, Real *crd, Real *vf, int *idiv, int *fast)
{
  int i, j, ii, jj,div;
  Real a, sum, viewint, viewint2, vf2, sumdvf;
  Real c1, c2;    /* Kiertokulman kosinin yl�- ja alaraja */
  Real _r1, _r2, _r3, _r4, _z1, _z2, _z3, _z4;
  Real ds1,ds2,dp1,dz1,dz2,dr1,dr2,err,maxerr;
  Real epsilon = 1.0e-5;

  nsurf = *n;
  coord = crd;
  surfEltop = surf; 
  
  div = *idiv;
  compact = *fast;


  if(!compact) {
    printf("Using original set of boundary elements for shading\n");
    nsurfShade = nsurf;
    surfEltopShade = surfEltop;
  }

  if(compact) {
    int *nodehits,*nodetable;
    Real _r0, _z0;
    int ind0,ind1,ind2;

    printf("Combining original boundary elements for shading\n");

    int maxind = 0;
    for (i=0; i<2*nsurf; i++) 
      if(maxind < surfEltop[i]) maxind = surfEltop[i];
    // printf("Maximum node index is %d\n",maxind);

    // nodehits = new int[maxind+1];
    nodehits = (int*) malloc((maxind+1)*sizeof(int));
    for(i=0;i<=maxind;i++)
      nodehits[i] = 0;
    for (i=0; i<2*nsurf; i++)  
      nodehits[surfEltop[i]] += 1; 

    int maxnodehits = 0;
    for (i=0; i<=maxind; i++) 
      if(nodehits[i] > maxnodehits) maxnodehits = nodehits[i];
    // printf("Maximum node hits is %d\n",maxnodehits);
    
    int tablesize = (maxind+1)*maxnodehits;
    // nodetable = new int[tablesize];
    nodetable = (int*) malloc(tablesize*sizeof(int));
    for (i=0; i< tablesize; i++) 
      nodetable[i] = 0;

    for(i=0;i<=maxind;i++)
      nodehits[i] = 0;
    for (i=0; i<nsurf; i++) {
      ind1 = surfEltop[2*i+1];
      ind2 = surfEltop[2*i+0];
      nodetable[maxnodehits*ind1 + nodehits[ind1]] = i;
      nodetable[maxnodehits*ind2 + nodehits[ind2]] = i;
      nodehits[ind1] += 1;
      nodehits[ind2] += 1;
    }
 
    // surfEltopShade = new int[2*nsurf];
    surfEltopShade = (int*) malloc(2*nsurf*sizeof(int));
    for(i=0;i<2*nsurf;i++)
      surfEltopShade[i] = surfEltop[i];
 
    for (i=0; i<=maxind; i++) {
      int elem1,elem2;
      ind0 = i;
      
      if( nodehits[ind0] != 2) continue;

      elem1 = nodetable[maxnodehits*ind0+0];
      if( surfEltopShade[2*elem1+1] == ind0 ) 
	ind1 = surfEltopShade[2*elem1];
      else 
	ind1 = surfEltopShade[2*elem1+1];

      elem2 = nodetable[maxnodehits*ind0+1];
      if( surfEltopShade[2*elem2+1] == ind0 ) 
	ind2 = surfEltopShade[2*elem2];
      else 
	ind2 = surfEltopShade[2*elem2+1];

      _r0 = coord[2 * ind0];
      _r1 = coord[2 * ind1];
      _r2 = coord[2 * ind2];
      
      _z0 = coord[2 * ind0 + 1];
      _z1 = coord[2 * ind1 + 1];
      _z2 = coord[2 * ind2 + 1];

      dr1 = _r1 - _r0;
      dr2 = _r2 - _r0;
      dz1 = _z1 - _z0;
      dz2 = _z2 - _z0;
      
      dp1 = dr1 * dr2 + dz1 * dz2;
      ds1 = sqrt(dr1*dr1+dz1*dz1);
      ds2 = sqrt(dr2*dr2+dz2*dz2);
      
      dp1 /= (ds1*ds2);

      // Boundary elements mush be aligned
      if( dp1 > epsilon - 1. ) continue;

      // printf("Eliminating node %d\n",ind0);
      
      // Make the 1st element bigger 
      if( surfEltopShade[2*elem1] == ind0 ) 
	surfEltopShade[2*elem1] = ind2;
      else 
	surfEltopShade[2*elem1+1] = ind2;
      
      // Destroy the 2nd element 
      surfEltopShade[2*elem2] = 0;
      surfEltopShade[2*elem2+1] = 0;

      // Update the node information 
      nodehits[ind0] = 0;
      if( nodetable[maxnodehits*ind2] == elem2) 
	nodetable[maxnodehits*ind2] = elem1;
      else 
	nodetable[maxnodehits*ind2+1] = elem1;
    }

    // Free, not needed anymore
    free((char*)(nodetable));

    // Cannibalism of already used vector which does not need to be used again!
    shadeParent = nodehits;
    for (i=0; i<nsurf; i++) 
      shadeParent[i] = -1;
    
    j = 0;
    for (i=0; i<nsurf; i++) {
      if(surfEltopShade[2*i+1] || surfEltopShade[2*i+0]) {
	surfEltopShade[2*j+1] = surfEltopShade[2*i+1];
	surfEltopShade[2*j+0] = surfEltopShade[2*i+0];	
	j++;
      }
    }
    nsurfShade = j;
    printf("The combined set includes %d line segments (vs. %d)\n",nsurfShade,nsurf);


    // This is a dummy N^2 algorithm of finding the parent superelements
    // The info could also be inhereted in time of creating the superelements...
    
    int k, hit, parents = 0;
    maxerr = 0.;
    for (i=0; i<nsurfShade; i++) {
      
      ind1 = surfEltopShade[2*i];
      ind2 = surfEltopShade[2*i+1];
      
      _r1 = coord[2 * ind1];
      _r2 = coord[2 * ind2];
      
      _z1 = coord[2 * ind1 + 1];
      _z2 = coord[2 * ind2 + 1];
      
      dz1 = _z2 - _z1;
      dr1 = _r2 - _r1;
      ds1 = sqrt(dz1*dz1 + dr1*dr1);
      
      // Unit vector in direction of superelement
      dz1 /= ds1;
      dr1 /= ds1;

       
      for (j=0; j<nsurf; j++) {

	if( shadeParent[j] >= 0) continue;
	hit = 1;
	
	for(k=0;k<2;k++) {
	  
	  ind0 = surfEltop[2*j+k];
	  
	  // if node is joined, it is still a good candidate
	  // this check avoids also singularity at division
	  if(ind0 == ind1 || ind0 == ind1) continue;
	  
	  _r3 = coord[2 * ind0];
	  _z3 = coord[2 * ind0 + 1];
	  
	  dz2 = _z3 - _z1;
	  dr2 = _r3 - _r1;
	  ds2 = sqrt(dz2*dz2 + dr2*dr2);
	  
	  // Dot product of the superelement and the candidate-node element
	  dp1 = dz1*dz2 + dr1*dr2;
	  
	  // check that the node is on the line defined by the superelement
	  if( dp1 / ds2 < 1-epsilon ) {
	    hit = 0;
	    break;
	  }
	  
	  // check that node is within the segment of the superelement
	  if( dp1 / ds1 > 1+epsilon  || dp1 / ds1 < -epsilon) {
	    hit = 0;
	    break;
	  }
	}

	if(hit) {
	  shadeParent[j] = i;
	  parents++;
	}
      }
    }
    if(parents != nsurf) printf("Inconsistent number of parents found %d (vs. %d)\n",parents,nsurf);
    





    // delete [] nodetable;
    // delete [] nodehits;
  }


  // ************************************************************
  // The main N^2*M loop where M is the size of the shading table  
  for (i=0; i<nsurf; i++) {

    inode = i;    
    sum = 0.;
    sumdvf = 0.;

    _r3 = coord[2 * surfEltop[2*i+1]];
    _r4 = coord[2 * surfEltop[2*i+0]];
    
    _z3 = coord[2 * surfEltop[2*i+1]+1];
    _z4 = coord[2 * surfEltop[2*i+0]+1];

    
    a = Area(_r3, _r4, _z3, _z4);

    for (j=0; j<nsurf; j++) {

      jnode = j;
      _r1 = coord[2 * surfEltop[2*j+1]];
      _r2 = coord[2 * surfEltop[2*j+0]];
      
      _z1 = coord[2 * surfEltop[2*j+1]+1];
      _z2 = coord[2 * surfEltop[2*j+0]+1];
      
      vf[i*nsurf+j] = 0.;
      vf2 = 0.;
 
      if (a < eps) continue;

      for (ii=0; ii<div; ii++) {
	r3 = _r3 * (div - ii)/div + _r4 * ii/div;
	r4 = _r3 * (div - ii - 1)/div + _r4 * (ii + 1)/div;
	z3 = _z3 * (div - ii)/div + _z4 * ii/div;
	z4 = _z3 * (div - ii - 1)/div + _z4 * (ii + 1)/div;

	r34 = .5*(r3+r4);
	z34 = .5*(z3+z4);
	zd3=z3-z4;
	rd3=r3-r4;
	
	for (jj=0; jj<div; jj++) {
	  r1 = _r1 * (div - jj)/div + _r2 * jj/div;
	  r2 = _r1 * (div - jj - 1)/div + _r2 * (jj + 1)/div;
	  z1 = _z1 * (div - jj)/div + _z2 * jj/div;
	  z2 = _z1 * (div - jj - 1)/div + _z2 * (jj + 1)/div;
	  r12 = .5*(r1+r2);
	  if ( r12 < eps || r34 < eps) continue;

	  if (r1 < eps) r1 = eps;
	  if (r2 < eps) r2 = eps;
	  if (r3 < eps) r3 = eps;
	  if (r4 < eps) r4 = eps;
	  
	  zd1=z1-z2;
	  rd1=r1-r2; 
	  z12 = .5*(z1+z2);
	  zd = z12-z34;
	  
	  if (!InitialInterval(&c1, &c2)) continue;	  
	  viewint = ViewIntegral(c1, c2, 0);

	  
	  // Code for verification
	  if(verify) {
	    int *surfEltopTmp;
	    int nsurfTmp;
	    
	    surfEltopTmp = surfEltopShade;
	    nsurfTmp = nsurfShade;	    
	    
	    surfEltopShade = surfEltop;
	    nsurfShade = nsurf;
	    	    
	    if (!InitialInterval(&c1, &c2)) continue;	  
	    viewint2 = ViewIntegral(c1, c2, 0);	    
	    vf2 = vf2 + 4. * viewint2;

	    surfEltopShade = surfEltopTmp;
	    nsurfShade = nsurfTmp;
	  }	      	    

	  vf[i*nsurf+j] += 4. * viewint;
	  /* Kerroin 4 koostuu tekij�ist� 2 (peilisymmetria), 2pi */
	  /* (kiertosymmetria) ja 1/pi (integraalin lausekkeessa esiintyv� */
	  /* vakio) */
	}
      }

      vf[i*nsurf+j] /= a;
      sum += vf[i*nsurf+j];

      if(verify) {
	vf2 /= a;
	sumdvf += ( fabs(vf[i*nsurf+j]-vf2 ) / a);
      }
    }

    if(verify) {
      if(sumdvf > maxerr) maxerr = sumdvf;    
      printf("Line sum: %d %g %g %g\n", i, sum, sumdvf, maxerr);
    }
    else {
      //      printf("Line sum: %d %g\n", i, sum );      
    }
  }


  // Deallocate stuff
  if(compact) {
    free((char*)(surfEltopShade));
    free((char*)(shadeParent));
  }

}


BOOL InitialInterval(Real *c1, Real *c2)
{
  /* M��r�� rajat katseltavan pisteen kiertokulman kosinille ehdosta, ett� */
  /* yhdysjanan ja pintojen normaalien v�listen kulmien on oltava < pi/2.  */
  /* Palauta FALSE, jos ratkaisujoukko on tyhj� tai nollamittainen, */
  /* muutoin TRUE. */
  /* Funktio olettaa, ett� r12 ja r34 eiv�t ole nollia. */ 
  
  Real cc1, cc3; 
  
  *c1 = -1.; *c2 = 1.;
  if ( fabs(zd1) > eps ) {
    cc1 = (- zd * rd1 + r12 * zd1) / (r34 * zd1);        
    if ( fabs(zd3) > eps ) {
      cc3 = (zd * rd3 + r34 * zd3) / (r12 * zd3);        
      if (zd1 > 0.) {
	if (zd3 > 0.) *c1 = max(cc1, cc3);
	else { *c1 = cc1; *c2 = cc3; }
      } else {
	if (zd3 < 0.) *c2 = min(cc1, cc3);
	else { *c1 = cc3; *c2 = cc1; }
      }
    } else {
      if ( sgn(rd3) && sgn(rd3) == -sgn(zd) ) {
	if (zd1 > 0.) *c1 = cc1;
	else *c2 = cc1;
      } else { *c1 = 1.; *c2 = -1.; } /* Joukko tyhj� */
    }
  } else {
    if ( fabs(zd3) > eps ) {
      cc3 = (zd * rd3 + r34 * zd3) / (r12 * zd3);
      if ( sgn(rd1) && sgn(rd1) == sgn(zd) ) {
	if (zd3 > 0.) *c1 = cc3;
	else *c2 = cc3;
      } else { *c1 = 1.; *c2 = -1.; } /* Joukko tyhj� */
    } else {
      if ( !sgn(rd1) || sgn(rd1) != sgn(zd) || sgn(rd1) != -sgn(rd3) )
	{ *c1 = 1.; *c2 = -1.; }  /* Muutoin joukko = [-1, 1] */
    }
  }
  
  *c1 = max(-1.+eps, *c1); *c2 = min(1.-eps, *c2);
  /* Epsilonilla estet��n nollalla jako integroinnissa */
  if (*c2 - *c1 < eps) return FALSE;
  return TRUE;
}


Real ViewIntegral (Real c1, Real c2, int k)
{
  /*
    T�m� funktio laskee view factorin yhdelle elementtiparille.
    Integrointialuetta rajoitetaan tutkimalla kartiopintojen aiheuttama
    varjostus. Jos integrointialue jakautuu kahtia, suoritetaan rekursiivinen
    kutsu molemmille osille. Jos integrointialue kutistuu mit�tt�m�ksi
    tai tyhj�ksi, palautetaan nolla.
    Funktio olettaa globaalit muuttujat r12 ja r34 nollasta poikkeaviksi.
    */

  static Real r5, r6, z5, z6;    /* Varjostavan pinnan reunojen koordinaatit */
  static Real zd5, t1, tt1, t2, tt2, t0;
  Real cc1,cc2;
  rratio = r34/r12;

  while (k < nsurfShade) {
    
    r5 = coord[2 * surfEltopShade[2*k+1]];
    r6 = coord[2 * surfEltopShade[2*k+0]];
    
    z5 = coord[2 * surfEltopShade[2*k+1]+1];
    z6 = coord[2 * surfEltopShade[2*k+0]+1];

    k++;

    // either element cannot shade one another
    if(selfshading) {
      // Condition for superelements
      if( nsurf != nsurfShade ) {
	if( k-1 == shadeParent[inode] || k-1 == shadeParent[jnode]) continue;
      }
      else if( nsurf == nsurfShade ) {
	if( k-1 == inode || k-1 == jnode) continue;
      }
    }

    if (r5+r6 < eps) continue;


    zd5 = z5-z6;
    if ( fabs(zd5) < eps ) {
      /* Varjostava pinta on tasorengas */
      
      /* Tasorengas ei voi varjostaa itse��n */
      /* T�m� lis�ys korjaa alirutiinissa pitk��n ollen bugin (P.R. 23.4.2004) */
      if(nsurf == nsurfShade && inode == k-1) continue;

      if ( fabs(zd) < eps ) continue;

      t1 = (z12-z5)/zd; 
      tt1 = 1.-t1;
      if (t1 < eps || tt1 < eps) continue;


      t = rratio * t1/tt1;
      cc1 = .5*(r5*r5/(r12*r34*t1*tt1) - t - 1./t);
      cc2 = .5*(r6*r6/(r12*r34*t1*tt1) - t - 1./t);

      if (cc1 > cc2) { t = cc1; cc1 = cc2; cc2 = t; }
    } 
    else  {
      /* Varjostava pinta on kartio tai lieri�       */
      /* Laske yhdysjanasta varjoon j��v� v�li z-suunnassa  */

      if ( fabs(zd) < eps ) {
	if ( (z12-z5 < eps && z12-z6 > eps) ||
	     (z12-z5 > eps && z12-z6 < eps) )
	  { t1 = 0.; t2 = 1; }
	else continue;
      } else {
	t1 = (z12-z5)/zd; 
	t2 = t1 + zd5/zd;
	if (t1 > t2) { t = t1; t1 = t2; t2 = t; }
      }

      if (! IntervalIsect(0., 1., t1, t2, &t1, &t2)) continue;
      tt1 = 1.-t1; 
      tt2 = 1.-t2;
      
      /* Laske, mit� arvoja kiertokulman kosini saa v�lill� [t1, t2] */
      cc1 = 1.; cc2 = -1.;
      g1 = (r5 * (z12-z6) - r6 * (z12-z5)) / (r12 * zd5);
      g3 = (r5 * (z34-z6) - r6 * (z34-z5)) / (r34 * zd5);
      d1 = g1*g1 - 1; 
      d3 = g3*g3 - 1;  
      /* N�m� ilmaisevat, kummalla */
      /* puolen kartiota ovat katseleva ja katseltava piste */

      /* Tutki v�lin p��tepiste */
      ExaminePoint (t1, &cc1, &cc2);
      ExaminePoint (t2, &cc1, &cc2);

      /* Jos kumpikin piste kartion ulkopuolella, tutki derivaatan */
      /* nollakohta, mik�li se on v�lill� [t1, t2] */
      if (d1 <= -eps && d3 <= -eps) {
	t0 = 1. / (1. + sqrt(rratio * d3/d1));
	if (t0 - t1 > eps && t2 - t0 > eps) {
	  ExaminePoint(t0, &cc1, &cc2);
	}
      }
      if (cc1 > cc2) {
	cc1 = cc2; /* N�in voi k�yd� py�ristysvirheiden takia */
      }

    }


    if (IntervalIsect(c1, c2, cc1, cc2, &cc1, &cc2)) {

      if (cc1 - c1 < delta) {
	if (c2 - cc2 < delta) {
	  return 0.;
	}
	else {
	  c1 = cc2;
	}
      }
      else if (c2 - cc2 < delta) {
	c2 = cc1;
      }
      else {
	return ViewIntegral(c1, cc1, k) + ViewIntegral(cc2, c2, k);
      }
    }
  }
  return Integrate(c1, c2);
}


BOOL IntervalIsect(Real x1, Real x2, Real y1, Real y2, Real *z1, Real *z2)
{
  /* Laske v�lien [x1, x2] ja [y1, y2] leikkaus ja palauta FALSE, jos */
  /* t�m� on tyhj� tai mit�t�n. Input-parametrien j�rjestyksen on oltava */
  /* oikea.*/
  
  *z1 = x1; *z2 = x2;
  if (x2 - y1 < eps) return FALSE;
  if (y1 - x1 > eps) *z1 = y1;
  if (y2 - x1 < eps) return FALSE;
  if (x2 - y2 > eps) *z2 = y2;
  return (*z2 - *z1 >= eps);
}


void ExaminePoint (Real x, Real *mi, Real *ma)
{
  Real y;
  if (x > eps) {
    if (1.-x > eps) {
      t = rratio*x/(1.-x);
      y = .5*(d1/t + d3*t) + g1*g3;
    } else 
      if ( fabs(d3) < eps ) y = g1*g3;
      else y = sgn(d3);
  } else
    if ( fabs(d1) < eps ) y = g1*g3;
    else y = sgn(d1);
  if (y > *ma) *ma = y;
  if (y < *mi) *mi = y;
}


Real Integrate(Real c1, Real c2)
{
  /* c1 ja c2 ovat integrointiv�lin kulman kosinin rajat. */ 
  /* Integraali lasketaan ilman nimitt�j�n pi-tekij��. */
  
  /* Ensimm�inen ja viimeinen integrointipiste eiv�t saa olla tasan */
  /* 0 ja 1, jottei vierekk�isten elementtien tapauksessa tule */
  /* nollalla jakoa */
  /*	static const Real qp[] = { 1e-6, .25, .5, .75, 1.-1e-6 }, */
  /*    					w[] = { 1./12., 1./3., 1./6., 1./3., 1./12. }; */
  /*  static const Real qp[] = { 0.211324865, 0.788675134 }, */
  /*	    w[] = { .5, .5 }; */
  static const Real qp[] = { 0.112701665, 0.5, 0.887298334 },
			     w[] = { 0.277777777, 0.444444444, 0.277777777 };
  static const int nqp = 3;
    
  int i;
  Real c = zd1*zd1 + rd1*rd1;
  if (c < eps2) return 0.; /* Pinta kutistunut ympyr�nkaareksi; t�m� testi */
  /* tarvitaan nollalla jaon v�ltt�miseksi */
  
  Real z, r, h, hh1, hh2, g1, g2, gg1, gg2, value, integral;
  Real d1, d2, e1, e2, f1, f2;
  Real zrd = r2*z1-r1*z2;
  Real a1 = rd3*r1, a2 = rd3*r2;
  Real b1 = zd3*z1, b2 = zd3*z2; 
  Real s1 = sqrt(1. - c1*c1), s2 = sqrt(1. - c2*c2);
  /* kosineissa ja sineiss� indeksit 1 ja 2 toisin p�in kuin */
  /* muissa muuttujissa! */
  Real cs = (1.+c1)*(1.+c2), cd = (1.-c1)*(1.-c2);

  
  integral = 0.;
  for (i=0; i<nqp; i++) {
    z = z3 - qp[i] * zd3;  /* qp on integroimismuuttuja */
    r = r3 - qp[i] * rd3;
    e1 = (z1-z)*(z1-z) + r1*r1 + r*r;
    f1 = (z2-z)*(z2-z) + r2*r2 + r*r;
    hh1 = 2*r1*r;
    hh2 = 2*r2*r;
    g1 = - e1 / hh1;
    g2 = - f1 / hh2;
    e2 = e1 - c1*hh1;
    f2 = f1 - c1*hh2;
    e1 -= c2*hh1;
    f1 -= c2*hh2;
    h = zd3*z + rd3*r;
    gg1 = (g1+c2)*(g1+c1);
    gg2 = (g2+c2)*(g2+c1);
    
    /* Kaarien osuus: */
    value  = (-.5 * (a1 + (h-b1)*g1) / sqrt(g1*g1-1) ) *
      acos( .5 * ( (1.-g1) * sqrt(cd/gg1) - 
		   (1.+g1) * sqrt(cs/gg1) ) );
    value -= (-.5 * (a2 + (h-b2)*g2) / sqrt(g2*g2-1) ) *
      acos( .5 * ( (1.-g2) * sqrt(cd/gg2) - 
		   (1.+g2) * sqrt(cs/gg2) ) );
    value += .25 * (b1-b2) * acos(c1*c2 + s1*s2);

    /* Suorien sivujen osuus: */
    gg1 = e1+f1-c; gg2 = e2+f2-c;
    hh1 = 4*e1*f1; hh2 = 4*e2*f2;
    d1 = hh1 - gg1*gg1; d2 = hh2 - gg2*gg2;
    h = r * (rd1*h + zrd*zd3);
    value -= h * (s1 / sqrt(d2)) * acos( gg2 / sqrt(hh2) );
    value += h * (s2 / sqrt(d1)) * acos( gg1 / sqrt(hh1) );

    integral += w[i] * value;
  }

  return integral;
}

Real Area(Real r1, Real r2, Real z1, Real z2)
{
    return 3.1415926535 * (r1+r2) *
        sqrt( (z1-z2)*(z1-z2) + (r1-r2)*(r1-r2) );
}


