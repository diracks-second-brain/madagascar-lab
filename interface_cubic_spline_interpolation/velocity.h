#include <rsf.h>

#ifndef __VELOCITY__H__
#define __VELOCITY__H__

void calcInterfacesZcoord(      float *zi, /* Interfaces depth coordinates */
                                int nint, /* Number of interfaces */
                                float xs, /* x coordinate */
                                int si, /* Spline index */
                                float **coef /* Cubic spline coefficients */)
/*< Calculate depth coordinates of the interfaces
 * Note: This function calculates interfaces depth coordinates and stores it
 * in the zi vector.
  >*/
{
        int i; // Loop counter

        for(i=0;i<nint;i++)
                zi[i] = coef[i][si*4+0]*xs*xs*xs+coef[i][si*4+1]*xs*xs+coef[i][si*4+2]*xs+coef[i][si*4+3];
}

void calculateSplineCoefficients(int n, /* Vectors (x,y) dimension */
                                float* x, /* x coordinates */
                                float* y, /* y coordinates */
                                float* coef /* Spline coefficients */)
/*< Function to calculate natural cubic spline coefficients

Note: It Receives n points and two vectors x and y with n dimension.
It returns a coefficients vector with 4 coefficients for each of the
n-1 natural cubic splines, coef[(n-1)*4].

IMPORTANT: The number of points must be equal or major than 3 (n>3)
and x vector must be in crescent order.

>*/
{

        float s2[n]; // Second derivatives matrix
        int i, ip1, ip2, im1, m; // Loop counter
        float hb, ha, deltaa, deltab, t; // temporary variables
        float e[n-2]; // hi's vector
        float dp[n-2]; // main diagonal

        /* Vectors dimension must be major than 3 */
        if(n<3) sf_error("Vectors dimension n must be major than 3\n");

        /* x vector must be in crescent order */
        for(i=1;i<n;i++){
                if(x[i-1]>x[i]) sf_error("Vector x should be in ascending order\n");
        }

        /* Simetric tridiagonal linear system build */
        ha = x[1]-x[0]; deltaa = (y[1]-y[0])/ha; m=n-2;
        for(i=0;i<m;i++){
                ip1 = i+1; ip2 = i+2;
                hb = x[ip2]-x[ip1];
                deltab = (y[ip2]-y[ip1])/hb;
                e[i] = hb; dp[i] = 2*(ha+hb);
                s2[ip1] = 6*(deltab-deltaa);
                ha=hb; deltaa=deltab;
        }

        /* Gauss elimination */
        for(i=1;i<m;i++){
                ip1=i+1; im1=i-1;
                t = e[im1]/dp[im1];
                dp[i] = dp[i]-t*e[im1];
                s2[ip1] = s2[ip1]-t*s2[i];
        }

        /* Retroactive substitutive solution */
        s2[m]=s2[m]/dp[m-1];
        for(i=m-1;i>0;i--){
                ip1=i+1; im1=i-1;
                s2[i]=(s2[i]-e[im1]*s2[ip1])/dp[im1];
        }
        s2[0]=0; s2[n-1]=0;
        /* Calculate spline coefficients */
        for(i=0;i<n-1;i++){
                ha = x[i+1]-x[i];
                coef[0+i*4] = (s2[i+1]-s2[i])/(6*ha);
                coef[1+i*4] = s2[i]/2;
                coef[2+i*4] = (y[i+1]-y[i])/ha-(s2[i+1]+2*s2[i])*(ha/6);
                coef[3+i*4] = y[i];
        }
}

void updateVelocityModel(
                           float *vel,
                           int *n, /* Velocity model dimension n1=n[0] n2=n[1] */
                           float *o, /* Velocity model axis origin o1=o[0] o2=o[1] */
                           float *d, /* Velocity model sampling d1=d[0] d2=d[1] */
                           float *sv, /* Velocity model disturbance */
                           int *nsv, /* Dimension of sv the vector */
                           float *osv,
                           float *dsv,
			   float *sz,
			   int *nsz,
			   float *osz,
			   float *dsz,
			   int itf,
			   float *vl)
/*< Velocity model update
Note: This function uses a sv (layers velocity) vector and sz (depth interfaces
coordinates) vector to build the depth velocity model. There is nsv constant
velocity layers in the model and nsv-1 interfaces separating them.
These interfaces are described with nsz control points (nodes) in the sz vector and
they are interpolated using natural cubic spline interpolation.
 >*/
{
        sf_eno2 map;
        float f2[2];
        int i, j, i1, i2; 
        float x, y;
        float *vprof;
	float *szz, *xs;
	float **coef;
	float zi[1];
	float zs;
	int l=0;
	float xx;
	int rect[2]={10,10}, nrep=2;

        vprof=sf_floatalloc(nsv[0]*nsv[1]);

     	xs = sf_floatalloc(nsz[0]);
	szz = sf_floatalloc(nsz[0]);

	for(i=0;i<nsz[0];i++){
		xs[i] = i*dsz[0]+osz[0];
		szz[i] = sz[(itf-1)*nsz[0]+i];
	}
		
        /* Calculate coefficients matrix (interfaces interpolation) */
	coef = sf_floatalloc2(4*(nsz[0]-1),1);
	calculateSplineCoefficients(nsz[0],xs,szz,coef[0]);

      	/* Calculate velocity function */
	for(j=0;j<nsv[1];j++){
		xx = dsv[1]*j+osv[1];
		l = (int) (xx-osz[0])/dsz[0];
		/* Calculate interfaces z coordinates */
		calcInterfacesZcoord(zi,1,xx-xs[l],l,coef);
		for(i=0;i<nsv[0];i++){
			zs = i*dsv[0]+osv[0];
			if(zs>zi[0]){
				vel[nsv[0]*j+i] = vl[itf];
			}else{
				vel[nsv[0]*j+i] = vl[itf-1];
			}
		} /* Loop over depth */
	} /* Loop over distance */

	free(coef);
	free(xs);
	free(szz);

}

#endif