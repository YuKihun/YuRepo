#include "CDirectMethod.h"


CDirectMethod::CDirectMethod()
{
	
}

CDirectMethod::~CDirectMethod()
{
	
}

/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/* Vincenty Direct Solution of Geodesics on the Ellipsoid (c) Chris Veness 2005-2011              */
/*                                                                                                */
/* from: Vincenty direct formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the  */
/*       Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975    */
/*       http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf                                             */
/* - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  */
/**
 * Calculates destination point given start point lat/long, bearing & distance, 
 * using Vincenty direct formula for ellipsoids
 *
 * e pointer to ellipsoid struct
 * lat1, lon1, brng in radians 
 * s in meters
 * lat2, lon2, revAz in radians
 */
int CDirectMethod::VincentyDirect( sEllipsoid * e, double lat1, double lon1, double alpha1, 
		double s, double *lat2, double *lon2, double * revAz) 
{

	lat1 = lat1 * D2R;
	lon1 = lon1 * D2R;
	alpha1 = alpha1 * D2R;
	
	double sinAlpha1 = sin(alpha1);
	double cosAlpha1 = cos(alpha1);

	double tanU1 = (1 - e->dFlatness) * tan(lat1);
	double cosU1 = 1 / sqrt((1 + tanU1 * tanU1)), sinU1 = tanU1 * cosU1;
	double sigma1 = atan2(tanU1, cosAlpha1);
	double sinAlpha = cosU1 * sinAlpha1;
	double cosSqAlpha = 1 - sinAlpha * sinAlpha;
	double uSq = cosSqAlpha * (e->dMajorAxis * e->dMajorAxis - e->dMinorAxis * e->dMinorAxis) / (e->dMinorAxis * e->dMinorAxis);
	double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
	double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));

	double sigma = s / (e->dMinorAxis * A), sigmaP = 2 * PI;
	double cos2SigmaM = NaN, sinSigma = NaN, cosSigma = NaN, deltaSigma = NaN;
	double eps = fabs(sigma - sigmaP);

	while( eps > 1e-12 ) 
	{
		cos2SigmaM = cos(2 * sigma1 + sigma);
		sinSigma = sin(sigma);
		cosSigma = cos(sigma);
		deltaSigma = B * sinSigma * 
				( cos2SigmaM + B / 4 * ( cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM) - B / 6 * cos2SigmaM * 
						(-3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM) ) );
		sigmaP = sigma;
		sigma = s / (e->dMinorAxis * A) + deltaSigma;
		eps = fabs(sigma - sigmaP);
	}

	double tmp = sinU1 * sinSigma - cosU1 * cosSigma * cosAlpha1;
	*lat2 = atan2(sinU1 * cosSigma + cosU1 * sinSigma * cosAlpha1,
			(1 - e->dFlatness) * sqrt(sinAlpha * sinAlpha + tmp * tmp)) * R2D; //Degree
	double lambda = atan2(sinSigma * sinAlpha1,
			cosU1 * cosSigma - sinU1 * sinSigma * cosAlpha1);
	double C = e->dFlatness / 16 * cosSqAlpha * (4 + e->dFlatness * (4 - 3 * cosSqAlpha));
	double L = lambda - (1 - C) * e->dFlatness * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM
													+ C * cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM)));

	*lon2 = ( fmod((lon1 + L + 3 * PI), (2 * PI)) - PI ) * R2D; ; // normalise to -180...+180
	*revAz = ( atan2(sinAlpha, -tmp) ) * R2D; // final bearing, if required

	return true;
}
