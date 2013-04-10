#include "CInverseMethod.h"

//===========================================================================
// Vincenty 알고리즘
// 두위치의거리와방위각을산출
//===========================================================================
CInverseMethod::CInverseMethod()
{
	m_dDistance = 0.0;
	m_dFwdAz = 0.0;
	m_dRevAz = 0.0;
}

CInverseMethod::~CInverseMethod( )
{
	
}
/**
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  
 Vincenty Inverse Solution of Geodesics on the Ellipsoid (c) Chris Veness 2002-2011             
                                                                                                
 from: Vincenty inverse formula - T Vincenty, "Direct and Inverse Solutions of Geodesics on the 
       Ellipsoid with application of nested equations", Survey Review, vol XXII no 176, 1975    
       http://www.ngs.noaa.gov/PUBS_LIB/inverse.pdf                                             
 - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - - -  

 * Calculates geodetic distance between two points specified by latitude/longitude using 
 * Vincenty inverse formula for ellipsoids
 *
 * source: http://www.movable-type.co.uk/scripts/latlong-vincenty.html
 *
 * e pointer to ellipsoid struct
 * lat, lon in radians, 
 * fwdAz, revAz in radians, 
 * s in meters 
 * 
*/
//int CInverseMethod::VincentyInverse( sEllipsoid *e, double lat1, double lon1, double lat2,
//		double lon2, double *s, double *fwdAz, double *revAz) 

//위경도만 타원체 정보만 입력하도록 수정
int CInverseMethod::VincentyInverse( sEllipsoid *e, double lat1, double lon1, double lat2, double lon2 )
{

	lat1 = lat1 * D2R;
	lon1 = lon1 * D2R;
	lat2 = lat2 * D2R;
	lon2 = lon2 * D2R;

	double L = (lon2 - lon1);
	double U1 = atan((1 - e->dFlatness) * tan(lat1));
	double U2 = atan((1 - e->dFlatness) * tan(lat2));
	double sinU1 = sin(U1), cosU1 = cos(U1);
	double sinU2 = sin(U2), cosU2 = cos(U2);

	double lambda = L, lambdaP, iterLimit = 100.0;
	double cosSqAlpha, cosSigma, sigma, cos2SigmaM, sinLambda, sinSigma, cosLambda, sinAlpha;
	double eps = 1000;
	do 
	{
		sinLambda = sin(lambda);
		cosLambda = cos(lambda);
		sinSigma = sqrt( (cosU2 * sinLambda) * (cosU2 * sinLambda) + (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda) * (cosU1 * sinU2 - sinU1 * cosU2 * cosLambda));
		if (sinSigma == 0) 
		{
			m_dDistance = 0.0;
			return true;
		} 

		cosSigma = sinU1 * sinU2 + cosU1 * cosU2 * cosLambda;
		sigma = atan2(sinSigma, cosSigma);
		sinAlpha = cosU1 * cosU2 * sinLambda / sinSigma;
		cosSqAlpha = 1 - sinAlpha * sinAlpha;
		cos2SigmaM = cosSigma - 2 * sinU1 * sinU2 / cosSqAlpha;

		double C = e->dFlatness / 16 * cosSqAlpha * (4 + e->dFlatness * (4 - 3 * cosSqAlpha));
		lambdaP = lambda;
		lambda = L + (1 - C) * e->dFlatness * sinAlpha * (sigma + C * sinSigma * (cos2SigmaM + C * cosSigma * ( -1 + 2 * cos2SigmaM * cos2SigmaM)));
		eps = fabs(lambda - lambdaP);
	} 
	while (eps > 1e-12 && --iterLimit > 0);

	if (iterLimit == 0)
		return false; // formula failed to converge

	double uSq = cosSqAlpha * (e->dMajorAxis * e->dMajorAxis - e->dMinorAxis * e->dMinorAxis) / (e->dMinorAxis * e->dMinorAxis);
	double A = 1 + uSq / 16384 * (4096 + uSq * (-768 + uSq * (320 - 175 * uSq)));
	double B = uSq / 1024 * (256 + uSq * (-128 + uSq * (74 - 47 * uSq)));
	double deltaSigma = B * sinSigma * ( cos2SigmaM + B / 4 * ( cosSigma * (-1 + 2 * cos2SigmaM * cos2SigmaM ) - B / 6 * cos2SigmaM * ( -3 + 4 * sinSigma * sinSigma) * (-3 + 4 * cos2SigmaM * cos2SigmaM)));
	double _s = e->dMinorAxis * A * (sigma - deltaSigma);

	double _fwdAz = atan2(cosU2 * sinLambda, cosU1 * sinU2 - sinU1 * cosU2 * cosLambda); //y,x
	double _revAz = atan2(cosU1 * sinLambda, -sinU1 * cosU2 + cosU1 * sinU2 * cosLambda);

	if( _fwdAz < 0.0 )
		_fwdAz = _fwdAz + PI + PI;

	if( _revAz < 0.0 )
		_revAz = _revAz + PI + PI;

	m_dDistance =  _s;
	m_dFwdAz = _fwdAz;
	m_dRevAz = _revAz;

	return true;
}
