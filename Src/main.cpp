#include "CInverseMethod.h"
#include "VincentyParam.h"

int main( )
{
	double dLatOfSource = 30.00; //시작 Lat
	double dLongOfSource = 120.00; //시작 Lat

	double dLatOfDest = 40.00; //목적지 Lat
	double dLongOfDest = 120.00; //목적지 Lat

	sEllipsoid ellipsoid; //타원체 정의
	ellipsoid.dMajorAxis = WGS84_MAJOR;
	ellipsoid.dMinorAxis = WGS84_MINOR;
	ellipsoid.dFlatness = WGS84_FLATNESS;

	CInverseMethod Inverse;

	for( int i = 0 ; i < 10 ; i++ )
	{
		Inverse.VincentyInverse( &ellipsoid, dLatOfSource, dLongOfSource, dLatOfDest, dLongOfDest );
		printf( "S(Lat : %f, Long : %f) -- ", dLatOfSource, dLongOfSource );
		printf( "D(Lat : %f, Long : %f) \n", dLatOfDest, dLongOfDest );
		printf( "Angle1 : %f, Angle2 : %f, Distance : %f m \n\n", Inverse.GetFwdAz( ), Inverse.GetRevAz( ), Inverse.GetDistance( ) );
		dLongOfDest = dLongOfDest + 1;
	}
	return 0;
}

