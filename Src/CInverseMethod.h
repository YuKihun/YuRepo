/**
 * @brief 두 지점 -> 거리, 방위각 산출
 */

#pragma once

#include "VincentyParam.h" 	 

class CInverseMethod
{
public:
	CInverseMethod( );
	~CInverseMethod( );

private:
	double m_dDistance;
	double m_dFwdAz;
	double m_dRevAz;
	
public:
//	int VincentyInverse( sEllipsoid *e, double lat1, double lon1, double lat2, double lon2, double *s, double *fwdAz, double *revAz);
	int VincentyInverse( sEllipsoid *e, double lat1, double lon1, double lat2, double lon2 );

	//두 좌표간의 거리..
	double GetDistance( ) const { return m_dDistance; }

	//두 좌표에서 서로 바라보는 방위각..
	double GetFwdAz( ) const { return m_dFwdAz * R2D; }
	double GetRevAz( ) const { return m_dRevAz * R2D; }
};
