/**
 * @brief �� ���� -> �Ÿ�, ������ ����
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

	//�� ��ǥ���� �Ÿ�
	double GetDistance( ){ return m_dDistance; }

	//�� ��ǥ���� ���� �ٶ󺸴� ������
	double GetFwdAz( ){ return m_dFwdAz * R2D; }
	double GetRevAz( ){ return m_dRevAz * R2D; }
};
