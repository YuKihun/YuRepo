/**
 * @brief �� ����, �Ÿ�, ������ -> ��ǥ ������ ��ǥ ���� 
 */

#pragma once

#include "VincentyParam.h"

class CDirectMethod
{
public:
	CDirectMethod();
	~CDirectMethod();

	int VincentyDirect( sEllipsoid* e, double lat1, double lon1, double alpha1, double s, double *lat2, double *lon2, double * revAz);
};
