/**
 * @brief 한 지점, 거리, 방위각 -> 목표 지점의 좌표 산출 
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
