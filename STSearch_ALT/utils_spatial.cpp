//
// Project:    Spatially Combined Text Searches
// Filename:   utils_spatial.cpp
// Created on: 17.07.20.
//
// Developer(s):
//    I) Artur Titkov
//       ArturTitkov@icloud.com
//       https://github.com/TitkovArtur
//    II) unknown
//
// Description:
//
//
// Copyright © 2020 Artur Titkov. All rights reserved.
//


#include "def.h"
#include "relation.h"




/// SPATIAL 
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




///MBRs
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

//Probe for intersection of two MBRs
inline bool intersect(int dim, double *mbrR, double *mbrS)
{
    for (auto i = 0; i < dim * 2; i += 2)
    {
        if ((mbrR[i] > mbrS[i + 1]) || (mbrS[i] > mbrR[i + 1]))
            return false;
    }

    return true;
}
//Probe for intersction and stores the insersection into mbr
inline bool intersect(int dim, double *mbr, double *mbrR, double *mbrS)
{
    for (auto i = 0; i < dim * 2; i += 2)
    {
        mbr[i] = max(mbrR[i], mbrS[i]);
        mbr[i + 1] = min(mbrR[i + 1], mbrS[i + 1]);
        if (mbr[i] > mbr[i + 1])
            return false;
    }

    return true;
}



inline bool inside(int dim, double *p, double *mbr)
{
    for (auto i = 0; i < dim; i++)
    {
        if ((p[i] < mbr[2 * i]) || (p[i] > mbr[2 * i + 1]))
            return false;
    }

    return true;
}


///Distances
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline float distance2N(const Record &r, const Record &s){
    return sqrt( (r.locx - s.locx) * (r.locx - s.locx) + (r.locy - s.locy) * (r.locy - s.locy) );
}
inline float distance1N(const Record &r, const Record &s){
    return  (r.locx - s.locx) * (r.locx - s.locx) + (r.locy - s.locy) * (r.locy - s.locy) ;
}

//return dist_sqrt
inline float MinMinDist(double *mbrR, double *mbrS, double threshold){
    //judge spatial threshold
    double v1 = max(mbrR[0], mbrS[0]);
    double v2 = min(mbrR[1], mbrS[1]);
    double d2 = ((v1 < v2) ? ((v1 - v2 - 2* threshold) * (v1 - v2 - 2*threshold)) : 0);
    if (d2 < 0)
        return -1;

    v1 = max(mbrR[2], mbrS[2]);
    v2 = min(mbrR[3], mbrS[3]);
    d2 += ((v1 < v2) ? ((v1 - v2 - 2* threshold) * (v1 - v2 - 2* threshold)) : 0);
    if (d2 < 0)
        return -1;

    return d2;
}


inline bool qualify2(Record* r, Record* s, double dthreshold_sqr){
    
    auto dist = (r->locx - s->locx) * (r->locx - s->locx); //compute sqr distance
    
    if (dist > dthreshold_sqr)
        return false;
    if (dist + (r->locy - s->locy) * (r->locy - s->locy) > dthreshold_sqr)
        return false;
    
//    cout << "\t" << sqrt(dist) << endl;
    
    return true;
}

inline bool qualify(const Record &r, const Record &s, const double dthreshold_sqr)
{
    auto dist = (r.locx - s.locx) * (r.locx - s.locx); //compute sqr distance
    
    if (dist > dthreshold_sqr)
        return false;
    if (dist + (r.locy - s.locy) * (r.locy - s.locy) > dthreshold_sqr)
        return false;
    
//    cout << "\t" << sqrt(dist) << endl;
    
    return true;
}


inline bool qualify(const Record &r, const double locx, const double locy, const double dthreshold_sqr)
{
    auto dist = (r.locx - locx) * (r.locx - locx);
    
    if (dist > dthreshold_sqr)
        return false;
    
    if (dist + (r.locy - locy) * (r.locy - locy) > dthreshold_sqr)
        return false;
    
    //    cout << "\t" << sqrt(dist) << endl;
    
    return true;
}


inline bool qualify(const double rlocx, const double rlocy, const double slocx, const double slocy, const double dthreshold_sqr)
{
    auto dist = (rlocx - slocx) * (rlocx - slocx);
    
    if (dist > dthreshold_sqr)
        return false;
    
    if (dist + (rlocy - slocy) * (rlocy - slocy) > dthreshold_sqr)
        return false;
    
    //    cout << "\t" << sqrt(dist) << endl;
    
    return true;
}


inline bool qualify(double *mbrR, double *mbrS, double dthreshold_sqr)
{
    //judge spatial threshold
    double v1 = max(mbrR[0], mbrS[0]);
    double v2 = min(mbrR[1], mbrS[1]);
    double d2 = ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
    if (d2 > dthreshold_sqr)
        return false;
    
    v1 = max(mbrR[2], mbrS[2]);
    v2 = min(mbrR[3], mbrS[3]);
    d2 += ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
    if (d2 > dthreshold_sqr)
        return false;
    
    return true;
}



//SORT
inline bool QualifyMinMinDistSqr(double *mbrR, double *mbrS, double threshold_sqr)
{
    //judge spatial threshold
    double v1 = max(mbrR[0], mbrS[0]);
    double v2 = min(mbrR[1], mbrS[1]);
    double d2 = ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
    if (d2 > threshold_sqr)
        return false;

    v1 = max(mbrR[2], mbrS[2]);
    v2 = min(mbrR[3], mbrS[3]);
    d2 += ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
    if (d2 > threshold_sqr)
        return false;

    return true;
}
