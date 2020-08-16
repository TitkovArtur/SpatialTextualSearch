//This file is included in every header in GSTOOLS
//It is useful to have some global control definitions
//by GS

#ifndef _COMMON_H_
#define _COMMON_H_

#pragma once

//control macro definition
//----------------------------------------
//For r*-tree implemented by berchtold
#define _USE_RSTAR_BERCHTOLD_
#define ZAEHLER //open this macro to activate page_access
#define CHECK_MBR

//#define _SHOW_SPLIT_INFO_
//----------------------------------------

//Some Constant Definitions
static const int defaultBufferMax=4096;
static const int defaultPageSize=4096; //4096
static const float defaultSigmaF=0.05f;
static const double defaultSigmaD=0.05;

#define rangeMin 0
#define rangeMax 1

#endif
