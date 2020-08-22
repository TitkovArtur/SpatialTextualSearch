#ifndef _DEF_H_
#define _DEF_H_


// Includes
#include "RStarTree.h"
#include "Common.h"
#include <cstdio>
#include <cstdlib>
#include <cstring>
#include <cmath>
#include <algorithm>
#include <sstream>
#include <iostream>
#include <fstream>
#include <limits>
#include <vector>
#include <list>
#include <map>
#include <queue>
#include <deque>
#include <tuple>
#include <chrono>
#include <unordered_map>
#include <unordered_set>
using namespace std;


typedef float Coordinates;
typedef size_t RecordId;
typedef long termSetSize;
typedef int* termSet;
typedef float Score;

typedef RSTNonLeafNode<TreeDataP,double> NonLeaf;
typedef RSTLeafNode<TreeDataP,double> Leaf;
typedef RSTNode<TreeDataP,double> Node;
typedef RStarTree<TreeDataP, double> RTree;
typedef TreeDataP<double> Data;
//typedef tuple<Record*, Node*> RecNodePair;


//#define poly
#define RetrieveResults



#define SPATIO_TEXTUAL
// Defines
#define hash                                 unordered_map
#define hashset                              unordered_set

#define EPS             1e-8
#define NUM_DIMENSIONS  2
#define PAGE_SIZE       4096

// Specify intersection strategy
#define INTERSECTION_VECTOR
//#define INTERSECTION_VECTOR_BINSEARCH
//#define INTERSECTION_VECTOR_SMART
//#define INTERSECTION_MAP
//#define INTERSECTION_HASHMAP


#define ASCENDING_KEYWORDS_ORDER

// Imports
class Record;
class Relation;
class InvertedListEntry;
class InvertedIndexVector;
class InvertedIndexVectorBinSearch;
class InvertedIndexVectorSmart;
class InvertedIndexMap;
class InvertedIndexHashMap;


#ifdef INTERSECTION_VECTOR
typedef InvertedIndexVector InvertedIndex;
#elif defined INTERSECTION_VECTOR_BINSEARCH
typedef InvertedIndexVectorBinSearch InvertedIndex;
#elif defined INTERSECTION_VECTOR_SMART
typedef InvertedIndexVectorSmart InvertedIndex;
#elif defined INTERSECTION_MAP
typedef InvertedIndexMap InvertedIndex;
#elif defined INTERSECTION_HASHMAP
typedef InvertedIndexHashMap InvertedIndex;
#endif






//bool inside(int dim, Coordinates *p, Coordinates *mbr)
//{
//    for (auto i = 0; i < dim; i++)
//    {
//        if ((p[i] < mbr[2 * i]) || (p[i] > mbr[2 * i + 1]))
//            return false;
//    }
//
//    return true;
//}
//
//bool intersect(int dim, Coordinates *mbr, Coordinates *mbrR, Coordinates *mbrS)
//{
//    for (auto i = 0; i < dim * 2; i += 2)
//    {
//        mbr[i] = max(mbrR[i], mbrS[i]);
//        mbr[i + 1] = min(mbrR[i + 1], mbrS[i + 1]);
//        if (mbr[i] > mbr[i + 1])
//            return false;
//    }
//
//    return true;
//}
//
//bool intersect(int dim, Coordinates *mbrR, Coordinates *mbrS)
//{
//    for (auto i = 0; i < dim * 2; i += 2)
//    {
//        if ((mbrR[i] > mbrS[i + 1]) || (mbrS[i] > mbrR[i + 1]))
//            return false;
//    }
//
//    return true;
//}
//
//
//
//struct mbr_sort_dim_x_less
//{
//    RSTNonLeafNode<TreeDataP, Coordinates> *nonleaf;
//
//    mbr_sort_dim_x_less(RSTNonLeafNode<TreeDataP, Coordinates> *node) : nonleaf(node) {}
//    bool operator()(const int id1, const int id2) const
//    {
//        return (nonleaf->entry_mbr[id1][0]<nonleaf->entry_mbr[id2][0]);
//    }
//};
//
//
//struct p_sort_dim_x_less
//{
//    RSTLeafNode<TreeDataP, Coordinates> *leaf;
//
//    p_sort_dim_x_less(RSTLeafNode<TreeDataP, Coordinates> *node) : leaf(node) {}
//    bool operator()(const int id1, const int id2) const
//    {
//        return (leaf->data[id1]->data[0]<leaf->data[id2]->data[0]);
//    }
//};
//
//
//struct p_sort_dim_y_less
//{
//    RSTLeafNode<TreeDataP, Coordinates> *leaf;
//
//    p_sort_dim_y_less(RSTLeafNode<TreeDataP, Coordinates> *node) : leaf(node) {}
//    bool operator()(const int id1, const int id2) const
//    {
//        return (leaf->data[id1]->data[1]<leaf->data[id2]->data[1]);
//    }
//};






//Shuyao
//extern Relation *R;
//extern Relation *S;
//extern bool *invalidRecR;
//extern bool *invalidRecS;
//extern bool *invalidR;
//extern bool *invalidS;
//extern InvertedIndex *iidxR;
//extern InvertedIndex *iidxS;
//extern InvertedIndex *IIR;
//extern InvertedIndex *IIS;

//tools.cpp
class Timer{
private:
    using Clock = std::chrono::high_resolution_clock;
    Clock::time_point start_time, stop_time;
    
public:
    Timer(){
        start();
    }
    
    void start(){
        start_time = Clock::now();
    }
    
    
    double getElapsedTimeInSeconds(){
        return std::chrono::duration<double>(stop_time - start_time).count();
    }
    
    
    double stop(){
        stop_time = Clock::now();
        double tmp = getElapsedTimeInSeconds();
        start();
        return tmp;
    }
};



void Test();
void ChooseKeywords();
bool PruneByKw(InvertedIndex *iidx, bool *invalid, int page_id, int numKw, int *kw);
bool QualifySpatial(Record *r, Record *s, double dthreshold_sqr);
void SpatialDistanceFilter(RSTNonLeafNode<TreeDataP, double> *nonleafR, RSTNonLeafNode<TreeDataP, double>* nonleafS, deque<pair<RSTNode<TreeDataP, double> *, RSTNode<TreeDataP, double> *> > &q, vector<int> &childR, vector<int> &childS, double sthreshold, double sthreshold_sqr);
unsigned int SpatialDistanceFilter(Relation *R, RSTLeafNode<TreeDataP, double> *leafR, Relation *S, RSTLeafNode<TreeDataP, double>* leafS, vector<int> &childR, vector<int> &childS, double sthreshold, double sthreshold_sqr);
unsigned int SpatialDistanceFilter(Relation *R, RSTLeafNode<TreeDataP, double> *leafR, Relation *S, RSTLeafNode<TreeDataP, double>* leafS, vector<int> &childR, vector<int> &childS, double threshold, double threshold_sqr, int RNumKeywords, int *RKeywords, int SNumKeywords, int *SKeywords);
bool ContainKeywords(Record *r, int numKeywords, int *keywords);
void CreateIRTreeNode(RSTNode<TreeDataP, double> *node, Relation *R, InvertedIndex *iidx);
void IRTreeFilter(RSTNonLeafNode<TreeDataP, double> *nonleafR, RSTNonLeafNode<TreeDataP, double>* nonleafS, deque<pair<RSTNode<TreeDataP, double> *, RSTNode<TreeDataP, double> *> > &q, vector<int> &childR, vector<int> &childS, double threshold, double threshold_sqr, int numKwR, int *kwR, int numKwS, int *kwS);
unsigned int IRTreeFilter(Relation *R, RSTLeafNode<TreeDataP, double> *leafR, Relation *S, RSTLeafNode<TreeDataP, double>* leafS, vector<int> &childR, vector<int> &childS, double threshold, double threshold_sqr, int numKwR, int *kwR, int numKwS, int *kwS);
RStarTree<TreeDataP, double>* BuildRTree(const char *datafile, Relation *RL, bool build_new);



//queries.cpp
unsigned int SetContainmentQuery(int numKeywords, int *keywords, InvertedIndex *iidx);
unsigned int TextualFirst(int RNumKeywords, int *RKeywords, int SNumKeywords, int *SKeywords, double threshold);
unsigned int TextualFirst(Relation *R, Relation *S, int RNumKeywords, int *RKeywords, int SNumKeywords, int *SKeywords, double threshold);
unsigned int SpatialFirst(Relation *R, RStarTree<TreeDataP, double> *rtR, Relation *S, RStarTree<TreeDataP, double> *rtS, double sthreshold, int RNumKeywords, int *RKeywords, int SNumKeywords, int *SKeywords);
unsigned int SpatialJoin(Relation *R, RStarTree<TreeDataP, double> *rtR, Relation *S, RStarTree<TreeDataP, double> *rtS, double sthreshold);
unsigned int SpatialJoinNaive(Relation *R, Relation *S, double threshold);


//BRQ_Queries.cpp

Relation* SetContainmentQueryWithResult(Relation& result, int numKeywords, int *keywords, InvertedIndex *iidx);


#endif
