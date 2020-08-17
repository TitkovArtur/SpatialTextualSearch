//
//  ResultSet.hpp
//  STSearch_ALT
//
//  Created by Artur Titkov on 11.06.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.
//

#ifndef _ResultSet_H_
#define _ResultSet_H_

#include "def.h"
#include "relation.h"

typedef InvertedIndexVector InvertedIndex;
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
typedef float Coordinates;
typedef size_t RecordId;
typedef long termSetSize;
typedef int* termSet;
typedef float Score;
typedef pair<int, float> WeightedTerm;

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
class PairRecord{
public:
    RecordId LID;
    Record* lRec;
    RecordId RID;
    Record* rRec;
    Score score;
    
    
    PairRecord(Record* lRec, Record* rRec, RecordId l, RecordId r, Score s);
    ~PairRecord();

    bool operator < (const PairRecord& rhs);
    bool operator >= (const PairRecord& rhs);
    void print(char c);
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct CDJResult{
    vector<PairRecord> result;
    size_t numResult;
    bool getResult;

    Coordinates theta;
    Coordinates theta_sqr;

    termSet termSetL;
    termSetSize termSetSizeL;

    termSet termSetR;
    termSetSize termSetSizeR;


    CDJResult(Coordinates t, termSet tl, termSetSize sl, termSet tr, termSetSize sr);
    ~CDJResult();
    inline void insert();
    inline void insert(PairRecord rec);
    void clear();
    void print(char c);


};


//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct BkQResult{
    vector<PairRecord> result;
    size_t numResult;
    size_t numInsertions;
    int k;
    priority_queue<float> dist;

    Coordinates theta;
    Coordinates theta_sqr;
    Coordinates tmp_theta;
    Coordinates tmp_theta_sqr;

    termSet termSetL;
    termSetSize termSetSizeL;

    termSet termSetR;
    termSetSize termSetSizeR;


    BkQResult(int t, termSet tl, termSetSize sl, termSet tr, termSetSize sr);
    ~BkQResult();
    inline void insert(Coordinates dist);
    inline void insert(PairRecord rec);
    void clear();
    void print(char c);
};

//-----TkRQ
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct TkRQResult{
    vector<PairRecord> result;
    size_t numResult;
    size_t numInsertions;
    int k;
    Score scoreTheshold = 0.0;
    Score scoreThesholdtmp = 0.0;

    Coordinates theta;
    Coordinates theta_sqr;


    termSet termSetL;
    termSetSize termSetSizeL;
    Score lThreshold;
    WeightedTerm* lTermVector;
    termSet termSetR;
    termSetSize termSetSizeR;
    Score rThreshold;
    WeightedTerm* rTermVector;


    TkRQResult(int numK, float distThreshold, termSet tl, termSetSize sl, termSet tr, termSetSize sr, InvertedIndex* ifl, InvertedIndex* ifr);
    ~TkRQResult();
    inline void insert(Score s);
    inline void insert(PairRecord rec);
    void clear();
    void print(char c);
};














#endif /* ResultSet_h */
