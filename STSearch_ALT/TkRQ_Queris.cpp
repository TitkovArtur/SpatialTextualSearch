//
//  TkRQ_Queris.cpp
//  TkRQ
//
//  Created by Artur Titkov on 16.07.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.
//

#include "def.h"
#include "inverted_index.h"
#include "ResultSet.h"
#include "relation.h"
//#include "utils_sim.cpp"
#include "utils_spatial.cpp"








Relation* textSearch(int numkeywords, WeightedTerm* terms, Relation* rel, InvertedIndex* iix);





struct compResult2 {
    bool operator()(const PairRecord& a, const PairRecord& b) const
    {
        return a.score > b.score;
    }
};







inline void TkRQResult::insert(Score s){
    numInsertions++;
    
    if(numInsertions > k){
        scoreTheshold = s;
    }else if (numInsertions == k){
        scoreTheshold = max(s, scoreThesholdtmp);
        theta_sqr = theta * theta;
    }else{
        scoreThesholdtmp = max(s, scoreThesholdtmp);
    }
};

inline void TkRQResult::insert(PairRecord rec){
    numInsertions++;
    result.push_back(rec);
    
    auto start = result.begin();
    auto end   = result.end();
    sort(start,end, compResult2());
    if(numInsertions > k){
        result.erase(result.end()-1);
        scoreTheshold = rec.score;
//    scoreThesholdtmp = (result.end()-1)->score;
    }
};







////////////////////////////// BASELINE x0
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void NestedLoopsDistanceJoin(TkRQResult* q, Relation* R, Relation* S){
    size_t result = 0;
    auto dthreshold_sqr = q->theta_sqr;
    Record* x;
    Record* y;
    WeightedRecord* wl;
    WeightedRecord* wr;


    for(int i = 0; i < R->numRecords; i++){
        wl = &R->weightedRecords[i];
        x  = wl->second;
        for(int j = 0; j < S->numRecords; j++){
            wr = &S->weightedRecords[j];
            y = wr->second;
//            cout << "x1 " <<x->locx << endl;
//            cout << "x2 " <<x->locy << endl;
//            cout << "y1 " <<y->locx << endl;
//            cout << "y2 " <<y->locy << endl;
//            cout << "s " << wl->first << endl;
//            cout << "s " << wr->first << endl;
            if(qualify(x->locx, x->locy, y->locx, y->locy, q->theta_sqr)){
//                cout << " SCORE " << (float)wl->first *  (float)wr->first;
                Score s = wl->first * wr->first;
//                cout << s;
                if(s > q->scoreTheshold ){
#ifdef RetrieveResults
            q->insert(PairRecord(x, y, x->id, y->id, s));
#else
            q->insert(s);
#endif
                }
            }
        }
    }
}



////////////////////////////// Setup 1 x1
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


inline void internalLoop(BRQResult* q, const SubRelationIterator rec, const SubRelationIterator firstFS, const SubRelationIterator lastFS){
    //unsigned long long result = 0;
    auto pivot = firstFS;

    // Sweep on X.
    while ((pivot < lastFS) && ((*rec)->locx >= (*pivot)->locx - q->theta)){
        // Verification on Y.
        if (((*rec)->locy <= (*pivot)->locy + q->theta) && ((*rec)->locy >= (*pivot)->locy - q->theta) && (qualify(**rec, **pivot, q->theta_sqr))){
            // Found a result pair.
            
#ifdef RetrieveResults
            q->insert(PairRecord(*rec, *pivot, (*rec)->id, (*pivot)->id, distance2N(**rec, **pivot)) );
#else
            q->insert();
#endif
            
        }
        pivot++;
    }
}


//void PlaneSweepDistanceJoin(BRQResult* q, Relation* R,Relation* S){
//    auto r = R->weightedRecords;
//    auto s = S->weightedRecords;
//    auto lastR = &(R->weightedRecords[R->numRecords]);
//    auto lastS = &(S->weightedRecords[S->numRecords]);
//
//    //Call Internal Loop
//    //With OR Without Result
//    while ((r < lastR) && (s < lastS)){
//        if ((*r)->locx < (*s)->locx){
////        if ((*r) < (*s)){
//            internalLoop(q, r, s, lastS);
//            r++;
//        }else{
//            internalLoop(q, s, r, lastR);
//            s++;
//        }
//    }
//}
