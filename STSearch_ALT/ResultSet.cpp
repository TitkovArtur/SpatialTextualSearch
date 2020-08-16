//
//  ResultSet.cpp
//  STSearch_ALT
//
//  Created by Artur Titkov on 11.06.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.
//


#include "ResultSet.h"
#include "inverted_index.h"

typedef InvertedIndexVector InvertedIndex;
class InvertedListEntry;



bool sortbysec2(const pair<int,float> &a, const pair<int,float> &b){
    return (a.second > b.second);
}



struct compResult {
    bool operator()(const PairRecord& a, const PairRecord& b) const
    {
        return a.score > b.score;
    }
};

//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
PairRecord::PairRecord(Record* lR, Record* rR, RecordId l, RecordId r, Score s){
    LID = l;
    lRec = lR;
    rRec = rR;
    RID = r;
    score = s;
}

bool PairRecord::operator < (const PairRecord& rhs){
    return this->score < rhs.score;
}

bool PairRecord::operator >= (const PairRecord& rhs){
    return !((this->score) < score);
}

void PairRecord::print(char c){
//    cout << c << " [(" << this->LID << ", "  << this->RID << "), ("<< this->lRec->id << " ," << this->rRec->id<< ") " << this->score << "]\n";
    cout << c << " [(" << this->LID << ", "  << this->RID << "), d = " << this->score << "]\n";

}

PairRecord::~PairRecord(){}




//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

BRQResult::BRQResult(Coordinates t, termSet tl, termSetSize sl, termSet tr, termSetSize sr){
//    getResult = result;
    numResult = 0;
    theta = t;
    theta_sqr = theta * theta;

    termSetL = tl;
    termSetSizeL = sl;
    termSetR = tr;
    termSetSizeR = sr;
}
BRQResult::~BRQResult(){}



inline void BRQResult::insert(){
    numResult++;
}

inline void BRQResult::insert(PairRecord rec){
    numResult++;
    result.push_back(rec);
}

void BRQResult::clear(){
    result.clear();
    numResult = 0;
}

void BRQResult::print(char c){
    cout << "RESULT " << c << " pairs " << result.size() <<"\n";
    if(result.size()==0){
        cout << "EMPTY " <<"\n";
    }else{
        for(int i = 0; i < result.size(); i++){
            result.at(i).print(c);
        }
        cout << "\n";
    }
}







//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
BkQResult::BkQResult(int num, termSet tl, termSetSize sl, termSet tr, termSetSize sr){
    numInsertions= 0;
    k = num;
    theta = numeric_limits<int>::max();
    tmp_theta = numeric_limits<int>::max();
    theta_sqr = theta * theta;
    termSetL = tl;
    termSetSizeL = sl;
    termSetR = tr;
    termSetSizeR = sr;
}



inline void BkQResult::insert(Coordinates dist){
    numInsertions++;
    
    if(numInsertions > k){
        theta = min(dist, theta);
        theta_sqr = theta * theta;
    }else if (numInsertions == k){
        theta = min(dist, tmp_theta);
        theta_sqr = theta * theta;
    }else{
        tmp_theta = min(dist, tmp_theta);
    }
}



inline void BkQResult::insert(PairRecord rec){
    numInsertions++;
    auto it = result.begin();
    for (; it != result.end(); it++){
        if( !(*it < rec) ){
            break;
        }
    }
    result.insert(it, rec);
    if(result.size() > k){
        result.erase(result.end()-1);
        theta = (result.end()-1)->score;
        theta_sqr = theta * theta;
    }
}

void BkQResult::clear(){
    result.clear();
    numInsertions = 0;
    numResult = 0;
    theta = numeric_limits<int>::max();
    tmp_theta = numeric_limits<int>::max();
    theta_sqr = theta * theta;
    while (!dist.empty()) {
        dist.pop();
    }
}


void BkQResult::print(char c){
    cout << "RESULT " << c << " pairs " << result.size() <<"\n";
    if(result.size()==0){
        cout << "EMPTY " <<"\n";
    }else{
        for(int i = 0; i < result.size(); i++){
            result.at(i).print(c);
        }
        cout << "\n";
    }
}









//-----TkRQ
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

TkRQResult::TkRQResult(int numK, float distThreshold, termSet tl, termSetSize sl, termSet tr, termSetSize sr, InvertedIndex* ifL, InvertedIndex* ifR){
    numResult = 0;
    numInsertions = 0;
    k = numK;
    theta = distThreshold;
    theta_sqr = distThreshold * distThreshold;
    scoreTheshold = 0;
    
    lTermVector = new WeightedTerm[6];
    rTermVector = new WeightedTerm[6];
    termSetL = tl;
    termSetSizeL = sl;
    lThreshold = 0;
    termSetR = tr;
    termSetSizeR = sr;
    rThreshold = 0;
    for(int i = 0; i < sl; i++){ //Computer weight L
        int termID = tl[i];
        float weigthedTerm = ifL->lists.at(termID).size();
        weigthedTerm = log2(weigthedTerm);
        
        lTermVector[i] = make_pair(termID, weigthedTerm);
    }
    for(int i = 0; i < sr; i++){ //Compute weight R
        int termID = tr[i];
        float weigthedTerm = ifR->lists.at(termID).size();
        weigthedTerm = log2(weigthedTerm);
        
        rTermVector[i] = make_pair(termID, weigthedTerm);
    }
    
    float norLenght = 0.0; //LEFT
    for(int i = 0; i < sl; i++){ //get vecor length
        norLenght += lTermVector[i].second * lTermVector[i].second;
    }
    norLenght = sqrt(norLenght);
    for(int i = 0; i < sl; i++){ //get vecor length
        lTermVector[i].second /= norLenght;
    }
    
    norLenght = 0.0; //RIGHT
    for(int i = 0; i < sr; i++){ //get vecor length
        norLenght += rTermVector[i].second * rTermVector[i].second;
    }
    norLenght = sqrt(norLenght);
    for(int i = 0; i < sr; i++){ //get vecor length
        rTermVector[i].second /= norLenght;
    }
    
    
    //SORT
    sort( lTermVector, &lTermVector[sl], sortbysec2 );
    sort( rTermVector, &rTermVector[sr], sortbysec2 );
    
//    for(int i = 0; i < sr; i++ ){
//        cout << "w: "<< lTermVector[i].second;
//    }
//    for(int i = 0; i < sr; i++ ){
//        cout << "w: "<< rTermVector[i].second;
//    }
    
    
    
    

};
//inline void BkQResult::insert(Coordinates dist){
//    numInsertions++;
//
//    if(numInsertions > k){
//        theta = min(dist, theta);
//        theta_sqr = theta * theta;
//    }else if (numInsertions == k){
//        theta = min(dist, tmp_theta);
//        theta_sqr = theta * theta;
//    }else{
//        tmp_theta = min(dist, tmp_theta);
//    }
//}
//
//
//
//inline void BkQResult::insert(PairRecord rec){
//    numInsertions++;
//    auto it = result.begin();
//    for (; it != result.end(); it++){
//        if( !(*it < rec) ){
//            break;
//        }
//    }
//    result.insert(it, rec);
//    if(result.size() > k){
//        result.erase(result.end()-1);
//        theta = (result.end()-1)->score;
//        theta_sqr = theta * theta;
//    }
//}

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
    auto start = result.begin();
    auto end   = result.end();
    result.push_back(rec);
    
    sort(start,end, compResult());
    result.erase(result.end()-1);
    scoreThesholdtmp = (result.end()-1)->score;
};

void TkRQResult::clear(){
    result.clear();
    lThreshold = 0;
    rThreshold = 0;
    numInsertions = 0;
    numResult = 0;
    scoreTheshold = 0;
    
    
    
};
void TkRQResult::print(char c){
    cout << "RESULT " << c << " pairs: " << result.size() <<"\n";
    if(result.size()==0){
        cout << "EMPTY " <<"\n";
    }else{
        for(int i = 0; i < result.size(); i++){
            result.at(i).print(c);
        }
        cout << "\n";
    }
};



















