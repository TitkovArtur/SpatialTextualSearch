//
//  CkJ_Queries.cpp
//  STSearch_ALT
//
//  Created by Artur Titkov on 17.06.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.
//

#include "def.h"
#include "inverted_index.h"
#include "ResultSet.h"
#include "relation.h"
#include "utils_boolean.cpp"
#include "utils_spatial.cpp"


#define RetrieveResults
int counter = 0;

struct compareScoredNodes;
struct compareScoredPair;

typedef tuple<Record*, Node*> RecNodePair;
typedef tuple<Node*, Node*> NodesPair;
typedef pair<Node*, float> ScoredNode;
typedef tuple<Node*, Node*, float> ScoredPair;
typedef priority_queue<ScoredNode, vector<ScoredNode>, compareScoredNodes> PqNode;
typedef priority_queue<ScoredPair, vector<ScoredPair>, compareScoredPair> PqNodes;




inline bool textSearchOnLeaf(int* keywords, int numKeywords, InvertedIndex* iidx, Node* leaf);
inline bool textProbeOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);
inline bool textSearchOnLeaf(int* keywords, int numKeywords, InvertedIndex* iidx, Node* leaf);
inline bool textSearchOnNodeBinSearch(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);
inline bool textVerificationOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);

inline bool textSearchOnLeaf(int* keywords, int numKeywords, InvertedIndex* iidx, Node* leaf);
inline bool textProbeOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);
inline bool textSearchOnLeaf(int* keywords, int numKeywords, InvertedIndex* iidx, Node* leaf);
inline bool textSearchOnNodeBinSearch(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);
inline bool textVerificationOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);




//////////////////////////////TOOLS Closest Pairs Tools
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
struct compareScoredNodes {
    bool operator()(const ScoredNode& x,const ScoredNode& y) const{
        return x.second > y.second;
    }
};
struct compareScoredPair {
    bool operator()(const ScoredPair& x,const ScoredPair& y) const{
        return get<2>(x) > get<2>(y);
    }
};

//calculate overlap of MBR
template <typename T>
T overlapRectangles(int dim, T *r1, T *r2)
{
    T sum=(T)1;
    for(int i=0;i<2*dim;i+=2)
    {
        T lower=max<T>(r1[i], r2[i]);
        T upper=min<T>(r1[i+1], r2[i+1]);
        if(lower>=upper) return (T)0; //do not overlap
        else sum*=upper-lower;
    }

    return sum;
}


//////////////////////////////TOOLS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

struct mbr_sort_dim_x_less
{
    RSTNonLeafNode<TreeDataP, double> *nonleaf;

    mbr_sort_dim_x_less(RSTNonLeafNode<TreeDataP, double> *node) : nonleaf(node) {}
    bool operator()(const int id1, const int id2) const
    {
        return (nonleaf->entry_mbr[id1][0]<nonleaf->entry_mbr[id2][0]);
    }
};


struct p_sort_dim_x_less
{
    RSTLeafNode<TreeDataP, double> *leaf;

    p_sort_dim_x_less(RSTLeafNode<TreeDataP, double> *node) : leaf(node) {}
    bool operator()(const int id1, const int id2) const
    {
        return (leaf->data[id1]->data[0]<leaf->data[id2]->data[0]);
    }
};


struct p_sort_dim_y_less
{
    RSTLeafNode<TreeDataP, double> *leaf;

    p_sort_dim_y_less(RSTLeafNode<TreeDataP, double> *node) : leaf(node) {}
    bool operator()(const int id1, const int id2) const
    {
        return (leaf->data[id1]->data[1]<leaf->data[id2]->data[1]);
    }
};

inline float QualifyMinMinDistSqr(float *mbrR, float *mbrS, float threshold_sqr)
{
    //judge spatial threshold
    float v1 = max(mbrR[0], mbrS[0]);
    float v2 = min(mbrR[1], mbrS[1]);
    float d2 = ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
    if (d2 > threshold_sqr)
        return -1;

    v1 = max(mbrR[2], mbrS[2]);
    v2 = min(mbrR[3], mbrS[3]);
    d2 += ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
    if (d2 > threshold_sqr)
        return -1;

    return 0;
}



bool qualify(const Record &r, const Record &s, const double dthreshold_sqr);
bool qualify(const Record &r, const double locx, const double locy, const double dthreshold_sqr);
bool qualify(const double rlocx, const double rlocy, const double slocx, const double slocy, const double dthreshold_sqr);
bool qualify(double *mbrR, double *mbrS, double dthreshold_sqr);
RStarTree<TreeDataP, double>* bulkload(int dim, int page_len, TreeDataP<double> **data, int data_num);
bool qualify(const Record &r, const Record &s, const double dthreshold_sqr);
bool qualify(const Record &r, const double locx, const double locy, const double dthreshold_sqr);
bool qualify(double *mbrR, double *mbrS, double dthreshold_sqr);
bool intersect(int dim, double *mbr, double *mbrR, double *mbrS);

//OLDINSERTION
//inline void CkJResult::insert(Coordinates d){
//    numInsertions++;
//
//
//    dist.push(d);
//
//    if(numInsertions > k){
////        dist.push(d);
//        float tmp = dist.top();
//        dist.pop();
//        theta = tmp;
//        theta_sqr = tmp * tmp;
//    }else if(numInsertions == k){
////        dist.push(d);
//        theta = dist.top();
//        theta_sqr = theta * theta;
//    }
//
//
//
////    if(numInsertions > k){
////        theta = min(dist, theta);
////        theta_sqr = theta * theta;
////    }else if (numInsertions == k){
////        theta = min(dist, tmp_theta);
////        theta_sqr = theta * theta;
////    }else{
////        tmp_theta = min(dist, tmp_theta);
////    }
//}



inline void CkJResult::insert(Coordinates d){
    numInsertions++;
    dist.push(d);
    
    if(numInsertions > k){
//        dist.push(d);
        dist.pop();
        float tmp = dist.top();
        theta = tmp;
        theta_sqr = tmp * tmp;
    }else if(numInsertions == k){
//        dist.push(d);
        theta = dist.top();
        theta_sqr = theta * theta;
    }

}

inline void CkJResult::insert(PairRecord rec){
        numInsertions++;
        result.push(rec);
        
        if(numInsertions > k){
    //        dist.push(d);
            result.pop();
            PairRecord tmprec =  result.top();
            float tmp = tmprec.score;
            theta = tmp;
            theta_sqr = tmp * tmp;
        }else if(numInsertions == k){
    //        dist.push(d);
            PairRecord tmprec = result.top();
            theta = tmprec.score;
            theta_sqr = theta * theta;
        }

    
    
    
    
    
    
    
    
    
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
}


//OLD
//inline void CkJResult::insert(PairRecord rec){
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

//class CompareKeywordsByFrequencyR
//{
//public:
//    bool operator() (const int& klhs, const int& krhs) const
//    {
//        return (iidxR->lists[klhs].size() < iidxR->lists[krhs].size());
//    }
//};
//
//Relation* SetContainmentQueryWithResult(Relation& T, int numKeywords, int *keywords, InvertedIndex *iidx){
//
//    std::vector<int>* result = new std::vector<int>;
//    std::vector<InvertedListEntry*> lists;
//    std::vector<int> actualPosition;
//    Relation* output = new Relation();
//    int listLengths[numKeywords];
//    int tmp[numKeywords];
//
//
//    //create tmp keyword set
//    for(int i = 0; i < numKeywords; i++){
//        tmp[i] = keywords[i];
//    }
//
//    //get postingList lengths for keywords
//    for(int i = 0; i < numKeywords; i++){
//        listLengths[i] = iidx->lists.at(tmp[i]).size();
//    }
//
//    //sort termSet regarding listLengths (bubble sort) // already sorted
////    for(int i = 0; i < terms; i++){
////        for(int j = 0; j < terms-1-i; j++){
////            if(listLengths[j] > listLengths[j+1]){
////                int t_tmp = tmp[j];
////                int l_tmp = listLengths[j];
////
////                tmp[j]     = tmp[j+1];
////                listLengths[j] = listLengths[j+1];
////                tmp[j+1]     = t_tmp;
////                listLengths[j+1] = l_tmp;
////            }
////        }
////    }
//
//    // get postingslists set positon vector
//    for(int i = 0; i < numKeywords ; i++){
//        lists.push_back( (iidx->lists.at(tmp[i]).data()) );
//        actualPosition.push_back(0);
//    }
//
//
//    // for each doc in first postingLists
//    for(int i = 0; i < listLengths[0]; i++){
//        int document = lists.at(0)[i].rec->id;
//
//        // for each postingsList from 1 ... numkeywords-1
//        for(int j = 1; j < numKeywords ; j++){
//            //for each doc in postinglist
//            for(int l = actualPosition[j]; l < listLengths[j]; l++){
//                if( document == lists.at(j)[l].rec->id ){
//                    actualPosition[j] = l;
//                    break;
//                } //lists.at(0)->rec->id
//            }
//        }
//        for(int j = 1; j < numKeywords ; j++){
//            if(document != lists.at(j)[actualPosition[j]].rec->id){
//                break;
//            }
//            if(j == numKeywords -1){
//                result->push_back(document);
//            }
//        }
//    }
//
//    vector<Record*>* m = new vector<Record*>;
//
//    for(int i = 0; i < result->size(); i++ ){
//        //(T.recs[(result->at(i))]).Print('s');
//        m->push_back(  &(T.recs[(result->at(i))]) );
//        output->numRecords++;
//    }
//    output->subTable = (m->data());
//    output->numRecords = m->size();
//    return output;
//};
//





////////////////////////////// SETUPS
////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////




////////////////////////////// BASELINE
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////x
void NestedLoopsDistanceJoin(CkJResult* q, Relation* R, Relation* S){
    size_t result = 0;
    auto dthreshold_sqr = q->theta_sqr;
    Record x,y;
        
        for(int r = 0; r < R->numRecords; r++){
            x = *R->beginSubRelation()[r];
            for(int s = 0; s < S->numRecords; s++){
                y = *S->beginSubRelation()[s];
                if (qualify(x, y, q->theta_sqr)){
#ifdef RetrieveResults
            q->insert(PairRecord(&x, &y, x.id, y.id, distance2N(x, y)));
#else
            q->insert(distance2N(x, y));
#endif
                }
            }
        }

}


void NestedLoopsDistanceJoinALT(CkJResult* q, Relation* R, Relation* S){
    size_t result = 0;
    auto dthreshold_sqr = q->theta_sqr;
    Record x,y;
        
        for(int r = 0; r < R->numRecords; r++){
            x = *R->beginSubRelation()[r];
            for(int s = 0; s < S->numRecords; s++){
                y = *S->beginSubRelation()[s];
                if (qualify(x, y, q->theta_sqr)){
                    q->insert(distance2N(x, y));
                }
            }
        }

}









////////////////////////////// SETUP 1 X1
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
inline void internalLoop(CkJResult* q, const SubRelationIterator rec, const SubRelationIterator firstFS, const SubRelationIterator lastFS){
    //unsigned long long result = 0;
    auto pivot = firstFS;

    // Sweep on X.
    while ((pivot < lastFS) && ((*rec)->locx >= (*pivot)->locx - q->theta)){

        // Verification on Y.                                                                               //TODO QUALIFY merge with insert
        if (((*rec)->locy <= (*pivot)->locy+ q->theta) && ((*rec)->locy >= (*pivot)->locy - q->theta) && (qualify(**rec, **pivot, q->theta_sqr))){
            // Found a result pair.
//        if(distance2N(**rec, **pivot) < q->theta ){
#ifdef RetrieveResults
            q->insert(PairRecord(*rec, *pivot, (*rec)->id, (*pivot)->id, distance2N(**rec, **pivot)) );
#else
            q->insert(distance2N(**rec, **pivot));
#endif
            
        }
        pivot++;
    }
}


void PlaneSweepDistanceJoin(CkJResult* q, Relation* R,Relation* S){
    auto r = R->beginSubRelation();
    auto s = S->beginSubRelation();
    auto lastR = R->endSubRelation();
    auto lastS = S->endSubRelation();

    //Call Internal Loop
    while ((r < lastR) && (s < lastS)){
        if ((*r)->locx < (*s)->locx){
            internalLoop(q, r, s, lastS);
            r++;
        }else{
            internalLoop(q, s, r, lastR);
            s++;
        }
    }
}





inline float QualifyMinMinDist(float *mbrR, float *mbrS, float threshold_sqr){
     //judge spatial threshold
     float v1 = max(mbrR[0], mbrS[0]);
     float v2 = min(mbrR[1], mbrS[1]);
     float d2 = ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
     if (d2 > threshold_sqr)
         return -1;

     v1 = max(mbrR[2], mbrS[2]);
     v2 = min(mbrR[3], mbrS[3]);
     d2 += ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
     if (d2 > threshold_sqr)
         return -1;

     return 0;
 }




////////////////////////////// SETUP 2___X2
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void ContainmentkNNQueryOnRTree(CkJResult& q, Record& l, Relation& R, const RStarTree<TreeDataP, double> *rtR){
    size_t result = 0;
    double qrange[4];
    double mbr[4], inter[4];
    float score;
    PqNode Q;
    counter = 0;
    
    qrange[0] = l.locx - q.theta;
    qrange[1] = l.locx + q.theta;
    qrange[2] = l.locy - q.theta;
    qrange[3] = l.locy + q.theta;
    
        
    
//    score = overlap(2, mbr, qrange);
//    if(score != 0)
//        Q.push(ScoredNode(rtR->root, overlap(2, mbr, qrange)));
    
    Q.push(make_pair(rtR->root, 1));
    stack<RSTNode<TreeDataP, double> *> S;
    S.push(rtR->root);
    

    
    while (!Q.empty()){
        if(q.numResult == q.k &&  Q.top().second > q.theta){
//        if( Q.top().second > q.theta){
            Q.pop();
            continue;
        }
        Node* cnode = Q.top().first;
//        cout << Q.top().second << "\n";
        Q.pop();

        if(cnode->is_leaf_node())
        {
            Leaf* leaf = static_cast<RSTLeafNode<TreeDataP, double> *>(cnode);

            leaf->get_mbr(mbr);
            mbr[0] -= q.theta;
            mbr[1] += q.theta;
            mbr[2] -= q.theta;
            mbr[3] += q.theta;

            if (!intersect(NUM_DIMENSIONS, inter, qrange, mbr))
                continue; //no intersection we return
            if (!(qualify(qrange, mbr, q.theta_sqr)))
                continue; //judge false hit
            
//            if(q.numResult >= q.k && )
            
            
            
            
            TreeDataP<double> **entry = leaf->data;
            vector<RecordId> child;
            for (auto i = 0; i < leaf->entry_num; i++)
                if ( (entry[i]->textScore != -1))
                    if(inside(NUM_DIMENSIONS, entry[i]->data, inter))
                            child.push_back(i);
            
            sort(child.begin(), child.end(), p_sort_dim_x_less(leaf));

            for (vector<RecordId>::iterator iter = child.begin(); iter != child.end(); ++iter){
                RecordId i = *iter;
                double locx = entry[i]->data[0];
                double locy = entry[i]->data[1];
                float* curScore = &(entry[i]->textScore);

                if (l.locx < (locx-q.theta))
                    continue;
                else if (l.locx > (locx+q.theta))
                    break;

                if ((l.locy+q.theta >= locy) && (l.locy - q.theta <= locy) && (qualify(l, locx, locy, q.theta_sqr)))

                    if( (entry[i]->textScore == 1) || textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[entry[i]->id], curScore) ){
#ifdef RetrieveResults
                        //TODO UPDATE PQ
                        q.insert(PairRecord(&l, R[entry[i]->id], l.id, entry[i]->id, distance2N(l, *R[entry[i]->id])));
#else
                        q.insert(distance2N(l, *R[entry[i]->id]));
#endif
                    }
                 }
        }else{
            NonLeaf* nonleaf = static_cast<RSTNonLeafNode<TreeDataP, double> *>(cnode);
            double mbr[4], inter[4];

            nonleaf->get_mbr(mbr);
            mbr[0] -= q.theta;
            mbr[1] += q.theta;
            mbr[2] -= q.theta;
            mbr[3] += q.theta;

            if (!intersect(NUM_DIMENSIONS, inter, qrange, mbr))
                continue; //no intersection we return
            if (!qualify(qrange, mbr, q.theta_sqr))
                continue; //judge false hit

            vector<RecordId> child;
            double **entry = nonleaf->entry_mbr;

            for (auto i = 0; i != nonleaf->entry_num; ++i)
            {
                if (intersect(NUM_DIMENSIONS, entry[i], inter))
                    child.push_back(i);
            }
            sort(child.begin(), child.end(), mbr_sort_dim_x_less(nonleaf));

            for (vector<RecordId>::iterator iter = child.begin(); iter != child.end(); ++iter)
            {
                RecordId i = *iter;
                if (l.locx+q.theta >= entry[i][0])
                {
                    if (l.locx-q.theta <= entry[i][1])
                    {
                        if ((l.locy+q.theta >= entry[i][2]) && (l.locy-q.theta <= (entry[i][3])))
                        {
                            Node* child = nonleaf->get_child(i);
//                            Q.push(ScoredNode(child, overlap(2, mbr, qrange))); // TODO
//                            cout << "dist " << MinMinDist(mbr, qrange, q.theta_sqr) << "\n";
                            float tmp = MinMinDist(mbr, qrange, q.theta);
//                            float tmp2 = sqrt(tmp);
//                            cout << tmp;
//                            q.insert(tmp2);
                            counter++;
//                            cout << "insertions\n";
                            Q.push(ScoredNode(child, sqrt(tmp) )); // TODO
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
    }
//    cout << "COUNTE " << counter << endl;
}


void secSetup(CkJResult* q, Relation& L, Relation& R, RStarTree<TreeDataP, double>* rtR ){

    auto l = L.beginSubRelation();

    while (l != L.endSubRelation()) {
        ContainmentkNNQueryOnRTree(*q, **l, R, rtR);
        l++;
    }
}




// 3. SETUP X3
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void RangeQueryWithTextProbe(CkJResult& q, Record& l, Relation& R, const RStarTree<TreeDataP, double> *rtR, InvertedIndex* iidx){
    size_t result = 0;
    double qrange[4];
    float score;
    PqNode Q;
    
    
    qrange[0] = l.locx - q.theta;
    qrange[1] = l.locx + q.theta;
    qrange[2] = l.locy - q.theta;
    qrange[3] = l.locy + q.theta;

    
    Q.push(make_pair(rtR->root, 1));
    while (!Q.empty())
    {
        
        if(q.numResult == q.k &&  Q.top().second > q.theta){
//        if( Q.top().second > q.theta){
            Q.pop();
            continue;
        }
        Node* cnode = Q.top().first;
//        cout << Q.top().second << "\n";
        Q.pop();
        
//        if(cnode->textScore == -1) // node do not satisfies the containment query
//            continue;
//
//        if(cnode->textScore == 0){ // node not visited as far
//            if( textSearchOnIRTree(q.termSetR, q.termSetSizeR, &iidx[cnode->page_id]) ){
//                cnode->textScore = 1;
//            }else{
//                cnode->textScore = -1;
//            }
//        }


        if(cnode->is_leaf_node() && cnode->textScore != -1){
            Leaf* leaf = static_cast<RSTLeafNode<TreeDataP, double> *>(cnode);
            double mbr[4], inter[4];
            
            if (leaf->verifiedRecords->size() == 0){
//                if(!textSearchOnNodeBinSearch(q.termSetR, q.termSetSizeR, &(iidx[leaf->page_id]), leaf)) //TBD CHANGE NOT WORKING
                if(!textSearchOnLeaf(q.termSetR, q.termSetSizeR, &(iidx[leaf->page_id]), leaf)) //TBD CHANGE NOT WORKING

                    continue;
            }
            
            leaf->get_mbr(mbr);
            mbr[0] -= q.theta;
            mbr[1] += q.theta;
            mbr[2] -= q.theta;
            mbr[3] += q.theta;

            if (!intersect(NUM_DIMENSIONS, inter, qrange, mbr))
                continue; //no intersection we return
            if (!(qualify(qrange, mbr, q.theta_sqr)))
                continue; //judge false hit

            TreeDataP<double> **entry = leaf->data;
            vector<RecordId> child;
            
            
//            for (auto i = 0; i < leaf->entry_num; i++)
//                if (inside(NUM_DIMENSIONS, entry[i]->data, inter))
//                    child.push_back(i);
            
            for (auto i = 0; i < leaf->verifiedRecords->size(); i++){
                if (inside(NUM_DIMENSIONS, entry[leaf->verifiedRecords->at(i)]->data, inter))
                    child.push_back(leaf->verifiedRecords->at(i));
            }
            
            
            
            
            
            sort(child.begin(), child.end(), p_sort_dim_x_less(leaf));

            for (vector<RecordId>::iterator iter = child.begin(); iter != child.end(); ++iter){
                RecordId i = *iter;
                double locx = entry[i]->data[0];
                double locy = entry[i]->data[1];

                if (l.locx < (locx-q.theta))
                    continue;
                else if (l.locx > (locx+q.theta))
                    break;

                if ((l.locy+q.theta >= locy) && (l.locy - q.theta <= locy) && (qualify(l, locx, locy, q.theta_sqr)))
//                    if( textProbe(q.termSetR, q.termSetSizeR, *R[entry[i]->id]) ){ //TODO i do not need this
#ifdef RetrieveResults
                        q.insert(PairRecord(&l, R[entry[i]->id], l.id, entry[i]->id, distance2N(l, *R[entry[i]->id])   ));
#else
                        q.insert( distance2N(l, *R[entry[i]->id]) );
#endif
//                    }
                 }
        }else{
            RSTNonLeafNode<TreeDataP, double> *nonleaf = static_cast<RSTNonLeafNode<TreeDataP, double> *>(cnode);
            double mbr[4], inter[4];

            nonleaf->get_mbr(mbr);
            mbr[0] -= q.theta;
            mbr[1] += q.theta;
            mbr[2] -= q.theta;
            mbr[3] += q.theta;

            if (!intersect(NUM_DIMENSIONS, inter, qrange, mbr))
                continue; //no intersection we return
            if (!qualify(qrange, mbr, q.theta_sqr))
                continue; //judge false hit

            vector<RecordId> child;
            double **entry = nonleaf->entry_mbr;

            for (auto i = 0; i != nonleaf->entry_num; ++i){
                if ( (nonleaf->children[i]->textScore > -1) )
                    if( (intersect(NUM_DIMENSIONS, entry[i], inter) ) )
                        child.push_back(i);
            }
            sort(child.begin(), child.end(), mbr_sort_dim_x_less(nonleaf));

            for (vector<RecordId>::iterator iter = child.begin(); iter != child.end(); ++iter)
            {
                RecordId i = *iter;
                if (l.locx+q.theta >= entry[i][0])
                {
                    if (l.locx-q.theta <= entry[i][1])
                    {
                        if ((l.locy+q.theta >= entry[i][2]) && (l.locy-q.theta <= (entry[i][3]))){
                            RSTNode<TreeDataP, double> *child = nonleaf->get_child(i);
                            if((child->textScore == 1) || textProbeOnNode(q.termSetR, q.termSetSizeR, &iidx[child->page_id], child) ){
//                                Q.push(ScoredNode(child, overlap(2, mbr, qrange)));
                            float tmp = MinMinDist(mbr, qrange, q.theta);
//                            float tmp2 = sqrt(tmp);
//                            cout << tmp;
//                            q.insert(tmp2);
                            Q.push(ScoredNode(child, sqrt(tmp) )); // TODO
                            }
                        }
                    }
                    else
                    {
                        break;
                    }
                }
            }
        }
    }
}






void thirdSetup(CkJResult& q, Relation& L, Relation& R, InvertedIndex* IFRT_R, RTree* RT_R){
    int RNumKeywords = q.termSetSizeR;
    int *RKeywords   = q.termSetR;
    double threshold = q.theta;

    deque<Node*>* treeIterator = new deque<Node*>;
    vector<Node*>* textCandidates = new vector<Node*>;
    
    NonLeaf* node = (NonLeaf*) RT_R->root;
    
//    treeIterator->push_back(node);
    
    //determineTextCandidates(&q, IFRT_R , node, 'R');

    Record** l = L.beginSubRelation();
    
    while (l != L.endSubRelation()) {
        RangeQueryWithTextProbe(q, **l, R, RT_R,IFRT_R);;
        l++;
    }


}



// 4. SETUP X4
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

void SpatialTextualDistanceFilterOnNonLeafs(CkJResult& q, NonLeaf* l, NonLeaf* r, InvertedIndex* irtL, InvertedIndex* irtR, PqNodes* Q){
    //according to T. Brinkhoff et al SIGMOD 93
    //first we restrict the search space, one difference is that we now expand the MBR -/+dist
    //we expanf the two nonleaf MBR and calculate the intersection I
    //we can prove that if any child MBR do not intersect I, then it will not be in the result set
    vector<int> childL;
    vector<int> childR;
    double mbrR[4], mbrS[4], inter[4];

    l->get_mbr(mbrR);
    mbrR[0] -= q.theta;
    mbrR[1] += q.theta; //expand x-axis
    mbrR[2] -= q.theta;
    mbrR[3] += q.theta; //expand y-axis
    r->get_mbr(mbrS);
    mbrS[0] -= q.theta;
    mbrS[1] += q.theta; //expand x-axis
    mbrS[2] -= q.theta;
    mbrS[3] += q.theta; //expand y-axis
    if (!intersect(2, inter, mbrR, mbrS))
        return; //no intersection we return
    if (!QualifyMinMinDistSqr(mbrR, mbrS, q.theta_sqr))
        return; //judge false hit

    //next we will pickup child which intersects with I
    //sort them according to the 1st dimension
    //since we are using distance intersection, wo we treat childL as normal
    //and childR as if they are all expanded by threshold
    double **entryR = l->entry_mbr, **entryS = r->entry_mbr;
    
    childL.clear();
    vector<RecordId>* lVerifiedRec = l->verifiedRecords;
    vector<RecordId>* rVerifiedRec = r->verifiedRecords;
    for (int i = 0; i < l->entry_num; i++){ // ORIGIN
        if (l->children[i]->textScore > -1 && intersect(2, entryR[i], inter)){
            childL.push_back(i);
        }
    }
    
//      for (int i = 0; i < lVerifiedRec->size(); i++) // MYVERSION
//          if (intersect(2, entryR[lVerifiedRec->at(i)], inter))
//              childL.push_back(lVerifiedRec->at(i));
    
    
    sort(childL.begin(), childL.end(), mbr_sort_dim_x_less(l));

    childR.clear();
    for (int i = 0; i < r->entry_num; i++) // ORIGIN
        if ( r->children[i]->textScore > -1 && intersect(2, entryS[i], inter))
            childR.push_back(i);
//    for (int i = 0; i < rVerifiedRec->size(); i++) // MYVERSION
//        if (intersect(2, entryS[rVerifiedRec->at(i)], inter))
//            childR.push_back(rVerifiedRec->at(i));
    
    sort(childR.begin(), childR.end(), mbr_sort_dim_x_less(r));

    //next we will use plane-sweep to test intersection
    vector<int>::iterator iter1 = childL.begin(), iter2 = childR.begin();
    vector<int>::iterator end1 = childL.end(), end2 = childR.end();
    while ((iter1 < end1) && (iter2 < end2))
    {
        if (entryR[*iter1][0] < (entryS[*iter2][0] - q.theta)) //now asymmetric, we use two hand-coded internal loops
        {
            //internal loop 1
            vector<int>::iterator iter = iter2;
            while ((iter < end2)
                    && ((entryS[*iter][0] - q.theta) < entryR[*iter1][1]))
            {
                if ((entryR[*iter1][2] > (entryS[*iter][3] + q.theta))
                    || ((entryS[*iter][2] - q.theta) > entryR[*iter1][3]))
                {
                    ++iter;
                    continue;
                }

                RSTNode<TreeDataP, double> *nodeR = l->get_child(*iter1);
                RSTNode<TreeDataP, double> *nodeS = r->get_child(*iter);
                nodeR->get_mbr(mbrR);
                nodeS->get_mbr(mbrS);
//                cout << "IAMHERNEWPAIR\n";
                if(nodeR->textScore == 1 || textProbeOnNode(q.termSetL, q.termSetSizeL, &irtL[nodeR->page_id], nodeR))
                    if(nodeS->textScore == 1 || textProbeOnNode(q.termSetR, q.termSetSizeR, &irtR[nodeS->page_id], nodeS)){
                        float tmp = MinMinDist(mbrR, mbrS, q.theta);
                        Q->push(make_tuple(nodeR, nodeS, tmp));

                    }
//                Q.push(make_tuple(nodeR, nodeS, overlap(2, mbrR, mbrS)));
                ++iter;
            }

            ++iter1;
        }
        else
        {
            //internal loop 2
            vector<int>::iterator iter = iter1;
            while ((iter < end1)
                    && (entryR[*iter][0] < (entryS[*iter2][1] + q.theta)))
            {
                if (((entryS[*iter2][2] - q.theta) > entryR[*iter][3])
                        || (entryR[*iter][2] > (entryS[*iter2][3] + q.theta)))
                {
                    ++iter;
                    continue;
                }

                RSTNode<TreeDataP, double> *nodeR = l->get_child(*iter);
                RSTNode<TreeDataP, double> *nodeS = r->get_child(*iter2);
                nodeR->get_mbr(mbrR);
                nodeS->get_mbr(mbrS);
//                cout << "IAMHERNEWPAIR\n";
                if(nodeR->textScore == 1 || textProbeOnNode(q.termSetL, q.termSetSizeL, &irtL[nodeR->page_id], nodeR))
                    if(nodeS->textScore == 1 || textProbeOnNode(q.termSetR, q.termSetSizeR, &irtR[nodeS->page_id], nodeS)){
                        float tmp = MinMinDist(mbrR, mbrS, q.theta);
                        Q->push(make_tuple(nodeR, nodeS, tmp));
                    }
//                        Q->push(make_tuple(nodeR, nodeS, Min(2, mbrR, mbrS))); //TODO

//                        Q->push(make_tuple(nodeR, nodeS, overlap(2, mbrR, mbrS))); //TODO

                ++iter;
            }

            ++iter2;
        }
    }
    
    
//    cout << "IN FILTER  " << Q.size() << "  \n\n";
}

void SpatialDistanceFilterONLeafAndNonLeaf(CkJResult& q, Leaf* l, NonLeaf* r, PqNodes* Q){
    vector<int> childR;
    double threshold = q.theta;
    double threshold_sqr = q.theta_sqr;
    double mbrL[4], mbrR[4], inter[4];
    
    l->get_mbr(mbrR);
    mbrL[0] -= threshold;
    mbrL[1] += threshold; //expand x-axis
    mbrL[2] -= threshold;
    mbrL[3] += threshold; //expand y-axis
    double **entryR = r->entry_mbr;
    
    
//    for (int i = 0; i < r->entry_num; i++){
//        if ( r->children[i]->textScore == 1 && intersect(2, entryR[i], mbrL))
////            stack->push(make_tuple((Node*)l, (Node*)entryR[i]));
//    }
}

void SpatialJoinOnLeafs(CkJResult& q, Leaf* lLeaf, Leaf* rLeaf, Relation& L, Relation& R){
    TreeDataP<double> **lEntry = lLeaf->data;
    TreeDataP<double> **rEntry = rLeaf->data;
    vector<RecordId> lChild;
    vector<RecordId> rChild;
    
    double mbrR[4], mbrS[4], inter[4];
    unsigned int numResults = 0;

    
    lLeaf->get_mbr(mbrR);
    mbrR[0] -= q.theta;
    mbrR[1] += q.theta; //expand x-axis
    mbrR[2] -= q.theta;
    mbrR[3] += q.theta; //expand y-axis
    rLeaf->get_mbr(mbrS);
    mbrS[0] -= q.theta;
    mbrS[1] += q.theta; //expand x-axis
    mbrS[2] -= q.theta;
    mbrS[3] += q.theta; //expand y-axis
    if (!intersect(2, inter, mbrR, mbrS))
        return ; //no intersection we return
    if (!QualifyMinMinDistSqr(mbrR, mbrS, q.theta_sqr))
        return ; //judge false hit

    //next we will pickup child which intersects with I
    //sort them according to the 1st dimension
    //since we are using distance intersection, wo we treat childR as normal
    //and childS as if they are all expanded by threshold
    TreeDataP<double> **entryL = lLeaf->data;
    lChild.clear();
    TreeDataP<double> **entryR = rLeaf->data;
    rChild.clear();
    vector<RecordId>* lVerifiedRec = lLeaf->verifiedRecords;
    vector<RecordId>* rVerifiedRec = rLeaf->verifiedRecords;
    
    
    for (int i = 0; i < lVerifiedRec->size(); i++) // ORIGIN
        if (inside(2, entryL[lVerifiedRec->at(i)]->data, inter)){
//            cout << "IAMHERE\n";

            lChild.push_back(lVerifiedRec->at(i));
        }
    
    for (int i = 0; i < rVerifiedRec->size(); i++) // ORIGIN
        if (inside(2, entryR[rVerifiedRec->at(i)]->data, inter))
            rChild.push_back(rVerifiedRec->at(i));
    
    
//    for (int i = 0; i < lLeaf->entry_num; i++){
//        if (inside(2, entryR[i]->data, inter))
//            lChild.push_back(i);
//    }
//
//
//    for (int i = 0; i < rLeaf->entry_num; i++){
//        if (inside(2, entryS[i]->data, inter))
//            rChild.push_back(i);
//    }
    
    
    
    
    sort(lChild.begin(), lChild.end(), p_sort_dim_x_less(lLeaf));
    sort(rChild.begin(), rChild.end(), p_sort_dim_x_less(rLeaf));




    vector<RecordId>::iterator iter1 = lChild.begin(), iter2 = rChild.begin();
    vector<RecordId>::iterator end1 = lChild.end(), end2 = rChild.end();
    while ((iter1 < end1) && (iter2 < end2))
    {
        if (entryL[*iter1]->data[0] < (entryR[*iter2]->data[0] - q.theta)) //now asymmetric, we use two hand-coded internal loops
        {
            //internal loop 1
            vector<RecordId>::iterator iter = iter2;
            while ((iter < end2) && ((entryR[*iter]->data[0] - q.theta) < entryL[*iter1]->data[0])){
                int id1 = entryL[*iter1]->id;
                int id2 = entryR[*iter]->id;

                if ((entryL[*iter1]->data[1] > (entryR[*iter]->data[1] + q.theta)) || ((entryR[*iter]->data[1] - q.theta) > entryL[*iter1]->data[1])){
                    ++iter;
                    continue;
                }

                if ((QualifySpatial(L[id1], R[id2], q.theta_sqr))) // && (QualifyTextual((*R)[id1], (*S)[id2], numKwR, kwR, numKwS, kwS)))
                {
#ifdef RetrieveResults
                    q.insert(PairRecord( L[id1], R[id2], id1, id2, distance2N(*L[id1], *R[id2])));
#else
                    q.insert();
#endif
                }

                ++iter;
            }

            ++iter1;
        }
        else
        {
            //internal loop 2
            vector<RecordId>::iterator iter = iter1;
            while ((iter < end1) && (entryL[*iter]->data[0] < (entryR[*iter2]->data[0] + q.theta)))
            {
                int id1 = entryL[*iter]->id;
                int id2 = entryR[*iter2]->id;

                if (((entryR[*iter2]->data[1] - q.theta) > entryL[*iter]->data[1]) || (entryL[*iter]->data[1] > (entryR[*iter2]->data[1] + q.theta))){
                    ++iter;
                    continue;
                }

                if ((QualifySpatial(L[id1], R[id2], q.theta_sqr))) // && (QualifyTextual((*R)[id1], (*S)[id2], numKwR, kwR,numKwS, kwS)))
                {
#ifdef RetrieveResults
                    q.insert(PairRecord( L[id1], R[id2], id1, id2, distance2N(*L[id1], *R[id2])));
#else
                    q.insert();
#endif
                }
                ++iter;
            }
            ++iter2;
        }
    }

    // NAIVE VERSION
//    for(int i = 0; i < lChild.size(); i++){
//        for(int j = 0; j < rChild.size(); j++){
//            if(distance2N(*L[lEntry[lChild.at(i)]->id], *R[rEntry[rChild.at(j)]->id]) < q.theta ){
//                q.insert(PairRecord(L[lEntry[lChild.at(i)]->id], R[rEntry[rChild.at(j)]->id], L[lEntry[lChild.at(i)]->id]->id, R[rEntry[rChild.at(j)]->id]->id, distance2N(*L[lEntry[lChild.at(i)]->id], *R[rEntry[rChild.at(j)]->id])   ));
//            }
//        }
//    }
}




void IRTreeJoin(CkJResult& q, Relation& L, Relation& R, RTree* rt_L, RTree* rt_R, InvertedIndex* iix_L, InvertedIndex* iix_R){
    int LNumKeywords = q.termSetSizeL;
    int* LKeywords   = q.termSetL;
    
    int RNumKeywords = q.termSetSizeR;
    int* RKeywords   = q.termSetR;

    double threshold = q.theta;
    
//    stack<NodesPair>* candidates = new stack<NodesPair>;
    PqNodes* Q = new PqNodes;
    
    
    stack<Node*>* lCandidates = new stack<Node*>;
    stack<Node*>* rCandidates = new stack<Node*>;

    
    ScoredPair* pair;
    Node* l = rt_L->root;
    Node* r = rt_R->root;
    Node* lChild;
    Node* rChild;
    NonLeaf* lNonLeaf;
    NonLeaf* rNonLeaf;
    Leaf* lLeaf;
    Leaf* rLeaf;
    l->is_leaf_node(); r->is_leaf_node();

//    candidates->push(make_tuple(l, r));
    (*Q).push(make_tuple(l, r, 1));
    
    
    
    while (!Q->empty()) {
//        pair = &(Q->top());
        float tmp = get<2>(Q->top());
        if(q.numResult == q.k &&  tmp > q.theta){
//        if( Q.top().second > q.theta){
            Q->pop();
            continue;
        }
        
        
        l = get<0>(Q->top());
        r = get<1>(Q->top());
        
//        l = get<0>(*pair);
//        r = get<1>(*pair);
        Q->pop();
//        cout << "IAMHER\n";

        if(r->is_leaf_node() && l->is_leaf_node()){ // if r leaf then l leaf -> compute pairs
            lLeaf = (Leaf*)l;
            rLeaf = (Leaf*)r;
            
            
            if(lLeaf->textScore == 0 || lLeaf->verifiedRecords->size() == 0)
                textSearchOnLeaf(LKeywords, LNumKeywords, &iix_L[lLeaf->page_id], (Node*)lLeaf);
            if(lLeaf->textScore == -1)
                continue;
            //textProbe on l
            if(rLeaf->textScore == 0 || rLeaf->verifiedRecords->size() == 0)
                textSearchOnLeaf(RKeywords, RNumKeywords, &iix_R[rLeaf->page_id], (Node*)rLeaf);
            if(rLeaf->textScore == -1)
                continue;
            SpatialJoinOnLeafs(q, lLeaf, rLeaf, L, R);


        }else if(l->is_leaf_node() && !r->is_leaf_node()){ // i only l leaf -> travers R
            cout << "R deph > L depht \n";
            lLeaf = (Leaf*)l;
            rNonLeaf = (NonLeaf*)r;
            //textProbe on l
            if(lLeaf->textScore == 0)
                textSearchOnLeaf(LKeywords, LNumKeywords, &iix_L[lLeaf->page_id], lLeaf);
            if(lLeaf->textScore == -1)
                continue;

            //textProbe on all children of r
            for(int j = 0; r->entry_num; j++){
                rChild = rNonLeaf->children[j];
                if(rChild->textScore == 0) // Containment Query on l if no did as far
                    textVerificationOnNode(RKeywords, RNumKeywords, &iix_R[rChild->page_id], rChild);

//                if(rChild->textScore == -1)
//                    continue;
            }
            SpatialDistanceFilterONLeafAndNonLeaf(q, lLeaf, rNonLeaf, Q);

        }else if ( !l->is_leaf_node() && !r->is_leaf_node() ){ //travers to next level in both trees
            lNonLeaf = (NonLeaf*)l;
            rNonLeaf = (NonLeaf*)r;

            textVerificationOnNode(q.termSetL, q.termSetSizeL, &iix_L[lNonLeaf->page_id], lNonLeaf);
            textVerificationOnNode(q.termSetR, q.termSetSizeR, &iix_R[rNonLeaf->page_id], rNonLeaf);
//            cout << "IAMHER\n";

            SpatialTextualDistanceFilterOnNonLeafs(q, lNonLeaf, rNonLeaf, iix_L, iix_R, Q);
//            cout << "AFTER FILTER  " << Q.size() << "  \n\n";

        }else{
            cout << "I SHOULD NEVER HERE: IR-TREE JOIN \n\n\n";
        }
    }

}















// 5. SETUP X5
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SpatialJoinOnLeafs2(CkJResult& q, Leaf& l, Leaf& r, Relation& L, Relation& R){
    //according to T. Brinkhoff et al SIGMOD 93
    //first we restrict the search space, one difference is that we now expand the MBR -/+dist
    //we expanf the two nonleaf MBR and calculate the intersection I
    //we can prove that if any child MBR do not intersect I, then it will not be in the result set
    double mbrR[4], mbrS[4], inter[4];
    unsigned int numResults = 0;

    Leaf* leafR = &l;
    Leaf* leafS = &r;
    vector<int> childR, childS;
    
    leafR->get_mbr(mbrR);
    mbrR[0] -= q.theta;
    mbrR[1] += q.theta; //expand x-axis
    mbrR[2] -= q.theta;
    mbrR[3] += q.theta; //expand y-axis
    leafS->get_mbr(mbrS);
    mbrS[0] -= q.theta;
    mbrS[1] += q.theta; //expand x-axis
    mbrS[2] -= q.theta;
    mbrS[3] += q.theta; //expand y-axis
    
    if (!intersect(2, inter, mbrR, mbrS))
        return; //no intersection we return
    if (!QualifyMinMinDistSqr(mbrR, mbrS, q.theta_sqr))
        return; //judge false hit

    //next we will pickup child which intersects with I
    //sort them according to the 1st dimension
    //since we are using distance intersection, wo we treat childR as normal
    //and childS as if they are all expanded by threshold
    TreeDataP<double> **entryR = leafR->data;
    childR.clear();
    TreeDataP<double> **entryS = leafS->data;
    childS.clear();

    for (int i = 0; i < leafR->entry_num; i++){
        if ( entryR[i]->textScore != -1){
            if(inside(2, entryR[i]->data, inter) ){ // +Filter data, which do not satisfies text probe
                childR.push_back(i);
            }
        }
    }
    sort(childR.begin(), childR.end(), p_sort_dim_x_less(leafR));
    for (int i = 0; i < leafS->entry_num; i++){
        if ( entryS[i]->textScore != -1){
            if( inside(2, entryS[i]->data, inter) ) // +Filter data, which do not satisfies text probe
                childS.push_back(i);
        }
    }
    sort(childS.begin(), childS.end(), p_sort_dim_x_less(leafS));

    
    
    vector<int>::iterator iter1 = childR.begin(), iter2 = childS.begin();
    vector<int>::iterator end1 = childR.end(), end2 = childS.end();
    
    
    
    
    while ((iter1 < end1) && (iter2 < end2))
    {
        if (entryR[*iter1]->data[0] < (entryS[*iter2]->data[0] - q.theta)) //now asymmetric, we use two hand-coded internal loops
        {
            //internal loop 1
            vector<int>::iterator iter = iter2;
            while ((iter < end2) && ((entryS[*iter]->data[0] - q.theta) < entryR[*iter1]->data[0]))
            {
                int id1 = entryR[*iter1]->id;
                int id2 = entryS[*iter]->id;
                if(
                    ((entryR[*iter1]->data[1] > (entryS[*iter]->data[1] + q.theta)) ||
                     ((entryS[*iter]->data[1] - q.theta) > entryR[*iter1]->data[1]))){
                    ++iter;
                    continue;
                }
                if(entryR[*iter1]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetL, q.termSetSizeL, L[id1], &(entryR[*iter1]->textScore)))
                    if(entryS[*iter]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[id2], &(entryS[*iter]->textScore))){
                       double d = distance2N(*L[id1], *R[id2]);
                      if (d < q.theta){

#ifdef RetrieveResults
                        q.insert(PairRecord( L[id1], R[id2], id1, id2, d));
#else
                        q.insert();
#endif
                      }
                }

                ++iter;
            }

            ++iter1;
        }
        else
        {
            //internal loop 2
            vector<int>::iterator iter = iter1;
            while ((iter < end1) && (entryR[*iter]->data[0] < (entryS[*iter2]->data[0] + q.theta))){

                int id1 = entryR[*iter]->id;
                int id2 = entryS[*iter2]->id;

                if (((entryS[*iter2]->data[1] - q.theta) > entryR[*iter]->data[1]) || (entryR[*iter]->data[1] > (entryS[*iter2]->data[1] + q.theta))){
                    ++iter;
                    continue;
                }

                if(entryR[*iter]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetL, q.termSetSizeL, L[id1], &(entryR[*iter]->textScore)))
                    if(entryS[*iter2]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[id2], &(entryS[*iter2]->textScore))){
                        double d = distance2N(*L[id1], *R[id2]);
                        if (d < q.theta){
                    
#ifdef RetrieveResults
                            q.insert(PairRecord( L[id1], R[id2], id1, id2, d));
#else
                            q.insert();
#endif
                     }
                            
                }

                ++iter;
            }

            ++iter2;
        }
    }
    
    
// NAIVE VERSION
//    for(int i = 0; i < leafR->entry_num; i++){
//        for(int j = 0; j < leafS->entry_num; j++){
//            int id1 = entryR[i]->id;
//            int id2 = entryS[j]->id;
//            if(
//               textProbe(q.termSetL, q.termSetSizeL, *L[id1]) &&
//               textProbe(q.termSetR, q.termSetSizeR, *R[id2]) &&
//               distance2N(*L[id1], *R[id2]) < threshold){
//                q.insert(PairRecord( L[id1], R[id2], id1, id2, distance2N(*L[id1], *R[id2])));
//            }
//        }
//    }


    
    
}


void SpatialDistanceFilterONLeafAndNonLeaf(CkJResult& q, Leaf& l, NonLeaf& r, stack<NodesPair>* stack){
    
}

void SpatialDistanceFilterOnNonLeafs(CkJResult& q, NonLeaf& l, NonLeaf& r, PqNodes* Q ){
    //according to T. Brinkhoff et al SIGMOD 93
    //first we restrict the search space, one difference is that we now expand the MBR -/+dist
    //we expanf the two nonleaf MBR and calculate the intersection I
    //we can prove that if any child MBR do not intersect I, then it will not be in the result set
    double mbrR[4], mbrS[4], inter[4];
    double threshold = q.theta;
    double threshold_sqr = q.theta_sqr;
    NonLeaf* nonleafR = &l;
    NonLeaf* nonleafS = &r;
    vector<int> childR, childS;

    nonleafR->get_mbr(mbrR);
    mbrR[0] -= threshold;
    mbrR[1] += threshold; //expand x-axis
    mbrR[2] -= threshold;
    mbrR[3] += threshold; //expand y-axis
    nonleafS->get_mbr(mbrS);
    mbrS[0] -= threshold;
    mbrS[1] += threshold; //expand x-axis
    mbrS[2] -= threshold;
    mbrS[3] += threshold; //expand y-axis
    if (!intersect(2, inter, mbrR, mbrS))
        return; //no intersection we return
    if (!QualifyMinMinDistSqr(mbrR, mbrS, threshold_sqr))
        return; //judge false hit

    //next we will pickup child which intersects with I
    //sort them according to the 1st dimension
    //since we are using distance intersection, wo we treat childR as normal
    //and childS as if they are all expanded by threshold
    double **entryR = nonleafR->entry_mbr, **entryS = nonleafS->entry_mbr;
    childR.clear();
    for (int i = 0; i < nonleafR->entry_num; i++)
        if (intersect(2, entryR[i], inter))
            childR.push_back(i);
    sort(childR.begin(), childR.end(), mbr_sort_dim_x_less(nonleafR));

    childS.clear();
    for (int i = 0; i < nonleafS->entry_num; i++)
        if (intersect(2, entryS[i], inter))
            childS.push_back(i);
    sort(childS.begin(), childS.end(), mbr_sort_dim_x_less(nonleafS));

    //next we will use plane-sweep to test intersection
    vector<int>::iterator iter1 = childR.begin(), iter2 = childS.begin();
    vector<int>::iterator end1 = childR.end(), end2 = childS.end();
    while ((iter1 < end1) && (iter2 < end2))
    {
        if (entryR[*iter1][0] < (entryS[*iter2][0] - threshold)) //now asymmetric, we use two hand-coded internal loops
        {
            //internal loop 1
            vector<int>::iterator iter = iter2;
            while ((iter < end2)
                    && ((entryS[*iter][0] - threshold) < entryR[*iter1][1]))
            {
                if ((entryR[*iter1][2] > (entryS[*iter][3] + threshold))
                        || ((entryS[*iter][2] - threshold) > entryR[*iter1][3]))
                {
                    ++iter;
                    continue;
                }

                RSTNode<TreeDataP, double> *nodeR = nonleafR->get_child(*iter1);
                RSTNode<TreeDataP, double> *nodeS = nonleafS->get_child(*iter);
                Q->push(make_tuple(nodeR, nodeS, overlap(2, nodeR->get_mbr(mbrR), nodeS->get_mbr(mbrS))));
                ++iter;
            }

            ++iter1;
        }
        else
        {
            //internal loop 2
            vector<int>::iterator iter = iter1;
            while ((iter < end1)
                    && (entryR[*iter][0] < (entryS[*iter2][1] + threshold)))
            {
                if (((entryS[*iter2][2] - threshold) > entryR[*iter][3])
                        || (entryR[*iter][2] > (entryS[*iter2][3] + threshold)))
                {
                    ++iter;
                    continue;
                }

                RSTNode<TreeDataP, double> *nodeR = nonleafR->get_child(*iter);
                RSTNode<TreeDataP, double> *nodeS = nonleafS->get_child(*iter2);
                Q->push(make_tuple(nodeR, nodeS, overlap(2, nodeR->get_mbr(mbrR), nodeS->get_mbr(mbrS))));
                ++iter;
            }

            ++iter2;
        }
    }

    
    
    
}


void RTreeJoinAndTextVerification(CkJResult& q, Relation& L, Relation& R, RTree* rt_L, RTree* rt_R){
    int LNumKeywords = q.termSetSizeL;
    int* LKeywords   = q.termSetL;
    int RNumKeywords = q.termSetSizeR;
    int* RKeywords   = q.termSetR;

    double threshold = q.theta;
    
    PqNodes* Q = new PqNodes;

    
    NodesPair* pair;
    Node* l = rt_L->root;
    Node* r = rt_R->root;
    Node* lChild;
    Node* rChild;
    NonLeaf* lNonLeaf;
    NonLeaf* rNonLeaf;
    Leaf* lLeaf;
    Leaf* rLeaf;
    
    Q->push(make_tuple(l, r, 1));
    
    while (!Q->empty()) {
        l = get<0>(Q->top());
        r = get<1>(Q->top());
        Q->pop();

        
        if(r->is_leaf_node() && l->is_leaf_node()){ // if r leaf then l leaf -> compute pairs
            lLeaf = (Leaf*)l;
            rLeaf = (Leaf*)r;

            SpatialJoinOnLeafs2(q, *lLeaf, *rLeaf, L, R);
            

        }else if(l->is_leaf_node() && !r->is_leaf_node()){ // i only l leaf -> travers R
            cout << "R deph > L depht \n";

            
            
//            SpatialDistanceFilterONLeafAndNonLeaf(q, lLeaf, rNonLeaf, Q);
                
        }else if ( !l->is_leaf_node() && !r->is_leaf_node() ){ //travers to next level in both trees
            lNonLeaf = (NonLeaf*)l;
            rNonLeaf = (NonLeaf*)r;
            
            
            SpatialDistanceFilterOnNonLeafs(q, *lNonLeaf, *rNonLeaf, Q);
        }else{
            cout << "I SHOULD NEVER HERE: IR-TREE JOIN \n\n\n";
        }
    }
}
