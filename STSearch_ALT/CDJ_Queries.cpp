//
//  CDJ_Queries.cpp
//  STSearch_
//
//  Created by Artur Titkov on 11.06.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.

#include "def.h"
#include "inverted_index.h"
#include "ResultSet.h"
#include "relation.h"
#include "utils_spatial.cpp"
#include "utils_boolean.cpp"


typedef tuple<Record*, Node*> RecNodePair;
typedef tuple<Node*, Node*> NodesPair;
typedef TreeDataP<double> Point;


//SPATIL
bool qualify(const Record &r, const Record &s, const double dthreshold_sqr);
bool qualify(const Record &r, const double locx, const double locy, const double dthreshold_sqr);
bool qualify(const double rlocx, const double rlocy, const double slocx, const double slocy, const double dthreshold_sqr);
bool qualify(double *mbrR, double *mbrS, double dthreshold_sqr);
RStarTree<TreeDataP, double>* bulkload(int dim, int page_len, TreeDataP<double> **data, int data_num);
bool qualify(const Record &r, const Record &s, const double dthreshold_sqr);
bool qualify(const Record &r, const double locx, const double locy, const double dthreshold_sqr);
bool qualify(double *mbrR, double *mbrS, double dthreshold_sqr);


//TEXTUAL
//Relation* ContainmentQueryWithBinSearchOfRelation(Relation& T, int numKeywords, int *keywords, InvertedIndex *iidx);

//inline bool textProbeOnDataPoint(termSet t, termSetSize n, Record* r, float* p);
//inline bool textProbeOnDataPointBinSearch(termSet t, termSetSize n, Record* r, float* p);


inline bool textSearchOnLeaf(int* keywords, int numKeywords, InvertedIndex* iidx, Node* leaf);
inline bool textProbeOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);
inline bool textSearchOnLeaf(int* keywords, int numKeywords, InvertedIndex* iidx, Node* leaf);
inline bool textSearchOnNodeBinSearch(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);
inline bool textVerificationOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);


//UTILS
//static void Print(vector<InvertedListEntry>* list){
//    cout << "Posting List size :" << list->size() << "\n";
//    for(int i = 0; i < list->size(); i++){
//        cout << list->at(i).rec->id << " ";
//    }
//    cout << "\n";
//
//
//}

inline void CDJResult::insert()
{
    numResult++;
};

inline void CDJResult::insert(PairRecord rec)
{
    numResult++;
    result.push_back(rec);
};


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




inline bool textProbeOnData(termSet t, termSetSize n, Record& r, Data* d){
    int i = 0;
    int j = 0;
    while( i < n && j < r.length){
        if(t[i] < r.keywords[j]){
            d->textScore = -1;
            return false;
        }else if(t[i] > r.keywords[j]){
            j++;
        }else if(t[i] == r.keywords[j]) {
            i++;
            j++;
        }
    }
    if(i == n){
        d->textScore = 1;
        return true;;
    }else{
        d->textScore = -1;
        return false;
    }
};




//bool intersect(int dim, double *mbr, double *mbrR, double *mbrS);

//inline float distance2N(const Record &r, const Record &s){
//    return sqrt( (r.locx - s.locx) * (r.locx - s.locx) + (r.locy - s.locy) * (r.locy - s.locy) );
//};
//inline float distance1N(const Record &r, const Record &s){
//    return  (r.locx - s.locx) * (r.locx - s.locx) + (r.locy - s.locy) * (r.locy - s.locy) ;
//};



//class CompareKeywordsByFrequencyR
//{
//public:
//    bool operator() (const int& klhs, const int& krhs) const
//    {
//        return (iidxR->lists[klhs].size() < iidxR->lists[krhs].size());
//    }
//};




//int binary_search_find_index(std::vector<int> v, int data) {
//    auto it = std::lower_bound(v.begin(), v.end(), data);
//    if (it == v.end() || *it != data) {
//        return -1;
//    } else {
//        std::size_t index = std::distance(v.begin(), it);
//        return index;
//    }
//}


//TBD REMOVE LATER










// Baseline x0
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
void NestedLoopsDistanceJoin(CDJResult* q, Relation* R, Relation* S){
    auto dthreshold_sqr = q->theta_sqr;
    Record x,y;
    auto r = R->beginSubRelation();
    auto s = S->beginSubRelation();
    auto lastR = R->endSubRelation();
    auto lastS = S->endSubRelation();
//    for(; r  != lastR; r++ ){
//        for(; s != lastS; s++){
//            if(distance2N(**r, **s) < q->theta){
//
        for(int r = 0; r < R->numRecords; r++){
            x = *R->beginSubRelation()[r];
            for(int s = 0; s < S->numRecords; s++){
                y = *S->beginSubRelation()[s];
                if (qualify(x, y, dthreshold_sqr)){
//                if(distance2N(x, y) < q->theta){
                    
#ifdef RetrieveResults
            q->insert(PairRecord(&x, &y, x.id, y.id, distance2N(x, y)));
#else
            q->insert();
#endif
                }
            }
        }

}


// 1. Setup x1
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////

inline void internalLoop(CDJResult* q, const SubRelationIterator rec, const SubRelationIterator firstFS, const SubRelationIterator lastFS){
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


void PlaneSweepDistanceJoin(CDJResult* q, Relation* R,Relation* S){
    auto r = R->beginSubRelation();
    auto s = S->beginSubRelation();
    auto lastR = R->endSubRelation();
    auto lastS = S->endSubRelation();

    //Call Internal Loop
    //With OR Without Result
    while ((r < lastR) && (s < lastS)){
        if ((*r)->locx < (*s)->locx){
//        if ((*r) < (*s)){
            internalLoop(q, r, s, lastS);
            r++;
        }else{
            internalLoop(q, s, r, lastR);
            s++;
        }
    }
}



// 2. Setup x2
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void ContainmentRangeQueryOnRTree(CDJResult& q, Record& l, Relation& R,const RStarTree<TreeDataP, double> *rtR){
    size_t result = 0;
    double qrange[4];

    qrange[0] = l.locx - q.theta; // [(x1,y1),(x2,y2)]
    qrange[1] = l.locx + q.theta;
    qrange[2] = l.locy - q.theta;
    qrange[3] = l.locy + q.theta;


    stack<RSTNode<TreeDataP, double> *> S;
    S.push(rtR->root);
    while (!S.empty()){
        RSTNode<TreeDataP, double> *cnode = S.top();
        S.pop();

        if(cnode->is_leaf_node()){
            Leaf* leaf = static_cast<RSTLeafNode<TreeDataP, double> *>(cnode);
            double mbr[4], inter[4];

            leaf->get_mbr(mbr);
            mbr[0] -= q.theta;
            mbr[1] += q.theta;
            mbr[2] -= q.theta;
            mbr[3] += q.theta;

            
            // spatial probe leaf mbr
            if (!intersect(NUM_DIMENSIONS, inter, qrange, mbr))
                continue; //no intersection we return
            if (!(qualify(qrange, mbr, q.theta_sqr)))
                continue; //judge false hit

            
            Point **entry = leaf->data;
            vector<RecordId> child;
            
//            for (auto i = 0; i < leaf->entry_num; i++)
//                cout << "Score " << entry[i]->textScore << "\n";
            
//            for (auto i = 0; i < leaf->entry_num; i++)
//                if ( (entry[i]->textScore != -1) && inside(NUM_DIMENSIONS, entry[i]->data, inter))
//                    child.push_back(i);
            
            // for each leaf child Probe text and spatial condition
            for (auto i = 0; i < leaf->entry_num; i++)
                if ( (entry[i]->textScore != -1))
                    if(inside(NUM_DIMENSIONS, entry[i]->data, inter))
//                        if( (entry[i]->textScore == 1) || textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[entry[i]->id], &(entry[i]->textScore)) )
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
                        q.insert(PairRecord(&l, R[entry[i]->id], l.id, entry[i]->id, distance2N(l, *R[entry[i]->id])));
#else
                        q.insert();
#endif
//                        q.insert(PairRecord(&l, R[entry[i]->id], l.id, i, distance2N(l, *R[entry[i]->id])   ));
                    }
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
                if (intersect(NUM_DIMENSIONS, entry[i], inter))
                    child.push_back(i);
            }
            
            sort(child.begin(), child.end(), mbr_sort_dim_x_less(nonleaf));

            for (vector<RecordId>::iterator iter = child.begin(); iter != child.end(); ++iter){
                RecordId i = *iter;
                
                if (l.locx+q.theta >= entry[i][0]){
                    if (l.locx-q.theta <= entry[i][1]){
                        if ((l.locy+q.theta >= entry[i][2]) && (l.locy-q.theta <= (entry[i][3]))){
                            
                            Node* child = nonleaf->get_child(i);
                            S.push(child);
                        
                        }
                    }else{
                        break;
                    }
                }
            }
        }
    }
}


void secSetup(CDJResult* q, Relation& L, Relation& R, RStarTree<TreeDataP, double>* rtR ){

    auto l = L.beginSubRelation();

    while (l != L.endSubRelation()) {
        ContainmentRangeQueryOnRTree(*q, **l, R, rtR);
        l++;
    }
}



// 3. SETUP X3
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////








void ContainmentRangeQueryOnIRTree(CDJResult& q, Record& l, Relation& R, const RStarTree<TreeDataP, double> *rtR, InvertedIndex* iidx){
    size_t result = 0;
    double qrange[4];
    
    qrange[0] = l.locx - q.theta;
    qrange[1] = l.locx + q.theta;
    qrange[2] = l.locy - q.theta;
    qrange[3] = l.locy + q.theta;


    stack<Node*> S;
    S.push(rtR->root);
    while (!S.empty())
    {
        Node* cnode = S.top();
        S.pop();
        
        
        
        // TEST cure NODE
//        if(cnode->textScore == -1) // node do not satisfies the containment query
//            continue;
//        if(cnode->textScore == 0){ // node not visited as far
//            textSearchOnLeaf(q.termSetR, q.termSetSizeR, &iidx[cnode->page_id], cnode);
//            if(cnode->textScore == -1)
//                continue;
//
//            }


        if(cnode->is_leaf_node()){ // ------- LEAF
            Leaf* leaf = static_cast<RSTLeafNode<TreeDataP, double> *>(cnode);
            double mbr[4], inter[4];
//            cout << "SCORE " << leaf->textScore << "\n";
//            cout << "SCORE " << leaf->textScore << "\n";
//            cout << "SIZE  " << leaf->verifiedRecords->size() << "\n";
            
            // TEST ENTRIES OF LEAF
            if (leaf->verifiedRecords->size() == 0){
//                if(!textSearchOnNodeBinSearch(q.termSetR, q.termSetSizeR, &(iidx[leaf->page_id]), leaf)) //TBD CHANGE NOT WORKING
                if(!textSearchOnLeaf(q.termSetR, q.termSetSizeR, &(iidx[leaf->page_id]), leaf)) //TBD CHANGE NOT WORKING

                    continue;
            }
//            if(leaf->textScore == -1)
//                continue;
            
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
//            cout << "BEFOREINSERT\n";
//            cout << "SIZE " << leaf->verifiedRecords->size() << "\n";
//            for (auto i = 0; i < leaf->entry_num; i++)
//                if (inside(NUM_DIMENSIONS, entry[i]->data, inter))
//                    child.push_back(i);
//            cout << "SIZE " <<leaf->verifiedRecords->size();
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
//                    if( textProbe(q.termSetR, q.termSetSizeR, *R[entry[i]->id]) ){
#ifdef RetrieveResults
                        q.insert(PairRecord(&l, R[entry[i]->id], l.id, entry[i]->id, distance2N(l, *R[entry[i]->id])   ));
                
#else
                        q.insert();
#endif
//                    }
                 }
        }else{ //----- inner NODE
            RSTNonLeafNode<TreeDataP, double> *nonleaf = static_cast<RSTNonLeafNode<TreeDataP, double> *>(cnode);
            double mbr[4], inter[4];

            nonleaf->get_mbr(mbr);
            mbr[0] -= q.theta;
            mbr[1] += q.theta;
            mbr[2] -= q.theta;
            mbr[3] += q.theta;
            
//            if(nonleaf->child[i]->textScore == 0){ // node not visited as far
//                if( textSearchOnIRTree(q.termSetR, q.termSetSizeR, &iidx[entry[i]->page_id]) ){
//                    nonleaf->textScore = 1;
//                }else{
//                    nonleaf->textScore = -1;
//                }
//            }
//            cout << "TEST " << q.theta << " " << q.theta << endl;

            if (!intersect(NUM_DIMENSIONS, inter, qrange, mbr))
                continue; //no intersection we return
            if (!qualify(qrange, mbr, q.theta_sqr))
                continue; //judge false hit

            vector<RecordId> child;
            double **entry = nonleaf->entry_mbr;
            //generate candiates list
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
                            if((child->textScore == 1) || textProbeOnNode(q.termSetR, q.termSetSizeR, &iidx[child->page_id], child) )
                                S.push(child);
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


void thirdSetup(CDJResult& q, Relation& L, Relation& R, InvertedIndex* IFRT_R, RTree* RT_R){
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
        ContainmentRangeQueryOnIRTree(q, **l, R, RT_R,IFRT_R);
        l++;
    }


}



// 4. SETUP X4
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////



void SpatialTextualDistanceFilterOnNonLeafs(CDJResult& q, NonLeaf* l, NonLeaf* r, InvertedIndex* irtL, InvertedIndex* irtR, stack<NodesPair>* stack){
    //according to T. Brinkhoff et al SIGMOD 93
    //first we restrict the search space, one difference is that we now expand the MBR -/+dist
    //we expanf the two nonleaf MBR and calculate the intersection I
    //we can prove that if any child MBR do not intersect I, then it will not be in the result set
    vector<int> childL;
    vector<int> childR;
    double threshold = q.theta;
    double threshold_sqr = q.theta_sqr;
    double mbrR[4], mbrS[4], inter[4];
    l->get_mbr(mbrR);
    mbrR[0] -= threshold;
    mbrR[1] += threshold; //expand x-axis
    mbrR[2] -= threshold;
    mbrR[3] += threshold; //expand y-axis
    r->get_mbr(mbrS);
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
//          if (intersect(2, entryR[lVerifiedRec->at(i)], inter)){
//              childL.push_back(lVerifiedRec->at(i));
//          }
    
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

                RSTNode<TreeDataP, double> *nodeR = l->get_child(*iter1);
                RSTNode<TreeDataP, double> *nodeS = r->get_child(*iter);
                
                if(nodeR->textScore == 1 || textProbeOnNode(q.termSetL, q.termSetSizeL, &irtL[nodeR->page_id], nodeR))
                    if(nodeS->textScore == 1 || textProbeOnNode(q.termSetR, q.termSetSizeR, &irtR[nodeS->page_id], nodeS))
                        stack->push(make_pair(nodeR, nodeS));

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

                RSTNode<TreeDataP, double> *nodeR = l->get_child(*iter);
                RSTNode<TreeDataP, double> *nodeS = r->get_child(*iter2);
                
                if(nodeR->textScore == 1 || textProbeOnNode(q.termSetL, q.termSetSizeL, &irtL[nodeR->page_id], nodeR))
                    if(nodeS->textScore == 1 || textProbeOnNode(q.termSetR, q.termSetSizeR, &irtR[nodeS->page_id], nodeS))
                        stack->push(make_tuple(nodeR, nodeS));

                ++iter;
            }

            ++iter2;
        }
    }
}

void SpatialDistanceFilterONLeafAndNonLeaf(CDJResult& q, Leaf* l, NonLeaf* r, stack<NodesPair>* stack){
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
    
    
    for (int i = 0; i < r->entry_num; i++){
        if ( r->children[i]->textScore == 1 && intersect(2, entryR[i], mbrL))
            stack->push(make_tuple((Node*)l, (Node*)entryR[i]));
    }
}

void SpatialTextualJoinOnLeafs(CDJResult& q, Leaf* lLeaf, Leaf* rLeaf, Relation& L, Relation& R){
    TreeDataP<double> **lEntry = lLeaf->data;
    TreeDataP<double> **rEntry = rLeaf->data;
    vector<RecordId> lChild;
    vector<RecordId> rChild;
    
    double mbrR[4], mbrS[4], inter[4];
    unsigned int numResults = 0;
    Coordinates theta = q.theta;
    Coordinates theta_sqr = q.theta_sqr;
//    cout << "TEST " << q.theta << " " << q.theta << endl;
    
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
    if (!QualifyMinMinDistSqr(mbrR, mbrS, theta_sqr))
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
        if (inside(2, entryL[lVerifiedRec->at(i)]->data, inter))
            lChild.push_back(lVerifiedRec->at(i));
    
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
//    cout << endl << "SIZE " << lChild.size() << " " << rChild.size() << endl;
//    cout << "T " << theta << " t2 " << theta_sqr << endl;



    vector<RecordId>::iterator iter1 = lChild.begin(), iter2 = rChild.begin();
    vector<RecordId>::iterator end1 = lChild.end(), end2 = rChild.end();
    while ((iter1 < end1) && (iter2 < end2))
    {
//        cout << "IAMAHEREINNERLOOP\n";
        if (entryL[*iter1]->data[0] < (entryR[*iter2]->data[0] - q.theta)) //now asymmetric, we use two hand-coded internal loops
        {
            //internal loop 1
//            cout << "BEFOREQUALI1\n";
            vector<RecordId>::iterator iter = iter2;
            while ((iter < end2) && ((entryR[*iter]->data[0] - q.theta) < entryL[*iter1]->data[0])){
                int id1 = entryL[*iter1]->id;
                int id2 = entryR[*iter]->id;
//                cout << "BEFOREQUALI\n";
                if ((entryL[*iter1]->data[1] > (entryR[*iter]->data[1] + q.theta)) || ((entryR[*iter]->data[1] - q.theta) > entryL[*iter1]->data[1])){
                    ++iter;
                    continue;
                }
//                cout << "IAMAHERE";
//                cout << "BEFOREQUALI\n";
                if ((QualifySpatial(L[id1], R[id2], q.theta_sqr))) // && (QualifyTextual((*R)[id1], (*S)[id2], numKwR, kwR, numKwS, kwS)))
//                if(distance2N(*L[id1], *R[id2]) < q.theta)
                {
//                    cout << "pp1\n";
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
//            cout << "BEFOREQUALI2\n";
            vector<RecordId>::iterator iter = iter1;
            while ((iter < end1) && (entryL[*iter]->data[0] < (entryR[*iter2]->data[0] + q.theta)))
            {
                int id1 = entryL[*iter]->id;
                int id2 = entryR[*iter2]->id;
//                cout << "BEFOREQUALI\n";
                if (((entryR[*iter2]->data[1] - q.theta) > entryL[*iter]->data[1]) || (entryL[*iter]->data[1] > (entryR[*iter2]->data[1] + q.theta))){
                    ++iter;
                    continue;
                }
//                cout << "BEFOREQUALI\n";
                if ((QualifySpatial(L[id1], R[id2], q.theta_sqr))) // && (QualifyTextual((*R)[id1], (*S)[id2], numKwR, kwR,numKwS, kwS)))
                {
//                    cout << "pp2\n";
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




void IRTreeJoin(CDJResult& q, Relation& L, Relation& R, RTree* rt_L, RTree* rt_R, InvertedIndex* iix_L, InvertedIndex* iix_R){
    int LNumKeywords = q.termSetSizeL;
    int* LKeywords   = q.termSetL;
    
    int RNumKeywords = q.termSetSizeR;
    int* RKeywords   = q.termSetR;

    double threshold = q.theta;
    
    stack<NodesPair>* candidates = new stack<NodesPair>;
    stack<Node*>* lCandidates = new stack<Node*>;
    stack<Node*>* rCandidates = new stack<Node*>;

    
    NodesPair* pair;
    Node* l = rt_L->root;
    Node* r = rt_R->root;
    Node* lChild;
    Node* rChild;
    NonLeaf* lNonLeaf;
    NonLeaf* rNonLeaf;
    Leaf* lLeaf;
    Leaf* rLeaf;
    
    candidates->push(make_tuple(l, r));
    
    while (!candidates->empty()) {
        pair = &candidates->top();
        candidates->pop();
        l = get<0>(*pair);
        r = get<1>(*pair);

        
        if(r->is_leaf_node() && l->is_leaf_node()){ // if r leaf then l leaf -> compute pairs
            lLeaf = (Leaf*)l;
            rLeaf = (Leaf*)r;
            //textProbe on l
            if(lLeaf->textScore == 0 || lLeaf->verifiedRecords->size() == 0)
                textSearchOnLeaf(LKeywords, LNumKeywords, &iix_L[lLeaf->page_id], (Node*)lLeaf);
            if(lLeaf->textScore == -1)
                continue;
            //textProbe on l
            if(rLeaf->textScore == 0 || rLeaf->verifiedRecords->size() == 0)
                textSearchOnLeaf(RKeywords, RNumKeywords, &iix_R[rLeaf->page_id], (Node*)rLeaf);
            if(rLeaf->textScore == -1)
                continue;
            SpatialTextualJoinOnLeafs(q, lLeaf, rLeaf, L, R);
            

        }else if(l->is_leaf_node() && !r->is_leaf_node()){ // i only l leaf -> travers R
            cout << "R deph > L depht \n";
            lLeaf = (Leaf*)l;
            rNonLeaf = (NonLeaf*)r;
            //textProbe on l
            if(lLeaf->textScore == 0)
                textSearchOnLeaf(LKeywords, LNumKeywords, &iix_L[lLeaf->page_id], (Node*)lLeaf);
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
            SpatialDistanceFilterONLeafAndNonLeaf(q, lLeaf, rNonLeaf, candidates);
                
        }else if ( !l->is_leaf_node() && !r->is_leaf_node() ){ //travers to next level in both trees
            lNonLeaf = (NonLeaf*)l;
            rNonLeaf = (NonLeaf*)r;
//            cout << "IAMINLEAF\n";
            
            
            textVerificationOnNode(q.termSetL, q.termSetSizeL, &iix_L[lNonLeaf->page_id], lNonLeaf);
            textVerificationOnNode(q.termSetR, q.termSetSizeR, &iix_R[rNonLeaf->page_id], rNonLeaf);
            
            
            // Text Probe for all children of l and r
//            for(int i = 0; i < lNonLeaf->entry_num; i++){
//                lChild = lNonLeaf->children[i];
//                if(lChild->textScore == 0) // Containment Query on l if no did as far
//                    textVerificationOnNode(LKeywords, LNumKeywords, &iix_L[lChild->page_id], lChild);
////                    textSearchOnLeaf(LKeywords, LNumKeywords, &iix_L[lChild->page_id], (Node*)lLeaf);
//
//            }
//            for(int j = 0; j < rNonLeaf->entry_num; j++){
//                rChild = rNonLeaf->children[j];
//                if(rChild->textScore == 0) // Containment Query on l if no did as far
//                    textVerificationOnNode(RKeywords, RNumKeywords, &iix_R[rChild->page_id], rChild);
////                    textSearchOnLeaf(RKeywords, RNumKeywords, &iix_R[rChild->page_id], (Node*)rLeaf);
//
//            }
//            cout << "AFTERTEXTPROBE\n";

            SpatialTextualDistanceFilterOnNonLeafs(q, lNonLeaf, rNonLeaf, iix_L, iix_R, candidates);
        }else{
            cout << "I SHOULD NEVER HERE: IR-TREE JOIN \n\n\n";
        }
    }
   
}















// 5. SETUP X5
//////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////
///////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////////


void SpatialJoinOnLeafs2(CDJResult& q, Leaf& l, Leaf& r, Relation& L, Relation& R){
    //according to T. Brinkhoff et al SIGMOD 93
    //first we restrict the search space, one difference is that we now expand the MBR -/+dist
    //we expanf the two nonleaf MBR and calculate the intersection I
    //we can prove that if any child MBR do not intersect I, then it will not be in the result set
    double mbrR[4], mbrS[4], inter[4];
    unsigned int numResults = 0;

    double threshold = q.theta;
    double threshold_sqr = q.theta_sqr;
    Leaf* leafR = &l;
    Leaf* leafS = &r;
    vector<int> childR, childS;

    leafR->get_mbr(mbrR);
    mbrR[0] -= threshold;
    mbrR[1] += threshold; //expand x-axis
    mbrR[2] -= threshold;
    mbrR[3] += threshold; //expand y-axis
    leafS->get_mbr(mbrS);
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
    TreeDataP<double> **entryR = leafR->data;
    childR.clear();
    for (int i = 0; i < leafR->entry_num; i++){
        if ( entryR[i]->textScore > -1){
            if(inside(2, entryR[i]->data, inter) ){ // +Filter data, which do not satisfies text probe
                childR.push_back(i);
            }
        }

    }
    sort(childR.begin(), childR.end(), p_sort_dim_x_less(leafR));

    TreeDataP<double> **entryS = leafS->data;
    childS.clear();
    for (int i = 0; i < leafS->entry_num; i++){
        if ( entryS[i]->textScore > -1){
            if( inside(2, entryS[i]->data, inter) ) // +Filter data, which do not satisfies text probe
                childS.push_back(i);
        }
    }
    sort(childS.begin(), childS.end(), p_sort_dim_x_less(leafS));


    vector<int>::iterator iter1 = childR.begin(), iter2 = childS.begin();
    vector<int>::iterator end1 = childR.end(), end2 = childS.end();
    while ((iter1 < end1) && (iter2 < end2)){
        if (entryR[*iter1]->data[0] < (entryS[*iter2]->data[0] - threshold)){ //now asymmetric, we use two hand-coded internal loops


            //internal loop 1
            vector<int>::iterator iter = iter2;
            while ((iter < end2) && ((entryS[*iter]->data[0] - threshold) < entryR[*iter1]->data[0]))
            {
                int id1 = entryR[*iter1]->id;
                int id2 = entryS[*iter]->id;

                if(
                    ((entryR[*iter1]->data[1] > (entryS[*iter]->data[1] + threshold)) ||
                     ((entryS[*iter]->data[1] - threshold) > entryR[*iter1]->data[1]))){
                    ++iter;
                    continue;
                }
                if(entryR[*iter1]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetL, q.termSetSizeL, L[id1], &(entryR[*iter1]->textScore)))
                   if(entryS[*iter]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[id2], &(entryS[*iter]->textScore)))
//                       if(entryR[*iter1]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetL, q.termSetSizeL, L[id1], &(entryR[*iter1]->textScore)))
//                          if(entryS[*iter]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[id2], &(entryS[*iter]->textScore)))
                      if (distance2N(*L[id1], *R[id2]) < threshold){

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
            vector<int>::iterator iter = iter1;
            while ((iter < end1) && (entryR[*iter]->data[0] < (entryS[*iter2]->data[0] + threshold))){

                int id1 = entryR[*iter]->id;
                int id2 = entryS[*iter2]->id;

                if (((entryS[*iter2]->data[1] - threshold) > entryR[*iter]->data[1]) || (entryR[*iter]->data[1] > (entryS[*iter2]->data[1] + threshold))){
                    ++iter;
                    continue;
                }
                //TBD OPTIMIZATION
                if(entryR[*iter]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetL, q.termSetSizeL, L[id1], &(entryR[*iter]->textScore)))
                    if(entryS[*iter2]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[id2], &(entryS[*iter2]->textScore)))
//                        if(entryR[*iter]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetL, q.termSetSizeL, L[id1], &(entryR[*iter]->textScore)))
//                            if(entryS[*iter2]->textScore == 1 || textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[id2], &(entryS[*iter2]->textScore)))
                        if (distance2N(*L[id1], *R[id2]) < threshold){
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
//    for(int i = 0; i < leafR->entry_num; i++){
//        for(int j = 0; j < leafS->entry_num; j++){
//            int id1 = entryR[i]->id;
//            int id2 = entryS[j]->id;
//            float score = 10;
////            if(id1 == 195 || id1 ==403 || id1 ==1008 || id1 ==2040 || id1 ==2926 || id1 ==5279 || id1 ==7167 || id1 ==7580)
////            if(id1 == 195){
////                cout << "COUT I AM HERE: result (" << q.termSetL[0] << ")\n";
////                //cout <<
////                textProbeOnDataPoint(q.termSetL, q.termSetSizeL, L[id1], &score);// << " )MYSCORE";
////            }
//            if(textProbeOnDataPointBinSearch(q.termSetL, q.termSetSizeL, L[id1], &score)){
////                if(id2 == 99 || id2 == 686 || id2 == 800 || id2 == 968 || id2 == 1020 || id2 == 1338 || id2 == 2358 || id2 == 3207 || id2 == 3320 || id2 == 3756 || id2 == 3910 || id2 == 4712 || id2 == 4786 || id2 == 5636 || id2 == 5719 || id2 == 6846 || id2 == 6903 || id2 == 7358 || id2 == 7626 || id2 == 8427 || id2 == 8427 || id2 == 8459 || id2 == 8979 || id2 == 9359 || id2 == 10800)
//               if(textProbeOnDataPointBinSearch(q.termSetR, q.termSetSizeR, R[id2], &score))
//                   if(distance2N(*L[id1], *R[id2]) < threshold){
//                    q.insert(PairRecord( L[id1], R[id2], id1, id2, distance2N(*L[id1], *R[id2])));
//                   }
//            }
//        }
//    }


    
    
}


void SpatialDistanceFilterONLeafAndNonLeaf(CDJResult& q, Leaf& l, NonLeaf& r, stack<NodesPair>* stack){
    
}

void SpatialDistanceFilterOnNonLeaf(CDJResult& q, NonLeaf& l, NonLeaf& r, stack<NodesPair>* stack){
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
                stack->push(make_pair(nodeR, nodeS));
//                cout << "PAIR " << nodeR->page_id << " " << nodeS->page_id << endl;

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
//                cout << "PAIR " << nodeR->page_id << " " << nodeS->page_id << endl;
                stack->push(make_pair(nodeR, nodeS));

                ++iter;
            }

            ++iter2;
        }
    }
    
    
}


void RTreeJoinAndTextVerification(CDJResult& q, Relation& L, Relation& R, RTree* rt_L, RTree* rt_R){
    int LNumKeywords = q.termSetSizeL;
    int* LKeywords   = q.termSetL;
    int RNumKeywords = q.termSetSizeR;
    int* RKeywords   = q.termSetR;

    double threshold = q.theta;
    
    stack<NodesPair>* candidates = new stack<NodesPair>;

    
    NodesPair* pair;
    Node* l = rt_L->root;
    Node* r = rt_R->root;
    Node* lChild;
    Node* rChild;
    NonLeaf* lNonLeaf;
    NonLeaf* rNonLeaf;
    Leaf* lLeaf;
    Leaf* rLeaf;
    
    candidates->push(make_tuple(l, r));
    
    while (!candidates->empty()) {
        pair = &candidates->top();
        candidates->pop();
        l = get<0>(*pair);
        r = get<1>(*pair);

        
        if(r->is_leaf_node() && l->is_leaf_node()){ // if r leaf then l leaf -> compute pairs
            lLeaf = (Leaf*)l;
            rLeaf = (Leaf*)r;
            SpatialJoinOnLeafs2(q, *lLeaf, *rLeaf, L, R);
            

        }else if(l->is_leaf_node() && !r->is_leaf_node()){ // i only l leaf -> travers R
            cout << "R deph > L depht \n";
            SpatialDistanceFilterONLeafAndNonLeaf(q, lLeaf, rNonLeaf, candidates);
                
        }else if ( !l->is_leaf_node() && !r->is_leaf_node() ){ //travers to next level in both trees
            lNonLeaf = (NonLeaf*)l;
            rNonLeaf = (NonLeaf*)r;
            SpatialDistanceFilterOnNonLeaf(q, *lNonLeaf, *rNonLeaf, candidates);
        }else{
            cout << "I SHOULD NEVER HERE: IR-TREE JOIN \n\n\n";
        }
    }
   
    
    
    
    
    
    

}
