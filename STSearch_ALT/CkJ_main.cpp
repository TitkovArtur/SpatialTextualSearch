//
//  BkQ_main.cpp
//  STSearch_ALT
//
//  Created by Artur Titkov on 17.06.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.
//

#include "def.h"
#include "relation.h"
#include "inverted_index.h"
#include "RStarTree.h"
#include "ResultSet.h"
#include <set>
#include "utils_boolean.cpp"

//#include "BkQ_Queries.cpp"

typedef pair<Node*, int> NodeValue;




bool *invalidL;
bool *invalidR;
bool *invalidRecR;
bool *invalidRecS;
Timer timerQuery;
Timer timerStep;
Timer timerClean;
Relation LSubRelation = Relation();
Relation RSubRelation = Relation();
char *dataL;
char *dataR;
Relation *L;
Relation *R;
RTree* rtL;
RTree* rtR;
InvertedIndex *ifL;
InvertedIndex *ifR;
InvertedIndex *irtL;
InvertedIndex *irtR;
BkQResult* q;
int runs = 1;
int rounds = 1;
float excludeTime = 0;
float runtimes[6] = {0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
int insertions[6] = {0,0,0,0,0,0};
int kwL[6] = {0,0,0,0,0,0};
int kwR[6] = {0,0,0,0,0,0};



Record* rec1;
Record* rec2;
int recid1;
int recid2;

float numContainmentQueryL = 0.0;
float numContainmentQueryR = 0.0;
float numInsertions = 0.0;
float defaultTheta;
float theta = 0;
const int defaultTerms = 3;
int k = 10;
float ks[5] = {1, 5, 10, 20, 50};






//IMPORTS
static bool Comparator(Record** l, Record** r){
    return (*l)->locx < (*r)->locx;
};
RStarTree<TreeDataP, double>* bulkload(int dim, int page_len, TreeDataP<double> **data, int data_num);
Relation* ContainmentQueryWithBinSearch(Relation& result, int numKeywords, int *keywords, InvertedIndex *iidx);
Relation* ContainmentQueryWithBinSearchOfRelation(Relation& result, int numKeywords, int *keywords, InvertedIndex *iidx);
size_t NestedLoopsDistanceJoin(BkQResult* q, Relation* R, Relation* S);
void PlaneSweepDistanceJoin(BkQResult* q, Relation* R,Relation* S);
void RangeQuery(BkQResult& q, const Record &l, Relation& R, const RStarTree<TreeDataP, double> *rtR);
void secSetup(BkQResult* q, Relation& L, Relation& R, RStarTree<TreeDataP, double>* rtR );
void thirdSetup(BkQResult& q, Relation& L, Relation& R, InvertedIndex* IFRT_R, RTree* RT_R);
void thirdSetupALT(BkQResult& q, Relation& L, InvertedIndex* IFRT_R, RTree* RT_R);
void IRTreeJoin(BkQResult& q, Relation& L, Relation& R, RTree* rt_L, RTree* rt_R, InvertedIndex* iix_L, InvertedIndex* iix_R);
Relation* ContainmentQuery(Relation& T, int numKeywords, int *keywords, InvertedIndex *iidx);
void RTreeJoinAndTextVerification(BkQResult& q, Relation& L, Relation& R, RTree* rt_L, RTree* rt_R);








inline void sortTermSets(Relation& L){
    auto i = L.begin();
    while(i != L.end()){
        sort(i->keywords, i->keywords + i->length);
        i++;
    }
};


static void Print(vector<InvertedListEntry> l){
    cout << "k X" << ":";
    for (int i = 0; i < l.size(); i++)
        cout << " <r" << l[i].rec->id << "," << l[i].position << ">";
    cout << endl;
}

void BuildIRTree(const char *datafile, Relation *R, RStarTree<TreeDataP, double> *&rt, InvertedIndex *&iidx, bool *&valid, bool build_new){

    rt = BuildRTree(datafile, R, build_new);

    iidx = new InvertedIndexVector[rt->num_nonleaf + rt->num_leaf];
    //valid = new bool[rt->num_nonleaf + rt->num_leaf];
    //for(int i = 0; i != rt->num_nonleaf + rt->num_leaf; ++i)
    //    valid[i] = false;

    rt->load_root();
    CreateIRTreeNode(rt->root, R, iidx);
    //iidx->Print();
    //iidx->ShowInfo();
    //    cout << datafile << ": InvertedIndex created." << endl;
}


struct com {
    bool operator()(const NodeValue& x,const NodeValue& y) const
    {
        // return "true" if "p1" is ordered
        // before "p2", for example:
        return x.second < y.second;
    }
};




bool benchmark_test();
void benchmark();
void benchmark_0();
void benchmark_1();
void benchmark_2();
void benchmark_3();
void benchmark_4();
void benchmark_5();

void PrintBenchmark(int terms, int k){
    cout << "===================================================================" << endl;
    cout << "================================BkQ================================" << endl;
    cout << "===================================================================" << endl << endl;
    cout << "Number terms " << terms << " k " << k << endl;
    cout << "AVG Contaiment Query L " << numContainmentQueryL/(float)rounds << endl;
    cout << "AVG Contaiment Query R " << numContainmentQueryR/(float)rounds << endl;
//    cout << "AVG num Insertions " << numInsertions/(float)rounds << endl;
    cout << " time SETUP 0 --> " << runtimes[0]/rounds << endl;
    cout << " time SETUP 1 --> " << runtimes[1]/rounds << endl;
    cout << " time SETUP 2 --> " << runtimes[2]/rounds << endl;
    cout << " time SETUP 3 --> " << runtimes[3]/rounds << endl;
    cout << " time SETUP 4 --> " << runtimes[4]/rounds << endl;
    cout << " time SETUP 5 --> " << runtimes[5]/rounds << endl << endl;
    
    cout << " insertions SETUP 0 --> " << insertions[0]/rounds << endl;
    cout << " insertions SETUP 1 --> " << insertions[1]/rounds << endl;
    cout << " insertions SETUP 2 --> " << insertions[2]/rounds << endl;
    cout << " insertions SETUP 3 --> " << insertions[3]/rounds << endl;
    cout << " insertions SETUP 4 --> " << insertions[4]/rounds << endl;
    cout << " insertions SETUP 5 --> " << insertions[5]/rounds << endl << endl;
    
    
    //CLEARING
    for(int i = 0; i < 7; i++){
        insertions[i] = 0;
        runtimes[i] = 0.0;
    }
    numContainmentQueryL = 0.0;
    numContainmentQueryR = 0.0;
    numInsertions = 0.0;
}
void findRandomTermSets(int numterms, int k){
    int pos;
    set<int> keywordsL; set<int>::iterator itrL;
    set<int> keywordsR; set<int>::iterator itrR;
    timerStep.stop();
    recid1 = 1;
    recid2 = 1;
    rec1 = (*L)[recid1];
    rec2 = (*R)[recid2];
    while(recid1 == recid2 && rec1->length < 5 && rec2->length < 5){
        while(rec2->length < 5){
            recid2 = rand() % R->numRecords;
            rec2 = (*R)[recid2];
        }
        while(rec1->length < 5){
            recid1 = rand() % L->numRecords;
            rec1 = (*L)[recid1];
        }
    }
    
    while (keywordsL.size() < numterms) {
        int pos = rand() % rec1->length;
        keywordsL.insert( rec1->keywords[pos]);
    }
    while (keywordsR.size() < numterms) {
        int pos = rand() % rec2->length;
        keywordsR.insert( rec2->keywords[pos]);
    }


    itrL = keywordsL.begin();
    itrR = keywordsR.begin();
//            cout << "KEYWORDS L\n";
    for(int j = 0; j < numterms; j++){
//                cout << *itrL << " ";
        kwL[j] = *itrL;
        itrL++;
    }
//            cout << "\nKEYWORDS R\n";
    for(int j = 0; j < numterms; j++){
//                cout << *itrR << " ";
        kwR[j] = *itrR;
        itrR++;
    }
//            cout << "\n";
    keywordsL.clear();
    keywordsR.clear();

    q = new BkQResult(k, kwL, numterms, kwR, numterms);

    
}


int main(int argc, char **argv){
    cout << "===================================================================" << endl;
    cout << "================================BkQ================================" << endl;
    cout << "===================================================================" << endl << endl;
        
    dataL = argv[1];
    dataR = argv[2];
    int modus = stoi(argv[3]); // if 1 var term size else distance

    
    #ifdef poly
    int numKwL = 2;
    int numKwR = 2;
    int kwL[2] = {15000000, 16419556};
    int kwR[2] = {16000000, 16419556};
    #else

    #endif
    
    int numKwL = 2;
    int numKwR = 2;
//    int kwL[3] = {326206, 326212, 326213};
//    int kwR[4] = {326205, 326211, 326214, 326215};
    int kwL[6] = {3, 153};
    int kwR[6] = {25, 71};
    
    
    
    cout << "READ FILES\n";
    L = new Relation(dataL);
    R = new Relation(dataR);
    
    sortTermSets(*L);
    sortTermSets(*R);
    cout << "CREATE IF\n";
    ifL = new InvertedIndex(L);
    ifR = new InvertedIndex(R);

    cout << "CREATE R-TREES\n";
    rtL = BuildRTree("", L, true);
    rtR = BuildRTree("", R, true);
    cout << "CREATE IR-TREES\n";
    BuildIRTree(dataL, L, rtL, irtL, invalidR, true);
    BuildIRTree(dataR, R, rtR, irtR, invalidR, true);
    

    BkQResult* q = new BkQResult(10, kwL, 2, kwR, 2);


    cout << "================INPUT================" << endl << endl;
    cout << "L (" << argv[1] << "): " << L->numRecords << " objects loaded, num keywords "<< L->numKeywords << " AVG terms " << L->avgRecordLength << " MAX terms " << L->maxRecordLength <<endl;
    cout << "R (" << argv[2] << "): " << R->numRecords << " objects loaded, num keywords "<< R->numKeywords << " AVG terms " << R->avgRecordLength << " MAX terms " << R->maxRecordLength <<endl;
    cout << "L in Range [" << L->minX << ", " << L->maxX << "] X [" << L->minY << ", " << L->maxY << "] diameter: " << L->computeDiameter()    << endl;
    cout << "R in Range [" << R->minX << ", " << R->maxX << "] X [" << R->minY << ", " << R->maxY << "] diameter: " << R->computeDiameter()    << endl;
    cout << "L R-Tree: num nodes " << rtL->num_nonleaf + rtL->num_leaf << "\t depth " << rtL->root->depth <<endl;
    cout << "R R-Tree: num nodes " << rtR->num_nonleaf + rtR->num_leaf << "\t depth " << rtR->root->depth <<endl;
    cout << "k = " << k << " " << endl << endl;
    if(modus == 0){
        cout << "vary term size, k = " << k << " ROUNDS " << rounds << endl << endl;
    }else{
        cout << "vary k, termsize  = " << 3 << " ROUNDS " << rounds << endl << endl;
    }
    
    
    //PICK RANDOM TERM SETS
    srand ((long)time(NULL));
//    srand (NULL);
    Record* rec1;
    Record* rec2;
    int recid1 = 1;
    int recid2 = 1;

    if(modus == 0){ // TERMS
        k = 10;
        for(int numKeywords = 1; numKeywords < 6; numKeywords++){
            for(int i = 0; i < rounds; i++){
                findRandomTermSets(numKeywords, k);
                if(!benchmark_test()){
                    q->clear();
                    i--;
                    continue;
                }
                benchmark();
            }
            PrintBenchmark(numKeywords, k);
        }
    }else{ // k
        const int numKeywords = 3;
        for(int d = 0; d < 5; d++){
            k = ks[d];
            for(int i = 0; i < rounds; i++){
                findRandomTermSets(numKeywords, k);
                if(!benchmark_test()){
                    q->clear();
                    i--;
                    continue;
                }
                benchmark();
            }
            PrintBenchmark(numKeywords, k);
        }
    }
    
    
    
    
    
    cout << "\n\n\n\n\n\n\n\n";
    return 0;
}






bool benchmark_test(){
    LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
    if(LSubRelation.numRecords == 0)
        return false;
//    LSubRelation.sortByX();
    RSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeR, q->termSetR, ifR);
    if(RSubRelation.numRecords == 0)
        return false;
//    RSubRelation.sortByX();
    return true;
    
    
    
//    PlaneSweepDistanceJoin(q, &LSubRelation, &RSubRelation);
////    cout << "1 " << LSubRelation.numRecords << endl;
////    cout << "2 " << RSubRelation.numRecords << endl;
////    cout << "3 " << q->numResult << endl << endl;
//    if(q->numResult == 0)
//        return false;
//    else
//        return true;
}
void benchmark(){
    
    //BASELINE
    q->clear();
    timerQuery.stop();
    LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
    RSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeR, q->termSetR, ifR);
    NestedLoopsDistanceJoin(q, &LSubRelation, &RSubRelation);
    runtimes[0] += timerQuery.stop();
    //CLEARING
    numContainmentQueryL += LSubRelation.numRecords;
    numContainmentQueryR += RSubRelation.numRecords;
    insertions[0] += q->numInsertions;
//    q->print('0');
    q->clear();
    timerQuery.stop();
    
    
    //SETUP 1
    LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
    LSubRelation.sortByX();
    RSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeR, q->termSetR, ifR);
    RSubRelation.sortByX();
    PlaneSweepDistanceJoin(q, &LSubRelation, &RSubRelation);
    runtimes[1] = timerQuery.stop();
    insertions[1] += q->numInsertions;
    //CLEARING
    q->print('1');
    q->clear();
    timerQuery.stop();
    
    rtR->clean();

    //SETUP 2
    LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
    secSetup(q, LSubRelation, *R, rtR);
    runtimes[2] = timerQuery.stop();
    insertions[2] += q->numInsertions;
    //CLEARING
//    q->print('2');
    q->clear();
    rtR->clean();
    timerQuery.stop();

    //SETUP 3
    LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
    LSubRelation.sortByX();
    thirdSetup(*q, LSubRelation, *R, irtR, rtR);
    runtimes[3] = timerQuery.stop();
    insertions[3] += q->numInsertions;
    //CLEARING
//    q->print('3');
    q->clear();
    rtR->clean();
    timerQuery.stop();
    
    
    //SETUP 4
    IRTreeJoin(*q, *L, *R, rtL, rtR, irtL, irtR);
    runtimes[4] = timerQuery.stop();
    insertions[4] += q->numInsertions;
    //CLEARING
//    q->print('4');
    q->clear();
    rtR->clean();
    rtL->clean();
    timerQuery.stop();
    
    //SETUP 5
    RTreeJoinAndTextVerification(*q, *L, *R, rtL, rtR);
    runtimes[5] = timerQuery.stop();
    insertions[5] += q->numInsertions;
    //CLEARING
//    q->print('5');
    q->clear();
    rtR->clean();
    rtL->clean();
}

void benchmark_0(){
    cout << "===================================================================" << endl;
    cout << "Test - Containment Query and Nested Loops Distance Join" << endl;
    cout << "===================================================================" << endl << endl;
//INIT
    timerQuery.stop();
    timerStep.stop();

//Containment Queries
    for(int i = 0; i < runs; i++)
        LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
    cout << "Containment Query with Bin Search L: retrieves " << LSubRelation.numRecords << "\tin " << timerStep.stop()/runs << "\tsec\n";
    for(int i = 0; i < runs; i++){
        RSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeR, q->termSetR, ifR);

    }
    cout << "Containment Query with Bin Search R: retrieves " << RSubRelation.numRecords << "\tin " << timerStep.stop()/runs << "\tsec\n";

//Nested Loops Distance Join
    for(int i = 0; i < runs; i++)
        NestedLoopsDistanceJoin(q, &LSubRelation, &RSubRelation);
    cout << "Nested Loops Distance Join: retrieves " << q->numResult/runs << "\tin " << timerStep.stop()/runs << "\tsec\n";
    cout << "\n-----Test in " << timerQuery.stop()/runs << "\tsec\n\n";
//    q->print('0');
}



void benchmark_1(){
    cout << "===================================================================" << endl;
    cout << "BkQ 1. SETUP: containment query and plane sweep" << endl;
    cout << "===================================================================" << endl << endl;;
//INIT
    timerQuery.stop();
    timerStep.stop();

    for(int i = 0; i < runs; i++){
        LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
        LSubRelation.sortByX();
    }
    cout << "Containment Query and Sort L: retrieves " << LSubRelation.numRecords << "\tin " << timerStep.stop()/runs << "\tsec\n";
    for(int i = 0; i < runs; i++){
        RSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeR, q->termSetR, ifR);
        RSubRelation.sortByX();
    }
    cout << "Containment Query and Sort R: retrieves " << RSubRelation.numRecords << "\tin " << timerStep.stop()/runs << "\tsec\n";

//Plane Sweep
    for(int i = 0; i < runs; i++){
        PlaneSweepDistanceJoin(q, &LSubRelation, &RSubRelation);
    }
//    q->print('1');
    cout << "Plane Sweep: retrieves " << q->numResult/runs <<  " in \t" << timerStep.stop()/runs << "\tsec\n";
    cout << "\n-----1. Setup in\t" << timerQuery.stop()/runs << "\tsec\n\n";
}

void benchmark_2(){
    cout << "===================================================================" << endl;
    cout << "BkQ 2. Setup - Containment Query and R-Tree" << endl;
    cout << "===================================================================" << endl << endl;;

//INIT
    float time = 0;
    excludeTime = 0;
    timerQuery.stop();
    timerStep.stop();

//Containment Query
    for(int i = 0; i < runs; i++){
        LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
    }
    cout << "Containment Query and Sort L: retrieves " << LSubRelation.numRecords << "\tin " << timerStep.stop()/runs << "\tsec\n";
    
    for(int i = 0; i < runs; i++){
        timerClean.stop();
        rtR->clean();
        excludeTime += timerClean.stop();
        secSetup(q, LSubRelation, *R, rtR);
//        q->print('2');
    }
    cout << " : retrieves " << q->numResult/runs <<  " in \t" << (timerStep.stop()-excludeTime)/runs << "\tsec\n";
    cout << "\n-----2. Setup in\t" << (timerQuery.stop()-excludeTime)/runs << "\tsec\n\n";
}
void benchmark_3(){
    cout << "===================================================================" << endl;
    cout << "BkQ 3. Setup - Containment Query and IR-Tree" << endl;
    cout << "===================================================================" << endl << endl;;
//INIT
    timerQuery.stop();
    timerStep.stop();
    excludeTime = 0;
    rtR->clean();
    
    
//Containment Query
    for(int i = 0; i < runs; i++){
        LSubRelation = *ContainmentQueryWithBinSearchOfRelation(*L, q->termSetSizeL, q->termSetL, ifL);
        LSubRelation.sortByX();
    }
    cout << "Containment Query and Sort L: retrieves " << LSubRelation.numRecords << "\tin " << timerStep.stop()/runs << "\tsec\n";
    
//IR-Tree text verification and Range Query
    for(int i = 0; i < runs; i++){
        timerClean.stop();
        rtR->clean();
        excludeTime += timerClean.stop();
//        q->clear(true, 0.1);
        thirdSetup(*q, LSubRelation, *R, irtR, rtR);
//        if(q->numResult <= 128){
//            cout << "SIZE " <<  q->numResult << "\n";
//            q->print('3');
//
//        }
        
    }
    cout << " Containment Query and IR-Tree: retrieves " << q->numResult/runs <<  " in \t" << (timerStep.stop()-excludeTime)/runs << "\tsec\n";
    cout << "\n-----3. Setup in\t" << (timerQuery.stop()-excludeTime)/runs << "\tsec\n\n";
}



void benchmark_4(){
    cout << "===================================================================" << endl;
    cout << "BkQ 4. Setup - IR-Tree Join" << endl;
    cout << "===================================================================" << endl << endl;
//INIT
    timerQuery.stop();
    timerStep.stop();
    excludeTime = 0;
        
// IR-Tree Join
    for(int i = 0; i < runs; i++){
        timerClean.stop();
        rtL->clean();
        rtR->clean();
        excludeTime += timerClean.stop();
        IRTreeJoin(*q, *L, *R, rtL, rtR, irtL, irtR);
    }
    cout << " Containment Query and IR-Tree: retrieves " << q->numResult/runs <<  " in \t" << (timerStep.stop()-excludeTime)/runs << "\tsec\n";
    cout << "\n-----4. Setup in\t" << (timerQuery.stop()-excludeTime)/runs << "\tsec\n\n";
}
void benchmark_5(){
    cout << "===================================================================" << endl;
    cout << "BkQ 5. Setup - R-Tree Join and textVerification" << endl;
    cout << "===================================================================" << endl << endl;
    timerQuery.stop();
    timerStep.stop();
    excludeTime = 0;



//R-Tree Join and text verification
    for(int i = 0; i < runs; i++){
        timerClean.stop();
        rtL->clean();
        rtR->clean();
        excludeTime += timerClean.stop();
        RTreeJoinAndTextVerification(*q, *L, *R, rtL, rtR);
    }
    cout << "R-Tree Join and text probe: retrieves " << q->numResult/runs <<  " in \t" << (timerStep.stop()-excludeTime)/runs << "\tsec\n";
    cout << "\n-----5. Setup in\t" << (timerQuery.stop()-excludeTime)/runs << "\tsec\n\n";
}
