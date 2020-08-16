#include "def.h"
#include "relation.h"
#include "inverted_index.h"
#include "RStarTree.h"
#include "timing.h"
#include "ResultSet.h"
//#include "BRQ_Queries.cpp"





Relation *L;
Relation *R;
Relation *S;
InvertedIndex *IFL;
InvertedIndex *IFR;
InvertedIndex *IIL;
InvertedIndex *IIR;
InvertedIndex *iidxR;
InvertedIndex *iidxS;
bool *invalidL;
bool *invalidR;
bool *invalidRecR;
bool *invalidRecS;

bool SELF_JOIN = true;
bool BUILD_NEW = false;

const char *dataL = "/Users/Artur/Desktop/tripadvisor_restaurants_london_replayced.txt";
const char *dataR = "/Users/Artur/Desktop/tripadvisor_restaurants_london_replayced.txt";
//const char *dataR = fileR.c_str();
//const char *dataS = fileS.c_str();
RStarTree<TreeDataP, double> *rtL, *rtR;
Timing tim;
unsigned long long total_duration = 0;



//INPORTS
static bool Comparator(Record** l, Record** r){
    return (*l)->locx < (*r)->locx;
};
RStarTree<TreeDataP, double>* bulkload(int dim, int page_len, TreeDataP<double> **data, int data_num);
Relation* SetContainmentQueryWithResult(Relation& result, int numKeywords, int *keywords, InvertedIndex *iidx);
size_t NestedLoopsDistanceJoin(BRQResult* q, Relation* R, Relation* S);
void PlaneSweepDistanceJoin(BRQResult* q, Relation* R,Relation* S);
void RangeQuery(BRQResult& q, const Record &l, Relation& R, const RStarTree<TreeDataP, double> *rtR);
void secSetup(BRQResult* q, Relation& L, Relation& R, RStarTree<TreeDataP, double>* rtR );
void thirdSetup(BRQResult& q, Relation& L, InvertedIndex* IFRT_R, RTree* RT_R);
void thirdSetupALT(BRQResult& q, Relation& L, InvertedIndex* IFRT_R, RTree* RT_R);
void IRTreeJoin(BRQResult& q, Relation& L, Relation& R, RTree* rt_L, RTree* rt_R, InvertedIndex* iix_L, InvertedIndex* iix_R);
void RTreeJoinAndTextVerification(BRQResult& q, Relation& L, Relation& R, RTree* rt_L, RTree* rt_R);





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
	//	valid[i] = false;

	rt->load_root();
	CreateIRTreeNode(rt->root, R, iidx);
    //iidx->Print();
    //iidx->ShowInfo();
//	cout << datafile << ": InvertedIndex created." << endl;
}




int main(int argc, char **argv) {
    cout << "===================================================================" << endl;
    cout << "================================BRQ================================" << endl;
    cout << "===================================================================" << endl << endl;
    
    
    int numKwL = 2;
    int numKwR = 2;
    //double theta = 0.01;
    int kwL[2] = {3, 153};
    int kwR[2] = {25, 71};
    double theta = 0.1;  //stof(argv[3]);
    BRQResult* q = new BRQResult(true, theta, kwL, 2, kwR, 2);

	L = new Relation(dataL);
	R = new Relation(dataR);
    sortTermSets(*L);
    sortTermSets(*R);
    IFL = new InvertedIndex(L);
    IFR = new InvertedIndex(R);
    Relation LSubRelation = Relation();
    Relation RSubRelation = Relation();

    
    //IFL->Print();
    
    //cout << "\n\n" << IFL->lists.at(33).size() << "\n\n";
    
    
    //IFL->Print();

    Timer timerQuery;
    Timer timerStep;
    
    cout << "================INPUT================" << endl << endl;
    cout << "------" << endl;
    cout << "L (" << argv[1] << "): " << L->numRecords << " objects loaded" << endl;
    cout << "R (" << argv[2] << "): " << R->numRecords << " objects loaded" << endl;
    cout << "theta = " << theta << " " << endl << endl;



    cout << "===================================================================" << endl;
    cout << "Test - Containment Query and Nested Loops Distance Join" << endl;
    cout << "===================================================================" << endl << endl;
//INIT
    timerQuery.stop();
    timerStep.stop();

//Containment Queries
    LSubRelation = *SetContainmentQueryWithResult(*L, numKwL, kwL, IFL);
    cout << "Containment Query L: retrieves " << LSubRelation.numRecords << "\tin " << timerStep.stop() << "\tsec\n";
    RSubRelation = *SetContainmentQueryWithResult(*R, numKwR, kwR, IFR);
    cout << "Containment Query R: retrieves " << RSubRelation.numRecords << "\tin " << timerStep.stop() << "\tsec\n";

//Nested Loops Distance Join
    NestedLoopsDistanceJoin(q, &LSubRelation, &RSubRelation);
    cout << "Nested Loops Distance Join: retrieves" << q->numResult << "\tin " << timerStep.stop() << "\tsec\n";
    cout << "\n-----Test in " << timerQuery.stop() << "\tsec\n\n";
//    LSubRelation.PrintSubRelation('L');
//    RSubRelation.PrintSubRelation('R');


// Test & Clear
//    RSubRelation.Print('l');
//    SSubRelation.Print('r');

    //q->print('0');







    //x1
    cout << "===================================================================" << endl;
    cout << "BRQ 1. SETUP: containment query and plane sweep" << endl;
    cout << "===================================================================" << endl << endl;;
//INIT
    q->clear(true, theta);
    timerQuery.stop();
    timerStep.stop();



//Containment Queries
    LSubRelation = *SetContainmentQueryWithResult(*L, numKwL, kwL, IFL);
    cout << "Containment Query L: retrieves" << LSubRelation.numRecords << "\tin " << timerStep.stop() << "\tsec\n";
    RSubRelation = *SetContainmentQueryWithResult(*R, numKwR, kwR, IFR);
    cout << "Containment Query R: retrieves" << RSubRelation.numRecords << "\tin " << timerStep.stop() << "\tsec\n";
    LSubRelation = *SetContainmentQueryWithBinSearch(*L, numKwL, kwL, IFL);
    cout << "Containment Query with bin L: retrieves" << LSubRelation.numRecords << "\tin " << timerStep.stop() << "\tsec\n";

    
    
    
    LSubRelation.sortByX();
    cout << "Sorting L on X: " <<  "\tin " << timerStep.stop() << "\tsec\n";
    RSubRelation.sortByX();
    cout << "Sorting R on X: " <<  "\tin " << timerStep.stop() << "\tsec\n";
//Plane Sweep
    PlaneSweepDistanceJoin(q, &LSubRelation, &RSubRelation);
    cout << "Plane Sweep: retrieves " << q->numResult <<  " in \t" << timerStep.stop() << "\tsec\n";
    cout << "\n-----1. Setup in\t" << timerQuery.stop() << "\tsec\n\n";


// Test
    //q->print('1');









    //x2
    cout << "===================================================================" << endl;
    cout << "BRQ 2. Setup - Containment Query and R-Tree" << endl;
    cout << "===================================================================" << endl << endl;;

//Build R-Tree
//        TreeDataP<Coordinates> **dataR; // Whats that
//        RStarTree<TreeDataP, Coordinates> *rtR;
//        dataR = new TreeDataP<Coordinates>*[R->numRecords];
//
//
//
//        for (const Record &r : *R){
//            dataR[r.id] = new TreeDataP<Coordinates>(NUM_DIMENSIONS);
//            dataR[r.id]->id = r.id;
//            dataR[r.id]->data[0] = r.locx;
//            dataR[r.id]->data[1] = r.locy;
//        }
//        rtR = bulkload(NUM_DIMENSIONS, PAGE_SIZE, dataR, R->numRecords);
//
//        RStarTree<TreeDataP, double> *rtR, *rtS;

    RTree *RTR = BuildRTree("", R, true);
//INIT
    q->clear(true, theta);
    timerQuery.stop();
    timerStep.stop();



//Containment Query
    LSubRelation = *SetContainmentQueryWithResult(*L, numKwL, kwL, IFL);
    cout << "Containment Query L: retrieves" << LSubRelation.numRecords << "\tin " << timerStep.stop() << "\tsec\n";
//
    secSetup(q, LSubRelation, *R, RTR);
    cout << " : retrieves " << q->numResult <<  " in \t" << timerStep.stop() << "\tsec\n";
    cout << "\n-----2. Setup in\t" << timerQuery.stop() << "\tsec\n\n";



//TEST
//    q->print('s');



//X4A
    cout << "===================================================================" << endl;
    cout << "BRQ 3. Setup - Containment Query and IR-Tree" << endl;
    cout << "===================================================================" << endl << endl;;



//BUILD IR-Tree
    //BuildIRTreeALT(dataR, R, RTR, IF_RT, invalidR, true);
    //R = new Relation(dataR);
    Relation* R_ALT;
    R_ALT = new Relation(dataR);
    InvertedIndex* IF_RT_ALT;
    RTree *RTR_ALT = BuildRTree("", R_ALT, true);
    BuildIRTree(dataR, R_ALT, RTR_ALT, IF_RT_ALT, invalidR, true);

//    IF_RT_ALT[0].Print();


//INIT
    q->clear(true, theta);
    timerQuery.stop();
    timerStep.stop();

//Containment Query
    LSubRelation = *SetContainmentQueryWithResult(*L, numKwL, kwL, IFL);
    cout << "Containment Query L: retrieves " << LSubRelation.numRecords << "\tin " << timerStep.stop() << "\tsec\n";
//IR-Tree text verification and Range Query
    thirdSetup(*q, LSubRelation, IF_RT_ALT, RTR_ALT);
    cout << " : retrieves " << q->numResult <<  " in \t" << timerStep.stop() << "\tsec\n";
    cout << "\n-----3. Setup in\t" << timerQuery.stop() << "\tsec\n\n";

// Test
//    q->print('1');



        
        
    
    cout << "===================================================================" << endl;
    cout << "BRQ 4. Setup - IR-Tree Join" << endl;
    cout << "===================================================================" << endl << endl;

        
//BUILD
    Relation* L_4 = new Relation(dataL);;
    RTree* RT_L4 = BuildRTree("", L_4, true);
    InvertedIndex* IF_L4 = new InvertedIndex(L_4);
    InvertedIndex* IF_RT_L4;
    BuildIRTree(dataL, L_4, RT_L4, IF_RT_L4, invalidL, true);

    Relation* R_4 = new Relation(dataR);;
    RTree* RT_R4 = BuildRTree("", R_4, true);
    InvertedIndex* IF_R4 = new InvertedIndex(R_4);
    InvertedIndex* IF_RT_R4;
    BuildIRTree(dataR, R_4, RT_R4, IF_RT_R4, invalidR, true);
    

//INIT
    q->clear(true, theta);
    timerQuery.stop();
    timerStep.stop();
    
    
// IR-Tree Join
    IRTreeJoin(*q, *L_4, *R_4, RT_L4, RT_R4, IF_RT_L4, IF_RT_R4);
    cout << "IR-Tree join: retrieves " << q->numResult <<  " in \t" << timerStep.stop() << "\tsec\n";
    cout << "\n-----4. Setup in\t" << timerQuery.stop() << "\tsec\n\n";
    
    
    
        
    
//TEST
//    q->print('M');


    cout << "===================================================================" << endl;
    cout << "BRQ 5. Setup - R-Tree Join and textVerificationn" << endl;
    cout << "===================================================================" << endl << endl;

//BUILD
//    Relation* L_5 = new Relation(dataL);;
    RTree* RT_L5 = BuildRTree("", L, true);
//    Relation* R_5 = new Relation(dataR);;
    RTree* RT_R5 = BuildRTree("", R, true);
    
//INIT
    q->clear(true, theta);
    timerQuery.stop();
    timerStep.stop();

//R-Tree Join and text verification
    RTreeJoinAndTextVerification(*q, *L, *R, RT_L5, RT_R5);
    
    cout << "R-Tree Join and text probe: retrieves " << q->numResult <<  " in \t" << timerStep.stop() << "\tsec\n";
    cout << "\n-----5. Setup in\t" << timerQuery.stop() << "\tsec\n\n";
    
    q->print('d');
    
    

//    L->Print('s');
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
    
    
    

    
    
    
    
    
    
    
    
    
    
    
    
    
	//kwR[0] = SELF_JOIN ? 4489 : 17440;	//4489;	//17440;
	//kwR[1] = SELF_JOIN ? 4494 : 17447;	//4494;	//17447;
	//kwS[0] = 4493;
	//kwS[1] = 4495;

	//ChooseKeywords();
	//Test();
	//IRTreeMethod();

//	delete iidxR;
//	delete iidxS;
//	//delete IIR;
//	//delete IIS;
//	delete rtR;
//	delete rtS;
//	delete R;
//	delete S;
//	delete kwR;
//	delete kwS;
    //delete RSubRelation;
    std::cout << "\n\n\n\n\n\n";
	return 0;
}
