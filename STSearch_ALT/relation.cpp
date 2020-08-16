#include "relation.h"
#include "def.h"







Record::Record(){
	this->id = this->length = -1;
	this->keywords = NULL;
	this->locx = this->locy = -1;
    this->textScore = 0;
}


Record::Record(const Record& other){
	this->id = other.id;
	this->length = other.length;
	this->locx = other.locx;
	this->locy = other.locy;
    this->textScore = 0;
    this->simScore  = -1;

	if (other.keywords != NULL){
		this->keywords = new int[this->length];
		memcpy(this->keywords, other.keywords, this->length*sizeof(int));
	}
}


Record& Record::operator =(const Record &other){
	if (this != &other){
		this->id = other.id;
		this->length = other.length;
		this->locx = other.locx;
		this->locy = other.locy;

		if (this->keywords != NULL)
			delete[] this->keywords;

		if (other.keywords != NULL){
			this->keywords = new int[this->length];
			memcpy(this->keywords, other.keywords, this->length*sizeof(int));
		}
	}

	return *this;
}

bool Record::operator < (Record* rhs) const{
    Record l = *this;
    Record r = *rhs;
    
    cout << l.locx << " " << r.locx << "\n";
    return l.locx < r.locx;
}

bool Record::operator >= (Record* rhs) const{
    return !((*this) < rhs);
}
void Record::Print(char c)
{
    cout << c << ":" << this->id << ",\t\t (" << this->locx << "," << this->locy << "),\t\t length = " << this->length << ":";
    for (int kid = 0; kid < this->length; kid++)
        cout << " " << this->keywords[kid];
    cout << endl;
}
void Record::PrintWithWeights(char c)
{
    cout << c << ":" << this->id << ",\t\t (" << this->locx << "," << this->locy << "),\t\t length = " << this->length << ":";
    for (int kid = 0; kid < this->length; kid++)
        cout << " (" << this->weightedTerms[kid].first << ", " << this->weightedTerms[kid].second << ") ";
    cout << endl;
}


Record::~Record(){
//	if (this->keywords != NULL)
		delete[] this->keywords;
}

Relation::Relation()
{
    this->numRecords = 0;
    recs = nullptr;
}

//FOR FORAMT x,y \t n \t terms
#ifdef poly
Relation::Relation(const char *filename){
    
    int rid = 0, pos = 0;
    unsigned long long totalLength = 0;
    float minX, maxX, minY, maxY;
    float x, y;
    termSetSize n;
    termSet t;
    long keywords;
    long documents;
    numRecords = 0;

    // Shen: parsing file by reading every character is sometimes very slow
    this->maxRecordLength = 0;
    this->minRecordLength = numeric_limits<int>::infinity();

    minX = numeric_limits<double>::infinity();
    minY = numeric_limits<double>::infinity();
    maxX = -numeric_limits<double>::infinity();
    maxY = -numeric_limits<double>::infinity();

    ifstream inp(filename);

    string line;
    cout << line;
    getline(inp, line);
//    cout << line << endl;

    //get number documents
    documents = stoi(line.substr(0, line.find("\t")));
    line.erase(0, line.find("\t")+1);
    //get number terms
    keywords = stoi(line);
//    cout << "IAMHERE\n";
    this->recs = new Record[documents];
    this->numKeywords = keywords;
//    while(getline(inp, line)){
//    cout << "BEFORE WHILE";

//    cout << "NUM REC "<< this->numRecords << "\n";
//    cout << "NUM docs "<< documents << "\n";

    while(numRecords < documents){
//        cout << "IAMINTHELOOP";
        getline(inp, line);
//        cout << line << endl;

        x = stof( line.substr(0, line.find(",")) );
        line.erase(0, line.find(",")+1);

        y = stof( line.substr(0, line.find("\t")) );
        line.erase(0, line.find("\t")+1);

        n = stoi(line.substr(0, line.find("\t")));
        line.erase(0, line.find("\t")+1);

        int* tmp = new int[n];
        WeightedTerm* weights = new WeightedTerm[n];
        int i = 0;
        while(i < n){
            tmp[i] = stoi( line.substr(0, line.find(",")));
            line.erase( 0, line.find(",")+1 );
            i++;
        }
        this->recs[numRecords].id = numRecords;
        this->recs[numRecords].locx = x;
        this->recs[numRecords].locy = y;
        this->recs[numRecords].length = n;
        this->recs[numRecords].keywords = new int[n];
        this->recs[numRecords].keywords = tmp;
//        this->recs[numRecords].WeightedTerms = weights;
        this->recs[numRecords].weightedTerms = new WeightedTerm[n];



        if (this->maxRecordLength < this->recs[numRecords].length)
            this->maxRecordLength = this->recs[numRecords].length;
        if (this->minRecordLength > this->recs[numRecords].length)
            this->minRecordLength = this->recs[numRecords].length;

        minX = min(minX, x);
        minY = min(minY, y);
        maxX = max(maxX, x);
        maxY = max(maxY, y);
        totalLength += n;

        this->numRecords++;

    }
    inp.close();

    this->avgRecordLength = (int)(totalLength / this->numRecords);
}


#else
Relation::Relation(const char *filename){
	char c;
	stringstream ss;
	FILE *fp = 0;
	int rid = 0, pos = 0;
	unsigned long long totalLength = 0;
	bool locationRead = false;
	float minXtmp, maxXtmp, minYtmp, maxYtmp;
    float x, y;
    termSetSize n;
    termSet t;
    long keywords;
    long documents;

	// Shen: parsing file by reading every character is sometimes very slow
	this->maxRecordLength = 0;
	this->minRecordLength = numeric_limits<int>::infinity();

	minXtmp = numeric_limits<double>::infinity();
	minYtmp = numeric_limits<double>::infinity();
	maxXtmp = -numeric_limits<double>::infinity();
	maxYtmp = -numeric_limits<double>::infinity();

    ifstream inp(filename);
    string line;
    getline(inp, line);
    documents = stoi(line.substr(0, line.find("\t")));
    line.erase(0, line.find("\t")+1);
    keywords = stoi(line);
    this->recs = new Record[documents];
    this->numKeywords = keywords;
//    while(getline(inp, line)){

    while(numRecords < documents){
        getline(inp, line);
        n = stoi(line.substr(0, line.find("\t")));
        line.erase(0, line.find("\t")+1);
        x = stof( line.substr(0, line.find(",")) );
        line.erase(0, line.find(",")+1);
        y = stof( line.substr(0, line.find("\t")) );
        line.erase(0, line.find("\t")+1);

        int* tmp = new int[n];
        int i = 0;
        while(line.size() > 0){
            tmp[i] = stoi( line.substr(0, line.find(" ")));
            line.erase( 0, line.find(" ")+1 );
            i++;
        }
        
        this->recs[numRecords].id = numRecords;
        this->recs[numRecords].locx = x;
        this->recs[numRecords].locy = y;
        this->recs[numRecords].length = n;
        this->recs[numRecords].keywords = new int[n];
        this->recs[numRecords].keywords = tmp;
        this->recs[numRecords].weightedTerms = new WeightedTerm[n];


        if (this->maxRecordLength < this->recs[numRecords].length)
            this->maxRecordLength = this->recs[numRecords].length;
        if (this->minRecordLength > this->recs[numRecords].length)
            this->minRecordLength = this->recs[numRecords].length;

        minXtmp = min(minXtmp, this->recs[numRecords].locx);
        minYtmp = min(minYtmp, this->recs[numRecords].locy);
        maxXtmp = max(maxXtmp, this->recs[numRecords].locx);
        maxYtmp = max(maxYtmp, this->recs[numRecords].locy);
        totalLength += n;
//        totalLength += this->recs[numRecords].length;

        this->numRecords++;

    }
    inp.close();
    this->minX = minXtmp;
    this->maxX = maxXtmp;
    this->minY = minYtmp;
    this->maxY = maxYtmp;
//    cout << "[" << minX << "," << maxX <<"]X[" << minY << "," << maxY << "]\n";

//	for (int rid = 0; rid < numRecords; rid++){
//		// Normalization
//		this->recs[rid].locx = (this->recs[rid].locx-minX)/(maxX-minX+EPS);
//		this->recs[rid].locy = (this->recs[rid].locy-minY)/(maxY-minY+EPS);
//	}

	this->avgRecordLength = (float)(totalLength / (float)this->numRecords);
}
#endif
void Relation::clean(){
    for(int i = 0; i<numRecords; i++){
        recs[i].textScore = 0;
    }
}

Record* Relation::operator[](int rid){
	return &(this->recs[rid]);
}

static bool Comparator(Record** l, Record** r){
    return (*l)->locx < (*r)->locx;
};
static bool Comparator2(Record* l, Record* r){
    return (*l).locx < (*r).locx;
};

void Relation::sortByX(){
    sort(this->beginSubRelation(), this->endSubRelation(), Comparator2);
}


float Relation::computeDiameter(){
    return sqrt((minX-maxX)*(minX-maxX) + (minY-maxY)*(minY-maxY));
}




void Relation::Print(int num){
	Record *r;
    if(num == 0)
        num = this->numRecords;
	for (int rid = 0; rid < num; rid++)
	{
		r = &(this->recs[rid]);
		r->Print('R');
	}
    cout << "\n";
}

void Relation::PrintWithWeights(int num){
    Record *r;
    if(num == 0)
        num = this->numRecords;
    for (int rid = 0; rid < num; rid++)
    {
        r = &(this->recs[rid]);
        r->PrintWithWeights('R');
    }
    cout << "\n";
}

void Relation::PrintSubRelation(char c){
    Record *r;
    
    for (int rid = 0; rid < this->numRecords; rid++)
    {
        r = subTable[rid];
        r->Print(c);
    }
    cout << "\n";
};
void PrintWeightedRecord(WeightedRecord* rec){
    cout << "[id: " << rec->second->id << ",\t score: " << rec->first << "]\n";
}

void Relation::PrintTextSearch(int num){
    WeightedRecord* rec;
    cout << "Result textSearch num: " << this->numRecords << endl;
    for (int rid = 0; rid < this->numRecords; rid++){
        rec = &weightedRecords[rid];
        PrintWeightedRecord(rec);
    }
    cout << "\n";
};



RelationIterator Relation::begin(){
    return recs;
}
RelationIterator Relation::end(){
    return &recs[numRecords];
}
SubRelationIterator Relation::beginSubRelation(){
    return subTable;
}
SubRelationIterator Relation::endSubRelation(){
    return &subTable[numRecords];
    
};

Relation::~Relation(){
//    delete[] this->subTable;
	delete[] this->recs;
}



