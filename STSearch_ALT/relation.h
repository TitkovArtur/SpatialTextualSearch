#ifndef _RELATION_H_
#define _RELATION_H_

#include "def.h"

typedef Record* RelationIterator;
typedef Record** SubRelationIterator;
typedef pair<int, float> WeightedTerm;
typedef pair<float, Record*> WeightedRecord;


class Record{
public:
	int id;
	int length;
	int *keywords;
	float locx;
	float locy;
    Score textScore;
    Score simScore;
    vector<WeightedTerm>* TMPweightedTerms = new vector<WeightedTerm>;
    WeightedTerm* weightedTerms;

	Record();
	~Record();
	//copy constructors
	Record(const Record& other);
	Record & operator = (const Record &other);
    bool operator < (Record* rhs) const;
    bool operator >= (Record* rhs) const;

	void Print(char c);
    void PrintWithWeights(char c);
};

class Relation
{
public:
	int numRecords;
	int numKeywords;
	float avgRecordLength;
	int maxRecordLength;
	int minRecordLength;
    float minX;
    float maxX;
    float minY;
    float maxY;

	Record *recs;
    Record** subTable;
    WeightedRecord* weightedRecords;
public:
    Relation();
    Relation(int k);
	Relation(const char *filename);
	~Relation();
    void clean();
    
    
    void sortByX();
    float computeDiameter();
    
    
    RelationIterator begin();
    RelationIterator end();
    SubRelationIterator beginSubRelation();
    SubRelationIterator endSubRelation();
    
    
public:
	Record* operator[](int rid);
	void Print(int num);
    void PrintSubRelation(char c);
    void PrintWithWeights(int num);
    void PrintTextSearch(int num);
};

#endif _RELATION_H_
