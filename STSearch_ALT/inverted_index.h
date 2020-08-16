#ifndef _INVERTED_INDEX_H_
#define _INVERTED_INDEX_H_

#include "def.h"
#include "relation.h"



//inverted list entry
struct InvertedListEntry
{
public:
	Record *rec;
	int position;
    float weight = 0;

	bool operator < (const InvertedListEntry &other) const
	{
		return (this->rec->id < other.rec->id);
	}


};


typedef vector<InvertedListEntry> InvertedList;


class InvertedIndexVector
{
public:
	hash<int, vector<InvertedListEntry> > lists;
	
	InvertedIndexVector();
	InvertedIndexVector(Relation *R);
	~InvertedIndexVector();

	bool ExistsKeyword(int kid);
	int EstimateCost(int kid, vector<int> &candidate);
	int EstimateCost(int kid, vector<InvertedListEntry> &candidate);

	void MoveOutNonleaf(int kid, vector<int> &candidate);
	void MakeIntersectionNonleaf(int kid, vector<int> &candidate, bool *invalid); //assume id is sorted in ascending order

	void MakeIntersection(int kid, vector<int> &candidate); //assume id is sorted in ascending order
	void MakeIntersection(int kid, vector<int> &candidate, int bound); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<int> &candidate); //make candidate sorted in ascending order
	void MoveOut(int kid, vector<int> &candidate, int bound); //make candidate sorted in ascending order

	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate); //assume id is sorted in ascending order
	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate, int bound); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate); //make candidate sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate, int bound); //make candidate sorted in ascending order
		
	void Print();
    static void Print(hash<int, vector<InvertedListEntry>> l);
	void ShowInfo();
};

////////////////////////////////////////////////////////////////////////////////
class InvertedIndexVectorBinSearch
{
public:
	hash<int, vector<InvertedListEntry> > lists;
	
	InvertedIndexVectorBinSearch();
	InvertedIndexVectorBinSearch(Relation *R);
	~InvertedIndexVectorBinSearch();

	bool ExistsKeyword(int kid);
	int EstimateCost(int kid, vector<int> &candidate);
	int EstimateCost(int kid, vector<InvertedListEntry> &candidate);

	void MakeIntersection(int kid, vector<int> &candidate); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<int> &candidate); //make candidate sorted in ascending order

	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate); //assume id is sorted in ascending order
	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate, int bound); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate); //make candidate sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate, int bound); //make candidate sorted in ascending order
		
	void Print();
	void ShowInfo();
};

////////////////////////////////////////////////////////////////////////////////
class InvertedIndexVectorSmart
{
public:
	hash<int, vector<InvertedListEntry> > lists;
	
	InvertedIndexVectorSmart();
	InvertedIndexVectorSmart(Relation *R);
	~InvertedIndexVectorSmart();

	bool ExistsKeyword(int kid);
	int EstimateCost(int kid, vector<int> &candidate);
	int EstimateCost(int kid, vector<InvertedListEntry> &candidate);

	void MakeIntersection(int kid, vector<int> &candidate); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<int> &candidate); //make candidate sorted in ascending order


	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate); //assume id is sorted in ascending order
	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate, int bound); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate); //make candidate sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate, int bound); //make candidate sorted in ascending order
		
	void Print();
	void ShowInfo();
};

////////////////////////////////////////////////////////////////////////////////
class InvertedIndexMap
{
public:
	hash<int, map<int, InvertedListEntry> > lists;
	
	InvertedIndexMap();
	InvertedIndexMap(Relation *R);
	~InvertedIndexMap();

	bool ExistsKeyword(int kid);
	int EstimateCost(int kid, vector<int> &candidate);
	int EstimateCost(int kid, vector<InvertedListEntry> &candidate);

	void MakeIntersection(int kid, vector<int> &candidate); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<int> &candidate); //make candidate sorted in ascending order

	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate); //assume id is sorted in ascending order
	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate, int bound); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate); //make candidate sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate, int bound); //make candidate sorted in ascending order
		
	void Print();
	void ShowInfo();
};

////////////////////////////////////////////////////////////////////////////////
class InvertedIndexHashMap
{
public:
	hash<int, hash<int, InvertedListEntry> > lists;
	
	InvertedIndexHashMap();
	InvertedIndexHashMap(Relation *R);
	~InvertedIndexHashMap();

	bool ExistsKeyword(int kid);
	int EstimateCost(int kid, vector<int> &candidate);
	int EstimateCost(int kid, vector<InvertedListEntry> &candidate);

	void MakeIntersection(int kid, vector<int> &candidate); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<int> &candidate); //make candidate sorted in ascending order

	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate); //assume id is sorted in ascending order
	void MakeIntersection(int kid, vector<InvertedListEntry> &candidate, int bound); //assume id is sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate); //make candidate sorted in ascending order
	void MoveOut(int kid, vector<InvertedListEntry> &candidate, int bound); //make candidate sorted in ascending order
		
	void Print();
	void ShowInfo();
};
#endif _INVERTED_INDEX_H_
