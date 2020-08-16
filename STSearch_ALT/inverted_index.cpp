#include "inverted_index.h"



// Imports
bool CompareInvertedListEntries(InvertedListEntry lhs, InvertedListEntry rhs) 
{
	#ifdef ASCENDING_KEYWORDS_ORDER
	return (lhs.rec->id < rhs.rec->id);
	#else
	return (lhs.rec->id > rhs.rec->id);
	#endif
}


InvertedIndexVector::InvertedIndexVector()
{
}

InvertedIndexVector::InvertedIndexVector(Relation *R)
{
	InvertedListEntry e;


	#ifdef ASCENDING_KEYWORDS_ORDER
	for (int rid = 0; rid < R->numRecords; rid++)
	#else
	for (int rid = R->numRecords-1; rid >= 0; rid--)
	#endif
	{
		Record *r = (*R)[rid];

		for (int k = 0; k < r->length; k++)
		{
			e.rec = r;
			e.position = k;
			this->lists[r->keywords[k]].push_back(e);
		}
	}
}

InvertedIndexVector::~InvertedIndexVector()
{
}

bool InvertedIndexVector::ExistsKeyword(int kid)
{
	return (lists.count(kid) > 0);
}

int InvertedIndexVector::EstimateCost(int kid, vector<int> &candidate)
{
	return (int)(lists[kid].size()+candidate.size());
}

int InvertedIndexVector::EstimateCost(int kid, InvertedList &candidate)
{
	return (int)(lists[kid].size()+candidate.size());
}

//assume id is sorted in ascending order
void InvertedIndexVector::MakeIntersection(int kid, vector<int> &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	vector<int> temp;

	InvertedList::iterator iter_x = kvec.begin(), iter_x_end = kvec.end();
	vector<int>::iterator iter_y = candidate.begin(), iter_y_end = candidate.end();

	while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if (iter_x->rec->id < *iter_y) 
			iter_x++;
		else if (iter_x->rec->id > *iter_y)
			iter_y++;
		else
		{
			iter_x++;
			temp.push_back(*iter_y++);
		}
		#else
		if (iter_x->rec->id > *iter_y) 
			iter_x++;
		else if (iter_x->rec->id < *iter_y)
			iter_y++;
		else
		{
			iter_x++;
			temp.push_back(*iter_y++);
		}
		#endif
	}

	//make results go back to candidate
	candidate.swap(temp);		
}


void InvertedIndexVector::MoveOutNonleaf(int kid, vector<int> &candidate)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for(InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
	{
		candidate.push_back(iter->position);
	}
	sort(candidate.begin(), candidate.end());
}


void InvertedIndexVector::MakeIntersectionNonleaf(int kid, vector<int> &candidate, bool *invalid)
{
	if(lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	InvertedList &kvec = lists[kid];
	vector<int> pos;

	for(InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
	{
		pos.push_back(iter->position);
	}
	sort(pos.begin(), pos.end());

	vector<int> temp;

	vector<int>::iterator iter_x = pos.begin(), iter_x_end = pos.end();
	vector<int>::iterator iter_y = candidate.begin(), iter_y_end = candidate.end();

	while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
	{
		if (*iter_x < *iter_y)
		{
			invalid[*iter_x] = true;
			iter_x++;
		}
		else if (*iter_x > *iter_y)
			iter_y++;
		else
		{
			iter_x++;
			//cout << *iter_y << " ";
			temp.push_back(*iter_y++);
		}
	}
	//cout << temp.size() << " ";
	candidate.swap(temp);
	sort(candidate.begin(), candidate.end());
	//cout << RList.size() << endl;
}


//assume id is sorted in ascending order
void InvertedIndexVector::MakeIntersection(int kid, vector<int> &candidate, int bound)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	vector<int> temp;

	InvertedList::iterator iter_x = kvec.begin(), iter_x_end = kvec.end();
	vector<int>::iterator iter_y = candidate.begin(), iter_y_end = candidate.end();

	while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if (iter_x->rec->id < *iter_y) 
			iter_x++;
		else if (iter_x->rec->id > *iter_y) 
			iter_y++;
		else
		{
			if (iter_x->rec->length-iter_x->position >= bound)
				temp.push_back(*iter_y);
			iter_x++;
			iter_y++;
		}
		#else
		if (iter_x->rec->id > *iter_y) 
			iter_x++;
		else if (iter_x->rec->id < *iter_y) 
			iter_y++;
		else
		{
			if (iter_x->rec->length-iter_x->position >= bound)
				temp.push_back(*iter_y);
			iter_x++;
			iter_y++;
		}
		#endif
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//make candidate sorted in ascending order
void InvertedIndexVector::MoveOut(int kid, vector<int> &candidate)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
		candidate.push_back(iter->rec->id);
}

void InvertedIndexVector::MoveOut(int kid, vector<int> &candidate, int bound)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
	{
		if (iter->rec->length-iter->position >= bound)
			candidate.push_back(iter->rec->id);
	}
}

//assume id is sorted in ascending order
void InvertedIndexVector::MakeIntersection(int kid, InvertedList &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	InvertedList temp;

	InvertedList::iterator iter_x = kvec.begin(), iter_x_end = kvec.end();
	InvertedList::iterator iter_y = candidate.begin(), iter_y_end = candidate.end();

	#ifdef ASCENDING_KEYWORDS_ORDER
	while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
	{
		if (iter_x->rec->id < iter_y->rec->id) 
			iter_x++;
		else if (iter_x->rec->id > iter_y->rec->id) 
			iter_y++;
		else
		{
			temp.push_back(*iter_x++);
			iter_y++;
		}
	}
	#else
	while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
	{
		if (iter_x->rec->id > iter_y->rec->id) 
			iter_x++;
		else if (iter_x->rec->id < iter_y->rec->id) 
			iter_y++;
		else
		{
			temp.push_back(*iter_x++);
			iter_y++;
		}
	}
	#endif

	//make results go back to candidate
	candidate.swap(temp);
}

//assume id is sorted in ascending order
void InvertedIndexVector::MakeIntersection(int kid, InvertedList &candidate, int bound)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	InvertedList temp;

	InvertedList::iterator iter_x = kvec.begin(), iter_x_end = kvec.end();
	InvertedList::iterator iter_y = candidate.begin(), iter_y_end = candidate.end();

	while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if (iter_x->rec->id < iter_y->rec->id) 
			iter_x++;
		else if (iter_x->rec->id > iter_y->rec->id) 
			iter_y++;
		else
		{
			if (iter_x->rec->length-iter_x->position >= bound)
				temp.push_back(*iter_x);
			iter_x++;
			iter_y++;
		}
		#else
		if (iter_x->rec->id > iter_y->rec->id) 
			iter_x++;
		else if (iter_x->rec->id < iter_y->rec->id) 
			iter_y++;
		else
		{
			if (iter_x->position+1 >= bound)
				temp.push_back(*iter_x);
			iter_x++;
			iter_y++;
		}
		#endif
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//make candidate sorted in ascending order
void InvertedIndexVector::MoveOut(int kid, InvertedList &candidate)
{
	candidate.clear();

	InvertedList kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
	{
		candidate.push_back(*iter);
	}
}

//make candidate sorted in ascending order
void InvertedIndexVector::MoveOut(int kid, InvertedList &candidate, int bound)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if (iter->rec->length-iter->position >= bound)
		#else
		if (iter->position+1 >= bound)
		#endif
			candidate.push_back(*iter);
	}
}

void InvertedIndexVector::Print()
{
    
	for (hash<int, InvertedList >::iterator iter = this->lists.begin(); iter != this->lists.end(); ++iter)
	{
		cout << "k" << iter->first << ":";
		for (InvertedList::iterator iterL = iter->second.begin(); iterL != iter->second.end(); ++iterL)
			cout << " <r" << iterL->rec->id << ", " << iterL->position << ", " << iterL->weight << ">";
		cout << endl;
	}
}





void InvertedIndexVector::ShowInfo()
{
	printf("Using vector implementation for intersection.\n");
}

////////////////////////////////////////////////////////////////////////////////
InvertedIndexVectorBinSearch::InvertedIndexVectorBinSearch()
{
}

InvertedIndexVectorBinSearch::InvertedIndexVectorBinSearch(Relation *R)
{
	InvertedListEntry e;


	#ifdef ASCENDING_KEYWORDS_ORDER
	for (int rid = 0; rid < R->numRecords; rid++)
	#else
	for (int rid = R->numRecords-1; rid >= 0; rid--)
	#endif
	{
		Record *r = (*R)[rid];

		for (int k = 0; k < r->length; k++)
		{
			e.rec = r;
			e.position = -1; //CHANGD
			this->lists[r->keywords[k]].push_back(e);
		}
	}
}

InvertedIndexVectorBinSearch::~InvertedIndexVectorBinSearch()
{
}

bool InvertedIndexVectorBinSearch::ExistsKeyword(int kid)
{
	return (lists.count(kid) > 0);
}

int InvertedIndexVectorBinSearch::EstimateCost(int kid, vector<int> &candidate)
{
	return (int)(log((double)lists[kid].size())/log(2.0)*candidate.size());
}

int InvertedIndexVectorBinSearch::EstimateCost(int kid, InvertedList &candidate)
{
	return (int)(log((double)lists[kid].size())/log(2.0)*candidate.size());
}

//assume id is sorted in ascending order
void InvertedIndexVectorBinSearch::MakeIntersection(int kid, vector<int> &candidate){
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	vector<int> temp;

	Record temp_rec;
	InvertedListEntry temp_entry;
	temp_entry.rec = &temp_rec;

	for (vector<int>::iterator iter = candidate.begin(); iter != candidate.end(); ++iter){
		temp_rec.id = *iter;

		if (binary_search(kvec.begin(), kvec.end(), temp_entry))
			temp.push_back(*iter);
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//make candidate sorted in ascending order
void InvertedIndexVectorBinSearch::MoveOut(int kid, vector<int> &candidate)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
		candidate.push_back(iter->rec->id);
}

//assume id is sorted in ascending order
void InvertedIndexVectorBinSearch::MakeIntersection(int kid, InvertedList &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	InvertedList temp;

	for (InvertedList::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
	{
		if (binary_search(kvec.begin(), kvec.end(), *iter, CompareInvertedListEntries))
			temp.push_back(*iter);
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//assume id is sorted in ascending order
void InvertedIndexVectorBinSearch::MakeIntersection(int kid, InvertedList &candidate, int bound)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	InvertedList temp;

	for (InvertedList::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if ((iter->rec->length-iter->position >= bound) && (binary_search(kvec.begin(), kvec.end(), *iter, CompareInvertedListEntries)))
		#else
		if ((iter->position+1 >= bound) && (binary_search(kvec.begin(), kvec.end(), *iter, CompareInvertedListEntries)))
		#endif
			temp.push_back(*iter);
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//make candidate sorted in ascending order
void InvertedIndexVectorBinSearch::MoveOut(int kid, InvertedList &candidate)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
		candidate.push_back(*iter);
}

//make candidate sorted in ascending order
void InvertedIndexVectorBinSearch::MoveOut(int kid, InvertedList &candidate, int bound)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if (iter->rec->length-iter->position >= bound)
		#else
		if (iter->position+1 >= bound)
		#endif
			candidate.push_back(*iter);
	}
}

void InvertedIndexVectorBinSearch::Print()
{
	for (hash<int, InvertedList >::iterator iter = this->lists.begin(); iter != this->lists.end(); ++iter)
	{
		cout << "k" << iter->first << ":";
		for (InvertedList::iterator iterL = iter->second.begin(); iterL != iter->second.end(); ++iterL)
			cout << " <r" << iterL->rec->id << "," << iterL->position << ">";
		cout << endl;
	}
}

void InvertedIndexVectorBinSearch::ShowInfo()
{
	printf("Using vector storage and binary search for intersection.\n");
}

////////////////////////////////////////////////////////////////////////////////
InvertedIndexVectorSmart::InvertedIndexVectorSmart()
{
}

InvertedIndexVectorSmart::InvertedIndexVectorSmart(Relation *R)
{
	InvertedListEntry e;


	#ifdef ASCENDING_KEYWORDS_ORDER
	for (int rid = 0; rid < R->numRecords; rid++)
	#else
	for (int rid = R->numRecords-1; rid >= 0; rid--)
	#endif
	{
		Record *r = (*R)[rid];

		for (int k = 0; k < r->length; k++)
		{
			e.rec = r;
			e.position = k;
			this->lists[r->keywords[k]].push_back(e);
		}
	}
}

InvertedIndexVectorSmart::~InvertedIndexVectorSmart()
{
}

bool InvertedIndexVectorSmart::ExistsKeyword(int kid)
{
	return (lists.count(kid) > 0);
}

int InvertedIndexVectorSmart::EstimateCost(int kid, vector<int> &candidate)
{
	int estimate_1 = (int)(lists[kid].size()+candidate.size());
	int estimate_2 = (int)(log((double)lists[kid].size())/log(2.0)*candidate.size());

	return min(estimate_1, estimate_2);
}

int InvertedIndexVectorSmart::EstimateCost(int kid, InvertedList &candidate)
{
	int estimate_1 = (int)(lists[kid].size()+candidate.size());
	int estimate_2 = (int)(log((double)lists[kid].size())/log(2.0)*candidate.size());

	return min(estimate_1, estimate_2);
}

//assume id is sorted in ascending order
void InvertedIndexVectorSmart::MakeIntersection(int kid, vector<int> &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	vector<int> temp;

	int size_il = (int)kvec.size();
	int size_cl = (int)candidate.size();
	double judge = pow(2.0, ((double)size_il/(double)size_cl+1.0)) - (double)size_il;

	if (judge > 0.0) //do binary search
	{
		Record temp_rec;
		InvertedListEntry temp_entry;
		temp_entry.rec = &temp_rec;

		for (vector<int>::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
		{
			temp_rec.id = *iter;

			if (binary_search(kvec.begin(), kvec.end(), temp_entry, CompareInvertedListEntries))
				temp.push_back(*iter);
		}
	}
	else //do merge join
	{
		InvertedList::iterator iter_x = kvec.begin(), iter_x_end = kvec.end();
		vector<int>::iterator iter_y = candidate.begin(), iter_y_end = candidate.end();

		#ifdef ASCENDING_KEYWORDS_ORDER
		while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
		{
			if (iter_x->rec->id < *iter_y) 
				iter_x++;
			else if (iter_x->rec->id > *iter_y)
				iter_y++;
			else
			{
				iter_x++;
				temp.push_back(*iter_y++);
			}
		}
		#else
		while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
		{
			if (iter_x->rec->id > *iter_y)
				iter_x++;
			else if (iter_x->rec->id < *iter_y)
				iter_y++;
			else
			{
				iter_x++;
				temp.push_back(*iter_y++);
			}
		}
		#endif
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//make candidate sorted in ascending order
void InvertedIndexVectorSmart::MoveOut(int kid, vector<int> &candidate)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
		candidate.push_back(iter->rec->id);
}

//assume id is sorted in ascending order
void InvertedIndexVectorSmart::MakeIntersection(int kid, InvertedList &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	InvertedList temp;

	int size_il = (int)kvec.size();
	int size_cl = (int)candidate.size();
	double judge = pow(2.0, ((double)size_il/(double)size_cl+1.0)) - (double)size_il;

	if (judge > 0.0) //do binary search
	{
		for (InvertedList::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
		{
			if (binary_search(kvec.begin(), kvec.end(), *iter, CompareInvertedListEntries))
				temp.push_back(*iter);
		}
	}
	else //do merge join
	{
		InvertedList::iterator iter_x = kvec.begin(), iter_x_end = kvec.end();
		InvertedList::iterator iter_y = candidate.begin(), iter_y_end = candidate.end();

		#ifdef ASCENDING_KEYWORDS_ORDER
		while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
		{
			if (iter_x->rec->id < iter_y->rec->id)
				iter_x++;
			else if (iter_x->rec->id > iter_y->rec->id)
				iter_y++;
			else
			{
				temp.push_back(*iter_x++);
				iter_y++;
			}
		}
		#else
		while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
		{
			if (iter_x->rec->id > iter_y->rec->id)
				iter_x++;
			else if (iter_x->rec->id < iter_y->rec->id)
				iter_y++;
			else
			{
				temp.push_back(*iter_x++);
				iter_y++;
			}
		}
		#endif
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//assume id is sorted in ascending order
void InvertedIndexVectorSmart::MakeIntersection(int kid, InvertedList &candidate, int bound)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we do merge join
	InvertedList &kvec = lists[kid];
	InvertedList temp;

	int size_il = (int)kvec.size();
	int size_cl = (int)candidate.size();
	double judge = pow(2.0, ((double)size_il/(double)size_cl+1.0)) - (double)size_il;

	if (judge > 0.0) //do binary search
	{
		for (InvertedList::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
		{
			#ifdef ASCENDING_KEYWORDS_ORDER
			if ((iter->rec->length-iter->position >= bound) && (binary_search(kvec.begin(), kvec.end(), *iter, CompareInvertedListEntries)))
			#else
			if ((iter->position+1 >= bound) && (binary_search(kvec.begin(), kvec.end(), *iter, CompareInvertedListEntries)))
			#endif
				temp.push_back(*iter);
		}
	}
	else //do merge join
	{
		InvertedList::iterator iter_x = kvec.begin(), iter_x_end = kvec.end();
		InvertedList::iterator iter_y = candidate.begin(), iter_y_end = candidate.end();

		#ifdef ASCENDING_KEYWORDS_ORDER
		while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
		{
			if (iter_x->rec->id < iter_y->rec->id)
				iter_x++;
			else if (iter_x->rec->id > iter_y->rec->id)
				iter_y++;
			else
			{
				if (iter_x->rec->length-iter_x->position >= bound)
					temp.push_back(*iter_x);
				iter_x++;
				iter_y++;
			}
		}
		#else
		while ((iter_x != iter_x_end) && (iter_y != iter_y_end))
		{
			if (iter_x->rec->id < iter_y->rec->id)
				iter_x++;
			else if (iter_x->rec->id > iter_y->rec->id)
				iter_y++;
			else
			{
				if (iter_x->position+1 >= bound)
					temp.push_back(*iter_x);
				iter_x++;
				iter_y++;
			}
		}
		#endif
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//make candidate sorted in ascending order
void InvertedIndexVectorSmart::MoveOut(int kid, InvertedList &candidate)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
		candidate.push_back(*iter);
}

//make candidate sorted in ascending order
void InvertedIndexVectorSmart::MoveOut(int kid, InvertedList &candidate, int bound)
{
	candidate.clear();

	InvertedList &kvec = lists[kid];
	for (InvertedList::iterator iter = kvec.begin(); iter != kvec.end(); ++iter)
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if (iter->rec->length-iter->position >= bound)
		#else
		if (iter->position+1 >= bound)
		#endif
			candidate.push_back(*iter);
	}
}

void InvertedIndexVectorSmart::Print()
{
	for (hash<int, InvertedList >::iterator iter = this->lists.begin(); iter != this->lists.end(); ++iter)
	{
		cout << "k" << iter->first << ":";
		for (InvertedList::iterator iterL = iter->second.begin(); iterL != iter->second.end(); ++iterL)
			cout << " <r" << iterL->rec->id << "," << iterL->position << ">";
		cout << endl;
	}
}

void InvertedIndexVectorSmart::ShowInfo()
{
	printf("Using vector storage and smart judgement for intersection.\n");
}

////////////////////////////////////////////////////////////////////////////////
InvertedIndexMap::InvertedIndexMap()
{
}

InvertedIndexMap::InvertedIndexMap(Relation *R)
{
	InvertedListEntry e;


	#ifdef ASCENDING_KEYWORDS_ORDER
	for (int rid = 0; rid < R->numRecords; rid++)
	#else
	for (int rid = R->numRecords-1; rid >= 0; rid--)
	#endif
	{
		Record *r = (*R)[rid];

		for (int k = 0; k < r->length; k++)
		{
			e.rec = r;
			e.position = k;
			this->lists[r->keywords[k]][r->id] = e;
		}
	}
}

InvertedIndexMap::~InvertedIndexMap()
{
}

bool InvertedIndexMap::ExistsKeyword(int kid)
{
	return (lists.count(kid) > 0);
}

int InvertedIndexMap::EstimateCost(int kid, vector<int> &candidate)
{
	return (int)(log((double)lists[kid].size())/log(2.0)*candidate.size());
}

int InvertedIndexMap::EstimateCost(int kid, InvertedList &candidate)
{
	return (int)(log((double)lists[kid].size())/log(2.0)*candidate.size());
}

//assume id is sorted in ascending order
void InvertedIndexMap::MakeIntersection(int kid, vector<int> &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we perform item probe
	map<int, InvertedListEntry> &kmap = lists[kid];
	vector<int> temp;

	for (vector<int>::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
	{
		if (kmap.count(*iter) > 0) 
			temp.push_back(*iter);
	}

	//make results go back to candidate
	candidate.swap(temp);
}

//make candidate sorted in ascending order
void InvertedIndexMap::MoveOut(int kid, vector<int> &candidate)
{
	candidate.clear();

	map<int, InvertedListEntry> &kmap = lists[kid];
	for (map<int, InvertedListEntry>::iterator iter = kmap.begin(); iter != kmap.end(); ++iter)
		candidate.push_back(iter->first);
}

//assume id is sorted in ascending order
void InvertedIndexMap::MakeIntersection(int kid, InvertedList &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we perform item probe
	map<int, InvertedListEntry> &kmap = lists[kid];
	InvertedList temp;
	
	for (InvertedList::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
	{
		map<int, InvertedListEntry>::iterator find_pos = kmap.find(iter->rec->id);
		
		if (find_pos != kmap.end()) temp.push_back(find_pos->second);
	}

	//make results go back to candidate
	candidate.swap(temp);
}

//assume id is sorted in ascending order
void InvertedIndexMap::MakeIntersection(int kid, InvertedList &candidate, int bound)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we perform item probe
	map<int, InvertedListEntry> &kmap = lists[kid];
	InvertedList temp;
	
	for (InvertedList::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
	{
		if (iter->rec->length-iter->position >= bound)
		{
			map<int, InvertedListEntry>::iterator find_pos = kmap.find(iter->rec->id);
		
			if (find_pos != kmap.end()) 
				temp.push_back(find_pos->second);
		}
	}

	//make results go back to candidate
	candidate.swap(temp);
}

//make candidate sorted in ascending order
void InvertedIndexMap::MoveOut(int kid, InvertedList &candidate)
{
	candidate.clear();

	map<int, InvertedListEntry> &kmap = lists[kid];
	for (map<int, InvertedListEntry>::iterator iter = kmap.begin(); iter != kmap.end(); ++iter)
		candidate.push_back(iter->second);
}

//make candidate sorted in ascending order
void InvertedIndexMap::MoveOut(int kid, InvertedList &candidate, int bound)
{
	candidate.clear();

	map<int, InvertedListEntry> &kmap = lists[kid];
	for (map<int, InvertedListEntry>::iterator iter = kmap.begin(); iter != kmap.end(); ++iter)
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if (iter->second.rec->length-iter->second.position >= bound)
		#else
		if (iter->second.position+1 >= bound)
		#endif
			candidate.push_back(iter->second);
	}
}

void InvertedIndexMap::Print()
{
	for (hash<int, map<int, InvertedListEntry> >::iterator iter = this->lists.begin(); iter != this->lists.end(); ++iter)
	{
		cout << "k" << iter->first << ":";
		for (map<int, InvertedListEntry>::iterator iterL = iter->second.begin(); iterL != iter->second.end(); ++iterL)
			cout << " <r" << iterL->second.rec->id << "," << iterL->second.position << ">";
		cout << endl;
	}
}

void InvertedIndexMap::ShowInfo()
{
	printf("Using map implementation for intersection.\n");
}

////////////////////////////////////////////////////////////////////////////////
InvertedIndexHashMap::InvertedIndexHashMap()
{
}

InvertedIndexHashMap::InvertedIndexHashMap(Relation *R)
{
	InvertedListEntry e;


	#ifdef ASCENDING_KEYWORDS_ORDER
	for (int rid = 0; rid < R->numRecords; rid++)
	#else
	for (int rid = R->numRecords-1; rid >= 0; rid--)
	#endif
	{
		Record *r = (*R)[rid];

		for (int k = 0; k < r->length; k++)
		{
			e.rec = r;
			e.position = k;
			this->lists[r->keywords[k]][rid] = e;
		}
	}
}

InvertedIndexHashMap::~InvertedIndexHashMap()
{
}

bool InvertedIndexHashMap::ExistsKeyword(int kid)
{
	return (lists.count(kid) > 0);
}

int InvertedIndexHashMap::EstimateCost(int kid, vector<int> &candidate)
{
	return (int)(candidate.size());
}

int InvertedIndexHashMap::EstimateCost(int kid, InvertedList &candidate)
{
	return (int)(candidate.size());
}

//assume id is sorted in ascending order
void InvertedIndexHashMap::MakeIntersection(int kid, vector<int> &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we perform item probe
	hash<int, InvertedListEntry> &khashmap = lists[kid];
	vector<int> temp;

	for (vector<int>::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
	{
		if (khashmap.count(*iter) > 0) temp.push_back(*iter);
	}

	//make results go back to candidate
	candidate.swap(temp);		
}

//make candidate sorted in ascending order
void InvertedIndexHashMap::MoveOut(int kid, vector<int> &candidate)
{
	candidate.clear();

	hash<int, InvertedListEntry> &khashmap = lists[kid];
	for (hash<int, InvertedListEntry>::iterator iter = khashmap.begin(); iter != khashmap.end(); ++iter)
		candidate.push_back(iter->first);

	sort(candidate.begin(), candidate.end());
}

//assume id is sorted in ascending order
void InvertedIndexHashMap::MakeIntersection(int kid, InvertedList &candidate)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we perform item probe
	hash<int, InvertedListEntry> &khashmap = lists[kid];
	InvertedList temp;
	
	for (InvertedList::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
	{
		hash<int, InvertedListEntry>::iterator find_pos = khashmap.find(iter->rec->id);
		
		if (find_pos != khashmap.end())
			temp.push_back(find_pos->second);
	}

	//make results go back to candidate
	candidate.swap(temp);
}

//assume id is sorted in ascending order
void InvertedIndexHashMap::MakeIntersection(int kid, InvertedList &candidate, int bound)
{
	//do not contain this keyword, result should be empty
	if (lists.count(kid) <= 0)
	{
		candidate.clear();
		return;
	}

	//or else we perform item probe
	hash<int, InvertedListEntry> &khashmap = lists[kid];
	InvertedList temp;
	
	for (InvertedList::iterator iter = candidate.begin(); iter != candidate.end(); ++iter)
	{
		if (iter->rec->length-iter->position >= bound)
		{
			hash<int, InvertedListEntry>::iterator find_pos = khashmap.find(iter->rec->id);
	
			if (find_pos != khashmap.end())
				temp.push_back(find_pos->second);
		}
	}

	//make results go back to candidate
	candidate.swap(temp);
}

//make candidate sorted in ascending order
void InvertedIndexHashMap::MoveOut(int kid, InvertedList &candidate)
{
	candidate.clear();

	hash<int, InvertedListEntry> &khashmap = lists[kid];
	for (hash<int, InvertedListEntry>::iterator iter = khashmap.begin(); iter != khashmap.end(); ++iter)
		candidate.push_back(iter->second);
}

//make candidate sorted in ascending order
void InvertedIndexHashMap::MoveOut(int kid, InvertedList &candidate, int bound)
{
	candidate.clear();

	hash<int, InvertedListEntry> &khashmap = lists[kid];
	for (hash<int, InvertedListEntry>::iterator iter = khashmap.begin(); iter != khashmap.end(); ++iter)
	{
		#ifdef ASCENDING_KEYWORDS_ORDER
		if (iter->second.rec->length-iter->second.position >= bound)
		#else
		if (iter->second.position+1 >= bound)
		#endif
			candidate.push_back(iter->second);
	}
}

void InvertedIndexHashMap::Print()
{
	for (hash<int, hash<int, InvertedListEntry> >::iterator iter = this->lists.begin(); iter != this->lists.end(); ++iter)
	{
		cout << "k" << iter->first << ":";
		for (hash<int, InvertedListEntry>::iterator iterL = iter->second.begin(); iterL != iter->second.end(); ++iterL)
			cout << " <r" << iterL->second.rec->id << "," << iterL->second.position << ">";
		cout << endl;
	}
}

void InvertedIndexHashMap::ShowInfo()
{
	printf("Using unordered_map implementation for intersection.\n");
}


