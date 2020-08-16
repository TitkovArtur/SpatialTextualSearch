#include "def.h"
#include "inverted_index.h"



//class CompareKeywordsByFrequencyR
//{
//public:
//	bool operator() (const int& klhs, const int& krhs) const
//	{
//		return (iidxR->lists[klhs].size() < iidxR->lists[krhs].size());
//	}
//};

//unsigned int SetContainmentQuery(int numKeywords, int *keywords, InvertedIndex *iidx)
//{
//	priority_queue<int, vector<int>, CompareKeywordsByFrequencyR> Q;
//	InvertedList cand_list;
//	int kid;
//	unsigned int numResults = 0;
//
//
//	// Examine keywords in query in ascending order w.r.t. their frequency
//	for (int j = 0; j < numKeywords; j++)
//		Q.push(keywords[j]);
//
//	kid = Q.top();
//	Q.pop();
//	//if (!idx->ExistsKeyword(kid))
//	//	return;
//	iidx->MoveOut(kid, cand_list);
//
//	while (!Q.empty())
//	{
//		kid = Q.top();
//		Q.pop();
//		iidx->MakeIntersection(kid, cand_list);
//	}
//
//	for (InvertedList::iterator iter = cand_list.begin(); iter != cand_list.end(); iter++)
//	{
//		//iter->rec->Print();
//		numResults++;
//	}
//
//
//	return numResults;
//}

//unsigned int TextualFirst(int RNumKeywords, int *RKeywords, int SNumKeywords, int *SKeywords, double threshold)
//{
//	InvertedList RList, SList;
//	int kid, numResults = 0;
//
//	kid = RKeywords[0];
//	iidxR->MoveOut(kid, RList);
//	for(int i = 1; i != RNumKeywords; ++i)
//		iidxR->MakeIntersection(RKeywords[i], RList);
//
//	kid = SKeywords[0];
//	iidxS->MoveOut(kid, SList);
//	for(int i = 1; i != SNumKeywords; ++i)
//		iidxS->MakeIntersection(SKeywords[i], SList);
//
//	for (InvertedList::iterator iterR = RList.begin(); iterR != RList.end(); iterR++)
//	{
//		for(InvertedList::iterator iterS = SList.begin(); iterS != SList.end(); iterS++)
//		{
//			if( QualifySpatial((*iterR).rec, (*iterS).rec, threshold * threshold) )
//				numResults++;
//		}
//	}
//
//	return numResults;
//}
//
//unsigned int TextualFirst(Relation *R, Relation *S, int RNumKeywords, int *RKeywords, int SNumKeywords, int *SKeywords, double threshold)
//{
//	sort(RKeywords, RKeywords + RNumKeywords);
//	sort(SKeywords, SKeywords + SNumKeywords);
//
//	vector<int> validR, validS;
//	for(int i = 0; i != R->numRecords; ++i)
//	{
//		if(ContainKeywords((*R)[i], RNumKeywords, RKeywords))
//		{
//			validR.push_back(i);
//		}
//	}
//	for(int i = 0; i != S->numRecords; ++i)
//	{
//		if(ContainKeywords((*S)[i], SNumKeywords, SKeywords))
//		{
//			validS.push_back(i);
//		}
//	}
//
//	//cout << validR.size() << " " << validS.size() << endl;
//
//	int numResults = 0;
//	for(vector<int>::iterator riter = validR.begin(); riter != validR.end(); ++riter)
//	{
//		for(vector<int>::iterator siter = validS.begin(); siter != validS.end(); ++siter)
//		{
//			if(QualifySpatial((*R)[*riter], (*S)[*siter], threshold*threshold) )
//				numResults++;
//		}
//	}
//
//	return numResults;
//}
//
//unsigned int SpatialFirst(
//                          Relation *R,
//                          RStarTree<TreeDataP, double> *rtR,
//                          Relation *S,
//                          RStarTree<TreeDataP, double> *rtS,
//                          double sthreshold,
//                          int RNumKeywords,
//                          int *RKeywords,
//                          int SNumKeywords,
//                          int *SKeywords){
//	vector<int> childR, childS;
//	deque<pair<RSTNode<TreeDataP, double> *, RSTNode<TreeDataP, double> *> > Q;
//	RSTNode<TreeDataP, double> *nodeR, *nodeS;
//	RSTNonLeafNode<TreeDataP, double> *nonleafR, *nonleafS;
//	RSTLeafNode<TreeDataP, double> *leafR, *leafS;
//	double sthreshold_sqr = sthreshold*sthreshold;
//	unsigned int numResults = 0;
//
//	childR.reserve(rtR->file->page_len);
//	childS.reserve(rtS->file->page_len);
//
//	Q.push_back(make_pair(rtR->root, rtS->root));
//	while(!Q.empty())
//	{
//		nodeR = Q.front().first;
//		nodeS = Q.front().second;
//		Q.pop_front();
//
//		if (!nodeR->is_leaf_node())
//		{
//			//Both internal nodes.
//			nonleafR = static_cast<RSTNonLeafNode<TreeDataP, double> *>(nodeR);
//			nonleafS = static_cast<RSTNonLeafNode<TreeDataP, double> *>(nodeS);
//
//			SpatialDistanceFilter(nonleafR, nonleafS, Q, childR, childS, sthreshold, sthreshold_sqr);
//		}
//		else
//		{
//			//Both leaf nodes.
//			leafR = static_cast<RSTLeafNode<TreeDataP, double> *>(nodeR);
//			leafS = static_cast<RSTLeafNode<TreeDataP, double> *>(nodeS);
//			//numResults += SpatialDistanceFilter(R, leafR, S, leafS, childR, childS, sthreshold, sthreshold_sqr);
//			numResults += SpatialDistanceFilter(R, leafR, S, leafS, childR, childS, sthreshold, sthreshold_sqr, RNumKeywords, RKeywords, SNumKeywords, SKeywords);
//		}
//	}
//
//	return numResults;
//}
//
//
//
//unsigned int SpatialJoin(Relation *R, RStarTree<TreeDataP, double> *rtR, Relation *S, RStarTree<TreeDataP, double> *rtS, double sthreshold)
//{
//	//temp working memory, avoid construction and destruction in function calls
//	vector<int> childR, childS;
//	deque<pair<RSTNode<TreeDataP, double> *, RSTNode<TreeDataP, double> *> > Q;
//	RSTNode<TreeDataP, double> *nodeR, *nodeS;
//	RSTNonLeafNode<TreeDataP, double> *nonleafR, *nonleafS;
//	RSTLeafNode<TreeDataP, double> *leafR, *leafS;
//	double sthreshold_sqr = sthreshold*sthreshold;
//	unsigned int numResults = 0;
//
//
//	childR.reserve(rtR->file->page_len);
//	childS.reserve(rtS->file->page_len);
//
//	Q.push_back(make_pair(rtR->root, rtS->root));
//	while(!Q.empty())
//	{
//		nodeR = Q.front().first;
//		nodeS = Q.front().second;
//		Q.pop_front();
//
//		if (!nodeR->is_leaf_node())
//		{
//			//Both internal nodes.
//			nonleafR = static_cast<RSTNonLeafNode<TreeDataP, double> *>(nodeR);
//			nonleafS = static_cast<RSTNonLeafNode<TreeDataP, double> *>(nodeS);
//
//			SpatialDistanceFilter(nonleafR, nonleafS, Q, childR, childS, sthreshold, sthreshold_sqr);
//		}
//		else
//		{
//			//Both leaf nodes.
//			leafR = static_cast<RSTLeafNode<TreeDataP, double> *>(nodeR);
//			leafS = static_cast<RSTLeafNode<TreeDataP, double> *>(nodeS);
//
//			numResults += SpatialDistanceFilter(R, leafR, S, leafS, childR, childS, sthreshold, sthreshold_sqr);
//		}
//	}
//
//	return numResults;
//}
//
//
//unsigned int SpatialJoinNaive(Relation *R, Relation *S, double threshold)
//{
//	Record *r, *s;
//	double threshold_sqr = threshold*threshold;
//	unsigned int numResults = 0;
//
//
//	for (int rid = 0; rid < R->numRecords; rid++)
//	{
//		r = (*R)[rid];
//		for (int sid = 0; sid < S->numRecords; sid++)
//		{
//			s = (*S)[sid];
//			if (QualifySpatial(r, s, threshold_sqr))
//			{
//				numResults++;
//
//				//cout << rid << " " << sid << endl;
//			}
//		}
//	}
//
//	return numResults;
//}
