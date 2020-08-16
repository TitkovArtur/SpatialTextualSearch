#include "def.h"
#include "relation.h"
#include "inverted_index.h"







//void ChooseKeywords()
//{
//	int max = 0;
//	int numKeywords = 2;
//	int *keywords = new int[numKeywords];
//	int *keywords_max = new int[numKeywords];
//	for (int i = S->numKeywords - 2; i >= S->numKeywords - 20; --i)
//	{
//		keywords[0] = i;
//		for (int j = i - 1; j >= i - 40; --j)
//		{
//			keywords[1] = j;
//			int count = SetContainmentQuery(numKeywords, keywords, iidxS);
//			if (count > max)
//			{
//				max = count;
//				keywords_max[0] = i;
//				keywords_max[1] = j;
//			}
//		}
//	}
//	cout << keywords_max[0] << ", " << keywords_max[1] << ": " << max << endl;
//
//}
//
//bool PruneByKw(InvertedIndex *iidx, bool *invalid, int page_id, int numKw, int *kw)
//{
//	if(invalid[page_id] == true)
//		return true;
//	int i;
//	for(i = 0; i != numKw; ++i)
//	{
//		if(!iidx[page_id].ExistsKeyword(kw[i]))
//			break;
//	}
//	if(i != numKw)
//	{
//		invalid[page_id] = true;
//		return true;
//	}
//
//	return false;
//}


RStarTree<TreeDataP, double>* BuildRTree(const char *datafile, Relation *RL, bool build_new){
//	cout << datafile << ": ";

	//RL = new Relation(datafile);
	RStarTree<TreeDataP, double>* rt;

    string str = string(datafile) + ".rt";

	if (build_new){
		TreeDataP<double> temp(2);
		rt = new RStarTree<TreeDataP, double>(str.c_str(), 2, 4096, "", PageFile::C_FULLMEM, 100, true);
//        rt = new RStarTree<TreeDataP, double>( 2, 4096, "", PageFile::C_FULLMEM, 100);

		for (int rid = 0; rid < RL->numRecords; rid++){
			temp.data[0] = (*RL)[rid]->locx;
			temp.data[1] = (*RL)[rid]->locy;
			temp.dist = 0;
			temp.id = rid;
			rt->insert(&temp);
		}
//		cout << "R*tree created..." << endl;
	}
	else
	{
//		rt = new RStarTree<TreeDataP, double>(str.c_str(), 2, 4096, "", PageFile::C_FULLMEM, 100, false);
//		cout << "R*tree loaded..." << endl;
	}

	return rt;
}


//bool intersect(int dim, double *mbr, double *mbrR, double *mbrS)
//{
//	for (int i = 0; i < dim * 2; i += 2)
//	{
//		mbr[i] = max(mbrR[i], mbrS[i]);
//		mbr[i + 1] = min(mbrR[i + 1], mbrS[i + 1]);
//		if (mbr[i] > mbr[i + 1])
//			return false;
//	}
//
//	return true;
//}
//
//bool intersect(int dim, double *mbrR, double *mbrS)
//{
//	for (int i = 0; i < dim * 2; i += 2)
//	{
//		if ((mbrR[i] > mbrS[i + 1]) || (mbrS[i] > mbrR[i + 1]))
//			return false;
//	}
//
//	return true;
//}
//
//bool inside(int dim, double *p, double *mbr)
//{
//	for (int i = 0; i < dim; i++)
//	{
//		if ((p[i] < mbr[2 * i]) || (p[i] > mbr[2 * i + 1]))
//			return false;
//	}
//
//	return true;
//}

bool QualifyMinMinDistSqr(double *mbrR, double *mbrS, double threshold_sqr)
{
	//judge spatial threshold
	double v1 = max(mbrR[0], mbrS[0]);
	double v2 = min(mbrR[1], mbrS[1]);
	double d2 = ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
	if (d2 > threshold_sqr)
		return false;

	v1 = max(mbrR[2], mbrS[2]);
	v2 = min(mbrR[3], mbrS[3]);
	d2 += ((v1 > v2) ? ((v1 - v2) * (v1 - v2)) : 0);
	if (d2 > threshold_sqr)
		return false;

	return true;
}

bool QualifySpatial(Record *r, Record *s, double dthreshold_sqr)
{
	double dist;

	if ((dist = (r->locx - s->locx) * (r->locx - s->locx)) > dthreshold_sqr)
		return false;

	if (dist + (r->locy - s->locy) * (r->locy - s->locy) > dthreshold_sqr)
		return false;

	return true;
}

//Shuyao, modified from STSJOIN
void CreateIRTreeNode(RSTNode<TreeDataP, double> *node, Relation *R, InvertedIndex *iidx){
	Record *r;
	InvertedListEntry entry;
	RSTLeafNode<TreeDataP, double> *leaf;
	RSTNonLeafNode<TreeDataP, double> *nonleaf;
	int k;

	if (node->is_leaf_node()){
            leaf = static_cast<RSTLeafNode<TreeDataP, double> *>(node);
        
            for (int eid = 0; eid < leaf->entry_num; eid++){ // for each record this leaf
                r = (*R)[leaf->data[eid]->id];
                for (int i = 0; i < r->length; i++){ // for all terms in record r
                    k = r->keywords[i];
//                    entry.rec = r;
                    entry.rec = new Record();
                    entry.rec->id = eid;
                    entry.position = i;

                    iidx[leaf->page_id].lists[k].push_back(entry);
                }
            }
        }else{
		nonleaf = static_cast<RSTNonLeafNode<TreeDataP, double> *>(node);
            
		for (int eid = 0; eid < nonleaf->entry_num; eid++){
			//for each entry, shuyao
			RSTNode<TreeDataP, double> *child = nonleaf->get_child(eid);

			CreateIRTreeNode(child, R, iidx);
            
			//for each lists of this child
			for (hash<int, vector<InvertedListEntry> >::iterator iter = iidx[child->page_id].lists.begin(); iter != iidx[child->page_id].lists.end(); ++iter){
				//keyword
				k = iter->first;
                
				//non-leaf node, so no record
				entry.rec = new Record();
                entry.rec->id = eid;
				//position records node
				entry.position = child->page_id;

				iidx[nonleaf->page_id].lists[k].push_back(entry);
			}
		}
	}
}




//Shuyao
bool ContainKeywords(Record *r, int numKeywords, int *keywords)
{
	if (numKeywords == 0)
		return true;

	sort(keywords, keywords + numKeywords);
	sort(r->keywords, r->keywords + r->length);

	int x = 0, y = 0, contain = 0;
	while (x != numKeywords && y != r->length)
	{
		if (keywords[x] < r->keywords[y])
			x++;
		else if (keywords[x] > r->keywords[y])
			y++;
		else
		{
			//cout << keywords[x] << " " << r->keywords[y] << endl;
			contain++;
			x++;
			y++;
		}
	}
	//cout << "---------------------" << endl;

	//cout << contain << " " << numKeywords << endl;
	if (contain == numKeywords)
		return true;
	else
		return false;
}

bool ContainKeywords(Record *r, bool *invalidRec, int numKeywords,
		int *keywords)
{

	if (numKeywords == 0)
		return true;

	if (invalidRec[r->id] == true)
	{
		//cout << "invalid" << endl;
		return false;
	}

	sort(keywords, keywords + numKeywords);
	sort(r->keywords, r->keywords + r->length);

	int x = 0, y = 0, contain = 0;
	while (x != numKeywords && y != r->length)
	{
		if (keywords[x] < r->keywords[y])
			x++;
		else if (keywords[x] > r->keywords[y])
			y++;
		else
		{
			//cout << keywords[x] << " " << r->keywords[y] << endl;
			contain++;
			x++;
			y++;
		}
	}
	//cout << "---------------------" << endl;

	//cout << contain << " " << numKeywords << endl;
	if (contain == numKeywords)
		return true;
	else
	{
		invalidRec[r->id] = true;
		return false;
	}
}

//Shuyao
//bool QualifyTextual(Record *r, Record *s, int numKwR, int *kwR, int numKwS,
//		int *kwS)
//{
//	if (ContainKeywords(r, invalidRecR, numKwR, kwR)
//			&& ContainKeywords(s, invalidRecS, numKwS, kwS))
//		return true;
//	else
//		return false;
//}

//void SpatialDistanceFilter(RSTNonLeafNode<TreeDataP, double> *nonleafR,
//		RSTNonLeafNode<TreeDataP, double>* nonleafS,
//		deque<pair<RSTNode<TreeDataP, double> *, RSTNode<TreeDataP, double> *> > &q,
//		vector<int> &childR, vector<int> &childS, double threshold,
//		double threshold_sqr)
//{
//	//according to T. Brinkhoff et al SIGMOD 93
//	//first we restrict the search space, one difference is that we now expand the MBR -/+dist
//	//we expanf the two nonleaf MBR and calculate the intersection I
//	//we can prove that if any child MBR do not intersect I, then it will not be in the result set
//	double mbrR[4], mbrS[4], inter[4];
//
//	nonleafR->get_mbr(mbrR);
//	mbrR[0] -= threshold;
//	mbrR[1] += threshold; //expand x-axis
//	mbrR[2] -= threshold;
//	mbrR[3] += threshold; //expand y-axis
//	nonleafS->get_mbr(mbrS);
//	mbrS[0] -= threshold;
//	mbrS[1] += threshold; //expand x-axis
//	mbrS[2] -= threshold;
//	mbrS[3] += threshold; //expand y-axis
//	if (!intersect(2, inter, mbrR, mbrS))
//		return; //no intersection we return
//	if (!QualifyMinMinDistSqr(mbrR, mbrS, threshold_sqr))
//		return; //judge false hit
//
//	//next we will pickup child which intersects with I
//	//sort them according to the 1st dimension
//	//since we are using distance intersection, wo we treat childR as normal
//	//and childS as if they are all expanded by threshold
//	double **entryR = nonleafR->entry_mbr, **entryS = nonleafS->entry_mbr;
//	childR.clear();
//	for (int i = 0; i < nonleafR->entry_num; i++)
//		if (intersect(2, entryR[i], inter))
//			childR.push_back(i);
//	sort(childR.begin(), childR.end(), mbr_sort_dim_x_less(nonleafR));
//
//	childS.clear();
//	for (int i = 0; i < nonleafS->entry_num; i++)
//		if (intersect(2, entryS[i], inter))
//			childS.push_back(i);
//	sort(childS.begin(), childS.end(), mbr_sort_dim_x_less(nonleafS));
//
//	//next we will use plane-sweep to test intersection
//	vector<int>::iterator iter1 = childR.begin(), iter2 = childS.begin();
//	vector<int>::iterator end1 = childR.end(), end2 = childS.end();
//	while ((iter1 < end1) && (iter2 < end2))
//	{
//		if (entryR[*iter1][0] < (entryS[*iter2][0] - threshold)) //now asymmetric, we use two hand-coded internal loops
//		{
//			//internal loop 1
//			vector<int>::iterator iter = iter2;
//			while ((iter < end2)
//					&& ((entryS[*iter][0] - threshold) < entryR[*iter1][1]))
//			{
//				if ((entryR[*iter1][2] > (entryS[*iter][3] + threshold))
//						|| ((entryS[*iter][2] - threshold) > entryR[*iter1][3]))
//				{
//					++iter;
//					continue;
//				}
//
//				RSTNode<TreeDataP, double> *nodeR = nonleafR->get_child(*iter1);
//				RSTNode<TreeDataP, double> *nodeS = nonleafS->get_child(*iter);
//				q.push_back(make_pair(nodeR, nodeS));
//
//				++iter;
//			}
//
//			++iter1;
//		}
//		else
//		{
//			//internal loop 2
//			vector<int>::iterator iter = iter1;
//			while ((iter < end1)
//					&& (entryR[*iter][0] < (entryS[*iter2][1] + threshold)))
//			{
//				if (((entryS[*iter2][2] - threshold) > entryR[*iter][3])
//						|| (entryR[*iter][2] > (entryS[*iter2][3] + threshold)))
//				{
//					++iter;
//					continue;
//				}
//
//				RSTNode<TreeDataP, double> *nodeR = nonleafR->get_child(*iter);
//				RSTNode<TreeDataP, double> *nodeS = nonleafS->get_child(*iter2);
//				q.push_back(make_pair(nodeR, nodeS));
//
//				++iter;
//			}
//
//			++iter2;
//		}
//	}
//}

//unsigned int SpatialDistanceFilter(Relation *R,
//		RSTLeafNode<TreeDataP, double> *leafR, Relation *S,
//		RSTLeafNode<TreeDataP, double>* leafS, vector<int> &childR,
//		vector<int> &childS, double threshold, double threshold_sqr,
//		int RNumKeywords, int *RKeywords, int SNumKeywords, int *SKeywords)
//{
//	//according to T. Brinkhoff et al SIGMOD 93
//	//first we restrict the search space, one difference is that we now expand the MBR -/+dist
//	//we expanf the two nonleaf MBR and calculate the intersection I
//	//we can prove that if any child MBR do not intersect I, then it will not be in the result set
//	double mbrR[4], mbrS[4], inter[4];
//	unsigned int numResults = 0;
//
//	leafR->get_mbr(mbrR);
//	mbrR[0] -= threshold;
//	mbrR[1] += threshold; //expand x-axis
//	mbrR[2] -= threshold;
//	mbrR[3] += threshold; //expand y-axis
//	leafS->get_mbr(mbrS);
//	mbrS[0] -= threshold;
//	mbrS[1] += threshold; //expand x-axis
//	mbrS[2] -= threshold;
//	mbrS[3] += threshold; //expand y-axis
//	if (!intersect(2, inter, mbrR, mbrS))
//		return 0; //no intersection we return
//	if (!QualifyMinMinDistSqr(mbrR, mbrS, threshold_sqr))
//		return 0; //judge false hit
//
//	//next we will pickup child which intersects with I
//	//sort them according to the 1st dimension
//	//since we are using distance intersection, wo we treat childR as normal
//	//and childS as if they are all expanded by threshold
//	TreeDataP<double> **entryR = leafR->data;
//	childR.clear();
//	for (int i = 0; i < leafR->entry_num; i++)
//		if (inside(2, entryR[i]->data, inter))
//			childR.push_back(i);
//	sort(childR.begin(), childR.end(), p_sort_dim_x_less(leafR));
//
//	TreeDataP<double> **entryS = leafS->data;
//	childS.clear();
//	for (int i = 0; i < leafS->entry_num; i++)
//		if (inside(2, entryS[i]->data, inter))
//			childS.push_back(i);
//	sort(childS.begin(), childS.end(), p_sort_dim_x_less(leafS));
//
//	vector<int>::iterator iter1 = childR.begin(), iter2 = childS.begin();
//	vector<int>::iterator end1 = childR.end(), end2 = childS.end();
//	while ((iter1 < end1) && (iter2 < end2))
//	{
//		if (entryR[*iter1]->data[0] < (entryS[*iter2]->data[0] - threshold)) //now asymmetric, we use two hand-coded internal loops
//		{
//			//internal loop 1
//			vector<int>::iterator iter = iter2;
//			while ((iter < end2)
//					&& ((entryS[*iter]->data[0] - threshold)
//							< entryR[*iter1]->data[0]))
//			{
//				int id1 = entryR[*iter1]->id;
//				int id2 = entryS[*iter]->id;
//
//				if ((entryR[*iter1]->data[1]
//						> (entryS[*iter]->data[1] + threshold))
//						|| ((entryS[*iter]->data[1] - threshold)
//								> entryR[*iter1]->data[1]))
//				{
//					++iter;
//					continue;
//				}
//
//#ifdef SPATIO_TEXTUAL
//				if ((QualifySpatial((*R)[id1], (*S)[id2], threshold_sqr))
//						&& (QualifyTextual((*R)[id1], (*S)[id2], RNumKeywords,
//								RKeywords, SNumKeywords, SKeywords)))
//#else
//				if (QualifySpatial((*R)[id1], (*S)[id2], threshold_sqr))
//#endif
//				{
//					numResults++;
//					//cout << id1 << " " << id2 << endl;
//				}
//
//				++iter;
//			}
//
//			++iter1;
//		}
//		else
//		{
//			//internal loop 2
//			vector<int>::iterator iter = iter1;
//			while ((iter < end1)
//					&& (entryR[*iter]->data[0]
//							< (entryS[*iter2]->data[0] + threshold)))
//			{
//				int id1 = entryR[*iter]->id;
//				int id2 = entryS[*iter2]->id;
//
//				if (((entryS[*iter2]->data[1] - threshold)
//						> entryR[*iter]->data[1])
//						|| (entryR[*iter]->data[1]
//								> (entryS[*iter2]->data[1] + threshold)))
//				{
//					++iter;
//					continue;
//				}
//
//#ifdef SPATIO_TEXTUAL
//				if ((QualifySpatial((*R)[id1], (*S)[id2], threshold_sqr))
//						&& (QualifyTextual((*R)[id1], (*S)[id2], RNumKeywords,
//								RKeywords, SNumKeywords, SKeywords)))
//#else
//				if (QualifySpatial((*R)[id1], (*S)[id2], threshold_sqr))
//#endif
//				{
//					numResults++;
//					//cout << id1 << " " << id2 << endl;
//				}
//
//				++iter;
//			}
//
//			++iter2;
//		}
//	}
//
//	return numResults;
//}
//
////Shuyao
//unsigned int SpatialDistanceFilter(Relation *R,
//		RSTLeafNode<TreeDataP, double> *leafR, Relation *S,
//		RSTLeafNode<TreeDataP, double>* leafS, vector<int> &childR,
//		vector<int> &childS, double threshold, double threshold_sqr)
//{
//	return SpatialDistanceFilter(R, leafR, S, leafS, childR, childS, threshold,
//			threshold_sqr, 0, NULL, 0, NULL);
//}
//
//void IRTreeFilter(RSTNonLeafNode<TreeDataP, double> *nonleafR,
//		RSTNonLeafNode<TreeDataP, double>* nonleafS,
//		deque<pair<RSTNode<TreeDataP, double> *, RSTNode<TreeDataP, double> *> > &q,
//		vector<int> &childR, vector<int> &childS, double threshold,
//		double threshold_sqr, int numKwR, int *kwR, int numKwS, int *kwS)
//{
//	vector<int> RList, SList;
//	IIR[nonleafR->page_id].MoveOutNonleaf(kwR[0], RList);
//	for(int i = 1; i != numKwR; ++i)
//		IIR[nonleafR->page_id].MakeIntersectionNonleaf(kwR[i], RList, invalidR);
//
//	IIS[nonleafS->page_id].MoveOutNonleaf(kwS[0], SList);
//	for(int i = 1; i != numKwS; ++i)
//		IIS[nonleafS->page_id].MakeIntersectionNonleaf(kwS[i], SList, invalidS);
//
//	double mbrR[4], mbrS[4], inter[4];
//
//	nonleafR->get_mbr(mbrR);
//	mbrR[0] -= threshold;
//	mbrR[1] += threshold; //expand x-axis
//	mbrR[2] -= threshold;
//	mbrR[3] += threshold; //expand y-axis
//	nonleafS->get_mbr(mbrS);
//	mbrS[0] -= threshold;
//	mbrS[1] += threshold; //expand x-axis
//	mbrS[2] -= threshold;
//	mbrS[3] += threshold; //expand y-axis
//	if (!intersect(2, inter, mbrR, mbrS))
//		return; //no intersection we return
//	if (!QualifyMinMinDistSqr(mbrR, mbrS, threshold_sqr))
//		return; //judge false hit
//
//	//next we will pickup child which intersects with I
//	//sort them according to the 1st dimension
//	//since we are using distance intersection, wo we treat childR as normal
//	//and childS as if they are all expanded by threshold
//	double **entryR = nonleafR->entry_mbr, **entryS = nonleafS->entry_mbr;
//	childR.clear();
//
//	//cout << nonleafR->entry_num << " " << RList.size() << endl;
//	//cout << nonleafS->entry_num << " " << SList.size() << endl;
//
//	for (int i = 0; i < nonleafR->entry_num; i++)
//	{
//		if(find(RList.begin(), RList.end(), nonleafR->get_child(i)->page_id) == RList.end())
//		{
//			continue;
//		}
//		if (intersect(2, entryR[i], inter))
//			childR.push_back(i);
//	}
//	sort(childR.begin(), childR.end(), mbr_sort_dim_x_less(nonleafR));
//
//	childS.clear();
//	for (int i = 0; i < nonleafS->entry_num; i++)
//	{
//		if(find(SList.begin(), SList.end(), nonleafS->get_child(i)->page_id) == SList.end())
//		{
//			continue;
//		}
//		if (intersect(2, entryS[i], inter))
//			childS.push_back(i);
//	}
//	sort(childS.begin(), childS.end(), mbr_sort_dim_x_less(nonleafS));
//
//	//next we will use plane-sweep to test intersection
//	vector<int>::iterator iter1 = childR.begin(), iter2 = childS.begin();
//	vector<int>::iterator end1 = childR.end(), end2 = childS.end();
//	while ((iter1 < end1) && (iter2 < end2))
//	{
//		if (entryR[*iter1][0] < (entryS[*iter2][0] - threshold)) //now asymmetric, we use two hand-coded internal loops
//		{
//			//internal loop 1
//			vector<int>::iterator iter = iter2;
//			while ((iter < end2)
//					&& ((entryS[*iter][0] - threshold) < entryR[*iter1][1]))
//			{
//				if ((entryR[*iter1][2] > (entryS[*iter][3] + threshold))
//						|| ((entryS[*iter][2] - threshold) > entryR[*iter1][3]))
//				{
//					++iter;
//					continue;
//				}
//
//				RSTNode<TreeDataP, double> *nodeR = nonleafR->get_child(*iter1);
//				RSTNode<TreeDataP, double> *nodeS = nonleafS->get_child(*iter);
//				q.push_back(make_pair(nodeR, nodeS));
//
//				++iter;
//			}
//
//			++iter1;
//		}
//		else
//		{
//			//internal loop 2
//			vector<int>::iterator iter = iter1;
//			while ((iter < end1)
//					&& (entryR[*iter][0] < (entryS[*iter2][1] + threshold)))
//			{
//				if (((entryS[*iter2][2] - threshold) > entryR[*iter][3])
//						|| (entryR[*iter][2] > (entryS[*iter2][3] + threshold)))
//				{
//					++iter;
//					continue;
//				}
//
//				RSTNode<TreeDataP, double> *nodeR = nonleafR->get_child(*iter);
//				RSTNode<TreeDataP, double> *nodeS = nonleafS->get_child(*iter2);
//				q.push_back(make_pair(nodeR, nodeS));
//
//				++iter;
//			}
//
//			++iter2;
//		}
//	}
//}
//
//unsigned int IRTreeFilter(Relation *R, RSTLeafNode<TreeDataP, double> *leafR,
//		Relation *S, RSTLeafNode<TreeDataP, double>* leafS, vector<int> &childR,
//		vector<int> &childS, double threshold, double threshold_sqr, int numKwR,
//		int *kwR, int numKwS, int *kwS)
//{
////	if(PruneByKw(IIR, invalidR, leafR->page_id, numKwR, kwR))
////		continue;
////	if(PruneByKw(IIS, invalidS, leafS->page_id, numKwS, kwS))
////		continue;
//
//	//according to T. Brinkhoff et al SIGMOD 93
//	//first we restrict the search space, one difference is that we now expand the MBR -/+dist
//	//we expanf the two nonleaf MBR and calculate the intersection I
//	//we can prove that if any child MBR do not intersect I, then it will not be in the result set
//	double mbrR[4], mbrS[4], inter[4];
//	unsigned int numResults = 0;
//
//	leafR->get_mbr(mbrR);
//	mbrR[0] -= threshold;
//	mbrR[1] += threshold; //expand x-axis
//	mbrR[2] -= threshold;
//	mbrR[3] += threshold; //expand y-axis
//	leafS->get_mbr(mbrS);
//	mbrS[0] -= threshold;
//	mbrS[1] += threshold; //expand x-axis
//	mbrS[2] -= threshold;
//	mbrS[3] += threshold; //expand y-axis
//	if (!intersect(2, inter, mbrR, mbrS))
//		return 0; //no intersection we return
//	if (!QualifyMinMinDistSqr(mbrR, mbrS, threshold_sqr))
//		return 0; //judge false hit
//
//	//next we will pickup child which intersects with I
//	//sort them according to the 1st dimension
//	//since we are using distance intersection, wo we treat childR as normal
//	//and childS as if they are all expanded by threshold
//	TreeDataP<double> **entryR = leafR->data;
//	childR.clear();
//	for (int i = 0; i < leafR->entry_num; i++)
//	{
//		if (inside(2, entryR[i]->data, inter))
//			childR.push_back(i);
//	}
//	sort(childR.begin(), childR.end(), p_sort_dim_x_less(leafR));
//
//	TreeDataP<double> **entryS = leafS->data;
//	childS.clear();
//	for (int i = 0; i < leafS->entry_num; i++)
//	{
//		if (inside(2, entryS[i]->data, inter))
//			childS.push_back(i);
//	}
//	sort(childS.begin(), childS.end(), p_sort_dim_x_less(leafS));
//
//	vector<int>::iterator iter1 = childR.begin(), iter2 = childS.begin();
//	vector<int>::iterator end1 = childR.end(), end2 = childS.end();
//	while ((iter1 < end1) && (iter2 < end2))
//	{
//		if (entryR[*iter1]->data[0] < (entryS[*iter2]->data[0] - threshold)) //now asymmetric, we use two hand-coded internal loops
//		{
//			//internal loop 1
//			vector<int>::iterator iter = iter2;
//			while ((iter < end2)
//					&& ((entryS[*iter]->data[0] - threshold)
//							< entryR[*iter1]->data[0]))
//			{
//				int id1 = entryR[*iter1]->id;
//				int id2 = entryS[*iter]->id;
//
//				if ((entryR[*iter1]->data[1]
//						> (entryS[*iter]->data[1] + threshold))
//						|| ((entryS[*iter]->data[1] - threshold)
//								> entryR[*iter1]->data[1]))
//				{
//					++iter;
//					continue;
//				}
//
//				if ((QualifySpatial((*R)[id1], (*S)[id2], threshold_sqr))
//						&& (QualifyTextual((*R)[id1], (*S)[id2], numKwR, kwR,
//								numKwS, kwS)))
//				{
//					numResults++;
//					//cout << id1 << " " << id2 << endl;
//				}
//
//				++iter;
//			}
//
//			++iter1;
//		}
//		else
//		{
//			//internal loop 2
//			vector<int>::iterator iter = iter1;
//			while ((iter < end1)
//					&& (entryR[*iter]->data[0]
//							< (entryS[*iter2]->data[0] + threshold)))
//			{
//				int id1 = entryR[*iter]->id;
//				int id2 = entryS[*iter2]->id;
//
//				if (((entryS[*iter2]->data[1] - threshold)
//						> entryR[*iter]->data[1])
//						|| (entryR[*iter]->data[1]
//								> (entryS[*iter2]->data[1] + threshold)))
//				{
//					++iter;
//					continue;
//				}
//
//				if ((QualifySpatial((*R)[id1], (*S)[id2], threshold_sqr))
//						&& (QualifyTextual((*R)[id1], (*S)[id2], numKwR, kwR,
//								numKwS, kwS)))
//				{
//					numResults++;
//					//cout << id1 << " " << id2 << endl;
//				}
//
//				++iter;
//			}
//
//			++iter2;
//		}
//	}
//
//	return numResults;
//}


RStarTree<TreeDataP, double>* bulkload(int dim, int page_len, TreeDataP<double> **data, int data_num)
{
    //the memory policy of page file
    PageFile::c_policy pols[] =
    { PageFile::C_FULLMEM, PageFile::C_LRU, PageFile::C_MRU,
        PageFile::C_NO_CACHE };
    int pagefile_cache_size = 0; //we use full memory
    RStarTree<TreeDataP, double> *tree = new RStarTree<TreeDataP, double>(dim, page_len, pols[0], pagefile_cache_size);

    //construct tree using STR bulkloading
    tree->bulkload_str(data, data_num, 0.7);
    //cout << "bulkload: " << (double) tim.stopClock() / CLOCKS_PER_SEC << " secs" << endl;

    return tree;
}
