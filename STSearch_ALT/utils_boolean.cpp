//
// Project:    Spatially Combined Text Searches
// Filename:   utils_spatial.cpp
// Created on: 06.07.20.
//
// Developer(s):
//    I) Artur Titkov
//       ArturTitkov@icloud.com
//       https://github.com/TitkovArtur
//    II) unknown
//
// Description:
//
//
// Copyright © 2020 Artur Titkov. All rights reserved.
//

#include "def.h"
#include "relation.h"
#include "inverted_index.h"
class Relation;
class InvertedListEntry;



//FUNCtION DECLARATIONS
inline bool bin_search(InvertedListEntry* l, int start, int end, int* pos, int x);



inline Relation* ContainmentQueryOfRelation(Relation& T, int numKeywords, int *keywords, InvertedIndex *iidx);
inline Relation* ContainmentQueryWithBinSearchOfRelation(Relation& T, int numKeywords, int *keywords, InvertedIndex *iidx);


inline bool textProbeOnRecord(termSet t, termSetSize n, Record* r); // NOT USED
inline bool textProbeOnDataPoint(termSet t, termSetSize n, Record* r, float* p); // NOT USED
inline bool textVerificationOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node); // NOT USED



inline bool textProbeOnDataPointBinSearch(termSet t, termSetSize n, Record* r, float* p); // TEST ENTRY OF ONE LEAF


inline bool textSearchOnLeaf(int* keywords, int numKeywords, InvertedIndex* iidx, Node* leaf); // NAIVE APPROACH
inline bool textSearchOnNodeBinSearch(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);
inline bool textProbeOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node);


static void Print(vector<InvertedListEntry>* list){
    cout << "Posting List size :" << list->size() << "\n";
    for(int i = 0; i < list->size(); i++){
        cout << list->at(i).rec->id << " ";
    }
    cout << "\n";
    
    
}


//NAIVE
inline bool textProbeOnRecord(termSet t, termSetSize n, Record* r){

    int i = 0;
    int j = 0;

    while( i < n && j < r->length){
        if(t[i] < r->keywords[j]){
            r->textScore = -1;
            return false;
        }else if(t[i] > r->keywords[j]){
            j++;
        }else if(t[i] == r->keywords[j]) {
            i++;
            j++;
        }
    }
    if(i == n){
        r->textScore = 1;
        return true;
    }else{
        r->textScore = -1;
        return false;
    }
};


//DO NOT KNOW 
inline bool textProbe(termSet t, termSetSize n, Record& r){
  
    int i = 0;
    int j = 0;

    while( i < n && j < r.length){
        if(t[i] < r.keywords[j]){
            return false;
        }else if(t[i] > r.keywords[j]){
            j++;
        }else if(t[i] == r.keywords[j]) {
            i++;
            j++;
        }
    }
    return i == n ? true : false ;
};




inline Relation* ContainmentQueryOfRelation(Relation& T, int numKeywords, int *keywords, InvertedIndex *iidx){
    
    std::vector<int>* result = new std::vector<int>;
    std::vector<InvertedListEntry*> lists;
    std::vector<int> actualPosition;
    Relation* output = new Relation();
    int listLengths[numKeywords];
    int tmp[numKeywords];
    
    
    //create tmp keyword set
    for(int i = 0; i < numKeywords; i++){
        tmp[i] = keywords[i];
    }

    //get postingList lengths for keywords
    for(int i = 0; i < numKeywords; i++){
        listLengths[i] = iidx->lists.at(tmp[i]).size();
    }
    
    //sort termSet regarding listLengths (bubble sort) // already sorted
//    for(int i = 0; i < terms; i++){
//        for(int j = 0; j < terms-1-i; j++){
//            if(listLengths[j] > listLengths[j+1]){
//                int t_tmp = tmp[j];
//                int l_tmp = listLengths[j];
//
//                tmp[j]     = tmp[j+1];
//                listLengths[j] = listLengths[j+1];
//                tmp[j+1]     = t_tmp;
//                listLengths[j+1] = l_tmp;
//            }
//        }
//    }
    
    // get postingslists set positon vector
    for(int i = 0; i < numKeywords ; i++){
        lists.push_back( (iidx->lists.at(tmp[i]).data()) );
        actualPosition.push_back(0);
    }

    
    // for each doc in first postingLists
    for(int i = 0; i < listLengths[0]; i++){
        int document = lists.at(0)[i].rec->id;
        
        // for each postingsList from 1 ... numkeywords-1
        for(int j = 1; j < numKeywords ; j++){
            //for each doc in postinglist
            for(int l = actualPosition[j]; l < listLengths[j]; l++){
                if( document == lists.at(j)[l].rec->id ){
                    actualPosition[j] = l;
                    break;
                } //lists.at(0)->rec->id
            }
        }
        for(int j = 1; j < numKeywords ; j++){
            if(document != lists.at(j)[actualPosition[j]].rec->id){
                break;
            }
            if(j == numKeywords -1){
                result->push_back(document);
            }
        }
    }

    vector<Record*>* m = new vector<Record*>;
    
    for(int i = 0; i < result->size(); i++ ){
        //(T.recs[(result->at(i))]).Print('s');
        m->push_back(  &(T.recs[(result->at(i))]) );
        output->numRecords++;
    }
    output->subTable = (m->data());
    output->numRecords = m->size();
    return output;
};

inline Relation* ContainmentQueryWithBinSearchOfRelation(Relation& T, int numKeywords, int *keywords, InvertedIndex *iidx){
    
    std::vector<int>* result = new std::vector<int>;
    std::vector<InvertedListEntry*> lists;
    std::vector<int> actualPosition;
    Relation* output = new Relation();
    int listLengths[numKeywords];
    int tmp[numKeywords];
    
    
    //create tmp keyword set
    for(int i = 0; i <= numKeywords; i++){
        tmp[i] = keywords[i];
    }

    //get postingList lengths for keywords
    for(int i = 0; i < numKeywords; i++){
        listLengths[i] = iidx->lists.at(tmp[i]).size();
//        Print(&iidx->lists.at(tmp[i]));
    }
    
    //sort termSet regarding listLengths (bubble sort) // already sorted
    for(int i = 0; i < numKeywords; i++){
        for(int j = 0; j <numKeywords-1-i; j++){
            if(listLengths[j] > listLengths[j+1]){
                int t_tmp = tmp[j];
                int l_tmp = listLengths[j];

                tmp[j]     = tmp[j+1];
                listLengths[j] = listLengths[j+1];
                tmp[j+1]     = t_tmp;
                listLengths[j+1] = l_tmp;
            }
        }
    }
 
    vector<int> pos;
    vector<int> end;
    
    for(int i = 0; i < numKeywords; i++){
        lists.push_back( (iidx->lists.at(tmp[i]).data()) );
//        Print( &iidx->lists.at(tmp[i]));
        pos.push_back(0);
        end.push_back(listLengths[i]-1);
    }
    
    
    
    int start;
    int last;
    int sth = 0;
    int* cur_pos = &sth;
    InvertedListEntry* l;
    int searchkey;
    
    if(numKeywords == 1){
        for(int i = 0; i < listLengths[0]; i++){
            result->push_back(lists.at(0)[i].rec->id);
        }
    }else{
        for(int i = 0; i <= end[0]; i++){
            searchkey = lists.at(0)[i].rec->id;
            
            for(int j = 1; j < numKeywords; j++){
                start = pos[j];
                last = end[j];
                l = lists[j];
    //            cout << "TEST____ ID" << searchkey << "\n";
    //            cout << "start " << start << ", end " << last << "\n";
                if(bin_search(l, start, last, cur_pos, searchkey)){
                    pos[j] = *cur_pos;
                }else{
                    pos[j] = *cur_pos;
                    break;
                }
                if(j == numKeywords -1){
                    result->push_back(searchkey);
                }
            }
        }
    }
    
    vector<Record*>* m = new vector<Record*>;
    
    for(int i = 0; i < result->size(); i++ ){
        //(T.recs[(result->at(i))]).Print('s');
        m->push_back(  &(T.recs[(result->at(i))]) );
        output->numRecords++;
    }
    output->subTable = (m->data());
    output->numRecords = m->size();
    return output;
};

//helper function for bin search containment query
inline bool bin_search(InvertedListEntry* l, int start, int end, int* pos, int x){
    int m;
      while (start <= end) {
          m = start + (end - start) / 2;

          // Check if x is present at mid
//          cout << "TEST: start "<< start << " end " << end << "m " << l[m].rec->id << "\n";
          if (l[m].rec->id == x){
              *pos = m;
              return true;
          }
//              return m;

          // If x greater, ignore left half
          if (l[m].rec->id < x)
              start = m + 1;

          // If x is smaller, ignore right half
          else
              end = m - 1;
      }

      // if we reach here, then element was
      // not present
     *pos = m;
    return false;
}

//bool bin_search(InvertedListEntry* l, int start, int end, int* pos, int x){
//    int m;
//      while (start <= end) {
//          m = start + (end - start) / 2;
//
//          // Check if x is present at mid
////          cout << "TEST: start "<< start << " end " << end << "m " << l[m].rec->id << "\n";
//          if (l[m].rec->id == x){
//              *pos = m;
//              return true;
//          }
////              return m;
//
//          // If x greater, ignore left half
//          if (l[m].rec->id < x)
//              start = m + 1;
//
//          // If x is smaller, ignore right half
//          else
//              end = m - 1;
//      }
//
//      // if we reach here, then element was
//      // not present
//     *pos = m;
//    return false;
//}
inline bool textProbeOnDataPoint(termSet t, termSetSize n, Record* r, float* p){
    int i = 0;
    int j = 0;

    while( i < n && j < r->length){
        if(t[i] < r->keywords[j]){
            *p = -1;

            return false;
        }else if(t[i] > r->keywords[j]){
            j++;
        }else if(t[i] == r->keywords[j]) {
            i++;
            j++;
        }
    }
    if(i == n){
        *p = 1;

        return true;
    }else{
        *p = -1;
        return false;
    }
}

inline bool textProbeOnDataPointBinSearch(termSet t, termSetSize n, Record* r, float* p){

    int start = 0;
    int tmp = 0;
    int end = r->length - 1;
    int* list = r->keywords;
    int m;
    for(int i = 0; i < n; i++){ //for each keyword
        end = r->length - 1;
        while (start <= end) {
            
            m = start + (end - start) / 2;
            
            if (list[m] == t[i]){
                start = m+1;
                break;
            }
            if (list[m] < t[i])
                start = m + 1;
            else
                end = m - 1;
        }
//        start = tmp;
        if (list[m] != t[i]){
            *p = -1;

//            r->Print('dontpass');
            return false;
        }
    }
    *p = 1;
//        r->Print('9');
    return true;
}

inline bool textProbeOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node){
    for(int i = 0; i < numKeywords; i++){
        if(!iidx->ExistsKeyword(keywords[i])){
            node->textScore = -1;
            return false;
        }
    }
    node->textScore = 1;
    return true;
}




inline bool textSearchOnNodeBinSearch(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node){
            
    std::vector<int>* result = new std::vector<int>;
    std::vector<InvertedListEntry*> lists;
    std::vector<int> actualPosition;
    Relation* output = new Relation();
    int listLengths[numKeywords];
    int tmp[numKeywords];
    
    //create tmp keyword set
    for(int i = 0; i < numKeywords; i++){
        if(!iidx->ExistsKeyword(keywords[i])){
            node->textScore = -1;
            return false;
        }
        tmp[i] = keywords[i];
    }
    
//    if(node->is_leaf_node()){
//        Leaf* leaf = (Leaf*)node;
//        TreeDataP<double> **entry = leaf->data;
//        for(int i = 0; i < leaf->entry_num; i++){
//            if(entry[i]->id == 3756){
//                for(int j = 0; j < numKeywords; j++)
//                    Print(&iidx->lists.at(keywords[j]));
//            }
//        }
//    }
        
        


    //get postingList lengths for keywords
    for(int i = 0; i < numKeywords; i++){
        listLengths[i] = iidx->lists.at(tmp[i]).size();
        
//        Print(&iidx->lists.at(tmp[i]));
    }
    
    //sort termSet regarding listLengths (bubble sort) // already sorted
    for(int i = 0; i < numKeywords; i++){
        for(int j = 0; j <numKeywords-1-i; j++){
            if(listLengths[j] > listLengths[j+1]){
                int t_tmp = tmp[j];
                int l_tmp = listLengths[j];

                tmp[j]     = tmp[j+1];
                listLengths[j] = listLengths[j+1];
                tmp[j+1]     = t_tmp;
                listLengths[j+1] = l_tmp;
            }
        }
    }
 
    vector<int> pos;
    vector<int> end;
    
    for(int i = 0; i < numKeywords; i++){
        lists.push_back( (iidx->lists.at(tmp[i]).data()) );
//        Print( &iidx->lists.at(tmp[i]));
        pos.push_back(0);
        end.push_back(listLengths[i]-1);
    }
    
    
    
    int start;
    int last;
    int sth = 0;
    int* cur_pos = &sth;
    InvertedListEntry* l;
    int searchkey;
    for(int i = 0; i <= end[0]; i++){
        searchkey = lists.at(0)[i].rec->id;
        
        for(int j = 1; j < numKeywords; j++){
            start = pos[j];
            last = end[j];
            l = lists[j];
//            cout << "TEST____ ID" << searchkey << "\n";
//            cout << "start " << start << ", end " << last << "\n";
            if(bin_search(l, start, last, cur_pos, searchkey)){
                pos[j] = *cur_pos;
            }else{
                pos[j] = *cur_pos;
                break;
            }
            if(j == numKeywords -1){
                node->verifiedRecords->push_back(searchkey);
            }
        }

    }
    if(node->verifiedRecords->size() == 0){
        node->textScore = -1;
        return false;
    }else{
        node->textScore = 1;
        return true;
    }

}

inline bool textVerificationOnNode(int* keywords, int numKeywords, InvertedIndex* iidx, Node* node){
    if(node->is_leaf_node()){
        return textSearchOnNodeBinSearch(keywords, numKeywords, iidx, (Node*)node);
//        return;
    }
    for(int i = 0; i < numKeywords; i++){
        if(!iidx->ExistsKeyword(keywords[i])){
            node->textScore = - 1;
            return false;
        }
    }
    node->textScore = 1;
    return true;
}


inline bool textSearchOnLeaf(int* keywords, int numKeywords, InvertedIndex* iidx, Node* leaf){
    std::vector<InvertedListEntry*> lists;
    std::vector<int> actualPosition;
    int listLengths[numKeywords];
    int tmp[numKeywords];
    vector<RecordId>* output = new vector<RecordId>;

    
    
    //create tmp keyword set
    for(int i = 0; i < numKeywords; i++){
        if(!iidx->ExistsKeyword(keywords[i])){
            leaf->textScore = -1;
            return false;
        }
        tmp[i] = keywords[i];
    }
    
    
//    if(leaf->is_leaf_node()){
//        Leaf* leaf2 = (Leaf*)leaf;
//        TreeDataP<double> **entry = leaf2->data;
//        for(int i = 0; i < leaf2->entry_num; i++){
//            if(entry[i]->id == 3756){
//                cout << leaf2->page_id << "id\n\n\n";
//                for(int j = 0; j < numKeywords; j++)
//                    Print(&iidx->lists.at(keywords[j]));
//            }
//        }
//    }
    
    //get postingList lengths for keywords
    for(int i = 0; i < numKeywords; i++){
        listLengths[i] = iidx->lists.at(tmp[i]).size();
    }
    
    
    //sort termSet regarding listLengths (bubble sort) // already sorted
    for(int i = 0; i < numKeywords; i++){
        for(int j = 0; j < numKeywords-1-i; j++){
            if(listLengths[j] > listLengths[j+1]){
                int t_tmp = tmp[j];
                int l_tmp = listLengths[j];
                tmp[j]     = tmp[j+1];
                listLengths[j] = listLengths[j+1];
                tmp[j+1]     = t_tmp;
                listLengths[j+1] = l_tmp;
            }
        }
    }
    
    // get postingslists set positon vector
    for(int i = 0; i < numKeywords ; i++){
        lists.push_back( (iidx->lists.at(tmp[i]).data()) );
        actualPosition.push_back(0);
    }

    // for each doc in first postingLists
    
    
    if(numKeywords == 1){
        for(int i = 0; i < listLengths[0]; i++){
            leaf->verifiedRecords->push_back(lists.at(0)[i].rec->id);
        }
    }else{
        for(int i = 0; i < listLengths[0]; i++){
            int document = lists.at(0)[i].rec->id; //XXX
            
            // for each postingsList from 1 ... numkeywords-1
            for(int j = 1; j < numKeywords ; j++){
                //for each doc in postinglist
                for(int l = actualPosition[j]; l < listLengths[j]; l++){
                    if( document == lists.at(j)[l].rec->id ){
                        actualPosition[j] = l;
                        break;
                    } //lists.at(0)->rec->id
                }
            }
            for(int j = 1; j < numKeywords ; j++){
                if(document != lists.at(j)[actualPosition[j]].rec->id){
                    break;
                }
                if(j == numKeywords - 1){
                    leaf->verifiedRecords->push_back(document);
                }
            }
        }
    }
    
    if(leaf->verifiedRecords->size() == 0){
        leaf->textScore = -1;
        return false;
    }else{
        leaf->textScore = 1;
        return true;
    }
}
