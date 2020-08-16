//
//  utils_sim.cpp
//  TkRQ
//
//  Created by Artur Titkov on 17.07.20.
//  Copyright © 2020 Artur Titkov. All rights reserved.
//

#include "def.h"
#include "relation.h"
#include "inverted_index.h"

class Relation;
class InvertedListEntry;


typedef pair<float, float> boundedDoc;

typedef tuple<vector<int>, float, float> test;





Relation* textSearch(int numkeywords, WeightedTerm* terms, Relation* rel, InvertedIndex* iix){
    
    std::vector<InvertedListEntry*> lists;
    Relation* output = new Relation();
    int actualPosition[numkeywords];
    int listLengths[numkeywords];

    
    for(int i = 0; i < numkeywords; i++){ //INIT vectors
        int termId = terms[i].first;
        actualPosition[i] = 0;
        listLengths[i] = iix->lists.at(termId).size();
        lists.push_back( (iix->lists.at(termId).data()) );
    }
    
    map<int, float> weightedRecords;
    int recid;
    for(int i = 0; i < numkeywords; i++){
        InvertedListEntry* list = lists[i];
        float qTermWeight = terms[i].second;
        
        for(int j = 0; j < listLengths[i]; j++){
            
            recid = list[j].rec->id;
            float w2 = (float)list[j].weight;
            auto entry = weightedRecords.find(recid);
            if(entry == weightedRecords.end()){ // if not in map
                weightedRecords.insert(make_pair(recid, qTermWeight * w2));
            }else{
                float oldW = entry->second;
//                cout << "F " << qTermWeight*w2;
                entry->second = oldW + qTermWeight*w2;
            }
            
        }
        
    }
    WeightedRecord* tmpRel = new WeightedRecord[weightedRecords.size()];
    int counter = 0;
    for (std::map<int,float>::iterator it=weightedRecords.begin(); it!=weightedRecords.end(); ++it){
        int recid = it->first;
        Record* rec = (*rel)[recid];
        tmpRel[counter] = make_pair(it->second, rec);
        counter++;
    }
    output->weightedRecords = tmpRel;
    output->numRecords = weightedRecords.size();
    
//    for(int i = 0; i < output->numRecords; i++){
//        cout << "score " << output->weightedRecords[i].first;
//    }
    
    return output;
}

struct termDoc{
    int docID;
    float lowerScore;
    float upperScore;
    vector<int> terms;
    
    
    termDoc(int id,float l, float u, int term){
        docID = id;
        lowerScore = l;
        upperScore = u;
        terms.push_back(term);
    }
    
    
    
};
static bool Comparator2(termDoc* l, termDoc* r){
    return (*l).lowerScore > (*r).lowerScore;
};

Relation* topkTextSearch(int k, int numkeywords, int* kw, WeightedTerm* t, Relation* rel, InvertedIndex* iix){
    
    Relation* output = new Relation();
    vector<InvertedList::iterator> actualPosition;
    vector<InvertedList::iterator> last;
    float* thresholds = new float[numkeywords];
    float totalThreshold   = 0.0;
    float totalLower = 0.0;
    float totalUpper = 0.0;
    //print inverted lists
    int* keywords = new int[numkeywords];
    bool endConditions = true;
    int round = 0;
    
    //COPY terms
    WeightedTerm* Wterms = new WeightedTerm[numkeywords];
    for(int i = 0; i < numkeywords; i++){
        WeightedTerm* tmp = &t[i];
        Wterms[i] = make_pair(tmp->first, tmp->second);
        keywords[i] = kw[i];
    }
    
//    for(int i = 0; i < numkeywords; i++){
//        cout << "k " << keywords[i] <<" : ";
//        for (InvertedList::iterator iterL = iix->lists.at(keywords[i]).begin(); iterL !=  iix->lists.at(keywords[i]).end(); ++iterL)
//            cout << " <r" << iterL->rec->id << ", " << iterL->position << ", " << iterL->weight << ">";
//        cout << endl;
//    }
    
    
    for(int i = 0; i < numkeywords; i++){ // get iterator for inverted lists
        actualPosition.push_back( iix->lists.at(keywords[i]).begin()) ;
        last.push_back( iix->lists.at(keywords[i]).end());
        thresholds[i] = actualPosition[i]->weight * Wterms[i].second;
        totalThreshold += thresholds[i];
        iix->lists.at(Wterms[i].second);
        
    }
    
    map<int, termDoc> candidates;
    
    vector<termDoc*> sortCand;
    
    int nextdoc = 0;
    float nextdocImpact = 0;
    int recid;

    while (endConditions) {
        round++;
//        sortCand.clear();
        
        
        //find term with max impact
        nextdoc = 0;
        nextdocImpact = 0;
        int termI = 0;
        for(int i = 0; i < numkeywords; i++){ // find max impact
            if( thresholds[i] == max( thresholds[i], nextdocImpact ) ){
                nextdoc = actualPosition[i]->rec->id;
                nextdocImpact = thresholds[i];
                termI = i;
            }
        }
        int termID = keywords[termI];
        
//        cout << "ITERQATION\n";
        //UPDATE thresholds
        actualPosition[termI] += 1;
        if(actualPosition[termI] == last[termI]){ // no elements in list remaining
            //swap terms
            Wterms[termI] = Wterms[numkeywords-1];
            //swap thresholds
            thresholds[termI] = thresholds[numkeywords-1];
            //swap actualPos
            actualPosition[termI] = actualPosition[numkeywords-1];
            //swap last
            last[termI] = last[numkeywords-1];
            //UPDATE threshold
            totalThreshold -= nextdocImpact;
            //decrement numkeywords
            numkeywords -= 1;
            
            //UPDATE UPPERBOUNDS
            map<int, termDoc>::iterator iter = candidates.begin();
            for(; iter != candidates.end(); iter++){
                termDoc* curTermDoc = &(*iter).second;
                vector<int>* curterms = &curTermDoc->terms;
                bool found = false;
                for(int t = 0; t < curterms->size(); t++){
                    if(Wterms[termI].first == curterms->at(t)){
                        found = true;
                        break;
                    }
                }
                if(!found){
                    curTermDoc->upperScore -= nextdocImpact;
                }
            }
        }else{ // else only update thresholds and all candidates
            float a = Wterms[termI].second;
            auto b = actualPosition[termI]->weight;
            thresholds[termI] =  a * b;
            totalThreshold -= nextdocImpact;
            totalThreshold += thresholds[termI];
            
            map<int, termDoc>::iterator iter = candidates.begin();
            for(; iter != candidates.end(); iter++){ // update all candidates
                termDoc* curTermDoc = &(*iter).second;
                vector<int>* curterms = &curTermDoc->terms;
                bool found = false;
                for(int t = 0; t < curterms->size(); t++){
                    if(Wterms[termI].first == termID){
                        found = true;
                        break;
                    }
                }
                if(!found){
                    curTermDoc->upperScore -= nextdocImpact;
                    curTermDoc->upperScore += thresholds[termI];
                }

            }
        }
        
        auto entry = candidates.find(nextdoc); // find doc in candidates
        
        if(entry == candidates.end()){ // if not in candidates
            float lower = nextdocImpact;
            float upper = lower;
            for(int j = 0; j < numkeywords; j++){
//                if(j != termI){
                upper += thresholds[j];
//                }
            }
            candidates.insert({nextdoc, termDoc(nextdoc, lower, upper, termID)});
            termDoc* tmp = &candidates.find(nextdoc)->second;
            sortCand.push_back( tmp );
        }else{ // already in list update scores
            termDoc* doc     = &entry->second;
            doc->terms.push_back(termID);
            doc->lowerScore += nextdocImpact;
        }
        
//        //creatSortedList
//        map<int, termDoc>::iterator iter = candidates.begin();
//        map<int, termDoc>::iterator iterEnd = candidates.end();
//        for(; iter != iterEnd; iter++ ){
//            termDoc* tmp = &iter->second;
//            sortCand.push_back(tmp);
//        }
//        if(round == 100){
//        cout << "thres " << totalThreshold << endl;
//        cout << "Round " << round << endl;
        if( round > 10*k && sortCand.size() > k  ){
            cout << "BREAK?\n";
            round = 0;
            sort(sortCand.begin(), sortCand.end(), Comparator2);

            if(sortCand[k]->lowerScore >= totalThreshold ){
                cout << "\n1.step BREAK\n";
                for(int i = 0; i < k; i++){
                    if(sortCand[i]->lowerScore < sortCand[i+1]->upperScore)
                        break;
                }
                float kthScore = sortCand[k]->lowerScore;
                for(int i = k + 1; i < sortCand.size(); i++){
                    if( kthScore < sortCand[i]->upperScore )
                        break;
                }
                cout << "\nBREAK\n";
                endConditions = false;
            }
        }
//            round = 0;
//        }
        if(numkeywords == 0 )
            break;
//        if()
    }
    
    //OUTPUT
    vector<WeightedRecord>* mm = new vector<WeightedRecord>;
    
    for(int i = 0; i < k; i++){
        float s = sortCand[i]->lowerScore;
        Record* r = &(rel->recs[ sortCand[i]->docID ]);
        mm->push_back( make_pair(s,r ) );
//        mm.push_back(make_pair( iter->first, &rel[ iter->first ] ));
//        tmpRel[recid] = make_pair( iter->first, &rel[ iter->first ] );
    
    }
    output->numRecords = mm->size();
    output->weightedRecords = mm->data();
    return output;
    
//    //OUTPUT
//    vector<WeightedRecord>* mm = new vector<WeightedRecord>;
//    map<int, termDoc>::iterator iter = candidates.begin();
//    for(; iter != candidates.end(); iter++){
//        float s = iter->second.lowerScore;
//        Record* r = &(rel->recs[ iter->first ]);
//        mm->push_back( make_pair(s,r ) );
////        mm.push_back(make_pair( iter->first, &rel[ iter->first ] ));
////        tmpRel[recid] = make_pair( iter->first, &rel[ iter->first ] );
//
//    }
//    output->numRecords = mm->size();
//    output->weightedRecords = mm->data();
//    return output;
}

//Relation* topkTextSearch(int k, int numkeywords, WeightedTerm* terms, Relation* rel, InvertedIndex* iix){
//
//    std::vector<InvertedListEntry*> lists;
//    InvertedListEntry* l = new InvertedListEntry[numkeywords];
//    Relation* output = new Relation();
//    int actualPosition[numkeywords];
//    int listLengths[numkeywords];
//    float* thresholds = new float[numkeywords];
//    float threshold   = 0.0;
//    for(int i = 0; i < numkeywords; i++){ //INIT vectors
//        int termId = terms[i].first;
//        actualPosition[i] = 0;
//        listLengths[i] = iix->lists.at(termId).size();
//        lists.push_back( (iix->lists.at(termId).data()) );
////        l[i] = &lists.at(i);
//        thresholds[i] = lists[i][0].weight * terms[i].second;
//        threshold += thresholds[i];
//    }
//
//
//
//    map<int, termDoc> candidates;
//
//
//
//
//    int nextdoc = 0;
//    float nextdocImpact = 0;
//    int recid;
//
//    while (true) {
//        //find term with max impact
//        nextdoc = 0;
//        nextdocImpact = 0;
//        int i = 0;
//        for(; i < numkeywords; i++){ // find max impact
//            if( thresholds[i] == max( thresholds[i], nextdocImpact ) ){
//                nextdoc = actualPosition[i];
//                nextdocImpact = thresholds[i];
//            }
//        }
//        auto entry = candidates.find(nextdoc); // find doc in candidates
//
//        if(entry == candidates.end()){ // if not in candidates
//            float lower = nextdocImpact;
//            float upper = lower;
//            for(int j = 0; j < numkeywords; j++){
//                if(j != i){
//                    upper += thresholds[i];
//                }
//            }
//            candidates.insert({nextdoc, termDoc(lower, upper, terms[i].second)});
//        }else{ // already in list update scores
//            termDoc* doc = &entry->second;
//            doc->terms.push_back(i);
//            doc->lowerScore += nextdocImpact;
////            doc->upperScore += nextdocImpact;
////            doc->upperScore -= thresholds[i];
//        }
//        //UPDATE THRESHOLDS
//        actualPosition[i] += 1;
//        InvertedListEntry* tmpList = lists.at(i);
//        thresholds[i] = lists[i][actualPosition[i]].weight * terms[i].second;
//
//
//
//        //UPDATE UPPERSCORES;
//        map<int, termDoc>::iterator iter = candidates.begin();
//        for(; iter != candidates.end(); iter++){
//            termDoc* curTermDoc = &(*iter).second;
//            vector<int>* curterms = &curTermDoc->terms;
//            bool found = false;
//            for(int t = 0; t < curterms->size(); t++){
//                if(terms[i].first == curterms->at(t)){
//                    found = true;
//                    break;
//                }
//            }
//            if(!found){
//                curTermDoc->upperScore -= nextdocImpact;
//                curTermDoc->upperScore += thresholds[i];
//            }
//        }
//
//
////        actualPosition[nextdoc] += 1;
////        thresholds[nextdoc] =
////            lists[nextdoc][actualPosition[nextdocImpact]].weight
////            * terms[nextdoc].second;
////
////        auto entry = weightedRecords.find(recid);
////        if(entry == weightedRecords.end()){ // if not in map
////
////        }else{
////
////        }
//
//
//
//
//    }
////
//
//
//
////
////
////    for(int i = 0; i < numkeywords; i++){
////        termMaxImpactId = 0;
////        termMaxImpact   = 0;
////        for(int j = 0; j < numkeywords; i++){
//////            float curW = lists[i][actualPosition[j]].weight;
////
//////            if( curList-> == max(  ) )
////
////
////        }
////
////
////
////        InvertedListEntry* list = lists[i];
////        float qTermWeight = terms[i].second;
////
////        for(int j = 0; j < listLengths[i]; j++){
////
////            recid = list[j].rec->id;
////            float w2 = (float)list[j].weight;
////            auto entry = weightedRecords.find(recid);
////            if(entry == weightedRecords.end()){ // if not in map
////                weightedRecords.insert(make_pair(recid, qTermWeight * w2));
////            }else{
////                float oldW = entry->second;
//////                cout << "F " << qTermWeight*w2;
////                entry->second = oldW + qTermWeight*w2;
////            }
////        }
////    }
////
////
////
////
////    WeightedRecord* tmpRel = new WeightedRecord[weightedRecords.size()];
////    int counter = 0;
////    for (std::map<int,float>::iterator it=weightedRecords.begin(); it!=weightedRecords.end(); ++it){
////        int recid = it->first;
////        Record* rec = (*rel)[recid];
////        tmpRel[counter] = make_pair(it->second, rec);
////        counter++;
////    }
////    output->weightedRecords = tmpRel;
////    output->numRecords = weightedRecords.size();
////
//////    for(int i = 0; i < output->numRecords; i++){
//////        cout << "score " << output->weightedRecords[i].first;
//////    }
////
//    return output;
//}

