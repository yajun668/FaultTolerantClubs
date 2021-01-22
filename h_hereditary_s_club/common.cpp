#include <stdio.h>
#include "common.h"
#include <algorithm>
//This file is for common functions
//function: convert long to string
string itos (long a)
{
    std::ostringstream str;
    str<<a;
    return str.str();
}

// function: remove an item if it exists in S.
void RemoveVectorItem(vector<long> &S, long item)

{
    vector<long>::iterator it = find(S.begin(), S.end(), item);
    if (it!=S.end())
        S.erase(it);
}

// find common vertices between two sorted subsets
vector<long> findCommonV(vector<long> &S1, vector<long> &S2)
{
    /* Must sort 2 subsets before call this function */
    sort(S1.begin(),S1.end());
    sort(S2.begin(),S2.end());
    vector<long> commonV;
    if (S1.empty() || S2.empty()) {
        return commonV;
    }

    long pos1=0;
    long pos2 =0;
    long p1,p2;
    long end1 = S1.size(), end2 = S2.size();
    while(pos1<end1 && pos2<end2)
    {
        p1 = S1[pos1];
        p2 = S2[pos2];
        if (p1 == p2) {
            commonV.push_back(p1);
            pos1++;
            pos2++;
        }
        else if (p1>p2)
            pos2++;
        else
            pos1++;
    }
    return commonV;
}

// find common vertices between 2 sorted subsets; set isInCommon to true if it is in common V
vector<long> findCommonV(const vector<long> &S1, const vector<long> &S2, vector<bool> &isInCommon)
{
    /* Must sort 2 subsets before call this function */
    vector<long> commonV;
    if (S1.empty() || S2.empty()) {
        return commonV;
    }

    long pos1=0;
    long pos2 =0;
    long p1,p2;
    long end1 = S1.size(), end2 = S2.size();
    while(pos1<end1 && pos2<end2)
    {
        p1 = S1[pos1];
        p2 = S2[pos2];
        if (p1 == p2) {
            commonV.push_back(p1);
            isInCommon[p1] = true;
            pos1++;
            pos2++;
        }
        else if (p1>p2)
            pos2++;
        else
            pos1++;
    }
    return commonV;
}

// find union vertices between two subsets. Return a sorted set
vector<long> findUnionV(const vector<long> &S1, const vector<long> &S2)
{
    vector<long> unionV;
    for (long i = 0; i < S1.size(); i++)
        unionV.push_back(S1[i]);
    for (long j = 0; j < S2.size(); j++) {
        if (!IsInVector(unionV, S2[j]))
            unionV.push_back(S2[j]);
    }
    sort(unionV.begin(), unionV.end());
    return unionV;
}

//check if a item in an vector S
bool IsInVector(const vector<long> &S, long a)
{
    if(find(S.begin(),S.end(),a)!=S.end())
        return true;
    else
        return false;
}