#ifndef common_h
#define common_h
#endif /* common_h */
#include <vector>
#include <string>
#include <iostream>
#include <sstream>
using namespace std;
//This file is for common functions
string itos (long a);
void RemoveVectorItem(vector<long> &S, long item);
vector<long> findCommonV(vector<long> &S1, vector<long> &S2);
vector<long> findCommonV(const vector<long> &S1, const vector<long> &S2, vector<bool> &isInCommon);
vector<long> findUnionV(const vector<long> &S1, const vector<long> &S2);
bool IsInVector(const vector<long> &S, long a);

