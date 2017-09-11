//
//  utils.hpp
//  GRASP
//
//  Created by David Sembroiz on 03/05/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#ifndef utils_hpp
#define utils_hpp

#include <fstream>
#include <vector>
#include <set>
#include <map>
#include "definitions.hpp"
using namespace std;

void read_file(const string name, ifstream & f);
void readDistanceMatrix(ifstream & file, md & distance);
void readTransmissionLevels(ifstream & file, int levels, vd & transEnergy, vd & transRange);
mi fillDistanceMatrix(md & distance, int res);
void saveSolutions(map<set<int>, pair<int, int> > & solutions, int resId, int nodes, int critical, long time);
void saveConnectivity(const vi & conn);
void saveDistanceToBS(const md & distance);
void saveTest(int resultId, const mi & coords, const md & dist);
int getConnectivity(int node, const md & dist);
bool fileExists(const string & filename);

#endif /* utils_hpp */
