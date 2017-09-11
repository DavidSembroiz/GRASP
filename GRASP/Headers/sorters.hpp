//
//  sorters.hpp
//  GRASP
//
//  Created by David Sembroiz on 03/05/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#ifndef sorters_hpp
#define sorters_hpp
#include <vector>
#include "definitions.hpp"
using namespace std;

class ConnSorter {
    vi _conn;
    
public:
    ConnSorter(md distance);
    
    bool operator()(int n1, int n2) const;
};

class DistanceSorter {
    md _distance;
    
public:
    DistanceSorter(md distance);
    
    bool operator()(int n1, int n2) const;
};

class RelaySorter {
    vi _relay;
    
public:
    
    
    RelaySorter(int id);
    
    bool operator()(int n1, int n2) const;
    
private:
    
    vi fillRelay(int id);

};

vi getCriticalNodeSubset(int nVertex, const md & distance, const int RESULT_ID);

#endif /* sorters_hpp */
