//
//  grasp.hpp
//  GRASP
//
//  Created by David Sembroiz on 03/05/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#ifndef grasp_hpp
#define grasp_hpp
#include <thread>
#include "definitions.hpp"
using namespace std;

class Solution {
    
public:
    Solution();
    Solution(double res);
    double _best;
    vi _bestSol;
    vi _bestKey;
};

set<Solution> graspIteration(int nVertex,int numThreads,int qtt,const mi & criticalNodesPerms,const md & distance,const vd & transRange,const vd & transEnergy,thread *t);
#endif /* grasp_hpp */
