    //
//  grasp.cpp
//  GRASP
//
//  Created by David Sembroiz on 04/05/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#include <iostream>
#include <vector>
#include <thread>
#include <set>
#include <map>
#include "Headers/grasp.hpp"
#include "Headers/definitions.hpp"
using namespace std;

/* AUXILIARY FUNCTIONS */

Solution::Solution(double res) {
    _best = res;
    for (int i = 0; i < _bestKey.size(); ++i) _bestKey[i] = -1;
}

Solution::Solution() {
    _best = 0;
}

static bool operator<(const Solution& a1, const Solution& a2) {
    return a1._bestKey < a2._bestKey;
}




int getTransmissionLevel(const double dist, const vd & transRange) {
    for (int i = 0; i < transRange.size(); ++i) {
        if (transRange[i] >= dist) return i;
    }
    return -1;
}

bool inRange(const int i) {
    return i >= 0 && i <= 25;
}

void printSolution(const vi & parent, const vi & key, const double & cost) {
    printf("Edge   Level\n");
    for (int i = 1; i < parent.size(); i++) printf("%d - %d    %d \n", i, parent[i], key[i]);
    cout << "Cost: " << cost << endl;
}


/* ----------------------------- */


/* SOLUTION CONSTRUCTION FUNCTIONS */


double minCost(const vd & cost, const vi & nodeSet) {
    double min = INTMAX_MAX;
    for (int i = 0; i < cost.size(); ++i) if (!nodeSet[i] && cost[i] > 0 && cost[i] < min) min = cost[i];
    return min;
}

double maxCost(const vd & cost, const vi & nodeSet) {
    double max = INTMAX_MIN;
    for (int i = 0; i < cost.size(); ++i) if (!nodeSet[i] && cost[i] > 0 && cost[i] > max) max = cost[i];
    return max;
}

/**
 
 For calculating the cost of adding an element to the solution, the range and the energy is checked. Result vectors
 contain the closest neighbour since the energy for sending the data would be the lowest.
 
 The returned result depends on the objective function under use:
 - 0: return 1 since only one hop is added.
 - 1: return the energy consumed for sending data at the calculated level of transmission.
 
 **/

double costOfAddingNodeToSolution(int node, vi & sol, vi & key, const vi & nodeSet, const md & distance, const vd & transRange, const vd & transEnergy) {
    for (int i = 0; i < sol.size(); ++i) {
        if (nodeSet[i] == 1) {
            int level = getTransmissionLevel(distance[node][i], transRange);
            if (level != -1 && (key[node] == -1 || transEnergy[level] < transEnergy[key[node]])) {
                sol[node] = i;
                key[node] = level;
            }
        }
    }
    if (key[node] == -1) return -1;
    if (OBJECTIVE == 0) return 1;
    return transEnergy[key[node]];
}


/**
 
 Construct solution by constructing a tree in which every new node is connected to a node from the tree
 to which the energy needed to transmit data is minimum.
 
 This solution might not be the best one since a node might be able to transmit more efficiently to a node that is
 already outside the tree.
 
 **/

vi greedyRandomizedSolution(const double alpha, const vi & criticalNodes, const md & distance, const vd & transRange, const vd & transEnergy, vi & key) {
    size_t numNodes = distance.size();
    size_t numCriticalNodes = criticalNodes.size();
    vi sol (numNodes);
    vd cost (numNodes, 0);
    vi nodeSet (numNodes, 0);
    vi rcl;
    sol[0] = -1;
    key[0] = -1;
    nodeSet[0] = 1;
    for (int i = 0; i < numCriticalNodes; ++i) {
        nodeSet[criticalNodes[i]] = -1;
    }
    int setNodes = 1;
    double minC, maxC;
    while (setNodes < numNodes - numCriticalNodes) {
        for (int i = 0; i < numNodes; ++i) {
            if (!nodeSet[i]) {
                cost[i] = costOfAddingNodeToSolution(i, sol, key, nodeSet, distance, transRange, transEnergy);
            }
        }
        minC = minCost(cost, nodeSet);
        maxC = maxCost(cost, nodeSet);
        rcl.clear();
        for (int i = 0; i < numNodes; ++i) {
            if (!nodeSet[i] && cost[i] > 0 && cost[i] <= minC + alpha * (maxC - minC)) {
                rcl.push_back(i);
            }
        }
        if (rcl.size() > 0) {
            int selected = rcl[rand() % rcl.size()];
            nodeSet[selected] = true;
            ++setNodes;
        }
        else {
            cout << "UNABLE TO FIND SOLUTION " << setNodes << endl;
            exit(-1);
        }
        
    }
    return sol;
}

/* ----------------------------- */


/* LOCAL SEARCH FUNCTIONS */

bool isCriticalNode(int node, const vi & criticalNodes) {
    return find(criticalNodes.begin(), criticalNodes.end(), node) != criticalNodes.end();
}

int getRandomNode(size_t size, const vi & criticalNodes) {
    int res = 0;
    do {
        res = rand() % size;
    } while (res == 0 && isCriticalNode(res, criticalNodes));
    return res;
}

/**
 
 Computes the cost of the path from node to Base Station, taking into account the current objective function under use.
 
 **/

double calculatePathCost(int node, const vi & sol, const vi & key, const vd & transEnergy) {
    if (node == 0) return 0;
    double res = 0;
    do {
        res += OBJECTIVE == 0 ? 1 : transEnergy[key[node]];
        node = sol[node];
    } while (node != 0);
    return res;
}

double solutionCost(const vi & sol, const vi & key, const vd & transEnergy, const vi & criticalNodes) {
    double sum = 0;
    for (int i = 1; i < sol.size(); ++i) {
        if (!isCriticalNode(i, criticalNodes)) sum += calculatePathCost(i, sol, key, transEnergy);
    }
    return sum;
}

vi getNeighbors(const vi & sol, int node, const vd & transRange, const vi & criticalNodes, const md & distance) {
    vi neighbors;
    for (int i = 0; i < sol.size(); ++i) {
        if (!isCriticalNode(i, criticalNodes) && i != node && i != sol[node]) {
            int level = getTransmissionLevel(distance[node][i], transRange);
            if (level != -1) neighbors.push_back(i);
            if (neighbors.size() == LOCAL_SEARCH_MAX_NEIGHBORS) return neighbors;
        }
    }
    return neighbors;
}

/**
 
 For improving the solution, firstly a random node is selected (different from Base Station). For all the neighbours in range (for any kind of transmission level),
 the cost of the new path is computed. If the cost is less than current path, it is updated.
 Local search stops when no improvement is present for N subsequent iterations.
 
 **/

void localSearch(vi & sol, vi & key, const vd & transEnergy, const vd & transRange, const md & distance, const vi & criticalNodes) {
    int noImprovement = 0;
    while (noImprovement < LOCAL_SEARCH_NO_IMPROV_THRESHOLD) {
        int node = getRandomNode(sol.size(), criticalNodes);
        double bestCost = calculatePathCost(node, sol, key, transEnergy);
        vi neighbors = getNeighbors(sol, node, transRange, criticalNodes, distance);
        for (int i = 0; i < neighbors.size(); ++i) {
            int level = getTransmissionLevel(distance[node][neighbors[i]], transRange);
            double pathCost = OBJECTIVE == 0 ? 1 : transEnergy[level];
            pathCost += calculatePathCost(neighbors[i], sol, key, transEnergy);
            if (pathCost < bestCost) {
                noImprovement = 0;
                bestCost = pathCost;
                sol[node] = neighbors[i];
                key[node] = level;
            }
        }
        noImprovement++;
    }
}


/* ----------------------------- */

void graspThread(int tid, const int qtt, const vector<vi> & criticalNodes, const int & nVertex, const double & alpha, const md & distance, const vd & transRange, const vd & transEnergy, double *res, vi *sols, vi *keys) {
    
    /**
     
     GRASP START
     
     **/
    
    for (int i = tid; i < tid + qtt; ++i) {
        Solution sol(INTMAX_MAX);
        int iterations = GRASP_ITERATIONS;
        double curCost;
        while(iterations-- > 0) {
            vi key(nVertex, -1);
            vi s = greedyRandomizedSolution(alpha, criticalNodes[i], distance, transRange, transEnergy, key);
            localSearch(s, key, transEnergy, transRange, distance, criticalNodes[i]);
            curCost = solutionCost(s, key, transEnergy, criticalNodes[i]);
            if (curCost < sol._best) {
                sol._best = curCost;
                sol._bestSol = s;
                sol._bestKey = key;
            }
        }
        res[i] = sol._best;
        sols[i] = sol._bestSol;
        keys[i] = sol._bestKey;
    }
}

set<Solution> graspIteration(int nVertex, int numThreads, int qtt, const mi & criticalNodesPerms, const md & distance, const vd & transRange, const vd & transEnergy, thread *t) {
    
    int size = static_cast<int>(criticalNodesPerms.size());
    double *res = new double[size];
    vi *sols = new vi[size];
    vi *keys = new vi[size];
    int resto = size - (numThreads * qtt);
    int start = 0;
    for (int i = 0; i < numThreads; ++i) {
        if (resto-- > 0) {
            t[i] = thread(graspThread,start,qtt + 1,criticalNodesPerms,nVertex,GRASP_ALPHA,distance,transRange,transEnergy, res, sols, keys);
            start += qtt + 1;
        }
        else {
            t[i] = thread(graspThread,start,qtt,criticalNodesPerms,nVertex,GRASP_ALPHA,distance,transRange,transEnergy, res, sols, keys);
            start += qtt;
        }
    }
    
    for (int i = 0; i < numThreads; ++i) {
        t[i].join();
    }
    
    
    /**
     
     Variables to store the best result, nodes with key -1 are Critical Nodes (which are initially removed).
     
     **/
    
    set<Solution> solret;
    
    /**
     
     From all the solutions for all the permutations, the ones with highest latency difference
     are considered the best.
     
     **/
    Solution sol(INTMAX_MIN);
    for (int i = 0; i < size; ++i) {
        if (res[i] > sol._best) solret.clear();
        if (res[i] > sol._best || res[i] == sol._best) {
            sol._best = res[i];
            sol._bestSol = sols[i];
            sol._bestKey = keys[i];
            solret.insert(sol);
        }
    }
    return solret;
}



