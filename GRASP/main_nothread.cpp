//
//  main.cpp
//  GRASP
//
//  Created by David Sembroiz on 27/03/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#include <iostream>
#include <algorithm>
#include <ctime>
#include <fstream>
#include <vector>
#include <set>
#include <random>
using namespace std;


#define vi vector<int>
#define mi vector<v>
#define vd vector<double>
#define md vector<vd>
#define vb vector<bool>
#define mb vector<vb>


/* INPUT DATA FUNCTIONS */

ifstream read_file(const char* name) {
    const char* filename  = name;
    ifstream file(filename);
    if (!file) {
        cerr << "ERROR: could not open file '" << filename << "' for reading" << endl;
        throw(-1);
    }
    return file;
}

void readDistanceMatrix(ifstream & file, md & distance) {
    size_t nVertex = distance.size();
    for (int i = 0; i < nVertex; ++i) {
        for (int j = 0; j < nVertex; ++j) file >> distance[i][j];
    }
}

void readTransmissionLevels(ifstream & file, int levels, vd & transEnergy, vd & transRange) {
    for (int i = 0; i < levels; ++i) file >> transEnergy[i] >> transRange[i];
}

/* ----------------------------- */


/* AUXILIARY FUNCTIONS */

void generateDistanceMatrix() {
    ofstream file("/Users/Sembroiz/Dropbox/DAC/workspace/GRASP/GRASP/distances.txt");
    int numNodes = 100;
    int min = 15, max = 82;
    md distances(numNodes, vd (numNodes, 0));
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < i; ++j) {
            double dist = (max - min) * ((double)rand()/(double) RAND_MAX) + min;
            distances[i][j] = dist;
            distances[j][i] = dist;
        }
    }
    file << numNodes << endl;
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            file << distances[i][j] << " ";
        }
        file << endl;
    }
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

/* ----------------------------- */


/* SOLUTION CONSTRUCTION FUNCTIONS */

void printSolution(const vi & parent, const vi & key, const double & cost) {
    printf("Edge   Level\n");
    for (int i = 1; i < parent.size(); i++) printf("%d - %d    %d \n", i, parent[i], key[i]);
    cout << "Cost: " << cost << endl;
}

/*
int minKey(const vi & key, const vb & nodeSet) {
    int min = INT_MAX, minIndex = 0;
    for (int i = 0; i < key.size(); ++i) {
        if (!nodeSet[i] && key[i] < min) {
            min = key[i];
            minIndex = i;
        }
    }
    return minIndex;
}


vi constructSolution(const vi & criticalNodes, const md & distance, const double range, const vd & transRange, const vd & transEnergy, vi & key) {
    size_t nVertex = distance.size();
    vi parent(nVertex, -1);
    vector<bool> nodeSet(nVertex, false);
    key[0] = -1;
    for (int i = 0; i < nVertex - 1; ++i) {
        int id = minKey(key, nodeSet);
        nodeSet[id] = true;
        for (int j = 0; j < nVertex; ++j) {
            int level = getTransmissionLevel(distance[id][j], transRange);
            if (!nodeSet[j] && level != -1 && transEnergy[level] < transEnergy[key[j]]) {
                parent[j] = id;
                key[j] = level;
            }
        }
    }
    return parent;
}
*/


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
            if (!nodeSet[i] && cost[i] <= minC + alpha * (maxC - minC)) {
                rcl.push_back(i);
            }
        }
        int selected = rcl[rand() % rcl.size()];
        nodeSet[selected] = true;
        ++setNodes;
        
    }
    return sol;
}

/* ----------------------------- */


/* LOCAL SEARCH FUNCTIONS */

int getRandomNode(size_t size) {
    int res = 0;
    do {
        res = rand() % size;
    } while (res == 0);
    return res;
}

double calculatePathCost(int node, const vi & sol, const vi & key, const vd & transEnergy) {
    if (node == 0) return 0;
    double res = 0;
    do {
        res += transEnergy[key[node]];
        node = sol[node];
    } while (node != 0);
    return res;
}

double solutionCost(const vi & sol, const vi & key, const vd & transEnergy) {
    double sum = 0;
    for (int i = 1; i < sol.size(); ++i) sum += calculatePathCost(i, sol, key, transEnergy);
    return sum;
}

/**
 For improving the solution, firstly a node is selected (different from Base Station). For all the neighbours in range (for any kind of transmission level),
 the cost of the new path is computed. If the cost is less than current path, it is updated.
 Local search stops when no improvement is present for N subsequent iterations.
 **/
void localSearch(vi & sol, vi & key, const double & alpha, const vd & transEnergy, const vd & transRange, const md & distance) {
    int noImprovement = 0;
    while (noImprovement < 50) {
        int node = getRandomNode(sol.size());
        int sendTo = sol[node];
        double bestCost = calculatePathCost(node, sol, key, transEnergy);
        for (int i = 0; i < sol.size(); ++i) {
            if (i != node && i != sendTo) {
                int level = getTransmissionLevel(distance[node][i], transRange);
                if (level != -1) {
                    double pathCost = calculatePathCost(i, sol, key, transEnergy) + transEnergy[level];
                    if (pathCost < bestCost) {
                        noImprovement = 0;
                        bestCost = pathCost;
                        sol[node] = i;
                        key[node] = level;
                    }
                }
            }
        }
        noImprovement++;
    }
}


/* ----------------------------- */


int main(int argc, const char * argv[]) {
    srand(time(NULL));
    generateDistanceMatrix();
    double alpha = 0.5;
    int nVertex, transmissionLevels;
    
    /* Distance initialization */
    ifstream file = read_file("/Users/Sembroiz/Dropbox/DAC/workspace/GRASP/GRASP/distances.txt");
    file >> nVertex;
    vb perms(nVertex - 1);
    int numCriticalNodes = 2;
    vi criticalNodes;
    fill(perms.begin(), perms.begin() + numCriticalNodes, true);
    md distance(nVertex, vd (nVertex));
    readDistanceMatrix(file, distance);
    
    /* Transmission levels initialization */
    file = read_file("/Users/Sembroiz/Dropbox/DAC/workspace/GRASP/GRASP/transmissions.txt");
    file >> transmissionLevels;
    vd transEnergy(transmissionLevels), transRange(transmissionLevels);
    readTransmissionLevels(file, transmissionLevels, transEnergy, transRange);
    
    /**
     
     Variables to store the best result, nodes with key -1 are Critical Nodes (which are initially removed).
     
     **/
    
    vi bestGlobalSol, bestGlobalKey;
    double bestGlobal = INTMAX_MAX;
    
    /**
     
     CRITICAL NODES PERMUTATIONS
     
     We create permutations of size numCriticalNodes to remove them from the graph in order to calculate best solution.
     
     **/
    
    
    do {
        criticalNodes.clear();
        for (int i = 0; i < perms.size(); ++i) {
            if (perms[i]) {
                criticalNodes.push_back(i + 1);
                cout << i + 1 << " ";
            }
        }
        cout << endl;
        /**
         
         GRASP START
         
         **/
        
        vi bestSol, bestKey;
        int iterations = 10;
        double best = INTMAX_MAX, curCost;
        while(iterations-- > 0) {
            vi key(nVertex, -1);
            vi sol = greedyRandomizedSolution(alpha, criticalNodes, distance, transRange, transEnergy, key);
            localSearch(sol, key, alpha, transEnergy, transRange, distance);
            curCost = solutionCost(sol, key, transEnergy);
            if (curCost < best) {
                bestSol = sol;
                bestKey = key;
                best = curCost;
            }
        }
        //printSolution(bestSol, bestKey, best);
        if (best < bestGlobal) {
            bestGlobal = best;
            bestGlobalSol = bestSol;
            bestGlobalKey = bestKey;
        }
    } while (prev_permutation(perms.begin(), perms.end()));
    return 0;
}
