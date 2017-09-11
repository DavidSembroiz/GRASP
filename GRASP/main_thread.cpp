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
#include <thread>
using namespace std;


#define vi vector<int>
#define mi vector<vi>
#define vd vector<double>
#define md vector<vd>
#define vb vector<bool>
#define mb vector<vb>

/**
 
 Parameter to decide the Objective function to apply:
  - 0: latency (hop count).
  - 1: lifetime (energy consumed).
 
 **/

const int OBJECTIVE = 1;


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

void fillDistanceMatrix(md & distance) {
    ifstream file = read_file("/Users/Sembroiz/Dropbox/DAC/workspace/GRASP/GRASP/coordinates.txt");
    int numNodes;
    file >> numNodes;
    mi coords(numNodes, vi(2));
    for (int i = 0; i < numNodes; ++i) file >> coords[i][0] >> coords[i][1];
    for (int i = 0; i < numNodes; ++i) {
        for (int j = 0; j < numNodes; ++j) {
            distance[i][j] = sqrt((coords[i][0] - coords[j][0])*(coords[i][0] - coords[j][0]) +
                                  (coords[i][1] - coords[j][1])*(coords[i][1] - coords[j][1]));
        }
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

void printSolution(const vi & parent, const vi & key, const double & cost) {
    printf("Edge   Level\n");
    for (int i = 1; i < parent.size(); i++) printf("%d - %d    %d \n", i, parent[i], key[i]);
    cout << "Cost: " << cost << endl;
}

void printCriticalNodes(const vi & keys) {
    for (int i = 1; i < keys.size(); ++i) {
        if (keys[i] == -1) cout << i << " ";
    }
    cout << endl;
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

double solutionCost(const vi & sol, const vi & key, const vd & transEnergy) {
    double sum = 0;
    for (int i = 1; i < sol.size(); ++i) sum += calculatePathCost(i, sol, key, transEnergy);
    return sum;
}

/**
 
 For improving the solution, firstly a random node is selected (different from Base Station). For all the neighbours in range (for any kind of transmission level),
 the cost of the new path is computed. If the cost is less than current path, it is updated.
 Local search stops when no improvement is present for N subsequent iterations.
 
 **/

void localSearch(vi & sol, vi & key, const double & alpha, const vd & transEnergy, const vd & transRange, const md & distance, const vi & criticalNodes) {
    int noImprovement = 0;
    while (noImprovement < 10) {
        int node = getRandomNode(sol.size(), criticalNodes);
        int sendTo = sol[node];
        double bestCost = calculatePathCost(node, sol, key, transEnergy);
        for (int i = 0; i < sol.size(); ++i) {
            if (!isCriticalNode(i, criticalNodes) && i != node && i != sendTo) {
                int level = getTransmissionLevel(distance[node][i], transRange);
                if (level != -1) {
                    double pathCost = OBJECTIVE == 0 ? 1 : transEnergy[level];
                    pathCost += calculatePathCost(i, sol, key, transEnergy);
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

void callFromThread(int tid, const int qtt, const vector<vi> & criticalNodes, const int & nVertex, const double & alpha, const md & distance, const vd & transRange, const vd & transEnergy,
                    double* res, vi* sols, vi* keys) {
    /**
     
     GRASP START
     
     **/
    for (int i = tid; i < tid + qtt; ++i) {
        vi bestSol, bestKey;
        int iterations = 10;
        double best = INTMAX_MAX, curCost;
        while(iterations-- > 0) {
            vi key(nVertex, -1);
            vi sol = greedyRandomizedSolution(alpha, criticalNodes[i], distance, transRange, transEnergy, key);
            localSearch(sol, key, alpha, transEnergy, transRange, distance, criticalNodes[i]);
            curCost = solutionCost(sol, key, transEnergy);
            if (curCost < best) {
                bestSol = sol;
                bestKey = key;
                best = curCost;
            }
        }
        res[i] = best;
        sols[i] = bestSol;
        keys[i] = bestKey;
    }
}

class DistanceSorter {
    md _distance;
    
public:
    DistanceSorter(md distance) {
        _distance = distance;
    }
    
    bool operator()(int n1, int n2) const {
        return _distance[n1][0] < _distance[n2][0];
    }
};

class ConnSorter {
    md _distance;
    
public:
    ConnSorter(md distance) {
        _distance = distance;
    }
    
    int getConnectivity(int node, const md & distance) const {
        int res = 0;
        for (int i = 0; i < distance[node].size(); ++i) {
            if (distance[node][i] <= 82.92) ++res;
        }
        return res;
    }
    
    bool operator()(int n1, int n2) const {
        return getConnectivity(n1, _distance) > getConnectivity(n2, _distance);
    }
};

class RelaySorter {
    vi _relay;
    
public:
    
    vi fillRelay() {
        ifstream file = read_file("/Users/Sembroiz/Dropbox/DAC/workspace/GRASP/GRASP/100nodeSolution.txt");
        int nVertex, node, pkts;
        file >> nVertex;
        vi res(nVertex);
        for (int i = 0; i < nVertex; ++i) {
            file >> node >> pkts;
            res[node] = pkts;
        }
        return res;
    }
    
    RelaySorter() {
        _relay = fillRelay();
    }
    
    bool operator()(int n1, int n2) const {
        return _relay[n1] > _relay[n2];
    }
};



int main(int argc, const char * argv[]) {
    srand(time(NULL));
    //generateDistanceMatrix();
    double alpha = 0;
    int nVertex, transmissionLevels;
    
    /* Distance initialization */
    ifstream file = read_file("/Users/Sembroiz/Dropbox/DAC/workspace/GRASP/GRASP/distances.txt");
    file >> nVertex;
    int numCriticalNodes = 1;
    vi criticalNodes;
    md distance(nVertex, vd (nVertex));
    fillDistanceMatrix(distance);
    //readDistanceMatrix(file, distance);
    
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
    vector<vi> criticalNodesPerms;
    
    /**
     
     CRITICAL NODES SUBSET
     
     According to the following features, only nodes fulfilling the conditions are tested for being critical. The list of features is:
      - Distance to Base Station in meters.
      - Connectivity (i.e. number of nodes in range).
      - Relay (i.e. number of packets to transmit: own + parent nodes packets).
     
     **/
    
    vi distanceRank(nVertex - 1), connRank(nVertex - 1), relayRank(nVertex - 1);
    iota(distanceRank.begin(), distanceRank.end(), 1);
    iota(connRank.begin(), connRank.end(), 1);
    iota(relayRank.begin(), relayRank.end(), 1);
    
    sort(distanceRank.begin(), distanceRank.end(), DistanceSorter(distance));
    sort(connRank.begin(), connRank.end(), ConnSorter(distance));
    sort(relayRank.begin(), relayRank.end(), RelaySorter());
    
    /*for (int i = 0; i < distanceRank.size(); ++i) {
        for (int j = 0; j < distanceRank.size(); ++j) {
            for (int k = 0; k< distanceRank.size(); ++k) {
                if (distanceRank[i] == connRank[j] && distanceRank[i] == relayRank[k]) {
                    cout << "Node " << distanceRank[i] << " : " << i << " " << j  << " " << k << endl;
                }
            }
        }
    }*/
    
    vi criticalNodeSubset;
    int threshold = 15;
    for (int i = 0; i < threshold; ++i) {
        int n = distanceRank[i];
        ptrdiff_t connPos = find(connRank.begin(), connRank.end(), n) - connRank.begin();
        ptrdiff_t relayPos = find(relayRank.begin(), relayRank.end(), n) - relayRank.begin();
        if (connPos < threshold || relayPos < threshold) {
            criticalNodeSubset.push_back(n);
        }
    }
    
    vb perms(criticalNodeSubset.size());
    fill(perms.begin(), perms.begin() + min(numCriticalNodes, (int) criticalNodeSubset.size()), true);
    
    
    
    /**
     
     CRITICAL NODES PERMUTATIONS
     
     We create permutations of size numCriticalNodes to remove them from the graph in order to calculate best solution.
     
     **/
    
    
    do {
        criticalNodes.clear();
        for (int i = 0; i < perms.size(); ++i) {
            if (perms[i]) {
                criticalNodes.push_back(criticalNodeSubset[i]);
            }
        }
        criticalNodesPerms.push_back(criticalNodes);
    } while (prev_permutation(perms.begin(), perms.end()));
    
    int numThreads = min(256, (int) criticalNodesPerms.size());
    thread *t = new thread[numThreads];
    double *res = new double[numThreads];
    vi *sols = new vi[numThreads];
    vi *keys = new vi[numThreads];
    int qtt = (int) criticalNodesPerms.size() / numThreads;
    int resto = (int) criticalNodesPerms.size() - numThreads * qtt;
    
    for (int i = 0; i < numThreads; ++i) {
        if (resto > 0) --resto, t[i] = thread(callFromThread, i, qtt + 1, criticalNodesPerms, nVertex, alpha, distance, transRange, transEnergy, res, sols, keys);
        else t[i] = thread(callFromThread, i, qtt, criticalNodesPerms, nVertex, alpha, distance, transRange, transEnergy, res, sols, keys);
    }
    
    for (int i = 0; i < numThreads; ++i) {
        t[i].join();
    }
    
    for (int i = 0; i < numThreads; ++i) {
        if (res[i] < bestGlobal) {
            bestGlobal = res[i];
            bestGlobalSol = sols[i];
            bestGlobalKey = keys[i];
        }
    }
    printCriticalNodes(bestGlobalKey);
    return 0;
}


