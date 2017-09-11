//
//  main.cpp
//  GRASP
//
//  Created by David Sembroiz on 27/03/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#include <iostream>
#include <ctime>
#include <set>
#include <map>
#include <random>
#include <thread>
#include <chrono>
#include "Headers/sorters.hpp"
#include "Headers/utils.hpp"
#include "Headers/generator.hpp"
#include "Headers/grasp.hpp"
using namespace std;


int RESULT_ID = 62;

/* Number of critical nodes to extract from the initial topology */

int NUM_CRITICAL_NODES = 0;

int NUM_NODES = 2000;

int NET_RADIUS = 900;


void printCriticalNodes(const vi & keys) {
    for (int i = 1; i < keys.size(); ++i) {
        if (keys[i] == -1) cout << i << " ";
    }
    cout << endl;
}

double calculateRelay(int node, const vi & sol) {
    if (node == 0) return 0;
    double res = 0;
    do {
        ++res;
        node = sol[node];
    } while (node != 0);
    return res;
}

void saveRelay(set<int> dominantCritNodes, map<set<int>, Solution> relays, const md & distance) {
    /**
     
     Get the vector bestSol. For every node,
     
     **/
    string name = GRASP_PATH + "Inputs/solution_" + to_string(RESULT_ID) + ".txt";
    if (fileExists(name)) return;
    map<set<int>, Solution>::iterator it = relays.find(dominantCritNodes);
    if (it != relays.end()) {
        vi rels = it->second._bestSol;
        vi keys = it->second._bestKey;
        vi relCount(rels.size(), 0);
        for (int i = 1; i < rels.size(); ++i) {
            cout << i << " to " << rels[i] << " " << distance[i][rels[i]] << " level " << keys[i] << endl;
            if (keys[i] != -1) {
                int node = i;
                do {
                    ++relCount[node];
                    node = rels[node];
                } while (node != 0);
            }
        }
        int total = 0;
        for (int i = 0; i < relCount.size(); ++i) {
            total += relCount[i];
        }
        ofstream file(name);
        file << relCount.size() << endl;
        file << "1 0 0" << endl;
        for (int i = 1; i < relCount.size(); ++i) {
            file << i + 1 << " " << rels[i] << " " << relCount[i] << endl;
        }
    }
}

set<int> getDominantCriticalNodeSet(map<set<int>, pair<int, int> > & solutions) {
    set<int> dominants;
    int times = 0;
    for (map<set<int>, pair<int, int> >::iterator it = solutions.begin(); it != solutions.end(); ++it) {
        if (it->second.first > times) {
            times = it->second.first;
            dominants = it->first;
        }
    }
    return dominants;
}


int main(int argc, const char * argv[]) {
    
    if (argc == 3) {
        RESULT_ID = atoi(argv[1]);
        NUM_CRITICAL_NODES = atoi(argv[2]);
        
    }
    cout << "Running GRASP for Result ID " << RESULT_ID << " and " << NUM_CRITICAL_NODES << " Critical Nodes" << endl << endl;
    
    srand(static_cast<unsigned int>(time(NULL)));
    if (GEN_NEW_NETWORK) generateCoordinates(RESULT_ID, NUM_NODES, NET_RADIUS);
    
    
    /* Distance initialization */
    
    md distance;
    mi coords = fillDistanceMatrix(distance, RESULT_ID);
    
    /* Transmission levels initialization */
    int nVertex = static_cast<int>(coords.size()), transmissionLevels;
    string name = GRASP_PATH + "Resources/transmissions.txt";
    ifstream file;
    read_file(name, file);
    file >> transmissionLevels;
    vd transEnergy(transmissionLevels), transRange(transmissionLevels);
    readTransmissionLevels(file, transmissionLevels, transEnergy, transRange);
    
    if (NUM_CRITICAL_NODES == 0) saveTest(RESULT_ID, coords, distance);
    
    vi criticalNodeSubset;
    if (NUM_CRITICAL_NODES != 0) criticalNodeSubset = getCriticalNodeSubset(nVertex, distance, RESULT_ID);
    
    vb perms(criticalNodeSubset.size());
    fill(perms.begin(), perms.begin() + min(NUM_CRITICAL_NODES, (int) criticalNodeSubset.size()), true);
    
    vi criticalNodes;
    vi occurrences(nVertex);
    
    
    /**
     
     CRITICAL NODES PERMUTATIONS
     
     We create permutations of size numCriticalNodes to remove them from the graph in order to calculate best solution.
     
     **/
    
    mi criticalNodesPerms;
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
    
    map<set<int>, pair<int, int> > solutions;
    map<set<int>, Solution> relays;
    
    int qtt = (int) criticalNodesPerms.size() / numThreads;
    auto elapsed = 0;
    for (int i = 0; i < GLOBAL_GRASP_ITERATIONS; ++i) {
        auto start = chrono::steady_clock::now();
        cout << "GRASP: " << i << endl;
        set<Solution> s = graspIteration(nVertex, numThreads, qtt, criticalNodesPerms, distance, transRange, transEnergy, t);
        int cost = 0;
        for (set<Solution>::iterator it = s.begin(); it != s.end(); ++it) {
            set<int> sol;
            for (int i = 1; i < it->_bestKey.size(); ++i) {
                if (it->_bestKey[i] == -1) sol.insert(i);
            }
            map<set<int>, pair<int, int> >::iterator i = solutions.find(sol);
            if (i != solutions.end()) {
                i->second.first++;
                if (it->_best < i->second.second) {
                    i->second.second = it->_best;
                    relays[sol] = *it;
                }
            }
            else {
                solutions.insert(make_pair(sol, make_pair(1, it->_best)));
                relays.insert(make_pair(sol, *it));
            }
            printCriticalNodes(it->_bestKey);
            cost = it->_best;
        }
        cout << "Cost: " << cost << endl;
        auto end = chrono::steady_clock::now();
        elapsed = chrono::duration_cast<chrono::seconds>(end - start).count();
        cout << "Elapsed time: " << elapsed << " s"<< endl << endl;
    }
    
    if (NUM_CRITICAL_NODES == 0) saveRelay(getDominantCriticalNodeSet(solutions), relays, distance);
    saveSolutions(solutions, RESULT_ID, nVertex, (int) criticalNodes.size(), elapsed);
    
    return 0;
}


