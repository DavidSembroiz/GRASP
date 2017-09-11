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

struct node {
    /* Node identifier */
    int id;
    /* Number of input connections */
    int input_connections;
    /* Transmission level used to send data */
    int transmission_level;
    /* Identifier of the node to send data */
    int send_to;
};

static bool operator<(const node& a1, const node& a2) {
    return a1.id < a2.id;
}

struct transmission {
    /* Level identifier */
    double level;
    /* Ratio between energy consumed and covering range */
    double ratio;
    /* Covering distance */
    double range;
    /* Energy consumed */
    double energy;
};

#define v vector<int>
#define m vector<v>
#define vd vector<double>
#define md vector<vd>
#define vb vector<bool>
#define mb vector<vb>
#define nodes set<node>
#define vtrans vector<transmission>


int minKey(const vd & key, const vb & nodeSet) {
    int min = INT_MAX, min_index = 0;
    for (int i = 0; i < key.size(); ++i) {
        if (!nodeSet[i] && key[i] < min) {
            min = key[i];
            min_index = i;
        }
    }
    return min_index;
}

void printMST(const v & parent, const md & distance) {
    printf("Edge   Weight\n");
    for (int i = 0; i < parent.size(); i++)
        printf("%d - %d    %f \n", parent[i], i, distance[i][parent[i]]);
}

node getNode(int id, const nodes & graph) {
    node res;
    res.id = -1;
    for (node n : graph) if (n.id == id) return n;
    return res;
}


v primMST(const md & distance, const double range) {
    size_t nVertex = distance.size();
    v parent(nVertex, -1);
    vd key(nVertex, INT_MAX);
    vector<bool> nodeSet(nVertex, false);
    key[0] = 0;
    for (int i = 0; i < nVertex - 1; ++i) {
        int id = minKey(key, nodeSet);
        nodeSet[id] = true;
        for (int j = 0; j < nVertex; ++j) {
            if (!nodeSet[j] && distance[id][j] <= range && distance[id][j] < key[j]) {
                parent[j] = id;
                key[j] = distance[id][j];
            }
        }
    }
    return parent;
}


ifstream read_file(const char* name) {
    const char* filename  = name;
    ifstream file(filename);
    if (!file) {
        cerr << "ERROR: could not open file '" << filename << "' for reading" << endl;
        throw(-1);
    }
    return file;
}


int getTransmissionLevel(int threshold, const v & range) {
    for (int i = 0; i < range.size(); ++i) {
        if (threshold <= range[i]) return i;
    }
    return -1;
}


void readDistanceMatrix(ifstream & file, md & distance) {
    size_t nVertex = distance.size();
    for (int i = 0; i < nVertex; ++i) {
        for (int j = 0; j < nVertex; ++j) file >> distance[i][j];
    }
}

/* Reads transmission levels. Ratio is calculated as energy consumed per distance unit.
 That means that LESS ratio is better.
 */

void readTransmissionLevels(ifstream & file, int levels, vtrans & transmissions) {
    transmission trans;
    for (int i = 0; i < levels; ++i) {
        trans.level = i;
        file >> trans.energy >> trans.range;
        trans.ratio = trans.energy/trans.range;
        transmissions.push_back(trans);
    }
}

int getBestTransmissionLevel(const double dist, const vtrans & transmissions) {
    int best = -1;
    for (transmission t : transmissions) {
        if (t.range >= dist && (best == -1 || t.energy < transmissions[best].energy)) best = t.level;
    }
    return best;
}


nodes updateNodes(const nodes & solution, nodes & graph, const vtrans & transmissions, const md & distances) {
    nodes res;
    for (nodes::iterator it = graph.begin(); it != graph.end(); ++it) {
        int best = -1;
        int send_to = -1;
        for (node s : solution) {
            double dist = distances[(*it).id][s.id];
            int cur = getBestTransmissionLevel(dist, transmissions);
            if (cur != -1 && (best == -1 || transmissions[best].energy > transmissions[cur].energy)) {
                best = cur;
                send_to = s.id;
            }
        }
        if (best != -1) {
            node copy = *it;
            copy.transmission_level = best;
            copy.send_to = send_to;
            res.insert(copy);
        }
    }
    return res;
}

void printGraph(const nodes & graph, const double cost) {
    for (node n : graph) {
        cout << n.id << " - " << n.send_to << " Level: " << n.transmission_level << endl;
    }
    cout << "Cost " << cost << endl;
}

double minCost(const nodes & graph, const vtrans & transmissions) {
    double min = INTMAX_MAX;
    for (node n : graph) {
        if (n.transmission_level > -1 && transmissions[n.transmission_level].energy < min) {
            min = transmissions[n.transmission_level].energy;
        }
    }
    return min;
}

double maxCost(const nodes & graph, const vtrans & transmissions) {
    double max = INTMAX_MIN;
    for (node n : graph) {
        if (n.transmission_level > -1 && transmissions[n.transmission_level].energy > max) {
            max = transmissions[n.transmission_level].energy;
        }
    }
    return max;
}

nodes greedyRandomizedSolution(const double alpha, nodes graph, const vtrans & transmissions, const md & distances) {
    /**
     Every candidate has cost: energy[necessaryRange]+nodeEnergyConsumption.
     We choose the good elements (i.e. ratio is as low as possible).
     With this we ensure good nerwork balancing and less energy consumption.
     The ratio also introduces dinamicity in the choosing process.
     **/
    nodes solution, rcl;
    node bs = getNode(0, graph);
    if (bs.id == -1) return solution;
    /* Add base station to solution */
    
    solution.insert(bs);
    graph.erase(bs);
    while (graph.size() > 0) {
        rcl.clear();
        nodes updated = updateNodes(solution, graph, transmissions, distances);
        double min = minCost(updated, transmissions);
        double max = maxCost(updated, transmissions);
        for (node n : updated) {
            if (n.transmission_level > -1 && transmissions[n.transmission_level].energy <= min + alpha * (max - min)) {
                rcl.insert(n);
            }
        }
        if (rcl.size() > 0) {
            nodes::iterator it(rcl.begin());
            advance(it, rand() % rcl.size());
            solution.insert(*it);
            graph.erase(*it);
        }
        else {
            cout << "Unable to find a solution" << endl;
            exit(-1);
        }
    }
    return solution;
}

void fillNodes(nodes & graph, int nVertex) {
    for (int i = 0; i < nVertex; ++i) {
        node n;
        n.id = i;
        n.input_connections = 0;
        n.send_to = -1;
        n.transmission_level = -1;
        graph.insert(n);
    }
}

double calculateSolutionValue(const nodes & sol, const vtrans & transmissions) {
    double res = 0;
    for (node n : sol) {
        if (n.transmission_level > -1) {
            res += transmissions[n.transmission_level].energy;
        }
    }
    return res;
}


void improvePath(nodes & sol, node & random) {
    
}

nodes localSearch(nodes & sol, const vtrans & transmissions, const md & distance) {
    int adv = rand() % sol.size();
    nodes::iterator it(sol.begin());
    advance(it, adv);
    node random = *it;
    improvePath(sol, random);
    return sol;
}


int main(int argc, const char * argv[]) {
    srand(time(NULL));
    ifstream file = read_file("/Users/Sembroiz/Dropbox/DAC/workspace/GRASP/GRASP/input.txt");
    int iterations = 10;
    double alpha = 0.5;
    int nVertex, transmissionLevels;
    file >> nVertex;
    nodes graph;
    fillNodes(graph, nVertex);
    vtrans transmissions;
    md distance(nVertex, vd (nVertex));
    readDistanceMatrix(file, distance);
    file >> transmissionLevels;
    readTransmissionLevels(file, transmissionLevels, transmissions);
    
    /**
     
     GRASP START
     
     **/
    
    /**
     Construct initial solution based on the minimum spanning tree using the maximum
     possible transmission range (i.e. 82,92m).
     **/
    v init = primMST(distance, transmissions[transmissionLevels - 1].range);
    double best = INTMAX_MAX;
    nodes best_sol, sol;
    while(iterations-- > 0) {
        sol = greedyRandomizedSolution(alpha, graph, transmissions, distance);
        //sol = localSearch(sol, transmissions, distance);
        double cost = calculateSolutionValue(sol, transmissions);
        if (cost < best) {
            best = cost;
            best_sol = sol;
        }
        printGraph(best_sol, best);
    }
    
    return 0;
}
