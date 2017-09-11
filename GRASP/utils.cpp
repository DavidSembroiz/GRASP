//
//  utils.cpp
//  GRASP
//
//  Created by David Sembroiz on 03/05/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#include <iostream>
#include <cmath>
#include <string>
#include <sys/stat.h>
#include "Headers/utils.hpp"
using namespace std;

void read_file(const string name, ifstream & f) {
    f.open(name);
    if (!f) {
        cerr << "ERROR: could not open file '" << name << "' for reading" << endl;
        throw(-1);
    }
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


mi fillDistanceMatrix(md & distance, int res) {
    string name = GRASP_PATH + "Inputs/coordinates_" + to_string(res) + ".txt";
    ifstream file;
    read_file(name, file);
    int numNodes;
    file >> numNodes;
    mi coords(numNodes, vi(2));
    for (int i = 0; i < numNodes; ++i) file >> coords[i][0] >> coords[i][1];
    for (int i = 0; i < numNodes; ++i) {
        vd dists(numNodes);
        for (int j = 0; j < numNodes; ++j) {
            dists[j] = sqrt((coords[i][0] - coords[j][0])*(coords[i][0] - coords[j][0]) +
                            (coords[i][1] - coords[j][1])*(coords[i][1] - coords[j][1]));
        }
        distance.push_back(dists);
    }
    return coords;
}

void saveSolutions(map<set<int>, pair<int, int> > & solutions, int resId, int nodes, int critical, long time) {
    string name = GRASP_PATH + "Results/latency_" + to_string(resId) + "_" + to_string(nodes) + "_" + to_string(critical) + ".txt";
    ofstream file(name);
    file << "Time " << time << "s" << endl << endl;
    for (map<set<int>, pair<int, int> >::iterator it = solutions.begin(); it != solutions.end(); ++it) {
        set<int> cns = (*it).first;
        file << "Critical Nodes ";
        for (set<int>::iterator it = cns.begin(); it != cns.end(); ++it) {
            file << (*it) << " ";
        }
        file << endl;
        file << "Occurrences " << (*it).second.first << endl;
        file << "Cost " << (*it).second.second << endl << endl;
    }
}

bool fileExists(const string & filename) {
    struct stat buf;
    if (stat(filename.c_str(), &buf) != -1) return true;
    return false;
}

int getConnectivity(int node, const md & distance) {
    int res = 0;
    for (int i = 0; i < distance[node].size(); ++i) {
        if (distance[node][i] <= 82.92) ++res;
    }
    return res;
}

void saveTest(int resultId, const mi & coords, const md & dist) {
    string name = GRASP_PATH + "Results/test_" + to_string(resultId) + "_" + to_string(coords.size()) + ".txt";
    if (fileExists(name)) return;
    ofstream file(name);
    for (int i = 0; i < coords.size(); ++i) {
        file << i + 1 << " " << coords[i][0] << " " << coords[i][1] << " ";
        file << i + 1 << " " << getConnectivity(i, dist) << " ";
        file << i + 1 << " " << dist[i][0] << endl;
    }
}
