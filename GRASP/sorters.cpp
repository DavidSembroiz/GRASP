//
//  connSorter.cpp
//  GRASP
//
//  Created by David Sembroiz on 03/05/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#include <iostream>
#include <numeric>
#include <vector>
#include "Headers/sorters.hpp"
#include "Headers/utils.hpp"
#include "Headers/definitions.hpp"
using namespace std;

ConnSorter::ConnSorter(md distance) {
    for (int i = 0; i < distance.size(); ++i) _conn.push_back(getConnectivity(i, distance));
}

    
bool ConnSorter::operator()(int n1, int n2) const {
    return _conn[n1] > _conn[n2];
}


DistanceSorter::DistanceSorter(md distance) {
    _distance = distance;
}
    
bool DistanceSorter::operator()(int n1, int n2) const {
    return _distance[n1][0] < _distance[n2][0];
}


    
vi RelaySorter::fillRelay(int id) {
    string name = GRASP_PATH + "Inputs/solution_" + to_string(id) + ".txt";
    ifstream file;
    read_file(name, file);
    if (!file) {
        
    }
    int nVertex, node, send, pkts;
    file >> nVertex;
    vi res(nVertex);
    for (int i = 0; i < nVertex; ++i) {
        file >> node >> send >> pkts;
        res[node - 1] = pkts;
    }
    return res;
}
    
RelaySorter::RelaySorter(int id) {
    _relay = fillRelay(id);
}
    
bool RelaySorter::operator()(int n1, int n2) const {
    return _relay[n1] > _relay[n2];
}

/**
 
 CRITICAL NODES SUBSET
 
 According to the following features, only nodes fulfilling the conditions are tested for being critical. The list of features is:
 - Distance to Base Station in meters.
 - Connectivity (i.e. number of nodes in range).
 - Relay (i.e. number of packets to transmit: own + parent nodes packets).
 
 **/

vi getCriticalNodeSubset(int nVertex, const md & distance, const int RESULT_ID) {
    vi distanceRank(nVertex - 1), connRank(nVertex - 1), relayRank(nVertex - 1);
    iota(distanceRank.begin(), distanceRank.end(), 1);
    iota(connRank.begin(), connRank.end(), 1);
    iota(relayRank.begin(), relayRank.end(), 1);
    
    sort(distanceRank.begin(), distanceRank.end(), DistanceSorter(distance));
    sort(connRank.begin(), connRank.end(), ConnSorter(distance));
    sort(relayRank.begin(), relayRank.end(), RelaySorter(RESULT_ID));
    
    vi criticalNodeSubset;
    for (int i = 0; i < distanceRank.size() && distance[distanceRank[i]][0] <= 82.92*DISTANCE_HC_THRESH; ++i) {
        int n = distanceRank[i];
        ptrdiff_t connPos = find(connRank.begin(), connRank.end(), n) - connRank.begin();
        ptrdiff_t relayPos = find(relayRank.begin(), relayRank.end(), n) - relayRank.begin();
        if (n != 0 && (connPos < CONN_RANK_THRESH || relayPos < RELAY_RANK_THRESH)) {
            criticalNodeSubset.push_back(n);
        }
    }
    for (int i = 0; i < CONN_RANK_THRESH; ++i) {
        int n = connRank[i];
        if (n != 0 && (find(relayRank.begin(), relayRank.end(), n) - relayRank.begin() < RELAY_RANK_THRESH)) {
            if (find(criticalNodeSubset.begin(), criticalNodeSubset.end(), n) == criticalNodeSubset.end())
            criticalNodeSubset.push_back(n);
        }
    }
    return criticalNodeSubset;
}


