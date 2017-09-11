//
//  generator.cpp
//  GRASP
//
//  Created by David Sembroiz on 04/05/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#include <iostream>
#include <fstream>
#include <random>
#include <set>
#include "Headers/definitions.hpp"
#include "Headers/generator.hpp"
using namespace std;


void generateCoordinates(int fileId, int numNodes, int rad) {
    set<pair<int, int> > coords;
    int i = 1;
    while (coords.size() < numNodes) {
        int x = rand() % (rad * 2) - rad;
        int y = rand() % (rad * 2) - rad;
        coords.insert(pair<int, int>(x, y));
        ++i;
    }
    string name = GRASP_PATH + "Inputs/coordinates_" + to_string(fileId) + ".txt";
    ofstream file(name);
    file << numNodes << endl;
    file << "0 0" << endl;
    for (set<pair<int, int> >::iterator it = coords.begin(); it != coords.end(); ++it) {
        file << (*it).first << " " << (*it).second << endl;
    }
}

