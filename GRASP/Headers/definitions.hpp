//
//  definitions.hpp
//  GRASP
//
//  Created by David Sembroiz on 03/05/17.
//  Copyright © 2017 David Sembroiz. All rights reserved.
//

#ifndef definitions_hpp
#define definitions_hpp

#define vi vector<int>
#define mi vector<vi>
#define vd vector<double>
#define md vector<vd>
#define vb vector<bool>
#define mb vector<vb>
using namespace std;


const string GRASP_PATH = "/Users/Sembroiz/Dropbox/DAC/workspace/GRASP/GRASP/";

/**
 
 Parameters to control the creation of a new network.
 - GEN_NEW_NETWORK: specifies whether a new network needs to be created or not.
 - NUM_NODES: number of nodes to create.
 - NET_RADIUS: radius in meters of the node placement region.
 
 **/

const bool GEN_NEW_NETWORK = false;


/* Values to restrict the rank to which critical nodes are searched using each feature */

const double DISTANCE_HC_THRESH = 4.5;
const int CONN_RANK_THRESH = 10;
const int RELAY_RANK_THRESH = 10;

/**
 
 Parameter to decide the Objective function to apply:
 - 0: latency (hop count).
 - 1: lifetime (energy consumed).
 
 **/

const int OBJECTIVE = 0;


/* Iterations that the whole GRASP (all instances) is run to extract mean value */

const int GLOBAL_GRASP_ITERATIONS = 1;

/* Iterations for a single GRASP instance */

const int GRASP_ITERATIONS = 8;

/* GRASP alpha parameter */

const double GRASP_ALPHA = 0;

/* Neighbors to explore inside the local search */

const int LOCAL_SEARCH_MAX_NEIGHBORS = 50;

/* Iterations with no improvement until local search is stopped */

const int LOCAL_SEARCH_NO_IMPROV_THRESHOLD = 50;


#endif /* definitions_hpp */
