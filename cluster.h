#include <stdio.h>
#include <stdlib.h>
#include <time.h>
#include <math.h>

const static int ROUTELENGTH = 100;
const static int NUMBEROFROUTES = 50;
const static int NUMBEROFPOINTS = ROUTELENGTH*NUMBEROFROUTES;

const static int STEP_NUM = 100;
const static float PARTICLE_SIZE = 1.0f;
const static float RADIUS = 0.1f;
const static float OUTERRADIUS = 1.25f*RADIUS;
const static float ALPHA = 4.0f/pow(RADIUS,2);
const static float BETA = 4.0f/pow(OUTERRADIUS,2);

const static int NUMOFSEGMENTSTODRAWCIRCLE = 20;

static float allRoutesCoOrdinatesX[NUMBEROFPOINTS];
static float allRoutesCoOrdinatesY[NUMBEROFPOINTS];
static float allRoutesCoOrdinatesZ[NUMBEROFPOINTS];


const static int NUMBEROFCLUSTERS = 6;

static float clusterWidth[NUMBEROFCLUSTERS];
static float clusterCenterX[NUMBEROFCLUSTERS];
static float clusterCenterY[NUMBEROFCLUSTERS];


void populateRoutes();

