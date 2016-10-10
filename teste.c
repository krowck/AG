void populateRoutes(){

	for(int count=0; count<maxClusters; count++){
		clusterCenters[count] = -1;
		innerClusterCircleFlag[count] = 0;
		outerClusterCircleFlag[count] = 0;
	}

	for(int count=0; count<NUMBEROFPOINTS; count++){
		allRoutesCoOrdinatesZ[count] = 0.0f;
	}

	int pointFinishedIndex = 0;

	/// 10% random points
	srand(time(NULL));
	for(int count=0; count<NUMBEROFPOINTS/10; count++){
		allRoutesCoOrdinatesX[count] = (float)rand()/(float)RAND_MAX;
		allRoutesCoOrdinatesY[count] = (float)rand()/(float)RAND_MAX;
		pointFinishedIndex++;
	}

	// generating cluster centers
	srand(time(NULL));
	for(int count=0; count<NUMBEROFCLUSTERS; count++){
		clusterWidth[count] = ((float)rand()/(float)RAND_MAX)*0.6f+0.2f;
		clusterCenterX[count] = ((float)rand()/(float)RAND_MAX)*(1.0f-clusterWidth[count])+(clusterWidth[count]*0.5f);
		clusterCenterY[count] = ((float)rand()/(float)RAND_MAX)*(1.0f-clusterWidth[count])+(clusterWidth[count]*0.5f);
	}

	// maximum points in multiple of NUMBEROFCLUSTERS
	int remainingPoints = NUMBEROFPOINTS - pointFinishedIndex;
	int numOfPointsPerCluster = remainingPoints/NUMBEROFCLUSTERS;

	// fill clusters
	int temp = pointFinishedIndex;
	for(int count=0; count<NUMBEROFCLUSTERS; count++){
		float width = clusterWidth[count];
		srand(time(NULL));
		for(int count1=0; count1<numOfPointsPerCluster; count1++){
			int index = temp+(count*numOfPointsPerCluster)+count1;
			if(count1<numOfPointsPerCluster/5){
				allRoutesCoOrdinatesX[index]=
					((float)rand()/(float)RAND_MAX)*(width)+
					(clusterCenterX[count]-(width*0.5f));
				allRoutesCoOrdinatesY[index]=
					((float)rand()/(float)RAND_MAX)*(width)+
					(clusterCenterY[count]-(width*0.5f));
			}else if(count1<numOfPointsPerCluster/4){
				float temp_width = width*(float)rand()/(float)RAND_MAX;
				allRoutesCoOrdinatesX[index]=
					((float)rand()/(float)RAND_MAX)*(temp_width)+
					(clusterCenterX[count]-(temp_width*0.5f));
				allRoutesCoOrdinatesY[index]=
					((float)rand()/(float)RAND_MAX)*(temp_width)+
					(clusterCenterY[count]-(temp_width*0.5f));
			}else if(count1<numOfPointsPerCluster/3){
				float temp_width = width*(float)rand()/(float)RAND_MAX;
				allRoutesCoOrdinatesX[index]=
					((float)rand()/(float)RAND_MAX)*(temp_width)+
					(clusterCenterX[count]-(temp_width*0.5f));
				allRoutesCoOrdinatesY[index]=
					((float)rand()/(float)RAND_MAX)*(temp_width)+
					(clusterCenterY[count]-(temp_width*0.5f));
			}else if(count1<numOfPointsPerCluster/2){
				float temp_width = width*(float)rand()/(float)RAND_MAX;
				allRoutesCoOrdinatesX[index]=
					((float)rand()/(float)RAND_MAX)*(temp_width)+
					(clusterCenterX[count]-(temp_width*0.5f));
				allRoutesCoOrdinatesY[index]=
					((float)rand()/(float)RAND_MAX)*(temp_width)+
					(clusterCenterY[count]-(temp_width*0.5f));
			}else{
				float temp_width = width*(float)rand()/(float)RAND_MAX;
				allRoutesCoOrdinatesX[index]=
					((float)rand()/(float)RAND_MAX)*(temp_width)+
					(clusterCenterX[count]-(temp_width*0.5f));
				allRoutesCoOrdinatesY[index]=
					((float)rand()/(float)RAND_MAX)*(temp_width)+
					(clusterCenterY[count]-(temp_width*0.5f));
			}
			if(allRoutesCoOrdinatesX[index]<=0.0f){
				allRoutesCoOrdinatesX[index]=0.00001f;
			}else if(allRoutesCoOrdinatesX[index]>=1.0f){
				allRoutesCoOrdinatesX[index]=0.99999f;
			}

			if(allRoutesCoOrdinatesY[index]<=0.0f){
				allRoutesCoOrdinatesY[index]=0.00001f;
			}else if(allRoutesCoOrdinatesY[index]>=1.0f){
				allRoutesCoOrdinatesY[index]=0.99999f;
			}
			pointFinishedIndex++;
		}
	}

	// fill remaining points
	srand(time(NULL));
	for(int count = pointFinishedIndex; count < NUMBEROFPOINTS; count++){
		allRoutesCoOrdinatesX[count] = (float)rand()/(float)RAND_MAX;
		allRoutesCoOrdinatesY[count] = (float)rand()/(float)RAND_MAX;
		pointFinishedIndex++;
	}

	// shuffle the points in array
    if (NUMBEROFPOINTS > 1) {
        int i;
		for (i = 0; i < NUMBEROFPOINTS - 1; i++) {
		  int j = i + rand() / (RAND_MAX / (NUMBEROFPOINTS - i) + 1);
		  float tx = allRoutesCoOrdinatesX[j];
		  float ty = allRoutesCoOrdinatesY[j];
		  allRoutesCoOrdinatesX[j] = allRoutesCoOrdinatesX[i];
		  allRoutesCoOrdinatesY[j] = allRoutesCoOrdinatesY[i];
		  allRoutesCoOrdinatesX[i] = tx;
		  allRoutesCoOrdinatesY[i] = ty;
		}
    }
}