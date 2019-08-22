#include<stdio.h>
#include<stdlib.h>
#include<string.h>
#include<math.h>
#include<time.h>

#define DIMENSIONS 2  
#define MAX_LINE_LENGTH 128  

FILE* checkArgs( int, char**, int, char const* const );
int** getData( FILE*, int* );
int* kmeans( int, int**, int );
void writeOutput( int**, int*, FILE*, int );

float getDistance( float*, int* );
float** getInitialCentroids( int, int**, int );
int* getNearestCentroids( int, float**, int**, int );
float** getNewCentroids( int, int*, int**, int );
int compareNearestCentroids( int*, int*, int );

int main( int argc, char** argv ) {
	int k, **data, numTuples[1], *clusters;
	FILE* inputFile, *outputFile = NULL;
	char const* const fileName = argv[2];
	clock_t start, end;
	double runtime;
	k = atoi(argv[1]);

	inputFile = checkArgs( argc, argv, k, fileName );
	data = getData( inputFile, numTuples );
	if( *numTuples <= k ) {
		fprintf(stderr, "K must be less than the number of data tuples.\n");
		return 1;
	}
	clusters =(int*) calloc( *numTuples, sizeof(int) );
	printf("Starting k-means with k=%d on %d tuples.\n", k, *numTuples);
	start = clock();
	clusters = kmeans( k, data, *numTuples );
	end = clock();
	runtime = (double)(end - start) / CLOCKS_PER_SEC;
	printf("Finished k-means in %.3lf seconds.\n", runtime);
	writeOutput( data, clusters, outputFile, *numTuples );
	return 0;
}

FILE* checkArgs( int argc, char** argv, int k, char const* const fileName ) {
	FILE* file;
	if(argc != 3) {
		fprintf(stderr, "Usage: %s [k] [file name]\n", argv[0]);
		exit(1);
	}

	if(k < 0) {
		fprintf(stderr, "k must be an integer greater than 0.\n");
		exit(1);
	}

	if( (file = fopen(fileName, "r")) == NULL ) {
		fprintf(stderr, "Error opening input file.\n");
		exit(1);
           }
	return file;
}

int** getData( FILE* inputFile, int* numTuples ) {
	int* data_values;
	int **data; 
	char line[MAX_LINE_LENGTH];
	char* token = NULL;
	int i, j;
           int rows = 0;
	
    while( fgets( line, MAX_LINE_LENGTH, inputFile) != NULL ) {
		++rows;
    }
	rewind(inputFile);
	data_values = calloc(DIMENSIONS*rows, sizeof(int));
	data = malloc(rows*sizeof(int*));
	for(i=0; i<rows; ++i) {
        data[i] = data_values + i*DIMENSIONS;
	}
	for(i=0; i<rows; ++i) {
		fgets( line, MAX_LINE_LENGTH, inputFile);
		token = strtok(line, " \n");
		j=0;
		while(token) {
			data[i][j] = atoi(token);
			token = strtok(NULL, " \n");
			++j;
		}
	}

	*numTuples = rows;
	return data;
}

float getDistance( float* a, int* b ) {
	float distance = 0.0;
	distance = pow( pow(b[0]-a[0], 2.0) + pow(b[1]-a[1], 2.0), 0.5 );
	return distance;
}

float** getInitialCentroids( int k, int** data, int numTuples ) {
	int i, j, lastClusterSize, clusterSize = numTuples/k; 
	float *new_centroid_values, **new_centroid; 
	float sum[2] = { 0.0, 0.0 };  
	int numComputed = 0;   

	new_centroid_values = calloc(DIMENSIONS*k, sizeof(float));
	new_centroid = malloc(k*sizeof(float*));
	for(i=0; i<k; ++i) {
        new_centroid[i] = new_centroid_values + i*DIMENSIONS;
	}

	for(i=0; i<numTuples; ++i) {
		for(j=0; j<DIMENSIONS; ++j) {
			sum[j] += data[i][j];
		}

		if( numComputed == k-1 ) {
			for(j=i; j<numTuples; j++) {
				for(i=0; i<DIMENSIONS; ++i) {
					sum[i] += data[j][i];
				}
			}
			lastClusterSize = numTuples - (clusterSize * (k-1));
			sum[0] = sum[0]/lastClusterSize;
			sum[1] = sum[1]/lastClusterSize;
			new_centroid[numComputed][0] = sum[0];
			new_centroid[numComputed][1] = sum[1];
			++numComputed;
			sum[0] = sum[1] = 0.0;
			break;
		}
		else if( (i+1) % clusterSize == 0 ) {
			sum[0] = sum[0]/clusterSize;
			sum[1] = sum[1]/clusterSize;
			new_centroid[numComputed][0] = sum[0];
			new_centroid[numComputed][1] = sum[1];
			++numComputed;
			sum[0] = sum[1] = 0.0;
		}
	}

	return new_centroid;
}

int* getNearestCentroids( int k, float** centroid, int** data, int numTuples ) {
	int i, j;
	int indexClosestMean, *newNearestMeans;
	newNearestMeans = malloc( numTuples*sizeof(int) );
	for(i=0; i<numTuples; ++i) {
		indexClosestMean = 0;
		for(j=1; j<k; ++j) {
			if( getDistance(centroid[j], data[i]) < getDistance(centroid[indexClosestMean], data[i]) )
			    indexClosestMean = j;
		}
		newNearestMeans[i] = indexClosestMean;
	}

	return newNearestMeans;
}

float** getNewCentroids( int k, int* nearestMeans, int** data, int numTuples ) {
	int i, j;
	float *new_centroid_values, **new_centroid;
	float sum[2] = { 0.0, 0.0 };  
	int clusterSize = 0;   /* size of each cluster */

	new_centroid_values = calloc(DIMENSIONS*k, sizeof(float));
	new_centroid = malloc(k*sizeof(float*));
	for(i=0; i<k; ++i) {
        new_centroid[i] = new_centroid_values + i*DIMENSIONS;
	}

	for(i=0; i<k; ++i) {
		for(j=0; j<numTuples; ++j) {
			if(nearestMeans[j] == i) {
				sum[0] += data[j][0];
				sum[1] += data[j][1];
				++clusterSize;
			}
		}
		sum[0] = sum[0]/clusterSize;
		sum[1] = sum[1]/clusterSize;
		new_centroid[i][0] = sum[0];
		new_centroid[i][1] = sum[1];
		sum[0] = sum[1] = 0.0;
		clusterSize = 0;
	}

	return new_centroid;
}

int compareNearestCentroids( int* last, int* current, int numTuples ) {
	int i, count = 0;
	if( last == NULL || current == NULL ) {
		return 0;
	}
	for(i=0; i<numTuples; ++i) {
		if( last[i] != current[i] )
			return 0;
	}

    return 1;
}

int* kmeans( int k, int** data, int numTuples ) {
	int *nearestMeans, *lastMeans; 
	int convergence, i; 
	float *centroid_values, **centroid;

	nearestMeans = calloc(numTuples, sizeof(int));
	lastMeans = calloc(numTuples, sizeof(int));
	centroid_values = calloc(DIMENSIONS*k, sizeof(float));
	centroid = malloc(k*sizeof(float*));
	for(i=0; i<k; ++i) 
	{
        centroid[i] = centroid_values + i*DIMENSIONS;
	}

	convergence = 0;
    centroid = getInitialCentroids( k, data, numTuples );
	do {
	    nearestMeans = getNearestCentroids( k, centroid, data, numTuples );
	    centroid = getNewCentroids( k, nearestMeans, data, numTuples );
		convergence = compareNearestCentroids(lastMeans, nearestMeans, numTuples);
		memcpy( lastMeans, nearestMeans, numTuples*sizeof(int) );
	} while( convergence != 1 );

	return nearestMeans;
}

void writeOutput( int** data, int* clusters, FILE* outputFile, int numTuples ) {
	int i;
	char* line = NULL;
    if( (outputFile = fopen("output.txt", "w")) != NULL ) {	
	    for(i=0; i<numTuples; ++i) {
	        fprintf(outputFile, "%d %d %d \n", data[i][0], data[i][1], clusters[i]+1 );
        }
		printf("Wrote results to output.txt\n");
	}
	else {
        fprintf(stderr, "Error opening output file; printing to stdout.\n");
		for(i=0; i<numTuples; ++i) {
		    printf("%d %d %d \n", data[i][0], data[i][1], clusters[i]+1 );
	    }
	}
}


