#ifndef PERSPECTIVE3POINT_H
#define PERSPECTIVE3POINT_H
#define MAX_NUM_POINTS 1200  // maximum number of points
#define ARRAY_SIZE 20  // 5*4
#define MAXIMUM 1200
#define MAX_CORRESPONDENCES 1000

#include<stdio.h>
#include <string.h> 
#include <stdbool.h>

typedef struct {
    double data[3][3];
} Matrix3x3;

typedef struct {
    double data[4][4];
} Matrix4x4;

typedef struct {
    double x, y, z;
} ThreeDPoint;

typedef struct {
    double x, y;
} TwoDPoint;

typedef struct {
    double focalX, focalY, principalX, principalY;
    double invFocalX, invFocalY, principalXFocalX, principalYFocalY;
} Perspective3Point;

// typedef enum {
//     CV_8U = 0,   // 8-bit unsigned integer
//     CV_8S = 1,   // 8-bit signed integer
//     CV_16U = 2,  // 16-bit unsigned integer
//     CV_16S = 3,  // 16-bit signed integer
//     CV_32S = 4,  // 32-bit signed integer
//     CV_32F = 5,  // 32-bit floating-point
//     CV_64F = 6   // 64-bit floating-point
// } MatrixDepth;

typedef struct {
    float data[3][3];
} CameraMatrixFloat;

typedef struct {
    double data[3];
} Vector3;

// RANSAC Structs 
// Updated Structs
typedef struct {
    double x, y;
} RANSAC_2D_Point;

typedef struct {
    double x, y, z;
} RANSAC_3D_Point;

typedef struct {
    RANSAC_2D_Point p2DPoints[MAXIMUM];
    TwoDPoint imagePoints[MAX_CORRESPONDENCES];
    float sigma2[MAXIMUM];
    RANSAC_3D_Point p3DPoints[MAXIMUM];
    ThreeDPoint objectPoints[MAX_CORRESPONDENCES];
    int allIndices[MAXIMUM];
    int index;
    int maxCorrespondences;
    double probability;
    int minInliers;
    int maxIterations;
    int currentIterations;
    float epsilon;
    int minSet;
    float maxError[MAXIMUM];
    int inliersCount;
    int correspondencesCount;
    int allIndices2D[MAX_NUM_POINTS];
    int bestInliers2D[MAX_NUM_POINTS];
    int inliers2D[MAX_NUM_POINTS];
    int keyPointIndices[MAX_NUM_POINTS];
    int mapPointMatches[MAX_NUM_POINTS];
    int numCorrespondences;
    int minSetRANSAC;
    int minInliersRANSAC;
    int maxIterationsRANSAC;
    int maxIterationsCurrent;
    int bestInliersCount;
    int currentInliersCount;
    double bestTransformation[4][4];
} RANSAC_Solver;

// Updated Structs
typedef struct {
    RANSAC_2D_Point keyPointsUn[MAXIMUM];
    float levelSigma2[MAXIMUM];
} FrameData; // From Frame.h Assumed

typedef struct {
    bool isBad;
    RANSAC_3D_Point worldPosition;
} MapPointData;

float worldPointsFloat[MAX_NUM_POINTS][3];  
double worldPointsDouble[MAX_NUM_POINTS][3]; 

int countValidImagePoints32Bit();
int countValidImagePoints64Bit();

double focalX, focalY, principalX, principalY;

void extract3DPoints(ThreeDPoint* worldPoints, TwoDPoint* imagePoints, double* pointsArray);

// Initialization functions
void initializeP3P(Perspective3Point* solver);
void initializeP3PWithParameters(Perspective3Point* solver, double focalX, double focalY, double principalX, double principalY);
void setCameraMatrix(Perspective3Point* solver, const Matrix3x3* intrinsicMatrix);
void setCameraMatrixFloat(Perspective3Point* solver, const Matrix3x3* intrinsicMatrixFloat);
void initializeAP3PWithIntrinsicMatrix(Perspective3Point* solver, const Matrix3x3* intrinsicMatrix, int matrixDepth);
void computeInverseIntrinsics(Perspective3Point* solver);

int solveAP3P(Perspective3Point* solver, Matrix3x3* rotation, Vector3* translation,
                double imageX0, double imageY0, double worldX0, double worldY0, double worldZ0,
                double imageX1, double imageY1, double worldX1, double worldY1, double worldZ1,
                double imageX2, double imageY2, double worldX2, double worldY2, double worldZ2,
                double imageX3, double imageY3, double worldX3, double worldY3, double worldZ3,
                bool isP4P);

bool solveSingleAP3P(Perspective3Point* solver, Matrix3x3* rotation, Vector3* translation, 
                double imageX0, double imageY0, double worldX0, double worldY0, double worldZ0,
                double imageX1, double imageY1, double worldX1, double worldY1, double worldZ1,
                double imageX2, double imageY2, double worldX2, double worldY2, double worldZ2,
                double imageX3, double imageY3, double worldX3, double worldY3, double worldZ3);

bool solveForAP3P(Perspective3Point* solver, Matrix3x3* rotation, Vector3* translation, const Matrix3x3* worldPointsMatrix, const Matrix3x3* imagePointsMatrix);
int solveMultipleAP3P(Perspective3Point* solver, Matrix3x3* rotations, Vector3* translations, double* worldPoints, double* imagePoints, double focalX, double principalX, double focalY, double principalY);

int computeMultiplePoses(const double featureVectors[3][4], const double worldPoints[3][4], double rotationSolutions[4][3][3],
                     double translationSolutions[4][3], bool isP4P);


// RANSAC //
void RANSAC_Init(RANSAC_Solver* solver, const FrameData* frame, MapPointData* mapPointMatches[], size_t matchSizes);
void SetRansacParameters(RANSAC_Solver* solver, double probability, int minInliers, int maxIterations, int minSet, float epsilon, float th2);
void IterateRANSAC(RANSAC_Solver* solver, Perspective3Point* ap3pObj, int nIterations, int* nInliers, int vbInliers[MAXIMUM]);
void CheckInliersRANSAC(RANSAC_Solver* solver, Perspective3Point* ap3pObj, Matrix3x3* R, Vector3* t);

#endif
