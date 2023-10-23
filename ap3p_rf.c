#include "ap3p_rf.h"

#include <math.h>
#include <stdbool.h>
#include <string.h>
#include <complex.h>
#include <stdint.h>

#if defined (_MSC_VER) && (_MSC_VER <= 1700)
static double cbrt(double x) { 
    
    return pow(x, 1.0/3.0); 
}
#endif


void RANSAC_Init(RANSAC_Solver* solver, const FrameData* frame, MapPointData* mapPointMatches[], size_t matchSizes) {
    solver->index = 0;
    void initializeAP3PWithIntrinsicMatrix(Perspective3Point* solver, const Matrix3x3* intrinsicMatrix, int matrixDepth);
    for (size_t i = 0; i < matchSizes; i++) {
        MapPointData* mapPoint = mapPointMatches[i];

        if (mapPoint && !mapPoint->isBad) {
            solver->p2DPoints[solver->index] = frame->keyPointsUn[i];
            solver->sigma2[solver->index] = 1.0f;
            solver->p3DPoints[solver->index] = mapPoint->worldPosition;
            solver->allIndices[solver->index] = solver->index;

            solver->index++;
        }
    }

    SetRansacParameters(solver, 0.99, 8, 300, 4, 0.4, 5.991);
}

void IterateRANSAC(RANSAC_Solver* solver, Perspective3Point* ap3pObj, int nIterations, int* nInliers, int vbInliers[MAXIMUM]) {
    int noMore = 0;
    *nInliers = 0;

    if (solver->correspondencesCount < solver->minInliersRANSAC) {
        noMore = 1;
        return;
    }

    int nCurrentIterations = 0;
    while (nCurrentIterations < nIterations && !noMore) {
        nCurrentIterations++;
        solver->maxIterationsCurrent++;
        

        for (int i = 0; i < solver->numCorrespondences; i++) {
            solver->allIndices2D[i] = i;
        }

        // 1. Randomly select 4 points
        double featureVectors[3][4];
        double worldPoints[3][4];
        for (short i = 0; i < solver->minSetRANSAC; ++i) {
            int randi = rand() % (solver->numCorrespondences - i);
            int idx = solver->allIndices2D[randi];

            worldPoints[0][i] = solver->p3DPoints[idx].x;
            worldPoints[1][i] = solver->p3DPoints[idx].y;
            worldPoints[2][i] = solver->p3DPoints[idx].z;

            featureVectors[0][i] = solver->p2DPoints[idx].x;
            featureVectors[1][i] = solver->p2DPoints[idx].y;
            
            

            solver->allIndices2D[randi] = solver->allIndices2D[solver->numCorrespondences - i - 1];
        }

        // 2. Compute camera pose
        double solutionsR[4][3][3];
        double solutionsT[4][3];
        int nbSolutions = computePoses(featureVectors, worldPoints, solutionsR, solutionsT, false);

        // 3. Check inliers for each solution
        for (int i = 0; i < nbSolutions; i++) {
            Matrix3x3 R;
            Vector3 t;
            for (int row = 0; row < 3; row++) {
                for (int col = 0; col < 3; col++) {
                    R.data[row][col] = solutionsR[i][row][col];
                }
                t.data[row] = solutionsT[i][row];
            }
            CheckInliersRANSAC(solver, ap3pObj, &R, &t);

            // 4. If the current solution has more inliers, then record this solution as the best solution
            if (solver->inliersCount > solver->bestInliersCount) {
                for (int j = 0; j < solver->numCorrespondences; j++) {
                    solver->bestInliers2D[j] = solver->inliers2D[j];
                }
                solver->bestInliersCount = solver->inliersCount;

                for (int row = 0; row < 3; row++) {
                    for (int col = 0; col < 3; col++) {
                        solver->bestTransformation[row][col] = R.data[row][col];
                    }
                    solver->bestTransformation[row][3] = t.data[row];
                }
                for (int col = 0; col < 3; col++) {
                    solver->bestTransformation[3][col] = 0.0;
                }
                solver->bestTransformation[3][3] = 1.0;
            }
        }
    }

    // Check the stopping condition
    if (nCurrentIterations >= nIterations) {
        noMore = 1;
        if (solver->bestInliersCount >= solver->minInliersRANSAC) {
            *nInliers = solver->bestInliersCount;
            for (int i = 0; i < solver->numCorrespondences; i++) {
                vbInliers[i] = solver->bestInliers2D[i];
            }
        }
    }
}

void SetRansacParameters(RANSAC_Solver* solver, double probability, int minInliers, int maxIterations, int minSet, float epsilon, float th2) {
    solver->probability = probability;
    solver->minInliers = minInliers;
    solver->maxIterations = maxIterations;
    solver->epsilon = epsilon;
    solver->minSet = minSet;

    // Adjust Parameters according to the number of correspondences
    int nMinInliers = solver->index * solver->epsilon;
    if (nMinInliers < solver->minInliers)
        nMinInliers = solver->minInliers;
    if (nMinInliers < minSet)
        nMinInliers = minSet;
    solver->minInliers = nMinInliers;

    if (solver->epsilon < (float)solver->minInliersRANSAC / solver->index)
        solver->epsilon = (float)solver->minInliersRANSAC / solver->index;

    // RANSAC Iterations
    int nIterations;

    if (solver->minInliersRANSAC == solver->index)
        nIterations = 1;
    else
        nIterations = (int)ceil(log(1 - solver->probability) / log(1 - pow(solver->epsilon, 3)));
    
    solver->maxIterations = nIterations < maxIterations ? nIterations : maxIterations;

    for (int i = 0; i < solver->index; i++) {
        solver->maxError[i] = solver->sigma2[i] * th2;
    }
}

void CheckInliersRANSAC(RANSAC_Solver* solver, Perspective3Point* ap3pObj, Matrix3x3* R, Vector3* t) {
    solver->inliersCount = 0;

    for (int i = 0; i < MAXIMUM; i++) {
        RANSAC_3D_Point p3Dw = solver->p3DPoints[i];
        RANSAC_2D_Point p2D = solver->p2DPoints[i];

        double Xc = R->data[0][0] * p3Dw.x + R->data[0][1] * p3Dw.y + R->data[0][2] * p3Dw.z + t->data[0];
        double Yc = R->data[1][0] * p3Dw.x + R->data[1][1] * p3Dw.y + R->data[1][2] * p3Dw.z + t->data[1];
        double invZc = 1.0 / (R->data[2][0] * p3Dw.x + R->data[2][1] * p3Dw.y + R->data[2][2] * p3Dw.z + t->data[2]);

        double ue = ap3pObj->principalX + ap3pObj->focalX * Xc * invZc;
        double ve = ap3pObj->principalY + ap3pObj->focalY * Yc * invZc;

        double distX = p2D.x - ue;
        double distY = p2D.y - ve;

        double error2 = distX * distX + distY * distY;

        if (error2 < solver->maxError[i]) {
            // This point is an inlier
            solver->inliersCount++;
        }
    } 
}

void solveQuartic(const double *coefficients, double *roots) {
    const double coeffA4 = coefficients[0];
    const double coeffA3 = coefficients[1];
    const double coeffA2 = coefficients[2];
    const double coeffA1 = coefficients[3];
    const double coeffA0 = coefficients[4];

    double coeffA4Squared = coeffA4 * coeffA4;
    double coeffA3Squared = coeffA3 * coeffA3;
    double coeffA4Cubed = coeffA4Squared * coeffA4;
    double coeffA2TimesA4 = coeffA2 * coeffA4;

    double quarticP = (8 * coeffA2TimesA4 - 3 * coeffA3Squared) / (8 * coeffA4Squared);
    double quarticQ = (coeffA3Squared * coeffA3 - 4 * coeffA2TimesA4 * coeffA3 + 8 * coeffA1 * coeffA4Squared) / (8 * coeffA4Cubed);
    double quarticR = (256 * coeffA0 * coeffA4Cubed - 3 * (coeffA3Squared * coeffA3Squared) - 64 * coeffA1 * coeffA3 * coeffA4Squared + 16 * coeffA2TimesA4 * coeffA3Squared) / (256 * (coeffA4Cubed * coeffA4));

    double cubicP = ((quarticP * quarticP) / 12 + quarticR) / -3;
    double cubicQ = (72 * quarticR * quarticP - 2 * quarticP * quarticP * quarticP - 27 * quarticQ * quarticQ) / 216;

    double tempValue;
    double complex cubicRoot;
    if (cubicQ >= 0)
        cubicRoot = -csqrt(cpow(cubicQ * cubicQ - cubicP * cubicP * cubicP, 0.5) - cubicQ);
    else
        cubicRoot = csqrt(cpow(cubicQ * cubicQ - cubicP * cubicP * cubicP, 0.5) - cubicQ);
    if (cimag(cubicRoot) == 0.0) {
        cubicRoot = cpow(cubicRoot, 1.0 / 3);
        tempValue = 2.0 * (creal(cubicRoot) + cubicP / creal(cubicRoot));
    } else {
        cubicRoot = cpow(cubicRoot, 1.0 / 3);
        tempValue = 4.0 * creal(cubicRoot);
    }

    double complex sqrtTerm = csqrt(-2 * quarticP / 3 + tempValue);
    double quarticShift = -coeffA3 / (4 * coeffA4);
    double complex complexTerm1 = 4 * quarticP / 3 + tempValue;

    double complex complexTerm2 = 2 * quarticQ / sqrtTerm;

    double halfSqrtTerm = creal(sqrtTerm) / 2;
    double root1 = sqrt(-(creal(complexTerm1) + creal(complexTerm2)) / 2);
    roots[0] = quarticShift + halfSqrtTerm + root1;
    roots[1] = quarticShift + halfSqrtTerm - root1;
    double root2 = sqrt(-(creal(complexTerm1) - creal(complexTerm2)) / 2);
    roots[2] = quarticShift - halfSqrtTerm + root2;
    roots[3] = quarticShift - halfSqrtTerm - root2;
}

void refineQuarticRoots(const double *coefficients, double *roots) {
    const int maxIterations = 2;
    int iter, rootIndex;
    for (iter = 0; iter < maxIterations; ++iter) {
        for (rootIndex = 0; rootIndex < 4; ++rootIndex) {
            double rootError = 
                (((coefficients[0] * roots[rootIndex] + coefficients[1]) * roots[rootIndex] + coefficients[2]) * roots[rootIndex] + coefficients[3]) * roots[rootIndex] + 
                coefficients[4];
            double rootDerivative = 
                ((4 * coefficients[0] * roots[rootIndex] + 3 * coefficients[1]) * roots[rootIndex] + 2 * coefficients[2]) * roots[rootIndex] + coefficients[3];
            roots[rootIndex] -= rootError / rootDerivative;
        }
    }
}


void initializeSolver(Perspective3Point* solver) {
    solver->focalX = 0;
    solver->focalY = 0;
    solver->principalX = 0;
    solver->principalY = 0;
    solver->invFocalX = 0;
    solver->invFocalY = 0;
    solver->principalXFocalX = 0;
    solver->principalYFocalY = 0;
}

void initializeSolverWithValues(Perspective3Point* solver, double focalX, double focalY, double principalX, double principalY) {
    solver->focalX = focalX;
    solver->focalY = focalY;
    solver->principalX = principalX;
    solver->principalY = principalY;
    computeInverseParameters(solver);
}

void setIntrinsicParameters(Perspective3Point* solver, const Matrix3x3* intrinsicMatrix) {
    solver->principalX = intrinsicMatrix->data[0][2];
    solver->principalY = intrinsicMatrix->data[1][2];
    solver->focalX = intrinsicMatrix->data[0][0];
    solver->focalY = intrinsicMatrix->data[1][1];
}

void setIntrinsicParametersFromFloat(Perspective3Point* solver, const Matrix3x3* intrinsicMatrixFloat) {
    solver->principalX = (double)intrinsicMatrixFloat->data[0][2];
    solver->principalY = (double)intrinsicMatrixFloat->data[1][2];
    solver->focalX = (double)intrinsicMatrixFloat->data[0][0];
    solver->focalY = (double)intrinsicMatrixFloat->data[1][1];
}

void computeInverseParameters(Perspective3Point* solver) {
    solver->invFocalX = 1.0 / solver->focalX;
    solver->invFocalY = 1.0 / solver->focalY;
    solver->principalXFocalX = solver->principalX / solver->focalX;
    solver->principalYFocalY = solver->principalY / solver->focalY;
}


void vect_cross(const double a[3], const double b[3], double result[3]) {
    result[0] = a[1] * b[2] - a[2] * b[1];
    result[1] = -(a[0] * b[2] - a[2] * b[0]);
    result[2] = a[0] * b[1] - a[1] * b[0];
}

double vect_dot(const double a[3], const double b[3]) {
    return a[0] * b[0] + a[1] * b[1] + a[2] * b[2];
}

double vect_norm(const double a[3]) {
    return sqrt(a[0] * a[0] + a[1] * a[1] + a[2] * a[2]);
}

void vect_sub(const double a[3], const double b[3], double result[3]) {
    result[0] = a[0] - b[0];
    result[1] = a[1] - b[1];
    result[2] = a[2] - b[2];
}

void vect_divide(const double a[3], const double d, double result[3]) {
    result[0] = a[0] / d;
    result[1] = a[1] / d;
    result[2] = a[2] / d;
}

void vect_scale(const double s, const double *a, double *result) {
    result[0] = a[0] * s;
    result[1] = a[1] * s;
    result[2] = a[2] * s;
}

void mat_mult(const double a[3][3], const double b[3][3], double result[3][3]) {
    result[0][0] = a[0][0] * b[0][0] + a[0][1] * b[1][0] + a[0][2] * b[2][0];
    result[0][1] = a[0][0] * b[0][1] + a[0][1] * b[1][1] + a[0][2] * b[2][1];
    result[0][2] = a[0][0] * b[0][2] + a[0][1] * b[1][2] + a[0][2] * b[2][2];

    result[1][0] = a[1][0] * b[0][0] + a[1][1] * b[1][0] + a[1][2] * b[2][0];
    result[1][1] = a[1][0] * b[0][1] + a[1][1] * b[1][1] + a[1][2] * b[2][1];
    result[1][2] = a[1][0] * b[0][2] + a[1][1] * b[1][2] + a[1][2] * b[2][2];

    result[2][0] = a[2][0] * b[0][0] + a[2][1] * b[1][0] + a[2][2] * b[2][0];
    result[2][1] = a[2][0] * b[0][1] + a[2][1] * b[1][1] + a[2][2] * b[2][1];
    result[2][2] = a[2][0] * b[0][2] + a[2][1] * b[1][2] + a[2][2] * b[2][2];
}

void swap_double(double* a, double* b) {
        double temp = *a;
        *a = *b;
        *b = temp;
    }

void swap_matrices(double A[3][3], double B[3][3]) {
        double temp[3][3];
        for (int i = 0; i < 3; i++) {
            for (int j = 0; j < 3; j++) {
                temp[i][j] = A[i][j];
                A[i][j] = B[i][j];
                B[i][j] = temp[i][j];
            }
        }
    }

void initializeAP3PWithIntrinsicMatrix(Perspective3Point* solver, const Matrix3x3* intrinsicMatrix, int matrixDepth) {
    
        setIntrinsicParametersFromFloat(solver, intrinsicMatrix);
     
        setIntrinsicParameters(solver, intrinsicMatrix);
        computeInverseParameters(solver);
}

int countValidImagePoints32Bit() {
    int count = 0;
    for (int i = 0; i < MAX_NUM_POINTS; i++) {
        if (worldPointsFloat[i][0] != 0.0f || worldPointsFloat[i][1] != 0.0f || worldPointsFloat[i][2] != 0.0f) {
            count++;
        }
    }
    return count;
}

int countValidImagePoints64Bit() {
    int count = 0;
    for (int i = 0; i < MAX_NUM_POINTS; i++) {
        // Assuming a point is invalid if all its coordinates are 0
        if (worldPointsDouble[i][0] != 0.0 || worldPointsDouble[i][1] != 0.0 || worldPointsDouble[i][2] != 0.0) {
            count++;
        }
    }
    return count;
}

void extract3DPoints(ThreeDPoint* worldPoints, TwoDPoint* imagePoints, double* pointsArray) {
    
    memset(pointsArray, 0, sizeof(double) * MAX_NUM_POINTS * 5);     
    int validPointsCount = (countValidImagePoints32Bit() > countValidImagePoints64Bit()) ? countValidImagePoints32Bit() : countValidImagePoints64Bit();
    
    for (int idx = 0; idx < validPointsCount; idx++) {
        pointsArray[idx * 5] = imagePoints[idx].x * focalX + principalX;
        pointsArray[idx * 5 + 1] = imagePoints[idx].y * focalY + principalY;
        pointsArray[idx * 5 + 2] = worldPoints[idx].x;
        pointsArray[idx * 5 + 3] = worldPoints[idx].y;
        pointsArray[idx * 5 + 4] = worldPoints[idx].z;
    }
    // Fill the array with zeros for the p3p case
    for (int idx = validPointsCount; idx < 4; idx++) {
        for (int j = 0; j < 5; j++) {
            pointsArray[idx * 5 + j] = 0;
        }
    }
}

int computeMultiplePoses(const double featureVectors[3][4], 
                         const double worldPoints[3][4], 
                         double rotationSolutions[4][3][3], 
                         double translationSolutions[4][3], 
                         bool isP4P) {
    //world point vectors
    double worldPoint1[3] = {worldPoints[0][0], worldPoints[1][0], worldPoints[2][0]};
    double worldPoint2[3] = {worldPoints[0][1], worldPoints[1][1], worldPoints[2][1]};
    double worldPoint3[3] = {worldPoints[0][2], worldPoints[1][2], worldPoints[2][2]};
    
    // vector subtraction
    double vectorDifference[3];
    vect_sub(worldPoint1, worldPoint2, vectorDifference);

    double vectorMagnitude = vect_norm(vectorDifference);
    double normalizedVector[3];
    vect_divide(vectorDifference, vectorMagnitude, normalizedVector);

    // Feature Vectors
    double featureVector1[3] = {featureVectors[0][0], featureVectors[1][0], featureVectors[2][0]};
    double featureVector2[3] = {featureVectors[0][1], featureVectors[1][1], featureVectors[2][1]};
    double featureVector3[3] = {featureVectors[0][2], featureVectors[1][2], featureVectors[2][2]};
    
    // crosss product of featureeVector1 and 2
    double crossProductVector[3];
    vect_cross(featureVector1, featureVector2, crossProductVector);

    double crossProductMagnitude = vect_norm(crossProductVector);
    vect_divide(crossProductVector, crossProductMagnitude, crossProductVector);

    double tzVector[3];
    vect_cross(featureVector1, crossProductVector, tzVector);

    // Cross Products for feature vectors
    double crossProduct1[3];
    vect_cross(featureVector1, featureVector3, crossProduct1);
    double crossProduct2[3];
    vect_cross(featureVector2, featureVector3, crossProduct2);

    // Vector subtraction between worldPoint1 and worldPoint3
    double vectorDifference2[3];
    vect_sub(worldPoint1, worldPoint3, vectorDifference2);

    // Coefficients for the polynomial equation
    double dotProduct1 = vect_dot(vectorDifference2, normalizedVector);
    double dotProduct2 = vect_dot(crossProductVector, featureVector3);

    // Coefficients for the first feature vector
    double coeffFeature1_1 = dotProduct2;
    double coeffFeature1_3 = vect_dot(crossProductVector, crossProduct1);
    double coeffFeature1_5 = -dotProduct1 * coeffFeature1_1;

     // Cross product of vectorDifference2 and normalizedVector
    double crossProductDiffNorm[3];
    vect_cross(vectorDifference2, normalizedVector, crossProductDiffNorm);

    // Magnitude of the cross product
    double crossProductMagnitude2 = vect_norm(crossProductDiffNorm);
    vect_divide(crossProductDiffNorm, crossProductMagnitude2, crossProductDiffNorm);

    // Adjust coefficients based on magnitude
    coeffFeature1_1 *= crossProductMagnitude2;
    coeffFeature1_3 *= crossProductMagnitude2;

    // Coefficients for the second feature vector
    double dotProductDifference = dotProduct1 - vectorMagnitude;
    double coeffFeature2_1 = vect_dot(tzVector, crossProduct2);
    double coeffFeature2_2 = crossProductMagnitude * dotProduct2;
    double coeffFeature2_3 = vect_dot(crossProductVector, crossProduct2);
    double coeffFeature2_4 = dotProductDifference * coeffFeature2_2;
    double coeffFeature2_5 = -dotProductDifference * coeffFeature2_1;

    // Adjust coefficients based on magnitude
    coeffFeature2_1 *= crossProductMagnitude2;
    coeffFeature2_2 *= crossProductMagnitude2;
    coeffFeature2_3 *= crossProductMagnitude2;

    // Coefficients for the polynomial equation
    double polyCoeff1 = coeffFeature1_3 * coeffFeature2_2;
    double polyCoeff2 = coeffFeature1_3 * coeffFeature2_5 - coeffFeature1_5 * coeffFeature2_3;
    double polyCoeff3 = coeffFeature1_1 * coeffFeature2_3 - coeffFeature1_3 * coeffFeature2_1;
    double polyCoeff4 = -coeffFeature1_3 * coeffFeature2_4;
    double polyCoeff5 = coeffFeature1_1 * coeffFeature2_2;
    double polyCoeff6 = coeffFeature1_1 * coeffFeature2_5 - coeffFeature1_5 * coeffFeature2_1;
    double polyCoeff7 = -coeffFeature1_5 * coeffFeature2_4;

    // Polynomial coefficients array
    double polynomialCoefficients[5] = {polyCoeff5 * polyCoeff5 + polyCoeff1 * polyCoeff1 + polyCoeff3 * polyCoeff3,
                                        2 * (polyCoeff5 * polyCoeff6 + polyCoeff1 * polyCoeff2 + polyCoeff3 * polyCoeff4),
                                        polyCoeff6 * polyCoeff6 + 2 * polyCoeff5 * polyCoeff7 + polyCoeff2 * polyCoeff2 + polyCoeff4 * polyCoeff4 - polyCoeff1 * polyCoeff1 - polyCoeff3 * polyCoeff3,
                                        2 * (polyCoeff6 * polyCoeff7 - polyCoeff1 * polyCoeff2 - polyCoeff3 * polyCoeff4),
                                        polyCoeff7 * polyCoeff7 - polyCoeff2 * polyCoeff2 - polyCoeff4 * polyCoeff4};

    // Solutions array for the polynomial equation
    double solutions[4];
    solveQuartic(polynomialCoefficients, solutions);
    refineQuarticRoots(polynomialCoefficients, solutions);

    // Cross product of normalizedVector and crossProductDiffNorm
    double crossProductNormDiff[3];
    vect_cross(normalizedVector, crossProductDiffNorm, crossProductNormDiff);

    // Matrix combining normalizedVector, crossProductDiffNorm, and crossProductNormDiff
    double matrixNormDiffCross[3][3] =
            {{normalizedVector[0], crossProductDiffNorm[0], crossProductNormDiff[0]},
             {normalizedVector[1], crossProductDiffNorm[1], crossProductNormDiff[1]},
             {normalizedVector[2], crossProductDiffNorm[2], crossProductNormDiff[2]}};

    // Matrix combining featureVector1, crossProductVector, and tzVector
    double matrixFeatureCrossTz[3][3] =
            {{featureVector1[0], featureVector1[1], featureVector1[2]},
             {crossProductVector[0], crossProductVector[1], crossProductVector[2]},
             {tzVector[0], tzVector[1], tzVector[2]}};

    // Scaled featureVector3 based on delta and dotProduct2
    double scaledFeatureVector3[3];
    vect_scale((crossProductMagnitude2 / dotProduct2), featureVector3, scaledFeatureVector3);

    // Fourth world point coordinates
    double worldPoint4X = worldPoints[0][3];
    double worldPoint4Y = worldPoints[1][3];
    double worldPoint4Z = worldPoints[2][3];

    // Fourth feature vector coordinates
    double featureVector4U = featureVectors[0][3];
    double featureVector4V = featureVectors[1][3];

    // Array to store reprojection errors
    double reprojectionErrors[4];

    int solutionCount = 0;
    for (int i = 0; i < 4; ++i) {
        double cosTheta1Prime = solutions[i];
        if (abs(cosTheta1Prime) > 1)
            continue;
        double sinTheta1Prime = sqrt(1 - cosTheta1Prime * cosTheta1Prime);
        sinTheta1Prime = (dotProduct2 > 0) ? sinTheta1Prime : -sinTheta1Prime;

        double cosTheta3 = polyCoeff1 * cosTheta1Prime + polyCoeff2;
        double sinTheta3 = polyCoeff3 * cosTheta1Prime + polyCoeff4;
        double normTheta3 = sinTheta1Prime / ((polyCoeff5 * cosTheta1Prime + polyCoeff6) * cosTheta1Prime + polyCoeff7);
        cosTheta3 *= normTheta3;
        sinTheta3 *= normTheta3;

        double rotationMatrix13[3][3] =
                {{cosTheta3,            0,         -sinTheta3},
                 {sinTheta1Prime * sinTheta3, cosTheta1Prime,  sinTheta1Prime * cosTheta3},
                 {cosTheta1Prime * sinTheta3, -sinTheta1Prime, cosTheta1Prime * cosTheta3}};

        double intermediateMatrix[3][3];
        double finalRotationMatrix[3][3];
        mat_mult(matrixNormDiffCross, rotationMatrix13, intermediateMatrix);
        mat_mult(intermediateMatrix, matrixFeatureCrossTz, finalRotationMatrix);

        // Final rotation matrix applied to worldPoint3
        double rotatedWorldPoint3[3] =
                {worldPoint3[0] * finalRotationMatrix[0][0] + worldPoint3[1] * finalRotationMatrix[1][0] + worldPoint3[2] * finalRotationMatrix[2][0],
                 worldPoint3[0] * finalRotationMatrix[0][1] + worldPoint3[1] * finalRotationMatrix[1][1] + worldPoint3[2] * finalRotationMatrix[2][1],
                 worldPoint3[0] * finalRotationMatrix[0][2] + worldPoint3[1] * finalRotationMatrix[1][2] + worldPoint3[2] * finalRotationMatrix[2][2]};

        double scaledFeatureVector3Sin[3];
        vect_scale(sinTheta1Prime, scaledFeatureVector3, scaledFeatureVector3Sin);

        vect_sub(scaledFeatureVector3Sin, rotatedWorldPoint3, translationSolutions[solutionCount]);

        for (int j = 0; j < 3; j++) {
            for (int k = 0; k < 3; k++) {
                rotationSolutions[solutionCount][j][k] = finalRotationMatrix[j][k];
            }
        }


            // Check if P4P algorithm should be applied
    if (isP4P) {
        // Calculate the transformed coordinates of the fourth world point
        double transformedX3 = rotationSolutions[solutionCount][0][0] * worldPoint4X + rotationSolutions[solutionCount][0][1] * worldPoint4Y + rotationSolutions[solutionCount][0][2] * worldPoint4Z + translationSolutions[solutionCount][0];
        double transformedY3 = rotationSolutions[solutionCount][1][0] * worldPoint4X + rotationSolutions[solutionCount][1][1] * worldPoint4Y + rotationSolutions[solutionCount][1][2] * worldPoint4Z + translationSolutions[solutionCount][1];
        double transformedZ3 = rotationSolutions[solutionCount][2][0] * worldPoint4X + rotationSolutions[solutionCount][2][1] * worldPoint4Y + rotationSolutions[solutionCount][2][2] * worldPoint4Z + translationSolutions[solutionCount][2];
        
        // Calculate the projected feature vector coordinates
        double projectedFeatureU = transformedX3 / transformedZ3;
        double projectedFeatureV = transformedY3 / transformedZ3;
        
        // Calculate the reprojection error
        reprojectionErrors[solutionCount] = (projectedFeatureU - featureVector4U) * (projectedFeatureU - featureVector4U) + (projectedFeatureV - featureVector4V) * (projectedFeatureV - featureVector4V);
    }

    solutionCount++;

    // Sort the solutions based on reprojection errors if P4P algorithm is applied
    if (isP4P) {
        for (int i = 1; i < solutionCount; i++) {
            for (int j = i; j > 0 && reprojectionErrors[j-1] > reprojectionErrors[j]; j--) {
                // Swap reprojection errors
                double tempError = reprojectionErrors[j];
                reprojectionErrors[j] = reprojectionErrors[j-1];
                reprojectionErrors[j-1] = tempError;

                // Swap rotation matrices
                double tempRotation[3][3];
                for (int x = 0; x < 3; x++) {
                    for (int y = 0; y < 3; y++) {
                        tempRotation[x][y] = rotationSolutions[j][x][y];
                        rotationSolutions[j][x][y] = rotationSolutions[j-1][x][y];
                        rotationSolutions[j-1][x][y] = tempRotation[x][y];
                    }
                }

                // Swap translation vectors
                double tempTranslation[3];
                for (int x = 0; x < 3; x++) {
                    tempTranslation[x] = translationSolutions[j][x];
                    translationSolutions[j][x] = translationSolutions[j-1][x];
                    translationSolutions[j-1][x] = tempTranslation[x];
                }
            }
        }
    }

    return solutionCount;
}}

bool solveForAP3P(Perspective3Point* solver, Matrix3x3* rotation, Vector3* translation, const Matrix3x3* worldPointsMatrix, const Matrix3x3* imagePointsMatrix) {
    Matrix3x3 computedRotation = {{0}};
    Vector3 computedTranslation = {0};
    double extractedPoints[20] = {0}; 

    // if (getMatrixDepth(worldPointsMatrix) == getMatrixDepth(imagePointsMatrix)) {
    //     if (getMatrixDepth(worldPointsMatrix) == countValidImagePoints32Bit())
    //         extract3DPoints((ThreeDPoint*)worldPointsMatrix, (TwoDPoint*)imagePointsMatrix, extractedPoints);
    //     else
    //     extract3DPoints((ThreeDPoint*)worldPointsMatrix, (TwoDPoint*)imagePointsMatrix, extractedPoints);
    // } else if (getMatrixDepth(worldPointsMatrix) == countValidImagePoints32Bit()) {
    //     extract3DPoints((ThreeDPoint*)worldPointsMatrix, (TwoDPoint*)imagePointsMatrix, extractedPoints);
    // } else {
    //     extract3DPoints((ThreeDPoint*)worldPointsMatrix, (TwoDPoint*)imagePointsMatrix, extractedPoints);
    // }
    extract3DPoints((ThreeDPoint*)worldPointsMatrix, (TwoDPoint*)imagePointsMatrix, extractedPoints);

    bool isSuccessful = solveSingleAP3P(solver, &computedRotation, &computedTranslation,
                        extractedPoints[0], extractedPoints[1], extractedPoints[2], extractedPoints[3], extractedPoints[4],
                        extractedPoints[5], extractedPoints[6], extractedPoints[7], extractedPoints[8], extractedPoints[9],
                        extractedPoints[10], extractedPoints[11], extractedPoints[12], extractedPoints[13], extractedPoints[14],
                        extractedPoints[15], extractedPoints[16], extractedPoints[17], extractedPoints[18], extractedPoints[19]);

    // Copying data to the output matrices
    for (int i = 0; i < 3; i++) {
         translation->data[i] = computedTranslation.data[i];
         for (int j = 0; j < 3; j++) {
             rotation->data[i][j] = computedRotation.data[i][j];
         }
    }

    return isSuccessful;
}


int solveMultipleAP3P(Perspective3Point* solver, Matrix3x3* rotations, Vector3* translations, double* worldPoints, double* imagePoints, double focalX, double principalX, double focalY, double principalY) {
    Matrix3x3 computedRotations[4];
    Vector3 computedTranslations[4];
    double extractedPoints[MAX_NUM_POINTS] = {0};

    if (countValidImagePoints32Bit() == countValidImagePoints64Bit()) {
        extract3DPoints((ThreeDPoint*)worldPoints, (TwoDPoint*)imagePoints, extractedPoints);
    } else if (countValidImagePoints32Bit() > countValidImagePoints64Bit()) {
        extract3DPoints((ThreeDPoint*)worldPoints, (TwoDPoint*)imagePoints, extractedPoints);
    } else {
        extract3DPoints((ThreeDPoint*)worldPoints, (TwoDPoint*)imagePoints, extractedPoints);
    }

    bool isP4PCase = (countValidImagePoints32Bit() == 4 || countValidImagePoints64Bit() == 4);
    int solutionCount = solveAP3P(solver, computedRotations, computedTranslations,
                                       extractedPoints[0], extractedPoints[1], extractedPoints[2], extractedPoints[3], extractedPoints[4],
                                       extractedPoints[5], extractedPoints[6], extractedPoints[7], extractedPoints[8], extractedPoints[9],
                                       extractedPoints[10], extractedPoints[11], extractedPoints[12], extractedPoints[13], extractedPoints[14],
                                       extractedPoints[15], extractedPoints[16], extractedPoints[17], extractedPoints[18], extractedPoints[19],
                                       isP4PCase);

    for (int i = 0; i < solutionCount; i++) {
        rotations[i] = computedRotations[i];
        translations[i] = computedTranslations[i];
    }

    return solutionCount;
}

bool solveSingleAP3P(Perspective3Point* solver, Matrix3x3* rotation, Vector3* translation, 
                double imageX0, double imageY0, double worldX0, double worldY0, double worldZ0,
                double imageX1, double imageY1, double worldX1, double worldY1, double worldZ1,
                double imageX2, double imageY2, double worldX2, double worldY2, double worldZ2,
                double imageX3, double imageY3, double worldX3, double worldY3, double worldZ3) {
    double Rs[4][3][3] = {{0}};
    double ts[4][3] = {{0}};

    bool p4p = true;
    
    int n = solveAP3P(solver, Rs, ts, 
                           imageX0, imageY0, worldX0, worldY0, worldZ0, 
                           imageX1, imageY1, worldX1, worldY1, worldZ1, 
                           imageX2, imageY2, worldX2, worldY2, worldZ2, 
                           imageX3, imageY3, worldX3, worldY3, worldZ3, p4p);
    if (n == 0)
        return false;

    for (int i = 0; i < 3; i++) {
        for (int j = 0; j < 3; j++)
            rotation->data[i][j] = Rs[0][i][j];
        translation->data[i] = ts[0][i];
    }

    return true;
}


                                          
int solveAP3P(Perspective3Point* solver, Matrix3x3 rotation[4], Vector3 translation[4],
                double imageX0, double imageY0, double worldX0, double worldY0, double worldZ0,
                double imageX1, double imageY1, double worldX1, double worldY1, double worldZ1,
                double imageX2, double imageY2, double worldX2, double worldY2, double worldZ2,
                double imageX3, double imageY3, double worldX3, double worldY3, double worldZ3,
                bool isP4P) {
    double normalizationFactor0, normalizationFactor1, normalizationFactor2;
    double norm;

    imageX0 = solver->invFocalX * imageX0 - solver->principalXFocalX;
    imageY0 = solver->invFocalY * imageY0 - solver->principalYFocalY;
    norm = sqrt(imageX0 * imageX0 + imageY0 * imageY0 + 1);
    normalizationFactor0 = 1.0 / norm;
    imageX0 *= normalizationFactor0;
    imageY0 *= normalizationFactor0;

    imageX1 = solver->invFocalX * imageX1 - solver->principalXFocalX;
    imageY1 = solver->invFocalY * imageY1 - solver->principalYFocalY;
    norm = sqrt(imageX1 * imageX1 + imageY1 * imageY1 + 1);
    normalizationFactor1 = 1.0 / norm;
    imageX1 *= normalizationFactor1;
    imageY1 *= normalizationFactor1;

    imageX2 = solver->invFocalX * imageX2 - solver->principalXFocalX;
    imageY2 = solver->invFocalY * imageY2 - solver->principalYFocalY;
    norm = sqrt(imageX2 * imageX2 + imageY2 * imageY2 + 1);
    normalizationFactor2 = 1.0 / norm;
    imageX2 *= normalizationFactor2;
    imageY2 *= normalizationFactor2;

    imageX3 = solver->invFocalX * imageX3 - solver->principalXFocalX;
    imageY3 = solver->invFocalY * imageY3 - solver->principalYFocalY;
    double normalizationFactor3 = 1.0; // Not used

    double featureVectors[3][4] = {{imageX0, imageX1, imageX2, imageX3},
                                   {imageY0, imageY1, imageY2, imageY3},
                                   {normalizationFactor0, normalizationFactor1, normalizationFactor2, normalizationFactor3}};
    double worldPoints[3][4] = {{worldX0, worldX1, worldX2, worldX3},
                                {worldY0, worldY1, worldY2, worldY3},
                                {worldZ0, worldZ1, worldZ2, worldZ3}};

    return computeMultiplePoses(featureVectors, worldPoints, rotation, translation, isP4P);
}