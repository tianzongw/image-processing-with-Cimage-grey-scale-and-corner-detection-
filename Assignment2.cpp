#include "stdAfx.h"
#include "Assignment2.h"
#include <vector>
#include <iostream>
#include <fstream>
#include <cmath>
#include <math.h> 
#include <algorithm>

using namespace std;
const double PI = 3.1415926535;	

/////////////////////////////
// CCorner Source File //
/////////////////////////////


// this function convert a given CImage to a GrayScale image
void CCorner::RGBToGrayScale(CImage* pIn, CImage* pOut)
{
	//
	// INPUT:
	//     CImage* pIn:		The input image with 24bit depth
	//
	// OUTPUT:
	//     CImage* pOut:	The output image. It has ALREADY been initialized
	//                      with the same dimension as the input image (pIN) and 
	//                      formatted to 8bit depth (256 gray levels). So, please
	//                      use 'SetIndex' instead of 'SetRGB'.
	//
	
	// Begin your code here: //
	int width = pIn->GetWidth();
	byte r, g, b;
	int height = pIn->GetHeight();
	byte avg;
	byte grey;
	for (int i = 0; i < width; i++) {
		for (int j = 0; j < height; j++) {
			pIn->GetRGB(i, j, &r, &g, &b);
			grey = 0.3*(int)r + 0.59*(int)g + 0.11*(int)b;
			pOut->SetIndex(i, j, grey);
		}
	}
}

// this function obtains the corners from a given GrayScale image
void CCorner::ObtainCorners(CImage* pIn, double sigma, double threshold, vector<C2DPoint*>* pResultCorners)
{
//
	// INPUT:
	//     CImage* pIn:		The input grayscale image, which is exactly the output of
	//                      RGBToGrayScale function.
	//     double sigma:    The sigma value for your convolve function
	//     double threhold: The minimum value for qualifying a corner. Please refers to the lecture 3's notes
	//
	// OUTPUT:
	//     vector<C2DPoint*>* pResultCorners:	
	//                      A std::vector object that holds all detected corners. Each corner is represented by
	//                      a C2DPoint structure. An example is given below, showing how a corner object is
	//                      initialized and stored in pResultCorners:
	//                      
	//                      	C2DPoint* pPnt = new C2DPoint(x, y);
	//                      	pResultCorners.push_back(pPnt);
	//
	//
	
	// Begin your code here: //

	////
	// Step 1: Compute a proper size for Gaussian filter
		cin >> sigma;
		int size = (int)sigma*(pow(2 * log(1000), 0.5));

	////
	// Step 2: Define the Gaussian filter and partial filter
		//define Guassian filter
			double *gaus = new double[size];
			//for (int i = 0; i<size; i++)
			//{
			//	gaus[i] = new double[size];   
			//}
		
			int center = size / 2;
			double sum = 0;

			for (int i = 0; i < size; i++){
				gaus[i] = exp(-(i - center)*(i - center) / (2 * sigma*sigma));					
			}
			//for (int i = 0; i < size; i++) {
			//	for (int j = 0; j < size; j++) {
			//		gaus[i][j] = (1 / (2 * PI*sigma*sigma))*exp(-((i - center)*(i - center) + (j - center)*(j - center)) / (2 * sigma*sigma));
			//		sum += gaus[i][j];
			//	}
			//}

			//or (int i = 0; i < size; i++) {
			//	for (int j = 0; j < size; j++) {
			//		gaus[i][j] /= sum;
			//	}
			//}

		//define partial filter
		//	double *partial = new double[3];
		//	partial[0] = 0.5;
		//	partial[1] = 0;
		//	partial[2] = 0.5;

	////
	// Step 3: Compute Ix, Iy
		int width = pIn->GetWidth();
		int height = pIn->GetHeight();

		//calculate Ix
			double **Ix = new double *[height];
			for (int i = 0; i<size; i++)
			{
				Ix[i] = new double[width];
			}

			for (int i = 0; i < height; i++) {
				for (int j = 0; j < width; j++) {
					if (j == 0) {
						Ix[i][j] = pIn->GetIndex(i, j + 1) - pIn->GetIndex(i, j);
					}
					else if (j == width - 1) {
						Ix[i][j] = pIn->GetIndex(i, j) - pIn->GetIndex(i, j - 1);
					}
					else {
						Ix[i][j] = 0.5*pIn->GetIndex(i, j - 1) - 0.5*pIn->GetIndex(i, j + 1);
					}
				}
			}

		//calculate Iy
			double **Iy = new double *[height];
			for (int i = 0; i<size; i++)
			{
				Iy[i] = new double[width];
			}

			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					if (i == 0) {
						Iy[i][j] = pIn->GetIndex(i + 1,j) - pIn->GetIndex(i, j);
					}
					else if (i == height - 1) {
						Ix[i][j] = pIn->GetIndex(i, j) - pIn->GetIndex(i - 1, j);
					}
					else {
						Ix[i][j] = 0.5*pIn->GetIndex(i - 1, j) - 0.5*pIn->GetIndex(i + 1, j);
					}
				}
			}


	////
	// Step 4: Compute Ixx, Iyy, Ixy
			double **Ixx = new double *[height];
			double **Iyy = new double *[height];
			double **Ixy = new double *[height];
			
			for (int i = 0; i< height; i++)
			{
				Ixx[i] = new double[width];
				Iyy[i] = new double[width];
				Ixy[i] = new double[width];
			}

			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					Ixx[i][j] = pow(Ix[i][j], 2);
					Iyy[i][j] = pow(Iy[i][j], 2);
					Ixy[i][j] = Ix[i][j] * Iy[i][j];
				}
			}

	////
	// Step 5: Smooth Ixx, Iyy, Ixy
			for (int i = 1; i < width; i++) {
				for (int j = 1; j < height; j++) {
					double temp1xx = 0.0;
					double temp1yy = 0.0;
					double temp1xy = 0.0;
					//1DConvolution wrt x
					for (int k = 0; k < size; k++) {
						if (j - center >= 0) {
							temp1xx += Ixx[i][j + (k - center)] * gaus[k];
							temp1yy += Iyy[i][j + (k - center)] * gaus[k];
							temp1xy += Ixy[i][j + (k - center)] * gaus[k];
						}
						Ixx[i][j] = temp1xx;
						Iyy[i][j] = temp1yy;
						Ixy[i][j] = temp1xy;
					}
				}
			}

			for (int j = 1; j < height; j++) {
				for (int i = 1; i < width; i++) {
					double temp2xx = 0.0;
					double temp2yy = 0.0;
					double temp2xy = 0.0;
					//1DConvolution wrt y
					for (int k = 0; k < size; k++) {
						if (i - center >= 0) {
							temp2xx += Ixx[i + (k - center)][j] * gaus[k];
							temp2yy += Iyy[i + (k - center)][j] * gaus[k];
							temp2xy += Ixy[i + (k - center)][j] * gaus[k];
						}
						Ixx[i][j] = temp2xx;
						Iyy[i][j] = temp2yy; 
						Ixy[i][j] = temp2xy;

					}
				}
			}
	////
	// Step 6: Compute R
			double **R = new double *[height];
			for (int i = 0; i < height; i++) {
				R[i] = new double[width];
			}

			double k = 0.04;	
			for (int i = 0; i < width; i++) {
				for (int j = 0; j < height; j++) {
					R[i][j] = (Ixx[i][j] * Iyy[i][j] - Ixy[i][j] * Ixy[i][j]) - 0.04*pow(Ixx[i][j]+Iyy[i][j],2);
				}
			}
	////
	// Step 7: Locate Maxima in R
			vector<pair<int, int>> cordinatesR;
			for (int i = 1; i < height; i++) {
				for (int j = 1; j < width; j++) {
					if (R[i][j] > max(R[i-1][j],R[i][j+1]) && R[i][j] > max(R[i+1][j], R[i][j-1])) {
						cordinatesR.push_back(make_pair(i, j));
					}
				}
			}
	////
	// Step 8: Compute corner candidates up to sub-pixel accuracy and interpolate R value for corner candidates.
			vector<pair<int, int>>::iterator itr;
			vector<pair<int, int>> cordinatesCimg;
			for (itr = cordinatesR.begin(); itr != cordinatesR.end(); itr++) {
				double a = 0.5*(R[itr->first][itr->second - 1] + R[itr->first][itr->second + 1] - 2 * R[itr->first][itr->second]);
				double b = 0.5*(R[itr->first - 1][itr->second] + R[itr->first + 1][itr->second] - 2 * R[itr->first][itr->second]);
				double c = 0.5*(R[itr->first][itr->second + 1] - R[itr->first][itr->second - 1]);
				double d = 0.5*(R[itr->first + 1][itr->second] - R[itr->first - 1][itr->second]);
				int x = -1 * c / (2 * a) + itr->first;
				int y = -1 * d / (2 * b) + itr->second;
				cordinatesCimg.push_back(make_pair(x, y));
			}

	////
	// Step 9: Use the threshold value to identify strong corners for output

}
	// OUTPUT:
	//     vector<C2DPoint*>* pResultCorners:	
	//                      A std::vector object that holds all detected corners. Each corner is represented by
	//                      a C2DPoint structure. An example is given below, showing how a corner object is
	//                      initialized and stored in pResultCorners:
	//                      
	//                      	C2DPoint* pPnt = new C2DPoint(x, y);
	//                      	pResultCorners.push_back(pPnt);

