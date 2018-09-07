#ifndef CURVE_H
#define CURVE_H

#include <string>
#include <iostream>
#include <vector>
#include <math.h>
#include <stdlib.h>
#include <assert.h>
#include <stdio.h>
#include <time.h>
#include "cv.h"
#include "highgui.h"

using namespace std;

bool polynomial_curve_fit(std::vector<cv::Point>& key_point, int n, cv::Mat& A);

#endif