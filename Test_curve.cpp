#include "curve.h"
#include "cv.h"
#include <time.h>
#include <windows.h>
#include"opencv2/opencv.hpp"
#include <iostream>
#include <fstream>
#include <iostream>
#include <sstream>


#include <vector>
#include <stack>

using namespace cv;
using namespace std;
void curveFit()
{
	for (int k = 0; k < 1000; k++)
	{
		//创建用于绘制的深蓝色背景图像
		cv::Mat image = cv::Mat::zeros(480, 640, CV_8UC3);
		image.setTo(cv::Scalar(100, 0, 0));
		//输入拟合点  
		std::vector<cv::Point2d> points;
		srand((unsigned)time(NULL));
		double a = 0, b = 100;
		for (int x = 0; x < 50; x++)
		{
			double x_ = static_cast<double>(x)*40.0 + ((double)rand() / RAND_MAX)*(b - a) + a;
			double y_ = -287.699263738501 + 5.917786427084218*x_ - 0.03171264283596321*x_*x_ +
				0.00005776328905868311*x_*x_*x_;
			y_ += ((double)rand() / RAND_MAX)*(b - a) + a;
			points.push_back(cv::Point2d(x_, y_));
		}

		//将拟合点绘制到空白图上  
		for (int i = 0; i < points.size(); i++)
		{
			cv::circle(image, points[i], 5, cv::Scalar(0, 0, 255), 2, 8, 0);
		}

		//绘制折线
		//cv::polylines(image, points, false, cv::Scalar(0, 255, 0), 1, 8, 0);

		cv::Mat A;
		clock_t time1, time2;

		time1 = clock();
		polynomial_curve_fit(points, 3, A);
		time2 = clock();
		cout << "Time curve fiting:" << (double)(time2 - time1) / CLOCKS_PER_SEC*1000.0 << "ms" << endl;

		//std::cout << "A = " << A << std::endl;

		std::vector<cv::Point> points_fitted;

		for (int x = 0; x < 400; x++)
		{
			double y = A.at<double>(0, 0) + A.at<double>(1, 0) * x +
				A.at<double>(2, 0)*std::pow(x, 2) + A.at<double>(3, 0)*std::pow(x, 3);

			points_fitted.push_back(cv::Point(x, y));
		}
		cv::polylines(image, points_fitted, false, cv::Scalar(0, 255, 255), 1, 8, 0);

		cv::imshow("image", image);
		//cv::waitKey(0);
		cv::waitKey(3);
	}
	cv::waitKey(0);
	return;

}
void PerformanceTesting()
{
	clock_t time1, time2;
	int size_ = 10000000;
	double*max_z = new double[size_]();
	double*min_z = new double[size_]();

	double a = 0, b = 100;
	while (true)
	{
		memset(max_z, -100, sizeof(double) * size_);
		memset(min_z, 100, sizeof(double) * size_);
		time1 = clock();

		for (int i = 0; i < size_; i++)
		{
			double de = ((double)rand() / RAND_MAX)*(b - a) + a;
			max_z[i] = max_z[i] < de ? de : max_z[i];
			min_z[i] = min_z[i] > de ? de : min_z[i];
		}
		time2 = clock();
		cout << "Time 1 :" << (double)(time2 - time1) / CLOCKS_PER_SEC*1000.0 << "ms" << endl << endl;;
		memset(max_z, -100, sizeof(double) * size_);
		memset(min_z, 100, sizeof(double) * size_);
		time1 = clock();

		for (int i = 0; i < size_; i++)
		{
			double de = ((double)rand() / RAND_MAX)*(b - a) + a;
			if (max_z[i] < de)
				max_z[i] = de;
			if (min_z[i] > de)
				min_z[i] = de;
		}
		time2 = clock();
		cout << "Time 2 :" << (double)(time2 - time1) / CLOCKS_PER_SEC*1000.0 << "ms" << endl;
	}

	delete[]max_z;
	delete[]min_z;
	cv::waitKey(0);
}

int minAreaRectTest()
{
	Mat img(500, 500, CV_8UC3);
	RNG rng = theRNG();//随机数类

	for (;;)
	{
		int  count = rng.uniform(1, 101);//uniform()返回指定范围的随机数
		vector<Point>points;
		for (int i = 0; i < count; i++)
		{
			Point pt;
			pt.x = rng.uniform(img.cols / 4, img.cols * 3 / 4);
			pt.y = rng.uniform(img.rows / 4, img.rows * 3 / 4);
			points.push_back(pt);
		}
		clock_t time1, time2;

		time1 = clock();
		RotatedRect box = minAreaRect(Mat(points));//点集的最小外接旋转矩形
		time2 = clock();
		cout << "Time box fiting:" << (double)(time2 - time1) / CLOCKS_PER_SEC*1000.0 << "ms" << endl;
		Point2f tr[4], center;
		float radius = 0;
		box.points(tr);
		minEnclosingCircle(Mat(points), center, radius);//点集的最小外接圆
		img = Scalar::all(0);
		for (int i = 0; i < count; i++)
		{
			circle(img, points[i], 3, Scalar(0, 0, 255), CV_FILLED, CV_AA);
		}
		for (int i = 0; i < 4; i++)
		{
			line(img, tr[i], tr[(i + 1) % 4], Scalar(0, 255, 0), 1, CV_AA);
		}
		//circle(img, center, cvRound(radius), Scalar(0, 255, 255), 1, CV_AA);
		imshow("j", img);
		cv::waitKey(300);

	}
	return 0;
}





/*
Input:
src: 待检测连通域的二值化图像
Output:
dst: 标记后的图像
featherList: 连通域特征的清单
return：
连通域数量。
*/
int bwLabel(Mat & src, Mat & dst, vector<Feather> & featherList)
{
	int rows = src.rows;
	int cols = src.cols;

	int labelValue = 0;
	Point seed, neighbor;
	stack<Point> pointStack;    // 堆栈

	int area = 0;               // 用于计算连通域的面积
	int leftBoundary = 0;       // 连通域的左边界，即外接最小矩形的左边框，横坐标值，依此类推
	int rightBoundary = 0;
	int topBoundary = 0;
	int bottomBoundary = 0;
	Rect box;                   // 外接矩形框
	Feather feather;

	featherList.clear();    // 清除数组

	dst.release();
	dst = src.clone();
	for (int i = 0; i < rows; i++)
	{
		uchar *pRow = dst.ptr<uchar>(i);
		for (int j = 0; j < cols; j++)
		{
			if (pRow[j] == 255)
			{
				area = 0;
				labelValue++;           // labelValue最大为254，最小为1.
				seed = Point(j, i);     // Point（横坐标，纵坐标）
				dst.at<uchar>(seed) = labelValue;
				pointStack.push(seed);

				area++;
				leftBoundary = seed.x;
				rightBoundary = seed.x;
				topBoundary = seed.y;
				bottomBoundary = seed.y;

				while (!pointStack.empty())
				{
					neighbor = Point(seed.x + 1, seed.y);
					if ((seed.x != (cols - 1)) && (dst.at<uchar>(neighbor) == 255))
					{
						dst.at<uchar>(neighbor) = labelValue;
						pointStack.push(neighbor);

						area++;
						if (rightBoundary < neighbor.x)
							rightBoundary = neighbor.x;
					}

					neighbor = Point(seed.x, seed.y + 1);
					if ((seed.y != (rows - 1)) && (dst.at<uchar>(neighbor) == 255))
					{
						dst.at<uchar>(neighbor) = labelValue;
						pointStack.push(neighbor);

						area++;
						if (bottomBoundary < neighbor.y)
							bottomBoundary = neighbor.y;

					}

					neighbor = Point(seed.x - 1, seed.y);
					if ((seed.x != 0) && (dst.at<uchar>(neighbor) == 255))
					{
						dst.at<uchar>(neighbor) = labelValue;
						pointStack.push(neighbor);

						area++;
						if (leftBoundary > neighbor.x)
							leftBoundary = neighbor.x;
					}

					neighbor = Point(seed.x, seed.y - 1);
					if ((seed.y != 0) && (dst.at<uchar>(neighbor) == 255))
					{
						dst.at<uchar>(neighbor) = labelValue;
						pointStack.push(neighbor);

						area++;
						if (topBoundary > neighbor.y)
							topBoundary = neighbor.y;
					}

					seed = pointStack.top();
					pointStack.pop();
				}
				box = Rect(leftBoundary, topBoundary, rightBoundary - leftBoundary, bottomBoundary - topBoundary);
				//rectangle(src, box, 255);
				feather.area = area;
				feather.boundingbox = box;
				feather.label = labelValue;
				featherList.push_back(feather);
			}
		}
	}
	return labelValue;
}

int ConnectedComponent(Mat& src)
{
	if (src.empty())
		exit(-1);
	threshold(src, src, 127, 255, THRESH_BINARY);   // 二值化图像
	vector<Feather> featherList;                    // 存放连通域特征
	Mat dst;
	cout << "连通域数量： " << bwLabel(src, dst, featherList) << endl;

	// 为了方便观察，可以将label“放大”
	for (int i = 0; i < dst.rows; i++)
	{
		uchar *p = dst.ptr<uchar>(i);
		for (int j = 0; j < dst.cols; j++)
		{
			p[j] = 30 * p[j];
		}
	}

	cout << "标号" << "\t" << "面积" << endl;
	for (vector<Feather>::iterator it = featherList.begin(); it < featherList.end(); it++)
	{
		cout << it->label << "\t" << it->area << endl;
		rectangle(dst, it->boundingbox, 255);
	}

	imshow("src", src);
	imshow("dst", dst);

	cv::waitKey(100);
	return 0;
}

void ImageBasedObjDector(int file)
{
	stringstream ss;
	ss << "../data/" << file << "_LiDAR.txt";

	std::ifstream in_file_(ss.str(), std::ifstream::in);
	if (!in_file_.is_open()) {
		std::cerr << "Cannot open input file: " << ss.str() << std::endl;
		return;
	}
	std::string line;
	vector<PointsXYZInt>PointCloud; PointCloud.reserve(50000);
	getline(in_file_, line);
	PointsXYZInt tmp;
	while (getline(in_file_, line))
	{
		std::istringstream iss(line);
		iss >> tmp.x;
		iss >> tmp.y;
		iss >> tmp.z;
		iss >> tmp.index;
		iss >> tmp.horizontal_angle;
		iss >> tmp.intensity;
		iss >> tmp.ClassNum;
		iss >> tmp.Distance;
		iss >> tmp.Range_meter;
		iss >> tmp.Laser_Line_Number;
		iss >> tmp.sensorCount;
		PointCloud.push_back(tmp);
	}
	const int pointNum = PointCloud.size();
	int *Pixnum_02m = new int[pointNum]();  //0.2m 分辨率 
	int imgx = 0, imgy = 0;
	std::vector<cv::Point2d> pointsIndex; pointsIndex.resize(pointNum);
	for (int i = 0; i < pointNum; i += 1){
		imgx = pointsIndex[i].x = PointCloud[i].x / 20 + 500;//以1分米为单位
		imgy = pointsIndex[i].y = PointCloud[i].y / 20 + 500;
		Pixnum_02m[i] = imgy * 1000 + imgx;
	}

	//创建用于绘制的深蓝色背景图像
	cv::Mat image = cv::Mat::zeros(1000, 1000, CV_8UC1);
	image.setTo(cv::Scalar(0));
	//将点云投影到空白图上  
	for (int i = 0; i < pointNum; i++)
	{
		image.at <uchar>(Pixnum_02m[i]) = 255;
	}
	cv::imshow("image", image);
	//
	int SizeOfKernel = 5;
	cv::Mat kernel = cv::getStructuringElement(cv::MORPH_RECT, Size(SizeOfKernel, SizeOfKernel));
	cv::Mat dst1, dis2;
	//先膨胀
	cv::dilate(image, dst1, kernel);
	cv::imshow("result1", dst1);
	//再腐蚀
	SizeOfKernel = 5;
	kernel = cv::getStructuringElement(cv::MORPH_RECT, Size(SizeOfKernel, SizeOfKernel));
	cv::erode(dst1, dis2, kernel);
	cv::imshow("result2", dis2);

	//ConnectedComponent(dst1);
	cv::waitKey(100);
	//cv::waitKey(0);
	PointCloud.clear();
	pointsIndex.clear();
	delete[]Pixnum_02m;
	return;

}

int main(int argc, char** argv[])
{
	for (int i = 1; i < 6733;i++)
		ImageBasedObjDector(i);
	return 0;
}