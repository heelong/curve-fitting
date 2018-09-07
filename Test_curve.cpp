#include "curve.h"
#include "cv.h"
#include <time.h>
#include <windows.h>
void main()
{
	for (int k = 0; k < 1000; k++)
	{
		//创建用于绘制的深蓝色背景图像
		cv::Mat image = cv::Mat::zeros(480, 640, CV_8UC3);
		image.setTo(cv::Scalar(100, 0, 0));
		//输入拟合点  
		std::vector<cv::Point> points;
		srand((unsigned)time(NULL));
		double a = 0, b = 100;
		for (int x = 0; x < 50; x++)
		{
			double x_ = static_cast<double>(x)*40.0 + ((double)rand() / RAND_MAX)*(b - a) + a;
			double y_ = -287.699263738501 + 5.917786427084218*x_ - 0.03171264283596321*x_*x_ +
				0.00005776328905868311*x_*x_*x_;
			y_ += ((double)rand() / RAND_MAX)*(b - a) + a;
			points.push_back(cv::Point(x_, y_));
		}
		
		//points.push_back(cv::Point(110., 60.));
		//points.push_back(cv::Point(150., 70.));
		//points.push_back(cv::Point(200., 90.));
		//points.push_back(cv::Point(220., 95.));
		//points.push_back(cv::Point(230., 105.));
		//points.push_back(cv::Point(252., 140.));
		//points.push_back(cv::Point(301., 200.));
		//points.push_back(cv::Point(300., 220.));
		//points.push_back(cv::Point(350., 400.));
		//points.push_back(cv::Point(110., 48.));
		//points.push_back(cv::Point(117., 68.));
		//points.push_back(cv::Point(140., 60.));
		//points.push_back(cv::Point(212., 98.));
		//points.push_back(cv::Point(243., 90.));
		//points.push_back(cv::Point(233., 115.));
		//points.push_back(cv::Point(272., 120.));
		//points.push_back(cv::Point(321., 240.));
		//points.push_back(cv::Point(330., 210.));
		//points.push_back(cv::Point(342., 370.));

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