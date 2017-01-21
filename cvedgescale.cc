#include <cv/image.h>

void edge(cv::Image<float, 1> & im, cv::Image<uchar, 1> & res, int width){
	cv::Image<float,1> Ix,Iy,Ixx,Iyy,Ixy, Ixxx,Ixxy,Ixyy,Iyyy;
	cv::Image<float,1> L;
	cv::Image<uchar,1> cmp1,cmp2;
	cmp1.reshape(res);
	cmp2.reshape(res);
	L.reshape(im);
	Ix.reshape(im);
	Iy.reshape(im);
	Ixx.reshape(im);
	Ixy.reshape(im);
	Iyy.reshape(im);
	Ixxx.reshape(im);
	Iyyy.reshape(im);
	Ixxy.reshape(im);
	Ixyy.reshape(im);

	cvSmooth(&im, &L, CV_GAUSSIAN, width, width, 0);
	cvSobel(&L, &Ix, 1, 0, 3);  Ix = Ix * (1.0/9);
	Ix.imagesc("im"); cvWaitKey(-1);
	cvSobel(&L, &Iy, 0, 1, 3);Iy = Iy * (1.0/9);
	Iy.imagesc("im"); cvWaitKey(-1);
	cvSobel(&Ix, &Ixx, 1, 0, 3);Ixx = Ixx * (1.0/9);
	Ixx.imagesc("im"); cvWaitKey(-1);
	cvSobel(&Iy, &Iyy, 0, 1, 3);Iy = Iy * (1.0/9);
	Iyy.imagesc("im"); cvWaitKey(-1);
	cvSobel(&Ix, &Ixy, 0, 1, 3);Ixy = Ixy * (1.0/9);
	Ixy.imagesc("im"); cvWaitKey(-1);
	cvSobel(&Ixx, &Ixxx, 1, 0, 3);Ixxx = Ixxx * (1.0/9);
	Ixxx.imagesc("im"); cvWaitKey(-1);
	cvSobel(&Iyy, &Iyyy, 0, 1, 3);Iyyy = Iyyy * (1.0/9);
	Iyyy.imagesc("im"); cvWaitKey(-1);
	cvSobel(&Ixx, &Ixxy, 0, 1, 3);Ixxy = Ixxy * (1.0/9);
	Ixxy.imagesc("im"); cvWaitKey(-1);
	cvSobel(&Iyy, &Ixyy, 1, 0, 3);Ixyy = Ixyy * (1.0/9);
	Ixyy.imagesc("im"); cvWaitKey(-1);

	L = (Ix*Ix*Ixx) + (Ix*Iy*Ixy*2);
	L = L + (Iy*Iy*Iyy);
	L.imagesc("edge"); cvWaitKey(-1);
	cvCmpS(&L, 0, &cmp1, CV_CMP_EQ);
	L = (Ix*Ix*Ix*Ixx) + (Ix*Ix*Iy*Ixxy*3);
	L = L + (Ix*Iy*Iy*Ixyy*3);
	L = L + (Iy*Iy*Iy*Iyyy);
	L.imagesc("edge"); cvWaitKey(-1);
	cvCmpS(&L, 0, &cmp2, CV_CMP_LT);

	cvAnd(&cmp1, &cmp2, &res);
}

int main(int argc, char ** argv){
	cv::Image<uchar, 3> im3;
	cv::Image<uchar, 1> res;
	cv::Image<float, 1> imf;
	
	im3.open(argv[1]);
	res.reshape(im3);
	imf = im3;
	cvNamedWindow("im", 0);
	cvNamedWindow("win", 0);
	cvNamedWindow("edge", 0);
	im3.show("win");
	for(int width=3; width<33; width+=2){
		edge(imf, res, width);
		res.imagesc("edge");
		cvWaitKey(-1);
	}
}
