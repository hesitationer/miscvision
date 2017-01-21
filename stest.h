#include <cv/image.h>

void stest(cv::Image<uchar, 3> &rgb, CvRect roi){
	cv::Image<float,1> imf;
	cv::Image<uchar,1> intensity;
	cv::Image<uchar,3> smrgb(roi.width, roi.height);

	// convert to greyscale 
	rgb.getSubImage(roi, smrgb);
	intensity.convert(smrgb);

	imf.convert(intensity);

	intensity.show("win");
	cvWaitKey(1000);

	double sum=0;
	for(int j=0; j<imf.height; j++){
		float * row = imf[j];
		for(int i=0; i<imf.width; i++){
			sum+=row[i];
		}
	}
	printf("sum=%f\n", sum);

}
void stest(uchar * rgb_template, int w, int h, int roi_x,int roi_y, int roi_w, int roi_h){
	cv::Image<uchar,3> rgb;

	rgb.setData(rgb_template, w, h);
	std::cout<<"2"<<std::endl;
	CvRect roi = cvRect(roi_x, roi_y, roi_w, roi_h);	
	
	stest(rgb, roi);
	
	rgb.setData(0);
}
