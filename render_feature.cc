#include <highgui.h>
#include <cv/dog_feature_detector.hh>
#include <cv/wincompat.h>

int main(int argc, char ** argv){
	
	DoGFeatureDetector detector(10.0, 20, 10.0, 1.6, 16, 1);
	cv::Image<unsigned char, 3> im3;
	cv::Image<float, 3> im3f;
	cv::Image<float, 3> dest;
	cv::Image<unsigned char> im;
	cv::Image<float, 3> mask;
	cv::Image<float> imf;
	cvNamedWindow("p", 1);
	for(int i=1; i<argc; i++){
		std::cout<<"Opening "<<argv[i]<<std::endl;
		im.open(argv[i]);
		im3.open(argv[i]);
		std::cout<<"Convert to float: "<<argv[i]<<std::endl;	
		imf = im;
		std::cout<<"Detecting features: "<<argv[i]<<std::endl;
		im3f = im3;
		cvScale(&im3f, &im3f, 1/255.);
		mask.reshape(im3);
		detector.detect(imf);	
		/*for(std::list<detector::feature_t>::const_iterator it=detector.begin(); it!=detector.end(); ++it){
			continue;
					Keypoint k;
				if(!filterEdge(dog.getImage(ar[i][j].level(), ar[i][j].subLevel()), ar[i][j].x(), ar[i][j].y(), -20.0)){
					std::cerr<<"Edge failed"<<std::endl;
					continue;
				}
				if(!subPixelLocalize(dog, ar[i][j], k)){
					std::cerr<<"localization failed\n"<<std::endl;
					continue;
				}
				if(!filterContrast(k)){
					std::cerr<<"contrast failed"<<std::endl;
					continue;
				}
				
				std::cout<<"["<<it->flX()<<","<<it->flY()<<std::endl;
				cvCircle(&im3, cvPoint(it->imX(), it->imY()), it->getRadius(), CV_RGB(255,0,0));
		}*/
		//im3.show("p");
		//cvWaitKey(1);

		//detector.filter();
		for(DoGFeatureDetector::const_iterator it=detector.begin(); it!=detector.end(); ++it){
			// draw filled circle on mask
			cvSet(&mask, cvScalarAll(192.0));
			cvCircle(&mask, cvPoint(it->imX(), it->imY()), (int) it->getRadius(), CV_RGB(255,255,255), -1, CV_AA);
			cvScale(&mask, &mask, 1/192.);
			im3f *= mask;
			cvThreshold( &im3f, &im3f, 255, 255, CV_THRESH_TRUNC );
			// copy contents of circle to dest
			// draw a circle around it
			break;
		}
		cvScale(&im3f, &im3f, 255.);
		im3 = im3f; 
		im3.show("p");
		cvWaitKey(-1);
															
			    
	/*	cvWaitKey(-1);
		for(int j=0; j<detector.getPyramid().getNumScales(); j++){
			for(int k=1; k<detector.getPyramid().getNumSubScales()-1; k++){
				im3 = detector.getPyramid().getImage(j,k).convert(0.5,128);
				for(std::list<detector::feature_t>::const_iterator it=detector.begin(); it!=detector.end(); ++it){
					if(it->level!=j || it->subLevel!=k) continue;
					cvCircle(&im3, cvPoint(it->x, it->y), 2, CV_RGB(0,255,0));
				}
				im3.show("p");
				cvWaitKey(1);
				cvWaitKey(-1);
			}
		}*/
		/*for(int i=0; i<dog.getScaleSpace().length(); i++){
		  dog.getScaleSpace().getImage(i).convert(1,0).show("w");
		  cvWaitKey(-1);
		  }
		  for(int i=0; i<dog.length(); i++){
		  dog.getImage(i).convert(.5,128).show("w");
		  cvWaitKey(-1);
		  }*/
	}
}
