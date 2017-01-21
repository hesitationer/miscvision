#include <stdlib.h>
#include <iostream>
#include <vector>
#include "sift_descriptor.hh"
#include "match_points.hh"

int main(){
	std::vector<SIFTDescriptor> array1(50);
	std::vector<SIFTDescriptor> array2(50);
	std::vector<int> matches(50);
	for(int i=0; i<50; i++){
		for(int j=0; j<array1[i].size(); j++){
			array1[i][j] = drand48()*50;
			array2[i][j] = drand48()*50;
		}
	}
	
	match_points<SIFTDescriptor, float, 128>::match(array1, array2, matches, 1023.0f);
	for(int i=0; i<50; i++){
		if(matches[i]!=-1){
			std::cout<<array2[i]<<" matches index "<<matches[i]<<" value "<<array1[matches[i]]<<std::endl;
		}
		else{
			std::cout<<array2[i]<<" has no match"<<std::endl;
		}
	}
		
}
