#include <vector>
#include <cv/sort_by.hpp>

int main(int argc, char ** argv){
	std::vector<int> indices( 25 );
	std::vector<float> weights( 25 );
	for(int i=0; i<25; i++){
		indices[i] = i;
		weights[i] = drand48();
		printf("%d - %f\n", indices[i], weights[i]);
	}
	sort_by( indices.begin(), indices.end(), weights.begin(), weights.end());
	for(int i=0; i<25; i++){
		printf("%d - %f\n", indices[i], weights[i]);
	}
}
