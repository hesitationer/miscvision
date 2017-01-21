include ../make/Makefile.base

all: libcvext.a 

# the opencv extension library
libcvext.a: cvext.o cvransac.o cvprosac.o cvmat.o cvmat5.o
	$(LIBLINK)

harris: harris.o libcvext.a
	$(LINK)

cvpca: cvpca.o libcvext.a
	$(LINK)

corner_detector: corner_detector.o libcvext.a
	$(LINK)

median: test/median.o libcvext.a
	$(LINK)

gramschmidt: test/gramschmidt.o libcvext.a
	$(LINK)

autocrop: autocrop.o libcvext.a
	$(LINK)

stereo_camera: stereo_camera.o
	$(LINK)

cvmatlabio: cvmatlabio.o libcvext.a
	$(LINK)

cvmread: cvmread.o libcvext.a 
	$(LINK)

shape: shape.o libcvext.a 
	$(LINK)

cholesky: cholesky.o libcvext.a 
	$(LINK)

matrix: matrix.o libcvext.a
	$(LINK)

eigen: eigen.o libcvext.a
	$(LINK)

matrixnd: matrixnd.o
	$(LINK)

arraybase: test/arraybase.o
	$(LINK)

rotdft: rotdft.o libcvext.a 
	$(LINK)

ridge_detector_extract: ridge_detector_extract.o libcvext.a 
	$(LINK)

ridge_detector: ridge_detector.o
	$(LINK)

circle_detector: circle_detector.o libcvext.a 
	$(LINK)

seq: seq.o
	$(LINK)
circle: circle.o libcvext.a 
	$(LINK)
blobdetect: test/blobdetect.o
	$(LINK)
rotate: rotate.o libcvext.a 
	$(LINK)

localmax: test/localmax.o libcvext.a
	$(LINK)
geometricblur: geometricblur.o
	$(LINK)
rgbe_sub: rgbe_sub.o
	$(LINK)
rgbe: rgbe.o
	$(LINK)
contour: contour.o
	$(LINK)

bullseye_detect: bullseye_detect.o
	$(LINK)
	
matlab_compat_test: matlab_compat_test.o
	$(LINK)
	
cap: cap.o
	$(LINK)

sp_tree: sp_tree.o
	$(LINK)

render_feature: render_feature.o cvmat.o
	$(LINK)

find_eyes_features: test/find_eyes_features.o libcvext.a 
	$(LINK)

find_eyes2: test/find_eyes2.o libcvext.a 
	$(LINK)

find_eyes: test/find_eyes.o libcvext.a 
	$(LINK)

color_map: test/color_map.o
	$(LINK)

face_detector: test/face_detector.o
	$(LINK)
match_points: match_points.o
	$(LINK)

scale_space: test/scale_space.o libcvext.a 
	$(LINK)

array_test: test/array_test.o
	$(LINK)

operator_test: test/operator_test.o libcvext.a 
	$(LINK)

image: test/image.o
	$(LINK)

mat: test/mat.o
	$(LINK)

dog_scale_space: test/dog_scale_space.o cvsift.o libcvext.a
	$(LINK)

sift_feature_track: test/sift_feature_track.o libcvext.a
	$(LINK)
point_correspond_lk: test/point_correspond_lk.o libcvext.a
	$(LINK)
point_correspond_sift: test/point_correspond_sift.o libcvext.a
	$(LINK)
sift_feature_detector: test/sift_feature_detector.o libcvext.a
	$(LINK)

dog_feature_detector2: test/dog_feature_detector2.o libcvext.a
	$(LINK)

dog_feature_detector: test/dog_feature_detector.o cvsift.o libcvext.a
	$(LINK)

time_dog: time_dog.o
	$(LINK)

int: int.o
	$(LINK)

smooth: test/smooth.o
	$(LINK)
	
scale_view: test/scale_view.o libcvext.a
	$(LINK)
	
piped_viewer: piped_viewer.o
	$(LINK)

mmap_writer: mmap_writer.o
	$(LINK)
mmap_viewer: im_mmap.o
	$(LINK)

imview1: test/imview1.o
	$(LINK)
imview: test/imview.o
	$(LINK)

sift_detect_and_describe.dll: mex/sift_detect_and_describe.o libcvext.a
	$(MEX_DLL)

sift_detect.dll: mex/sift_detect.o libcvext.a
	$(MEX_DLL)

sift_correspond.dll: mex/sift_correspond.o libcvext.a
	$(MEX_DLL)

mex: sift_detect.dll sift_detect_and_describe.dll sift_correspond.dll
	cp $^ "$(MATLAB_DIR)/work"
clean:
	rm -f $(CLEAN_FILES) ; cd test && rm -f $(CLEAN_FILES) 


