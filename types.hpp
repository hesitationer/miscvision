#include<cxtypes.h>
// stupid template tricks to get correct depth flag
template <typename T> struct CvDepthType { static int depth() { assert(0); return 0; }};
template <> struct CvDepthType<double> {    static int depth() { return CV_64F; }};
template <> struct CvDepthType<float> { static int depth() { return CV_32F; }};
template <> struct CvDepthType<int> {   static int depth() { return CV_32S; }};
template <> struct CvDepthType<short> { static int depth() { return CV_16S; }};
template <> struct CvDepthType<unsigned short> {    static int depth() { return CV_16S; }};
template <> struct CvDepthType<char> {  static int depth() { return CV_8S; }};
template <> struct CvDepthType<unsigned char> { static int depth() { return CV_8U; }};

