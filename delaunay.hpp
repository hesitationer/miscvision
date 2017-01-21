#ifndef __DELAUNAY_HPP
#define __DELAUNAY_HPP

#include <cxcore.h>
#include <cvtypes.h>
#include <cv.h>

// calculate delaunay subdivision of a set of points
class Delaunay {
protected:
	CvSubdiv2D * m_subdiv;
	CvMemStorage * m_storage;

public:	
	Delaunay(): m_subdiv(NULL), m_storage(NULL) {}

	void init(const CvRect & rect ){
		if(!m_storage){
			m_storage = cvCreateMemStorage(0);
		}
		else{
			cvClearMemStorage(m_storage);
		}
		m_subdiv = cvCreateSubdiv2D( CV_SEQ_KIND_SUBDIV2D, sizeof(*m_subdiv),
				sizeof(CvSubdiv2DPoint),
				sizeof(CvQuadEdge2D),
				m_storage );
		cvInitSubdivDelaunay2D( m_subdiv, rect );
	}
	void insert( CvMat * points ){
		CvMat vec2;

		// if a vector, reshape into a 2 column matrix
		if(points->cols==1){
			cvReshape( points, &vec2, 1, points->rows/2*CV_MAT_CN(points->type) );
			points = &vec2;
		}

		for(int i=0; i<points->rows; i++){
			CvPoint2D32f point = * (CvPoint2D32f *) cvPtr2D( points, i, 0 );
	        cvSubdivDelaunay2DInsert( m_subdiv, point );
		}
	}
	void insert( const CvPoint2D32f & point ){
        cvSubdivDelaunay2DInsert( m_subdiv, point );
	}
	void calc(){
		cvCalcSubdivVoronoi2D( m_subdiv );
	}
	void clear(){

	}
	CvSubdiv2DPoint * operator [] (int i){
		return (CvSubdiv2DPoint *) cvGetSetElem( (CvSet *) m_subdiv, i+3);
	}
	int nFaces() {
		return m_subdiv->quad_edges*2;
	}
	void draw(IplImage * im, CvScalar color, bool delaunay=true){
		if(delaunay){
			this->draw_delaunay( im, m_subdiv, color );
		}
		else{
			this->draw_voronoi( im, m_subdiv, color );
		}
	}

	void draw_subdiv_edge( IplImage* img, CvSubdiv2DEdge edge, CvScalar color )
	{
		CvSubdiv2DPoint* org_pt;
		CvSubdiv2DPoint* dst_pt;
		CvPoint2D32f org;
		CvPoint2D32f dst;
		CvPoint iorg, idst;

		org_pt = cvSubdiv2DEdgeOrg(edge);
		dst_pt = cvSubdiv2DEdgeDst(edge);

		if( org_pt && dst_pt )
		{
			org = org_pt->pt;
			dst = dst_pt->pt;

			iorg = cvPoint( cvRound( org.x ), cvRound( org.y ));
			idst = cvPoint( cvRound( dst.x ), cvRound( dst.y ));

			cvLine( img, iorg, idst, color, 1, CV_AA, 0 );
		}
	}


	void draw_delaunay( IplImage* img, CvSubdiv2D* subdiv, CvScalar color )
	{
		CvSeqReader  reader;
		int i, total = subdiv->edges->total;
		int elem_size = subdiv->edges->elem_size;

		cvStartReadSeq( (CvSeq*)(subdiv->edges), &reader, 0 );

		for( i = 0; i < total; i++ )
		{
			CvQuadEdge2D* edge = (CvQuadEdge2D*)(reader.ptr);

			if( CV_IS_SET_ELEM( edge ))
			{
				draw_subdiv_edge( img, (CvSubdiv2DEdge)edge, color );
			}

			CV_NEXT_SEQ_ELEM( elem_size, reader );
		}
	}

	void draw_voronoi( IplImage* img, CvSubdiv2D* subdiv, CvScalar color )
	{
		CvSeqReader  reader;
		int i, total = subdiv->edges->total;
		int elem_size = subdiv->edges->elem_size;

		cvStartReadSeq( (CvSeq*)(subdiv->edges), &reader, 0 );

		for( i = 0; i < total; i++ )
		{
			CvQuadEdge2D* edge = (CvQuadEdge2D*)(reader.ptr);

			if( CV_IS_SET_ELEM( edge ))
			{
				draw_subdiv_edge( img, (CvSubdiv2DEdge)edge + 1, color );
			}

			CV_NEXT_SEQ_ELEM( elem_size, reader );
		}
	}
	const CvPoint2D32f & getPoint(int i){
		return ((CvSubdiv2DPoint *) cvGetSetElem( (CvSet *) m_subdiv, i ))->pt;
	}

	// 
	const CvPoint2D32f & getVoronoiPoint(int i){
		CvQuadEdge2D * quadedge = (CvQuadEdge2D *) cvGetSeqElem( (CvSeq *) m_subdiv->edges, i+3 );
		if( quadedge && CV_IS_SET_ELEM( quadedge ) &&  quadedge->pt[1] ){
			return quadedge->pt[1]->pt;
		}
		return cvPoint2D32f(0,0);
	}
	static int CMP_SUBDIV_PT(const void* _a, const void* _b, void* userdata){
		CvPoint2D32f a = ((CvSubdiv2DPoint *)_a)->pt;
		CvPoint2D32f b = ((CvSubdiv2DPoint *)_b)->pt;
		return a.x < b.x ? -1 : a.x > b.x ? 1 : a.y < b.y ? -1 : a.y > b.y ? 1 : 0;
	}
	
	void getTriangle(int i, int idx[3]){
		CvQuadEdge2D* edge = CV_GET_SEQ_ELEM( CvQuadEdge2D, m_subdiv->edges, i/2);
		if( CV_IS_SET_ELEM( edge ))
        {
            CvSubdiv2DEdge e = cvSubdiv2DRotateEdge( (CvSubdiv2DEdge)edge, (i%2)*2 );
            CvSubdiv2DEdge t = e;
            int j, count = 3;

            // gather points
            for( j = 0; j < count; j++ )
            {
				CvSubdiv2DPoint* pt = cvSubdiv2DEdgeOrg( t );
				
				// find point in set
				if(!pt){
					idx[j] = -1;
				}

				//|| cvSeqSearch( (CvSeq *)m_subdiv, pt, Delaunay::CMP_SUBDIV_PT, 0, idx+j)==NULL ){
				idx[j] = cvSeqElemIdx( (CvSeq *)m_subdiv, pt);
				

				t = cvSubdiv2DGetEdge( t, CV_NEXT_AROUND_LEFT );
			}
		}
		else{
			idx[0] = idx[1] = idx[2] = -1;
		}
	}
	void getTriangle(int i, CvPoint2D32f tri[3]){
		CvQuadEdge2D* edge = CV_GET_SEQ_ELEM( CvQuadEdge2D, m_subdiv->edges, i/2);
		if( CV_IS_SET_ELEM( edge ))
		{
			CvSubdiv2DEdge e = cvSubdiv2DRotateEdge( (CvSubdiv2DEdge)edge, (i%2)*2 );
			CvSubdiv2DEdge t = e;
			int j, count = 0;

			// count number of edges in facet
			do
			{
				count++;
				t = cvSubdiv2DGetEdge( t, CV_NEXT_AROUND_LEFT );
			} while (t != e );

			assert(count==3);

			// gather points
			t = e;
			for( j = 0; j < count; j++ )
			{
				CvSubdiv2DPoint* pt = cvSubdiv2DEdgeOrg( t );
				if( !pt ) break;
				tri[j] = pt->pt;
				t = cvSubdiv2DGetEdge( t, CV_NEXT_AROUND_LEFT );
			}
		}
		else{
			tri[0] = cvPoint2D32f(0,0);
			tri[1] = cvPoint2D32f(0,0);
			tri[2] = cvPoint2D32f(0,0);
		}
	}
};
#endif // __DELAUNAY_HPP
