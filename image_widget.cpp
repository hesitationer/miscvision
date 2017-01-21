#include <cxcore.h>
#include <cv.h>
#include <highgui.h>
#include <gtk/gtk.h>

// Public API
typedef struct _CvImageWidget        CvImageWidget;
typedef struct _CvImageWidgetClass   CvImageWidgetClass;

struct _CvImageWidget {
	GtkWidget widget;
	CvMat * original_image;
	CvMat * scaled_image;
	int flags;
};

struct _CvImageWidgetClass
{
  GtkWidgetClass parent_class;
};


/** Allocate new image viewer widget */
GtkWidget*     CvImageWidgetNew      (int flags);

/** Set the image to display in the widget */
void           CvImageWidgetSetImage(CvImageWidget * widget, CvArr *arr);

// standard GTK object macros
#define CV_IMAGE_WIDGET(obj)          GTK_CHECK_CAST (obj, CvImageWidget_get_type (), CvImageWidget)
#define CV_IMAGE_WIDGET_CLASS(klass)  GTK_CHECK_CLASS_CAST (klass, CvImageWidget_get_type (), CvImageWidgetClass)
#define CV_IS_IMAGE_WIDGET(obj)       GTK_CHECK_TYPE (obj, CvImageWidget_get_type ())

// Private API
GtkType        CvImageWidget_get_type (void);

static GtkWidgetClass * parent_class = NULL;

void CvImageWidgetSetImage(CvImageWidget * widget, CvArr *arr){
	CvMat mat;
	int origin=0;
	if( CV_IS_IMAGE_HDR( arr ))
		origin = ((IplImage*)arr)->origin;

	cvGetMat(arr, &mat);
	if(widget->original_image && !CV_ARE_SIZES_EQ(&mat, widget->original_image)){
		cvReleaseMat( &widget->original_image );
	}
	if(!widget->original_image){
		widget->original_image = cvCreateMat( mat.rows, mat.cols, CV_8UC3 );
		gtk_widget_queue_resize( GTK_WIDGET( widget ) );
	}
	cvConvertImage( &mat, widget->original_image,
			                (origin != 0 ? CV_CVTIMG_FLIP : 0) + CV_CVTIMG_SWAP_RB );
	if(widget->scaled_image){
		cvResize( widget->original_image, widget->scaled_image, CV_INTER_AREA );
	}
}

GtkWidget*
CvImageWidgetNew (int flags)
{
  CvImageWidget *image_widget;

  image_widget = CV_IMAGE_WIDGET( gtk_type_new (CvImageWidget_get_type ()) );
  image_widget->original_image = 0;
  image_widget->scaled_image = 0;
  image_widget->flags = flags;

  return GTK_WIDGET (image_widget);
}

static void
CvImageWidget_realize (GtkWidget *widget)
{
  CvImageWidget *image_widget;
  GdkWindowAttr attributes;
  gint attributes_mask;

  g_return_if_fail (widget != NULL);
  g_return_if_fail (CV_IS_IMAGE_WIDGET (widget));

  GTK_WIDGET_SET_FLAGS (widget, GTK_REALIZED);
  image_widget = CV_IMAGE_WIDGET (widget);

  attributes.x = widget->allocation.x;
  attributes.y = widget->allocation.y;
  attributes.width = widget->allocation.width;
  attributes.height = widget->allocation.height;
  attributes.wclass = GDK_INPUT_OUTPUT;
  attributes.window_type = GDK_WINDOW_CHILD;
  attributes.event_mask = gtk_widget_get_events (widget) | 
    GDK_EXPOSURE_MASK | GDK_BUTTON_PRESS_MASK | 
    GDK_BUTTON_RELEASE_MASK | GDK_POINTER_MOTION_MASK |
    GDK_POINTER_MOTION_HINT_MASK;
  attributes.visual = gtk_widget_get_visual (widget);
  attributes.colormap = gtk_widget_get_colormap (widget);

  attributes_mask = GDK_WA_X | GDK_WA_Y | GDK_WA_VISUAL | GDK_WA_COLORMAP;
  widget->window = gdk_window_new (widget->parent->window, &attributes, attributes_mask);

  widget->style = gtk_style_attach (widget->style, widget->window);

  gdk_window_set_user_data (widget->window, widget);

  gtk_style_set_background (widget->style, widget->window, GTK_STATE_ACTIVE);
}

static void 
CvImageWidget_size_request (GtkWidget      *widget,
                       GtkRequisition *requisition)
{
	CvImageWidget * image_widget = CV_IMAGE_WIDGET( widget );
	printf("CvImageWidget_size_request\n");
	if( (image_widget->flags & CV_WINDOW_AUTOSIZE) && image_widget->original_image){
		requisition->width = image_widget->original_image->cols;
		requisition->height = image_widget->original_image->rows;
	}
	else{
		requisition->width = 100;
		requisition->height = 100;
	}
}
CvSize icvCalcImageSize( int im_width, int im_height, int max_width, int max_height ){
    float aspect = (float)im_width/(float)im_height;
    float max_aspect = (float)max_width/(float)max_height;
    if(aspect > max_aspect){
        return cvSize( max_width, cvRound(max_width/aspect) );
    }
    return cvSize( cvRound(max_height*aspect), max_height );
}

static void
CvImageWidget_size_allocate (GtkWidget     *widget,
                        GtkAllocation *allocation)
{
  CvImageWidget *image_widget;
	printf("CvImageWidget_size_allocate\n");

  g_return_if_fail (widget != NULL);
  g_return_if_fail (CV_IS_IMAGE_WIDGET (widget));
  g_return_if_fail (allocation != NULL);

  widget->allocation = *allocation;
  image_widget = CV_IMAGE_WIDGET (widget);
  
  // (re) allocated scaled image
  if( (image_widget->flags & CV_WINDOW_AUTOSIZE)==0 && image_widget->original_image ){
	  CvSize scaled_image_size = icvCalcImageSize( image_widget->original_image->cols, 
			                                       image_widget->original_image->rows,
			                                       allocation->width, allocation->height );
	  if( image_widget->scaled_image && 
			  ( image_widget->scaled_image->cols != scaled_image_size.width ||
				image_widget->scaled_image->rows != scaled_image_size.height ))
	  {
		  cvReleaseMat( &image_widget->scaled_image );
	  }
	  if( !image_widget->scaled_image ){
		  image_widget->scaled_image = cvCreateMat( scaled_image_size.height, scaled_image_size.width, CV_8UC3 );
	  }
	  cvResize( image_widget->original_image, image_widget->scaled_image, CV_INTER_AREA );
  }

  if (GTK_WIDGET_REALIZED (widget))
    {
      image_widget = CV_IMAGE_WIDGET (widget);

	  if( (image_widget->flags & CV_WINDOW_AUTOSIZE) && image_widget->original_image ){
		  widget->allocation.width = image_widget->original_image->cols;
		  widget->allocation.height = image_widget->original_image->rows;
		  gdk_window_move_resize( widget->window, allocation->x, allocation->y, 
				  image_widget->original_image->cols, image_widget->original_image->rows );
	  }
	  else{
		  gdk_window_move_resize (widget->window,
				  allocation->x, allocation->y,
				  allocation->width, allocation->height);
	  }

    }
}

static gboolean
CvImageWidget_expose( GtkWidget      *widget,
                 GdkEventExpose *event )
{
  CvImageWidget *image_widget;
  GdkPoint points[3];
  gdouble s,c;
  gdouble theta;
  gint xc, yc;
  gint tick_length;
  gint i;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (CV_IS_IMAGE_WIDGET (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  if (event->count > 0)
    return FALSE;
  
  image_widget = CV_IMAGE_WIDGET (widget);

  gdk_window_clear_area (widget->window,
                         0, 0,
                         widget->allocation.width,
                         widget->allocation.height);

  if( image_widget->scaled_image ){
	  gdk_draw_rgb_image( widget->window, widget->style->fg_gc[GTK_STATE_NORMAL],
		  0, 0, MIN(image_widget->scaled_image->cols, widget->allocation.width), 
		  MIN(image_widget->scaled_image->rows, widget->allocation.height),
		  GDK_RGB_DITHER_MAX, image_widget->scaled_image->data.ptr, image_widget->scaled_image->step );
  }
  else if( image_widget->original_image ){
	  gdk_draw_rgb_image( widget->window, widget->style->fg_gc[GTK_STATE_NORMAL],
		  0, 0, widget->allocation.width, widget->allocation.height,
		  GDK_RGB_DITHER_MAX, image_widget->original_image->data.ptr, image_widget->original_image->step );
	}
  return FALSE;
}

static gboolean
CvImageWidget_button_press( GtkWidget      *widget,
                       GdkEventButton *event )
{
  CvImageWidget *image_widget;
  gint dx, dy;
  double s, c;
  double d_parallel;
  double d_perpendicular;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (CV_IS_IMAGE_WIDGET (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  image_widget = CV_IMAGE_WIDGET (widget);

  /* Determine if button press was within pointer region - we 
     do this by computing the parallel and perpendicular distance of
     the point where the mouse was pressed from the line passing through
     the pointer */
  
  dx = event->x - widget->allocation.width / 2;
  dy = widget->allocation.height / 2 - event->y;
  
  return FALSE;
}

static gboolean
CvImageWidget_button_release( GtkWidget      *widget,
                         GdkEventButton *event )
{
  CvImageWidget *image_widget;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (CV_IS_IMAGE_WIDGET (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  image_widget = CV_IMAGE_WIDGET (widget);

  return FALSE;
}

static gboolean
CvImageWidget_motion_notify( GtkWidget      *widget,
                        GdkEventMotion *event )
{
  CvImageWidget *image_widget;
  GdkModifierType mods;
  gint x, y, mask;

  g_return_val_if_fail (widget != NULL, FALSE);
  g_return_val_if_fail (CV_IS_IMAGE_WIDGET (widget), FALSE);
  g_return_val_if_fail (event != NULL, FALSE);

  image_widget = CV_IMAGE_WIDGET (widget);

  return FALSE;
}

static void
CvImageWidget_update_mouse (CvImageWidget *image_widget, gint x, gint y)
{
  gint xc, yc;
  gfloat old_value;

  g_return_if_fail (image_widget != NULL);
  g_return_if_fail (CV_IS_IMAGE_WIDGET (image_widget));

  xc = GTK_WIDGET(image_widget)->allocation.width / 2;
  yc = GTK_WIDGET(image_widget)->allocation.height / 2;

}

static void
CvImageWidget_update (CvImageWidget *image_widget)
{
  gfloat new_value;
  
  g_return_if_fail (image_widget != NULL);
  g_return_if_fail (CV_IS_IMAGE_WIDGET (image_widget));

  gtk_widget_draw (GTK_WIDGET(image_widget), NULL);
}



static void
CvImageWidget_destroy (GtkObject *object)
{
  CvImageWidget *image_widget;

  g_return_if_fail (object != NULL);
  g_return_if_fail (CV_IS_IMAGE_WIDGET (object));

  image_widget = CV_IMAGE_WIDGET (object);

  if (GTK_OBJECT_CLASS (parent_class)->destroy)
    (* GTK_OBJECT_CLASS (parent_class)->destroy) (object);
}

static void CvImageWidget_class_init (CvImageWidgetClass * klass)
{
  GtkObjectClass *object_class;
  GtkWidgetClass *widget_class;

  object_class = (GtkObjectClass*) klass;
  widget_class = (GtkWidgetClass*) klass;
  
  parent_class = GTK_WIDGET_CLASS( gtk_type_class (gtk_widget_get_type ()) );

  object_class->destroy = CvImageWidget_destroy;

  widget_class->realize = CvImageWidget_realize;
  widget_class->expose_event = CvImageWidget_expose;
  widget_class->size_request = CvImageWidget_size_request;
  widget_class->size_allocate = CvImageWidget_size_allocate;
  widget_class->button_press_event = CvImageWidget_button_press;
  widget_class->button_release_event = CvImageWidget_button_release;
  widget_class->motion_notify_event = CvImageWidget_motion_notify;
}

static void
CvImageWidget_init (CvImageWidget *image_widget)
{
	image_widget->original_image=0;
	image_widget->scaled_image=0;
	image_widget->flags=0;
}

GtkType CvImageWidget_get_type (void){
  static GtkType image_type = 0;

  if (!image_type)
    {
      static const GtkTypeInfo image_info =
      {
        "CvImageWidget",
        sizeof (CvImageWidget),
        sizeof (CvImageWidgetClass),
        (GtkClassInitFunc) CvImageWidget_class_init,
        (GtkObjectInitFunc) CvImageWidget_init,
        /* reserved_1 */ NULL,
        /* reserved_1 */ NULL,
        (GtkClassInitFunc) NULL
      };

      image_type = gtk_type_unique (GTK_TYPE_WIDGET, &image_info);
    }

  return image_type;
}


int 
main (int argc, char *argv[])
{
  GtkWidget *window;
  GtkWidget *vbox;
  GtkWidget *image_widget;
  
  gtk_init (&argc, &argv);

  window = gtk_window_new (GTK_WINDOW_TOPLEVEL);
  
  gtk_window_set_title (GTK_WINDOW (window), "Aspect Frame");
  
  gtk_signal_connect (GTK_OBJECT (window), "destroy",
                      GTK_SIGNAL_FUNC (gtk_exit), NULL);
  
  gtk_container_border_width (GTK_CONTAINER (window), 0);

  IplImage * image = cvLoadImage( argv[1] );
  vbox = gtk_vbox_new( FALSE, 0 );
  image_widget = CvImageWidgetNew( CV_WINDOW_AUTOSIZE ); 
  //image_widget = CvImageWidgetNew( 0 ); 
  gtk_box_pack_end( GTK_BOX(vbox), image_widget, TRUE, TRUE, 0 );
  gtk_container_add (GTK_CONTAINER(window), vbox);
  gtk_widget_show( vbox );
  gtk_widget_show (image_widget);
  gtk_widget_show (window);
  
  CvImageWidgetSetImage( CV_IMAGE_WIDGET(image_widget), image );

  gtk_main ();
  
  return 0;
}
