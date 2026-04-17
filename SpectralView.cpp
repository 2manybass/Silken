#include "SpectralView.h"
#include "Params.h"
#include <cmath>
#include <iostream>

const int FREQ_DELTA = 100;
const int offX=50;
const int offY=30;
extern const double CORR_THRESH;
extern const double RMS_THRESH;

gboolean SpectralView::move(GtkWidget *widget, GdkEventMotion *event, gpointer user_data){
	auto* view = static_cast<SpectralView*>(user_data);
    auto* engine = view->engine;
	double x=event->x;
	double y=event->y;
	if(x > offX){
		view->pointerX=x;
	}
	if(y > offY){
		view->pointerY=y;
	}
	gtk_widget_queue_draw(widget);
	return true;
}
gboolean SpectralView::click(GtkWidget *widget, GdkEventMotion *event, gpointer user_data){
	auto* view = static_cast<SpectralView*>(user_data);
    auto* engine = view->engine;
	double x=event->x;
	double y=event->y;
	if(x > offX){
		view->selectX=x;
	}
	if(y > offY){
		view->selectY=y;
	}
	gtk_widget_queue_draw(widget);
	std::cout << "CLICKED" << std::endl;
	return true;
}
SpectralView::SpectralView(AnalysisEngine* e) : engine(e) {
    drawing_area = gtk_drawing_area_new();
	gtk_widget_add_events(drawing_area, GDK_POINTER_MOTION_MASK);
	gtk_widget_add_events(drawing_area, GDK_BUTTON_PRESS_MASK);
    g_signal_connect(G_OBJECT(drawing_area), "draw", G_CALLBACK(on_draw_event), this);
    g_signal_connect(G_OBJECT(drawing_area), "motion-notify-event", G_CALLBACK(move), this);
    g_signal_connect(G_OBJECT(drawing_area), "button-press-event", G_CALLBACK(click), this);
}
double SpectralView::windowRatio(double freq, double minFreq, double maxFreq, bool linear){
	double ratio;
	double max_ratio = std::log2(maxFreq / minFreq);
	if(linear){
		ratio = (freq - minFreq) / (double) maxFreq;
	}else{
		ratio = std::log2(freq / minFreq) / max_ratio;
	}
	return ratio;
}
gboolean SpectralView::on_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data) {
    auto* view = static_cast<SpectralView*>(user_data);
    auto* engine = view->engine;
    
    int width = gtk_widget_get_allocated_width(widget)-offX;
    int height = gtk_widget_get_allocated_height(widget)-offY;
    
    cairo_set_source_rgb(cr, 0.8, 0.8, 0.8);
	cairo_rectangle(cr, 0, 0, width+offX, height+offY);
    cairo_fill(cr);
    
    for(int sec=0; sec < engine->getClipDuration(); sec++){
    	int x = sec * width / engine->getClipDuration() + offX;
    	cairo_set_source_rgb(cr, 0, 0, 0);
    	cairo_rectangle(cr, x, 0, 1, offY);
    }
    
    double min_freq = (double)engine->getSampleRate() / (double)engine->getFrameSize() / 2.0;
    double max_freq = (double)engine->getSampleRate() / 2.0;
	double max_log = std::log2(max_freq/min_freq); // find log2 of highest frequency bin #
    bool linear=false;
    int lineNum=0;
    if(view -> mode == 0){ // frequency analysis mode

		for(int freq=0; freq<max_freq; freq+=FREQ_DELTA){ // TODO: replace this with a helper function
			int rectY;
			double ratio;
			double max_ratio = std::log2(max_freq / min_freq);
			if(linear){
				ratio = ((double)freq - min_freq) / (double) max_freq;
			}else{
				ratio = std::log2(freq / min_freq) / max_ratio;
			}
			double ratioDiff;
			if(linear){
				ratioDiff = ((double)(freq + FREQ_DELTA) - min_freq) / (double) max_freq - ratio;
			}else{
				ratioDiff = std::log2((freq + FREQ_DELTA) / min_freq) / max_ratio - ratio;
			}
			bool draw=true;
			bool restrict=false;
			if(ratioDiff < 0.03){	// lines too close together
				draw = lineNum%10 == 0;
				restrict = true;
			}
			rectY=height - (int)std::round(ratio * (double)height) + offY;
			if(draw){
				cairo_set_source_rgb(cr, 0, 0, 0);
				cairo_rectangle(cr, offX*3/4, rectY, offX/4, 1);			
				cairo_fill(cr);
			}
			lineNum++;
		}
		cairo_set_source_rgb(cr, 0, 0, 0);
		cairo_rectangle(cr, offX, offY, width, height);
		cairo_fill(cr);

		for(int x=0; x<width; x++){
			int frame=x * engine->getFrameCount() / width;
			std::vector<double> peakFreqs = engine -> getPeakFreqs(frame);
			std::vector<double> peakMags = engine -> getPeakMags(frame);
			for(int i=0; i < peakFreqs.size(); i++){
				int rectY;
				double ratio = windowRatio(peakFreqs[i], min_freq, max_freq, linear);
				double val = peakMags[i]/100.0; // lazy way to convert dB to RGB
				rectY=height - (int)std::round(ratio * (double)height);
				if(rectY > 0 && rectY < height){ // if the rect is within bounds
					cairo_set_source_rgb(cr, val, 0, 0);
					cairo_rectangle(cr, x + offX, rectY + offY, 2, 2);
					cairo_fill(cr);
				}
			}
			
			double corr=engine->getCorrelation(frame);
			double rms = engine -> getRMS(frame);
			double ratio = rms / engine->getMaxRMS();
			if(ratio > RMS_THRESH){
				cairo_set_source_rgb(cr, 0, corr, corr); // draw the fundamental in CYAN
				//cairo_set_source_rgb(cr, 0, 1, 1);
				double ratio=windowRatio(engine->getFund(frame), min_freq, max_freq, linear);
				int rectY=height - (int)std::round(ratio * (double)height);
				if(rectY > 0 && rectY < height){ // if the rect is within bounds
					cairo_rectangle(cr, x + offX, rectY + offY, 2, 2);
					cairo_fill(cr);
				}
			}
		}
		cairo_set_source_rgb(cr, 1, 1, 0);
		cairo_rectangle(cr, view->getPointerX(), offY, 1, height);
		cairo_fill(cr);
		cairo_rectangle(cr, offX, view->getPointerY(), width, 1);
		cairo_fill(cr);

		cairo_set_source_rgb(cr, 1, 1, 1);
		cairo_rectangle(cr, view->getSelectX(), offY, 1, height);
		cairo_fill(cr);
		cairo_rectangle(cr, offX, view->getSelectY(), width, 1);
		cairo_fill(cr);
	
	}else if(view -> mode == 1){
		double xscale = (double)width / (double)engine->getFrameCount();
		double yscale = (double)height / (double)NUM_HARMS;
		for(int f=0; f<engine->getFrameCount(); f++){
		std::vector<double> harmData = engine->getHarmonicData(f);
			for(int h=0; h < NUM_HARMS; h++){
				int y = (int)std::round(height - yscale * (double) h) + offY;
				double col = (harmData[h] + 50) / 50;
				if(col < 0){
					col = 0;
				}
				if(engine->getCorrelation(f) > CORR_THRESH){
					cairo_set_source_rgb(cr, col, 0, 0);
				}else{
					cairo_set_source_rgb(cr, col/3, col/3, col/3);				
				}
				cairo_rectangle(cr, (int)std::round((double)f * xscale) + offX, y, std::ceil(xscale), std::ceil(yscale));
				cairo_fill(cr);
			}
		}
	}
	//std::cout << "DREW SPECTAL STUFF" << std::endl;
    return FALSE;
}