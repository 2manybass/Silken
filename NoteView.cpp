#include "NoteView.h"
#include <cmath>
#include <iostream>
#include <vector>
#include <gtkmm.h>
#include "Params.h"

extern const double RMS_THRESH;
extern const double CORR_THRESH;
const int MAX_DRAG_DIST = 4; // maximum distance between mouse and frame to "drag" pitch

NoteView::NoteView(AnalysisEngine* e) : engine(e) {
    drawing_area = gtk_drawing_area_new();
	gtk_widget_add_events(drawing_area, GDK_BUTTON_PRESS_MASK | GDK_POINTER_MOTION_MASK);
    g_signal_connect(G_OBJECT(drawing_area), "draw", G_CALLBACK(on_draw_event), this);
    g_signal_connect(G_OBJECT(drawing_area), "motion-notify-event", G_CALLBACK(drag), this);
}
gboolean NoteView::drag(GtkWidget *widget, GdkEventMotion *event, gpointer user_data){
	auto* view = static_cast<NoteView*>(user_data);
    auto* engine = view->engine;
	
    int offX=50;
    int offY=0;
    
	if(event->state & GDK_BUTTON1_MASK){
		double x=event->x;
		double y=event->y;
		//std::cout << x << ", " << y << std::endl;
		int frameNum;
		double minDist=100000;
		int closest=0;
		std::vector<int> selected;
		double targetTime = (double)(x - offX) / (double)view->width * engine->getClipDuration();
		//std::cout << targetTime << std::endl;
		for(int f=0; f<engine->getFrameCount(); f++){
			double dist = std::abs(engine->getTime(f) - targetTime);
			double pixDist = std::abs((engine->getTime(f) / engine->getClipDuration() * (double)view->width + offX) - x);
			if(dist < minDist){
				closest=f;
				minDist=dist;
			}
			if(pixDist < MAX_DRAG_DIST){
				selected.push_back(f);
			}
		}
		//std::cout << "Frame index: " << closest << std::endl;
		double mouseRatio = 1 - (double)(y - offY) / (double)view->height;
		//std::cout << "Y: " << y << ", height: " << view->height << std::endl;
		double targetPitch = view->minPitch + mouseRatio * (double)(view->pitchRange) + 0.5;
		//std::cout << "Target pitch: " << targetPitch << std::endl;
		//std::cout << "Mouse ratio: " << mouseRatio << std::endl;
		//std::cout << "Min pitch: " << view->minPitch << std::endl;
		//std::cout << "Pitch at frame " << closest << ": " << engine->getTunedFrequency(closest) << std::endl;
//		bool att=engine->tuneFrame(closest, engine->MIDIToFreq(targetPitch));
		for(int i=0; i<selected.size(); i++){
			engine->tuneFrame(selected[i], engine->MIDIToFreq(targetPitch));
		}
		gtk_widget_queue_draw(widget);
	}
	return true;
}
gboolean NoteView::on_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data) {
    
    //std::cout << "Drewn." << std::endl;
    auto* view = static_cast<NoteView*>(user_data);
    auto* engine = view->engine;
    
    int offX=50;
    int offY=0;
    view->width = gtk_widget_get_allocated_width(widget)-offX;
   	view->height = gtk_widget_get_allocated_height(widget)-offY;
    //std::cout << "Height: " << view->height << std::endl;
    view->minPitch = engine->freqToMIDI(engine->getMinimumPitch(CORR_THRESH)-2);
    view->pitchRange=engine->freqToMIDI(engine->getMaximumPitch(CORR_THRESH)+2) - view->minPitch;
    //view->minPitch = 32;
    //view->pitchRange = 24;
    std::vector<bool> blacks = {false, true, false, true, false, false, true, false, true, false, true, false}; // i'm canceled
    
    cairo_set_source_rgb(cr, 0.8, 0.8, 0.8);
	cairo_rectangle(cr, 0, 0, view->width+offX, view->height+offY);
    cairo_fill(cr);
    cairo_set_source_rgb(cr, 1, 1, 1);
    cairo_rectangle(cr, 0, 0, offX, view->height+offY);
    cairo_fill(cr);
 
 	for(int p=view->minPitch; p<=view->minPitch+view->pitchRange; p++){
 		int y = view->height - (p - view->minPitch) * view->height / view->pitchRange;
 		bool black = blacks[p % 12];
 		bool nextBlack = blacks[(p+1) % 12];
 		if(black){
 			cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
 		}else{
 			cairo_set_source_rgb(cr, 0.8, 0.8, 0.8);
 		}
 		cairo_rectangle(cr, offX, y+offY, view->width, view->height/view->pitchRange);
 		cairo_fill(cr);
 		if(black){
 			cairo_set_source_rgb(cr, 0, 0, 0);
 			cairo_rectangle(cr, 0, y+offY, offX*2/3, view->height/view->pitchRange);
 			cairo_fill(cr);
			cairo_set_source_rgb(cr, 0, 0, 0);
	 		cairo_rectangle(cr, 0, y+offY + view->height/view->pitchRange / 2, offX, 1);
	 		cairo_fill(cr);
 		}
 		if(!black && !nextBlack){ // if stuck between 2 whities
	 		cairo_set_source_rgb(cr, 0, 0, 0);
	 		cairo_rectangle(cr, 0, y+offY, offX, 1);
	 		cairo_fill(cr); 		
 		}
 	}
 	cairo_set_source_rgb(cr, 0, 0, 0);
 	cairo_rectangle(cr, 0, 0, view->width+offX, 1); // top border
	cairo_fill(cr);
 	cairo_set_source_rgb(cr, 0, 0, 0);
 	cairo_rectangle(cr, 0, view->height+offY-1, view->width, 1); // bottom border
	cairo_fill(cr);
 	cairo_set_source_rgb(cr, 0, 0, 0);
 	cairo_rectangle(cr, offX, 0, 1, view->height+offY); // piano border
	cairo_fill(cr);
 	
	for(int f=0; f < engine->getFrameCount(); f++){
    	double fund = engine -> getFund(f);
    	double fund2 = engine -> getTunedFrequency(f);
    	int x = f * view->width / engine -> getFrameCount();
    	double midiPitch=engine->freqToMIDI(fund);
    	double midiPitch2=engine->freqToMIDI(fund2);
    	double ratio = (midiPitch - view->minPitch - 0.5) / (double)view->pitchRange;
    	double ratio2 = (midiPitch2 - view->minPitch - 0.5) / (double)view->pitchRange;
    	int y = view->height + offY - ratio*view->height - 2;
    	int y2 = view->height + offY - ratio2*view->height - 1;
    	double maxRMS = engine -> getMaxRMS();
		double corr = engine -> getCorrelation(f);
		double rms = engine -> getRMS(f);
		if(corr > CORR_THRESH && rms / maxRMS > RMS_THRESH){
			cairo_set_source_rgb(cr, 0.6, 0.6, 0.6);
			cairo_rectangle(cr, x + offX, y, 4, 4);
			cairo_fill(cr);
			cairo_set_source_rgb(cr, 0, 0, 0);
			cairo_rectangle(cr, x + offX, y2, 4, 2);
			cairo_fill(cr);
    	}
    }
    
    return FALSE;
}