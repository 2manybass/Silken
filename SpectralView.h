#pragma once
#include <gtk/gtk.h>
#include "AnalysisEngine.h"

class SpectralView {
public:
    SpectralView(AnalysisEngine* engine);
    
    GtkWidget* getWidget() { return drawing_area; }
    int getPointerX() const { return pointerX; }
    int getPointerY() const { return pointerY; }
    int getSelectX() const { return selectX; }
    int getSelectY() const { return selectY; }

private:
	static gboolean move(GtkWidget *widget, GdkEventMotion *event, gpointer user_data);
	static gboolean click(GtkWidget *widget, GdkEventMotion *event, gpointer user_data);
    static gboolean on_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data);
    static double windowRatio(double freq, double minFreq, double maxFreq, bool linear);
    
    int pointerX=-1;
    int pointerY=-1;
    int selectX=-1;
    int selectY=-1;
    GtkWidget* drawing_area;
    AnalysisEngine* engine;
};