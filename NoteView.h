#pragma once
#include <gtk/gtk.h>
#include "AnalysisEngine.h"

class NoteView {
public:
    NoteView(AnalysisEngine* engine);
    GtkWidget* getWidget() { return drawing_area; }

private:
    static gboolean on_draw_event(GtkWidget *widget, cairo_t *cr, gpointer user_data);
	static gboolean drag(GtkWidget *widget, GdkEventMotion *event, gpointer user_data);

    GtkWidget* drawing_area;
    AnalysisEngine* engine;
    int minPitch;
    int pitchRange;
    int width;
    int height;
};