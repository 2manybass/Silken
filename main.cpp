#include <gtk/gtk.h>
#include <sndfile.h>
#include <fftw3.h>
#include <vector>
#include <iostream>
#include <cmath>
#include <fstream>

#include "AnalysisEngine.h"
#include "SpectralView.h"
#include "NoteView.h"

// Structure to hold our application data
struct AppData {
	AnalysisEngine* engine;
	SpectralView* sv;
	NoteView* nv;
};

void quickTune(GtkWidget *widget, gpointer data){
	AppData* app = static_cast<AppData*>(data);
	app->engine->tune();
	app->engine->calculateHarmonicData();
	gtk_widget_queue_draw(app->nv->getWidget());
	gtk_widget_queue_draw(app->sv->getWidget());
}
void smoothPitch(GtkWidget *widget, gpointer data){
	AppData* app = static_cast<AppData*>(data);
	app->engine->smoothPitch();
	app->engine->calculateHarmonicData();
	gtk_widget_queue_draw(app->nv->getWidget());
	gtk_widget_queue_draw(app->sv->getWidget());
}
void exportAll(GtkWidget *widget, gpointer data){
	AppData* app = static_cast<AppData*>(data);
	app->engine->calculateHarmonicData();
    app->engine->exportPitchedAudio("outputfile-pitched.wav");
    app->engine->exportUnpitchedAudio("output-unpitched.wav");
}
void exportTuned(GtkWidget *widget, gpointer data){
	AppData* app = static_cast<AppData*>(data);
	app->engine->calculateHarmonicData();
    app->engine->exportTunedAudio("outputfile-tuned.wav");
}
void untune(GtkWidget *widget, gpointer data){
	AppData* app = static_cast<AppData*>(data);
	app->engine->untune(0.5);
	gtk_widget_queue_draw(app->nv->getWidget());
	app->engine->calculateHarmonicData();
}

// GTK Draw Callback: This is where the "Spectrogram" is rendered

int main(int argc, char *argv[]) {

    if (argc < 2) {
        g_print("Usage: %s <audio_file>\n", argv[0]);
        return 1;
    }
    
    AnalysisEngine engine(argv[1]);
    engine.processAudio();
    engine.smooth();
	engine.calculateHarmonicData();
	std::cout << "Clip duration: " << engine.getClipDuration() << std::endl;
	
    gtk_init(&argc, &argv);

    //AppData data;
  	//data.engine = &engine;
  	
  	SpectralView view(&engine);
  	NoteView noteView(&engine);
  	
  	AppData app;
  	app.engine = &engine;
  	app.sv = &view;
  	app.nv = &noteView;

    GtkWidget *window = gtk_window_new(GTK_WINDOW_TOPLEVEL);
	GtkWidget *vbox;
	
	GtkWidget *menubar;
	GtkWidget *fileMenu;
	GtkWidget *optionsMenu;
	GtkWidget *fileMI;
	GtkWidget *audioMI;
	GtkWidget *tuneMI;
	GtkWidget *quickTuneMI;
	GtkWidget *smartTuneMI;
	GtkWidget *smoothMI;
	GtkWidget *audioMenu;
	GtkWidget *tuneMenu;
	GtkWidget *exitMI;
	GtkWidget *optionsMI;
	GtkWidget *exportMI;
	GtkWidget *exportMenu;
	GtkWidget *exportUnpitchedMI;
	GtkWidget *exportPitchedMI;
	GtkWidget *exportAllMI;
	GtkWidget *exportTunedMI;
	GtkWidget *spectSettingsMI;
	GtkWidget *untuneMI;

    g_signal_connect(window, "destroy", G_CALLBACK(gtk_main_quit), NULL);

    gtk_window_set_default_size(GTK_WINDOW(window), 800, 400);
	vbox=gtk_box_new(GTK_ORIENTATION_VERTICAL, 0);
	gtk_container_add(GTK_CONTAINER(window), vbox);
	
	menubar=gtk_menu_bar_new();
	fileMenu=gtk_menu_new();
	audioMenu=gtk_menu_new();
	optionsMenu=gtk_menu_new();
	exportMenu=gtk_menu_new();
	tuneMenu=gtk_menu_new();
	fileMI=gtk_menu_item_new_with_label("File");
	audioMI=gtk_menu_item_new_with_label("Audio");
	tuneMI=gtk_menu_item_new_with_label("Tune");
	quickTuneMI=gtk_menu_item_new_with_label("Quick Tune");
	smartTuneMI=gtk_menu_item_new_with_label("Smart Tune");
	smoothMI=gtk_menu_item_new_with_label("Smooth pitch");
	untuneMI=gtk_menu_item_new_with_label("Un-tune");
	exitMI=gtk_menu_item_new_with_label("Exit");
	exportMI=gtk_menu_item_new_with_label("Export");
	exportUnpitchedMI=gtk_menu_item_new_with_label("Export unpitched audio...");
	exportPitchedMI=gtk_menu_item_new_with_label("Export pitched audio...");
	exportAllMI=gtk_menu_item_new_with_label("Export audio...");
	exportTunedMI=gtk_menu_item_new_with_label("Export tuned audio...");
	optionsMI=gtk_menu_item_new_with_label("Options");
	spectSettingsMI=gtk_menu_item_new_with_label("Spectrogram settings");
	
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(fileMI), fileMenu);
	gtk_menu_shell_append(GTK_MENU_SHELL(menubar), fileMI);

	gtk_menu_shell_append(GTK_MENU_SHELL(fileMenu), exportMI);
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(exportMI), exportMenu);
	gtk_menu_shell_append(GTK_MENU_SHELL(exportMenu), exportUnpitchedMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(exportMenu), exportPitchedMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(exportMenu), exportAllMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(exportMenu), exportTunedMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(fileMenu), exitMI);
	
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(audioMI), audioMenu);
	gtk_menu_shell_append(GTK_MENU_SHELL(audioMenu), tuneMI);
	gtk_menu_item_set_submenu(GTK_MENU_ITEM(tuneMI), tuneMenu);
	gtk_menu_shell_append(GTK_MENU_SHELL(tuneMenu), quickTuneMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(tuneMenu), smartTuneMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(tuneMenu), untuneMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(tuneMenu), smoothMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(menubar), audioMI);

	gtk_menu_item_set_submenu(GTK_MENU_ITEM(optionsMI), optionsMenu);
	gtk_menu_shell_append(GTK_MENU_SHELL(optionsMenu), spectSettingsMI);
	gtk_menu_shell_append(GTK_MENU_SHELL(menubar), optionsMI);

	//gtk_toolbar_insert(toolbar, playButton);

	gtk_box_pack_start(GTK_BOX(vbox), menubar, FALSE, FALSE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), view.getWidget(), TRUE, TRUE, 0);
	gtk_box_pack_start(GTK_BOX(vbox), noteView.getWidget(), TRUE, TRUE, 0);
	
	g_signal_connect(G_OBJECT(window), "destroy", G_CALLBACK(gtk_main_quit), NULL);
	
	g_signal_connect(G_OBJECT(exitMI), "activate", G_CALLBACK(gtk_main_quit), NULL);
	g_signal_connect(G_OBJECT(quickTuneMI), "activate", G_CALLBACK(quickTune), &app);
	g_signal_connect(G_OBJECT(smoothMI), "activate", G_CALLBACK(smoothPitch), &app);
	g_signal_connect(G_OBJECT(exportAllMI), "activate", G_CALLBACK(exportAll), &app);
	g_signal_connect(G_OBJECT(exportTunedMI), "activate", G_CALLBACK(exportTuned), &app);
	g_signal_connect(G_OBJECT(untuneMI), "activate", G_CALLBACK(untune), &app);
	
    gtk_widget_show_all(window);
    gtk_main();

    return 0;
}