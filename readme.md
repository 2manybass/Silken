Silken
v0.1
https://github.com/2manybass/Silken

Silken is an open source program for tuning vocals.  It uses a combination of different algorithms to analyze and correct pitch with the goal of giving the user complete control over how the pitch is altered.

Unlike a real-time autotune VST or LADSPA plugin, Silken looks at the audio file as a whole so that it can do statistical analysis on pitch and assist the user in implementing the best possible tuning transformation.

That being said, as of 2026-04-10, Silken is in its absolute infancy.  After many excruciating hours of doing math and programming loops within loops to process audio, I have just now achieved a minimum viable product that can process one audio channel of up to 30 seconds length at a time.  The tuning is kind of rough and glitchy, but IT WORKS.  Also, there's no dialog for loading / saving files (it only accesses a file called input.wav), you can't zoom in on sections of audio, the spectral axes aren't labeled, you can't change parameters... all that is TO DO.

Without further ado, the TO DO LIST:

 - Create file load / save dialogs
 - Allow Silken to launch without specifying an input file
 - Create an "options" window
 - Create a window for export options (algorithm, window size, etc)
 - Tuning algorithms: PSOLA (default), [vocode, separate, and add], [vocode and fade]
 - Allow zooming and transport
 - Create progress bars
 - Allow audio playback in real time
 - Create "formant shift" graph
 - Set a threshold for minimum amplitude to register pitch OR multiply correlation by amplitude
 - Select scales, keys, tuning calibration
 - Protect against seg faults
 
 So yeah, it's a lot.  I'm gonna publish this now because it's technically a minimum viable product, but good God, I got a long way to go.  Gonna take a nice long break for a while first.
 
Anyways, license is GPL 2.0.  I'm too lazy to research the proper way to declare this.  I made it, and it's free.  You can do what you want with it, except restrict access to it or sell it.  Don't be a dick.