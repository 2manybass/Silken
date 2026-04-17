#ifndef ANALYSIS_ENGINE_H
#define ANALYSIS_ENGINE_H

#include <vector>
#include <string>
#include "Params.h"

extern const double RMS_THRESH;
extern const double CORR_THRESH;
const int NUM_HARMS=64;
const int DEFAULT_SAMPLE_RATE = 48000;

// Structure shared across the project
struct SpectrogramFrame {
    std::vector<double> peakFreqs;
    std::vector<double> peakMags;
    std::vector<double> harmData;
    double fundamental;
    double correlation;
    double tunedFrequency;
    double time;
    double rms; // rms stored in a double value representing amplitude
    double localAverage;
    double averageCorrelation;
};

class AnalysisEngine {
public:
    AnalysisEngine(const std::string& filename);
    bool processAudio(); // Returns true if successful
    bool smooth();
    bool smoothPitch();
    bool tune();
    bool untune(double amount);
    
    // Getters for the UI to consume
    const std::vector<SpectrogramFrame>& getFrames() const { return frames; }
    int getSampleRate() const { return sampleRate; }
    int getFrameSize() const { return frameSize; }
    int getFrameCount() const { return frames.size(); }
    int getBinCount() const { return binCount; }
    double getClipDuration() const { return (double)frames[frames.size() - 1].time + (double)frameSize/(double)sampleRate; }
    double getFund(int frame) const { return frames[frame].fundamental; }
    double getCorrelation(int frame) const { return frames[frame].correlation; }
    double getTime(int frame) const { return frames[frame].time; }
    double getTunedFrequency(int frame) const { return frames[frame].tunedFrequency; }
    std::vector<double> getPeakFreqs(int frame) const { return frames[frame].peakFreqs; }
    std::vector<double> getPeakMags(int frame) const { return frames[frame].peakMags; }
    double getMinimumPitch(double thresh);
    double getMaximumPitch(double thresh);
    double freqToMIDI(double freq);
    double MIDIToFreq(double midi);
    bool exportUnpitchedAudio(const std::string& outf);
    bool exportPitchedAudio(const std::string& outf);
    bool exportTunedAudio(const std::string& outf);
    bool calculateHarmonicData();
    bool tuneFrame(int f, double pitch);
    double getLocalAverage(int frame) const { return frames[frame].localAverage; }
	double getAverageCorrelation(int frame) const { return frames[frame].averageCorrelation; }
	std::vector<double> getHarmonicData(int frame) const { return frames[frame].harmData; }
	double getMaxRMS() const {return maxRMS;}
	double getRMS(int frame) const { return frames[frame].rms; }

private:
    std::string filename;
    std::vector<SpectrogramFrame> frames;
    int sampleRate;
    int frameSize;
    int binCount;
    double maxRMS;
};

#endif