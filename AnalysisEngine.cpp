#include "AnalysisEngine.h"
#include <sndfile.h>
#include <fftw3.h>
#include <cmath>
#include <random>
#include <iostream>

AnalysisEngine::AnalysisEngine(const std::string& f) : filename(f) {}

const double RATIO_STEP = 1.002; // frequency ratio by which to multiply the fundamental in search of the correct pitch
const int MIN_FREQ=80;
const int MAX_FREQ=1600;
const int BUFFER_SIZE = 4096;
const int OVERSAMPLE = 4;
const int MAX_DB_DIST = 18;
const int ZERO_THRESH = 8; // max distance between the expected and actual end of a wave (for OLA)

bool AnalysisEngine::processAudio() {
    SF_INFO sfinfo;
    SNDFILE* infile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    sampleRate=sfinfo.samplerate;

    if (!infile) {
        std::cerr << "Could not open file!" << std::endl;
        return 1;
    }

	// Parameters
    int N = BUFFER_SIZE;
    int hop_size = N / 2 / OVERSAMPLE;
    binCount=hop_size; // TODO: determine if these are redundant
    frameSize=hop_size;
    
    // Buffers
    std::vector<double> full_buffer(N, 0.0);
    double* in = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

	std::vector<double> hannWindow(N);
    for (int i = 0; i < N; i++) {
		hannWindow[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (N - 1)));
    }

	sf_count_t read_count;
	int frame_count = 0;
	long frame_offset = 0;
	maxRMS = 0;

	while((read_count = sf_read_double(infile, full_buffer.data() + (N - hop_size), hop_size)) > 0){
		double rms;
		for(int i=0; i < N; i++) {
			in[i]=full_buffer[i] * hannWindow [i];
			rms+=in[i]*in[i];
		}
		rms = std::pow(rms / N, 0.5); // avg and root
		if(maxRMS < rms){
			maxRMS = rms;
		}
		fftw_execute(p);
		
		// 3. Process Magnitudes (convert to dB)
		double maxMag=0;
		double maxdB=-100;
		double avgdB;
		double maxBand=0;
		std::vector<double> peakFreq;
		std::vector<double> peakMag;
		//std::vector<unsigned char> mags(N/4);
		for (int i = 1; i <= N / 4; i++) {
			double mag=sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1])*i;
			double db = 20 * log10(mag + 1e-9); // 1e-9 prevents log0
			//mags[i] = (unsigned char)std::round(mag/10.0);
			if(db>maxdB){
				maxdB=db;
				maxMag=mag;
				maxBand=i;
			} // loop thru the data once to find the max dB
			avgdB+=db;
		}
		avgdB /= (N/4);
		double lastdB=-100;
		bool wayUp=false;
		for(int i=1; i <= N/4; i++){
			double mag=sqrt(out[i][0] * out[i][0] + out[i][1] * out[i][1])*i; // sum sin and cos components, adjusting according to the pink noise curve
			double db = 20 * log10(mag + 1e-9); // 1e-9 prevents log0
			if(db > lastdB){ // on the way up
				wayUp=true;
			}else{ // on the way down
				if(wayUp){ // if we WERE PREVIOUSLY on the way up
					if(maxdB - lastdB < MAX_DB_DIST  && i < N/4){ // if the peak is within [dyn]dB of the max AND the frequency is less than 1/4 the sample rate
						peakFreq.push_back((sfinfo.samplerate * i / N));
						peakMag.push_back(lastdB);
					}
				}
				wayUp=false;
			}
			lastdB=db;
		}
		//if(frame_count<8){
		double max_correlation=0;
		double best_fit=0;
		for(double fund=MIN_FREQ; fund < MAX_FREQ; fund*=RATIO_STEP){ // search for the best fitting fundamental
			double correlation=0;
			for(int i=0; i<peakFreq.size(); i++){
				double ratio=peakFreq[i]/fund;
				correlation += std::cos(ratio * M_PI * 2); // add a value that aligns the peaks of a cosine wave with whole numbers
			}
			correlation /= peakFreq.size(); // divide by the number of peaks
			if(correlation > max_correlation){ // if this new frequency fits better than any others
				max_correlation = correlation;
				best_fit=fund;
			}
		}
		double time=((double)frame_offset + (N / 2)) / (double)sfinfo.samplerate; // the "time" of the frame is in the middle of the FFT window
		SpectrogramFrame sf;
		sf.peakFreqs=peakFreq;
		sf.peakMags=peakMag;
		sf.fundamental=best_fit;
		sf.tunedFrequency=best_fit;
		sf.correlation=max_correlation;
		sf.time=time;
		//std::cout << "Reading frame: " << frame_count << ", time: " << time << std::endl;
		
		std::vector<double> harmData(NUM_HARMS);
		double maxHarm=-100; // find the maximum harmonic amplitude in dB
		for(int i=0; i<harmData.size(); i++){ // gather harmonic data
			double test_freq=best_fit*(i+1);
			int test_bin=std::round((double)N / (double)sfinfo.samplerate * test_freq);
			if(test_bin < N / 4){ // if the frequency is within bounds
				double mag=std::sqrt(out[test_bin][0] * out[test_bin][0] + out[test_bin][1] * out[test_bin][1]);
				double db= 20 * log10(mag + 1e-9);
				harmData[i]=db;
				if(db > maxHarm){
					maxHarm=db;
				}
			}else{
				harmData[i]=-100;  // set all unknown / unmeasurable harmonics to -100db
			}
		}
		
		for(int i=0; i<harmData.size(); i++){
			double currentHarm = harmData[i] -  maxHarm; // normalize harms so that the max value is 0db
			harmData[i] = currentHarm;
		}
		for(int i=0; i<N-hop_size; i++){
			full_buffer[i]=full_buffer[i+hop_size];
		}
		sf.harmData=harmData;
		sf.rms=rms;
		//if(frame_count % 10 == 1){ // let's fuck with it
			frames.push_back(sf);
		//}

		frame_count++;
		frame_offset += hop_size;
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    sf_close(infile);
    return true;
}

bool AnalysisEngine::calculateHarmonicData(){
    SF_INFO sfinfo;
    SNDFILE* infile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    sampleRate=sfinfo.samplerate;

    if (!infile) {
        std::cerr << "Could not open file!" << std::endl;
        return 1;
    }

	// Parameters
    int N = BUFFER_SIZE;
    int hop_size = N / 2 / OVERSAMPLE;
    // Buffers
    std::vector<double> full_buffer(N, 0.0);
    double* in = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N / 2 + 1));
    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);

    // Prepare Hann window pre-calculated
	std::vector<double> hannWindow(N);
    for (int i = 0; i < N; i++) {
		hannWindow[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (N - 1)));
    }
    
	sf_count_t read_count;
	int frame = 0;
	long frame_offset = 0;

	while((read_count = sf_read_double(infile, full_buffer.data() + (N - hop_size), hop_size)) > 0 && frame < frames.size()){
    	for(int i=0; i < N; i++) {
			in[i]=full_buffer[i] * hannWindow [i];
		}

		fftw_execute(p);
		double maxHarm = -100; // amplitude (dB) of max harmonic
		if(frame < frames.size()){
			std::vector<double> newHarmData (NUM_HARMS, 0.0);
			double time = (double)(frame_offset + hop_size) / (double)sfinfo.samplerate;
			// the OLD WAY: only test 1 FFT bin per harmonic
			// the NEW WAY: test ALL FFT bins and weight them according to correlation
			for(int i=0; i<N/2; i++){
				double fund=frames[frame].tunedFrequency;
				double mag=std::sqrt(std::pow(out[i][0], 2) + std::pow(out[i][1], 2));
				double freq=((double)sfinfo.samplerate * (double)i / (double)N);
				double ratio=freq/fund;
				int whichHarm=std::round(ratio) - 1;
				double correlation = std::cos(ratio * M_PI * 2)*0.5 + 0.5;
				if(whichHarm >= 0 && whichHarm < NUM_HARMS){
					newHarmData[whichHarm] += mag * correlation;
				}
			}
			for(int i=0; i<NUM_HARMS; i++){
				newHarmData[i] = 20 * log10(newHarmData[i] + 1e-9); // convert all to dB
				if(newHarmData[i] > maxHarm){
					maxHarm = newHarmData[i];
				}
			}
			for(int i=0; i<NUM_HARMS; i++){
				double currentHarm = newHarmData[i] - maxHarm; // normalize harms so that the max value is 0db
				frames[frame].harmData[i] = currentHarm;
			}
			//std::cout << "Finished analyzing harmonic data for frame " << frame << std::endl;
			int origFrame=frame;
			while(time > frames[frame].time && frame < frames.size()-1){ // if we have to skip a frame
				if(frame > origFrame){ // if we've already skipped a frame
					//std::cout << "Time: " << time << ", frame " << frame << " time: " << frames[frame].time << std::endl;
					//std::cout << "Copying harmonic data from frame " << origFrame << " to frame " << frame << std::endl;
					for(int h=0; h<NUM_HARMS; h++){
						frames[frame].harmData = frames[origFrame].harmData; // copy harmonic data from previous frame
					}
				}
				frame++;
			}
			for(int i=0; i<N-hop_size; i++){
				full_buffer[i]=full_buffer[i+hop_size];
			}
			frame_offset += hop_size;
		}else{
			std::cerr << "WARNING: input stream too long (check if file was externally modified)" << std::endl;
		}
    }

    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    sf_close(infile);
    std::cout << "Recalculated harmonic data..." << std::endl;
    return true;
		
}
bool AnalysisEngine::smooth(){ // use statistical analysis to get a more accurate pitch detection
	double avgFund;
	double numFunds;
	int windowSize = OVERSAMPLE * 4;
	std::vector<double> localAverage(getFrameCount(), 0.0);
	for(int f=0; f<getFrameCount(); f++){
		avgFund += getFund(f) * getCorrelation(f);
		numFunds += getCorrelation(f);
		int minFrame = f - windowSize;
		int maxFrame = f + windowSize + 1; // never actually "reaches" max frame but that's OK
		if(minFrame < 0){
			minFrame = 0;
		}
		if(maxFrame > getFrameCount()){
			maxFrame = getFrameCount();
		}
		double totalCorr = 0;
		for(int i=minFrame; i<maxFrame; i++){
			localAverage[f] += freqToMIDI(getFund(i)); // convert to MIDI to get a multiplicative average
			totalCorr += getCorrelation(i);
		}
		localAverage[f] /= (maxFrame - minFrame);
		frames[f].averageCorrelation = totalCorr / (maxFrame - minFrame);
		frames[f].localAverage=MIDIToFreq(localAverage[f]);
	}
	avgFund /= numFunds;
	std::cout << "Weighted average fundamental: " << avgFund << std::endl;
	
	double MAX_JUMP = 1.6;
	for(int f=0; f<getFrameCount(); f++){
		double max_correlation=0;
		double best_fit=0;
		std::vector<double> peakFreqs=frames[f].peakFreqs;
		std::vector<double> peakMags=frames[f].peakMags;
		double minSearch = avgFund / 2;
		double maxSearch = avgFund * 2;
		if(f > 0 && f < getFrameCount() -1){ // if this is not the first or last frame
			double expectedFund = MIDIToFreq(localAverage[f]);
			minSearch = expectedFund / MAX_JUMP;
			maxSearch = expectedFund * MAX_JUMP;
		}
		/*if(maxSearch > avgFund * 2 || minSearch < avgFund / 2){
			maxSearch = avgFund * 2;
			minSearch = avgFund / 2;
			std::cout << "Range restricted" << std::endl;
		}*/ // keep things within bounds if they get crazy
		for(double fund=minSearch; fund < maxSearch; fund*=RATIO_STEP){ // search for the best fitting fundamental, THIS TIME only searching within 2 octaves from the average
			double correlation=0;
			for(int i=0; i<peakFreqs.size(); i++){
				double ratio=peakFreqs[i]/fund;
				correlation += std::cos(ratio * M_PI * 2); // add a value that aligns the peaks of a cosine wave with whole numbers
			}
			correlation /= peakFreqs.size(); // divide by the number of peaks
			if(correlation > max_correlation){ // if this new frequency fits better than any others
				max_correlation = correlation;
				best_fit=fund;
			}
		}
		//std::cout << "Best fit: " << best_fit << std::endl;
		frames[f].fundamental=best_fit;
		frames[f].tunedFrequency=best_fit;
		//frames[f].fundamental=localAverage[f];
		frames[f].correlation=max_correlation;
		//std::cout << "New fundamental: " << frames[f].fundamental << std::endl;
		//std::cout << "New fundamental for frame (" << f << "): " << best_fit << std::endl;
	}
	return true;
}
bool AnalysisEngine::smoothPitch(){ // use a low-pass averaging algorithm to force the pitch tracking to be a little more human friendly
	int windowSize=8;
	std::vector<double> hannWindow(windowSize);
	std::vector<double> newPitches(frames.size());
    for (int i = 0; i < windowSize; i++) {
		hannWindow[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (windowSize - 1)));
    }
	for(int i=0; i<getFrameCount(); i++){
		double totalPitch=0;
		double totalCorr=0;
		for(int f=0; f<windowSize; f++){
			int off=f-(windowSize)/2;
			if(i+off>=0 && i+off < frames.size()){
				totalPitch += freqToMIDI(frames[i+off].tunedFrequency)*frames[i+off].correlation * hannWindow[f];
				totalCorr += frames[i+off].correlation * hannWindow[f];
			}
		}
		newPitches[i] = MIDIToFreq(totalPitch / totalCorr);
	}
	for(int i=0; i<frames.size(); i++){
		frames[i].tunedFrequency = newPitches[i];
	}
	return true;
}
bool AnalysisEngine::tune(){
	for(int i=0; i<getFrameCount(); i++){
		double fund=frames[i].tunedFrequency;
		double midi = std::round(freqToMIDI(fund));
		frames[i].tunedFrequency = MIDIToFreq(midi);
	}
	return true;
}
bool AnalysisEngine::untune(double amount){
	for(int i=0; i<getFrameCount(); i++){
		double fund=frames[i].tunedFrequency;
		double pitch=frames[i].fundamental;
		frames[i].tunedFrequency = (fund * amount) + (pitch * (1 - amount));
	}
	return true;
}
bool AnalysisEngine::exportUnpitchedAudio(const std::string& outf){
    SF_INFO sfinfo;
    SNDFILE* infile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    SNDFILE* outfile = sf_open(outf.c_str(), SFM_WRITE, &sfinfo);
    sampleRate=sfinfo.samplerate;
    //std::cout << "Sample rate: " << sfinfo.samplerate << std::endl;

    if (!infile) {
        std::cerr << "Could not open file!" << std::endl;
        return 1;
    }

	// Parameters
    int N = BUFFER_SIZE;
    int hop_size = N / 4;
    binCount=hop_size; // TODO: determine if these are redundant
    frameSize=hop_size;
    
    // Buffers
    std::vector<double> full_buffer(N, 0.0);
    std::vector<double> full_buffer_out(N, 0.0);
    double* in = (double*) fftw_malloc(sizeof(double) * N);
    double* outAudio = (double*) fftw_malloc(sizeof(double) * N);
    fftw_complex* out = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
	fftw_complex* out2 = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * (N));
    fftw_plan p = fftw_plan_dft_r2c_1d(N, in, out, FFTW_MEASURE);
    fftw_plan p2 = fftw_plan_dft_c2r_1d(N, out2, outAudio, FFTW_MEASURE);

    // Prepare Hann window pre-calculated
	std::vector<double> hannWindow(N);
    for (int i = 0; i < N; i++) {
		hannWindow[i] = 0.5 * (1.0 - cos(2.0 * M_PI * i / (N - 1)));
    }

	sf_count_t read_count;
	int frame = 0;
	long frame_offset = 0;
	std::default_random_engine generator;
	std::uniform_real_distribution<double> distribution(-M_PI, M_PI);
	std::vector<double> phase(BUFFER_SIZE, 0.0);

	// THE CORE LOOP
	// We read 'hop_size' new samples and shift the old ones back
	double fund=0.0;

	while((read_count = sf_read_double(infile, full_buffer.data() + (N - hop_size), hop_size)) > 0){
		double time = (double)frame_offset / (double)sfinfo.samplerate;
		while(time > frames[frame].time && frame < frames.size() -1){
			frame++;
		}
		double rms;
		for(int i=0; i < N; i++) {
			in[i]=full_buffer[i] * hannWindow [i];
			rms+=in[i]*in[i];
		}
		rms = std::pow(rms / N, 0.5); // avg and root
		fftw_execute(p);
		
		if(frame < getFrameCount()){ // if we have analysis data for this frame
			fund=frames[frame].fundamental; // TODO: interpolate this
		}
		if(fund != 0){ //filter the frequencies
			for(int i=0; i<N/2; i++){
				double freq = (double)sfinfo.samplerate / (double)N * (double) i;
				double ratio=freq/fund;
				double corr = std::cos(ratio * M_PI * 2)*0.5+0.5;
				double mag = std::sqrt(std::pow(out[i][0], 2) + std::pow(out[i][1], 2));
				mag *= std::pow(1 - corr, 4); // apply correlation "filtering" to magnitude
				double measured_phase = std::atan2(out[i][1], out[i][0]);
				if(frame == 64 && i < 200){
					std::cout << "Frame 64 fundamental: " << fund << std::endl;
					std::cout << "Bin: " << i << std::endl;
					std::cout << "Frequency of bin: " << freq << std::endl;
					std::cout << "Ratio: " << ratio << std::endl;
					std::cout << "Magnitude: " << mag << std::endl;
					std::cout << "Correlation: " << corr << std::endl;
				}
				double randomAng = distribution(generator);
				double amount = 0 * (double)i / (double) (N/2); // increase smearing with frequency
				phase[i] += randomAng * amount;
				if(phase[i] > M_PI){
					phase[i] -= M_PI;
				}
				if(phase[i] < -M_PI){
					phase[i] += M_PI;
				}
				out2[i][0] = mag * std::cos(measured_phase + phase[i]);
				out2[i][1] = mag * std::sin(measured_phase + phase[i]);
			}
		}
		fftw_execute(p2);

		for(int i=0; i< N; i++){
			double normalized = (outAudio[i] / N) * hannWindow[i] / 1.5;
			full_buffer_out[i] += normalized;
		}
		// write here?
		sf_write_double(outfile, full_buffer_out.data(), hop_size);
		
		for(int i=0; i<N-hop_size; i++){
			full_buffer[i]=full_buffer[i+hop_size];
			full_buffer_out[i]=full_buffer_out[i+hop_size];
		}
		for(int i=N-hop_size; i<N; i++){
			full_buffer_out[i] = 0.0;
		}
		frame_offset += hop_size;
    }
    fftw_destroy_plan(p);
    fftw_free(in);
    fftw_free(out);
    sf_close(infile);
    sf_close(outfile);
    return true;
}
bool AnalysisEngine::exportTunedAudio(const std::string& outf){
    SF_INFO sfinfo;
    SNDFILE* infile = sf_open(filename.c_str(), SFM_READ, &sfinfo);
    SNDFILE* outfile = sf_open(outf.c_str(), SFM_WRITE, &sfinfo);
    sampleRate=sfinfo.samplerate;
    
    if (!infile) {
        std::cerr << "Could not open file!" << std::endl;
        return 1;
    }

	// Parameters
    int N = BUFFER_SIZE;
    int hop_size = N/2;
    // Buffers
    std::vector<double> full_buffer(N, 0.0);
    std::vector<double> full_buffer_out(N, 0.0);
    std::vector<double> wave_period(N, 0.0);
	sf_count_t read_count;
	long frame_offset = 0;

	double fund = 0.0;
	double pitch = 0.0;
	double readStart = 0.0;
	double readEnd = 0.0;
	double writeStart = 0.0;
	double writeEnd = 0.0;
	int readOffset = 0; // the offset of the read RELATIVE to the write "heads" : read = write + off
	int targetDist = 0; // the write period length
	int readDist = 0; // read period length
	double phase = 0;
	int lastFullReadStart = 0;
	int lastFullWriteStart = 0;
	bool lastBufferTuned = false;
	int lastReadOffset = 0;
		
	while((read_count = sf_read_double(infile, full_buffer.data() + (N - hop_size), hop_size)) > 0){
		double time = (double)frame_offset / (double)sfinfo.samplerate;
		int frameLeft = 0;
		while(time > frames[frameLeft].time && frameLeft < frames.size() -1){
			frameLeft++; // find the correct frame
		}
		if(frameLeft > 0){
			frameLeft--;
		}
		double corr = frames[frameLeft].correlation;
		double rms = frames[frameLeft].rms;
		fund = frames[frameLeft].fundamental; // roughly calculate pitch-dependent variables
		pitch = frames[frameLeft].tunedFrequency;
		targetDist = sfinfo.samplerate * 2 / pitch; // the wavelength to aim for (2 periods)
		readDist = sfinfo.samplerate * 2 / fund;
		//std::cout << "Read offset at BEGINNING of buffer: " << readOffset << std::endl;
		if(corr > CORR_THRESH && rms / maxRMS > RMS_THRESH){
			bool finished = false;
			if(!lastBufferTuned || writeStart == 0){ // if read/write variables not initialized
				writeStart = N/4 - (targetDist / 2.0); // Set it up so that phase will be at 0.5; therefore we'll be in the middle of the Hann window
				readStart = writeStart + readOffset;
				phase = 0.5;
				//std::cout << "Time: " << time << ", phase-dependent data not set" << std::endl;
			}
			while(!finished){
				// Pitch and time dependent variables:  calculate them once per wave cycle
				int sampleIndex = frame_offset + writeStart;
				time = ((double)sampleIndex) / (double)sfinfo.samplerate;
				if(frameLeft > 0 && frames[frameLeft].time > time){ // if we're too far along
					frameLeft --;
				}
				if(frameLeft < frames.size() -1 && frames[frameLeft + 1].time < time){ // if the next frame's time is in the past
					frameLeft ++; // advance the frame
				}
				double interpolation = 0;
				if(frameLeft < frames.size() -1){ // if not last frame
					interpolation = (time - frames[frameLeft].time) / (frames[frameLeft + 1].time - frames[frameLeft].time);
				}
				if(time > frames[frameLeft].time && frameLeft < frames.size() - 1){ // if we're not at the beginning or end of the frame sequence
					fund = frames[frameLeft].fundamental + (frames[frameLeft + 1].fundamental - frames[frameLeft].fundamental)*interpolation;
					pitch = frames[frameLeft].tunedFrequency + (frames[frameLeft + 1].tunedFrequency - frames[frameLeft].tunedFrequency)*interpolation;
				}else{ // interpolate pitch between frames
					fund = frames[frameLeft].fundamental;
					pitch = frames[frameLeft].tunedFrequency;
				}
				targetDist = sfinfo.samplerate * 2 / pitch; // the wavelength to aim for (2 periods)
				readDist = sfinfo.samplerate * 2 / fund; // the wavelength of the data being read (2 periods)
				writeEnd = writeStart + targetDist;
				readEnd = readStart + readDist;
				int i = writeStart;
				while(i < writeEnd){
					phase = (double)(i - writeStart) / (double)targetDist; // Phase is determined based on WRITE position
					double fade = 0.5 - std::cos(phase * M_PI * 2.0) / 2.0; // hann window type thing for each wave cycle
					readOffset = readStart - writeStart; // offset tells you where to go from writeStart
					full_buffer_out[i] += full_buffer[i + readOffset] * fade; // read with offset
					if(writeStart > N * 3 / 4){ // if we're out of bounds
						finished = true; // exit the loop after finishing this write
					}
					if(i == N * 3 / 4){ // if we're EXACTLY at the window end
						//int progression = (readEnd - writeEnd) - (readStart - writeStart); // how much we expect the wave nodes to move over one cycle;
						lastReadOffset = readOffset;
						//std::cout << "Reached end of usable buffer" << std::endl;
						//std::cout << "Offset at end: " << lastReadOffset << std::endl;
					}
					i++;
				}
				if(i < N * 3 / 4){ // if we haven't touched the end yet
					lastFullReadStart = readStart;
					lastFullWriteStart = writeStart;
				}
				// use inner product to calculate how much to shift read window
				double minDiff = -1;
				int minDist = -1;
				int addOff = 0; // additional offset for moving forward or backward one wave cycle
				if(readOffset > readDist){ // if we're one whole wave cycle ahead
					addOff = -readDist; // move back one cycle
					//std::cout << "Reading too far forward, moving back one cycle at sample " << frame_offset + i << std::endl;
				}else if(readOffset < -readDist){ // if we're one cycle behind
					addOff = readDist; // move one cycle ahead
					//std::cout << "Reading too far behind, moving forward one cycle at sample " << frame_offset + i << std::endl;
				}
				for(int testDist = readDist * 0.4; testDist < readDist * 0.6; testDist ++){ // look for "best fit" around 1 wave period away
					double total = 0;
					int totalSamples = 0;
					for(int x = readStart; x < readEnd; x++){
						total += std::pow(full_buffer[x + testDist + addOff] - full_buffer[x], 2);
						totalSamples++;
					}
					double rms = std::sqrt(total / (double)totalSamples);
					if(minDist == -1){
						minDiff = rms;
						minDist = testDist;
					}
					if(rms < minDiff){
						minDiff = rms;
						minDist = testDist;
					}
				}
				readStart += minDist + addOff; // move the read head to the next wave period
				readEnd = readStart + readDist;
				writeStart += targetDist / 2; // advance the write head one calculated period
				writeEnd = writeStart + targetDist;
				phase = 0;
			} // "while(!finished)"
			readOffset = readStart - writeStart;
			writeStart = lastFullWriteStart - hop_size; // resume overlapping in the next buffer
			readStart = lastFullReadStart - hop_size;
			lastBufferTuned = true;
			readOffset = lastReadOffset;
		}else{ // if not enough pitch information
			for(int i=0; i < N; i++){
				double fade = 0;
				if(i < N / 4){ // fade in
					fade = (double) i / (double)N * 4.0;
				}else if(i > N * 3 / 4){ // fade out
					fade = 1.0 - (double) (i - (N * 3 / 4)) / (double)N * 4.0;
				}else{
					fade = 1;
				}
				full_buffer_out[i] += full_buffer[i + readOffset] * fade; // just copy the audio with offset
			}
			lastBufferTuned = false;
		} // end of correlation thresh "if"
		if(frame_offset > 0){ // if it's not the first loop and we have a full input buffer
			sf_write_double(outfile, full_buffer_out.data() + hop_size / 2, hop_size); // write the MIDDLE HALF of the buffer
		}
		frame_offset += hop_size;
		for(int i=0; i<N-hop_size; i++){
			full_buffer[i]=full_buffer[i+hop_size];
			full_buffer_out[i] = full_buffer_out[i + hop_size]; // shift 2nd half of output buffer to 1st half
		} // shift the data to the first half of array
		for(int i = 0; i < N; i++){
			full_buffer_out[i] = 0; // fill 2nd half of output buffer with 0s
		}
		//std::cout << "Read offset at END of buffer: " << readOffset << std::endl;
    }
    sf_close(infile);
    sf_close(outfile);
    return true;
}
bool AnalysisEngine::exportPitchedAudio(const std::string& outf){
    SF_INFO sfinfo;
	sfinfo.samplerate=DEFAULT_SAMPLE_RATE;
	sfinfo.channels=1;
	sfinfo.format=SF_FORMAT_WAV|SF_FORMAT_PCM_16;	
    SNDFILE* outfile = sf_open(outf.c_str(), SFM_WRITE, &sfinfo);
    sampleRate=sfinfo.samplerate;

    double amp=0;
    double theta=0;
    //vector<double> lastHarmData(32);
    //vector<double> nextHarmData(32);
    //vector<double> currentHarmData(32);
    double frame=0;
    double maxTime = getClipDuration();
    long total_samples=(long)std::round((double)maxTime*(double)sfinfo.samplerate);
	std::cout << "Clip duration: " << maxTime << std::endl;
    std::cout << "Number of samples: " << total_samples << std::endl;
	std::vector<double> buffer(BUFFER_SIZE, 0.0);
    
	double lastTime=0;
	double nextTime=0;
    for(int offset=0; offset<total_samples - BUFFER_SIZE; offset+=BUFFER_SIZE){ // buffer loop
		for(int i=0; i<BUFFER_SIZE; i++){ // loop for processing a single buffer
			double time=(double)(i+offset)/(double)sfinfo.samplerate;
			double freq;
			double amp;
			double corr;

			bool lastFrame=(frame == frames.size()-1);
			if(time >= nextTime){ // if we're onto the next frame
				if(!lastFrame){ // move to the next frame
					frame++;
					lastTime=frames[frame].time;
					nextTime=frames[frame+1].time;
				}else{ // final frame
					lastTime=frames[frame].time;
					nextTime = lastTime+0.0001;
				}
			}
			lastFrame=(frame == frames.size()-1); // re-check if it's the last frame
			
			double interpolation=0;
			if(!lastFrame && nextTime - lastTime !=0){ // avoid divide by 0
				interpolation = (time - lastTime) / (nextTime - lastTime);
			}
			double sample=0;
			for(int h=0; h<NUM_HARMS; h++){
				double thisHarm; // interpolate harmonic data
				if(lastFrame){
					thisHarm=frames[frame].harmData[h];
					freq=frames[frame].tunedFrequency;
					amp=frames[frame].rms;
					corr=frames[frame].correlation;
				}else{
					thisHarm=frames[frame].harmData[h] + (interpolation * (frames[frame+1].harmData[h]-frames[frame].harmData[h]));
					freq=frames[frame].tunedFrequency + (interpolation * (frames[frame+1].tunedFrequency - frames[frame].tunedFrequency));
					amp=frames[frame].rms + (interpolation * (frames[frame+1].rms - frames[frame].rms));
					corr=frames[frame].correlation + (interpolation * (frames[frame+1].correlation - frames[frame].correlation));
				}
				double harmAmp=std::pow(10, thisHarm/20); // convert from dB to ratio
				sample+=std::sin(theta * ((double)h+1.0)) * (amp * std::sqrt(2))*harmAmp*std::pow(corr, 1); // h+1 because the "zeroth" harmonic is not useful here, harm 0 = fundamental
			}
			buffer[i]=sample;
			theta+=freq / sfinfo.samplerate * M_PI * 2;
		}
		sf_count_t count = sf_write_double(outfile, buffer.data(), BUFFER_SIZE);
		
		if (count != BUFFER_SIZE) {
		    std::cerr << "Warning: Failed to write all samples." << std::endl;
		}
    }
    sf_close(outfile);
    return true;
}
bool AnalysisEngine::tuneFrame(int frame, double pitch){
	if(frame > 0 && frame < frames.size()){
		frames[frame].tunedFrequency = pitch;
		return true;
	}else{
		return false;
	}
}
double AnalysisEngine::getMinimumPitch(double thresh){
	double min = MAX_FREQ;
	for(int i=0; i<getFrameCount(); i++){
		if(getCorrelation(i)>thresh && getRMS(i) / maxRMS > RMS_THRESH && getFund(i)<min){
			min = getFund(i);
		}
	}
	return min;
}
double AnalysisEngine::getMaximumPitch(double thresh){
	double max = MIN_FREQ;
	for(int i=0; i<getFrameCount(); i++){
		if(getCorrelation(i)>thresh && getRMS(i) / maxRMS > RMS_THRESH && getFund(i)>max){
			max = getFund(i);
		}
	}
	return max;
}
double AnalysisEngine::freqToMIDI(double freq){
	return std::log2(freq/440.0)*12.0 + 69.0;
}
double AnalysisEngine::MIDIToFreq(double midi){
	return std::pow(2, (double)(midi - 69.0) / 12.0) * 440.0;
}