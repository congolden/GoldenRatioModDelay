#pragma once

#include "fxobjects.h"

struct ConnorModulatedDelayParameters
{
	ConnorModulatedDelayParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	ConnorModulatedDelayParameters& operator=(const ConnorModulatedDelayParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		algorithm = params.algorithm;
		lfoRate_Hz = params.lfoRate_Hz;
		lfoDepth_Pct = params.lfoDepth_Pct;
		feedback_Pct = params.feedback_Pct;
		width_Pct = params.width_Pct;
		return *this;
	}

	// --- individual parameters
	modDelaylgorithm algorithm = modDelaylgorithm::kFlanger; ///< mod delay algorithm
	double lfoRate_Hz = 0.0;	///< mod delay LFO rate in Hz
	double lfoDepth_Pct = 0.0;	///< mod delay LFO depth in %
	double feedback_Pct = 0.0;	///< feedback in %
	double width_Pct = 50.0;
};

/**
\class ModulatedDelay
\ingroup FX-Objects
\brief
The ModulatedDelay object implements the three basic algorithms: flanger, chorus, vibrato.

Audio I / O :
	-Processes mono input to mono OR stereo output.

Control I / F :
	-Use ModulatedDelayParameters structure to get / set object params.

\author Will Pirkle http ://www.willpirkle.com
\remark This object is included in Designing Audio Effects Plugins in C++ 2nd Ed.by Will Pirkle
\version Revision : 1.0
\date Date : 2018 / 09 / 7
*/
class ConnorModulatedDelay : public IAudioSignalProcessor
{
public:
	ConnorModulatedDelay() {
	}		/* C-TOR */
	~ConnorModulatedDelay() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- create new buffer, 100mSec long
		delay.reset(_sampleRate);
		delay.createDelayBuffers(_sampleRate, 100.0);

		// --- lfo
		lfo.reset(_sampleRate);
		OscillatorParameters params = lfo.getParameters();
		params.waveform = generatorWaveform::kTriangle;
		lfo.setParameters(params);

		return true;
	}

	/** process input sample */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		float input = xn;
		float output = 0.0;
		processAudioFrame(&input, &output, 1, 1);
		return output;
	}

	/** return true: this object can process frames */
	virtual bool canProcessAudioFrame() { return true; }

	/** process STEREO audio delay of frames */
	virtual bool processAudioFrame(const float* inputFrame,		/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
		float* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		// --- make sure we have input and outputs
		
		if (inputChannels == 0 || outputChannels == 0)
			return false;

		// --- render LFO
		SignalGenData lfoOutput = lfo.renderAudioOutput();

		// --- setup delay modulation
		AudioDelayParameters params = delay.getParameters();
		double minDelay_mSec = 0.0;
		double maxDepth_mSec = 0.0;

		// --- set delay times, wet/dry and feedback
		if (parameters.algorithm == modDelaylgorithm::kFlanger)
		{
			minDelay_mSec = 0.1;
			maxDepth_mSec = 7.0;
			params.wetLevel_dB = -3.0;
			params.dryLevel_dB = -3.0;
		}
		if (parameters.algorithm == modDelaylgorithm::kChorus)
		{
			minDelay_mSec = 10.0;
			maxDepth_mSec = 5.0;
			params.wetLevel_dB = -3.0;
			params.dryLevel_dB = -0.0;
			params.feedback_Pct = 0.0;
		}
		if (parameters.algorithm == modDelaylgorithm::kVibrato)// in my case a doubler
		{
			minDelay_mSec = 70.0;
			maxDepth_mSec = 5.0;
			params.wetLevel_dB = -3.0;
			params.dryLevel_dB = -3.0;
			params.feedback_Pct = 0.0;
		}

		// --- calc modulated delay times
		double depth = parameters.lfoDepth_Pct / 100.0;
		double modulationMin = minDelay_mSec;
		double modulationMax = minDelay_mSec + maxDepth_mSec;

		// --- flanger - unipolar
		if (parameters.algorithm == modDelaylgorithm::kFlanger)
			params.leftDelay_mSec = doUnipolarModulationFromMin(bipolarToUnipolar(depth * lfoOutput.normalOutput),
				modulationMin, modulationMax);
		else
			params.leftDelay_mSec = doBipolarModulation(depth * lfoOutput.normalOutput, modulationMin, modulationMax);


		// --- set right delay to match (*Hint Homework!)
		params.rightDelay_mSec = params.leftDelay_mSec + cookedWidth;

		// --- modulate the delay  
		delay.setParameters(params);

		// --- just call the function and pass our info in/out
		return delay.processAudioFrame(inputFrame, outputFrame, inputChannels, outputChannels);
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return ModulatedDelayParameters custom data structure
	*/
	ConnorModulatedDelayParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param ModulatedDelayParameters custom data structure
	*/
	void setParameters(ConnorModulatedDelayParameters _parameters)
	{
		// --- bulk copy
		parameters = _parameters;
		cookedWidth = (parameters.width_Pct / 100) * 6;

		OscillatorParameters lfoParams = lfo.getParameters();
		lfoParams.frequency_Hz = parameters.lfoRate_Hz;
		if (parameters.algorithm == modDelaylgorithm::kVibrato)
			lfoParams.waveform = generatorWaveform::kSin;
		else
			lfoParams.waveform = generatorWaveform::kTriangle;

		lfo.setParameters(lfoParams);

		AudioDelayParameters adParams = delay.getParameters();
		adParams.feedback_Pct = parameters.feedback_Pct;
		delay.setParameters(adParams);
	}

private:
	ConnorModulatedDelayParameters parameters; ///< object parameters
	AudioDelay delay;	///< the delay to modulate
	double cookedWidth;
	LFO lfo;			///< the modulator
};