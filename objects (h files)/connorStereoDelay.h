#pragma once

#include "fxobjects.h"
struct connorStereoDelayParameters
{
	connorStereoDelayParameters() {}
	/** all FXObjects parameter objects require overloaded= operator so remember to add new entries if you add new variables. */
	connorStereoDelayParameters& operator=(const connorStereoDelayParameters& params)	// need this override for collections to work
	{
		if (this == &params)
			return *this;

		algorithm = params.algorithm;
		wetLevel_dB = params.wetLevel_dB;
		dryLevel_dB = params.dryLevel_dB;
		
		feedback_Pct = params.feedback_Pct;

		updateType = params.updateType;
		leftDelay_mSec = params.leftDelay_mSec;
		rightDelay_mSec = params.rightDelay_mSec;
		delayRatio_Pct = params.delayRatio_Pct;
		tap1Level_db = params.tap1Level_db;
		tap2Level_db = params.tap2Level_db;
		tap3Level_db = params.tap3Level_db;
		tap4Level_db = params.tap4Level_db;
		delay2_mSec = params.delay2_mSec;
		delay3_mSec = params.delay3_mSec;
		delay4_mSec = params.delay4_mSec;

		return *this;
	}

	// --- individual parameters
	delayAlgorithm algorithm = delayAlgorithm::kNormal; ///< delay algorithm
	double wetLevel_dB = -3.0;	///< wet output level in dB
	
	double dryLevel_dB = -3.0;	///< dry output level in dB
	double feedback_Pct = 0.0;	///< feedback as a % value

	delayUpdateType updateType = delayUpdateType::kLeftAndRight;///< update algorithm
	double leftDelay_mSec = 0.0;	///< left delay time
	double rightDelay_mSec = 0.0;	///< right delay time
	double delayRatio_Pct = 100.0;	///< dela ratio: right length = (delayRatio)*(left length)
	double tap1Level_db = -3.0;
	double tap2Level_db = 0.0;
	double tap3Level_db = 0.0;
	double tap4Level_db = 0.0;
	double delay2_mSec = 0.0;
	double delay3_mSec = 0.0;
	double delay4_mSec = 0.0;
};

class ConnorStereoDelay : public IAudioSignalProcessor
{
public:
	ConnorStereoDelay() {}		/* C-TOR */
	~ConnorStereoDelay() {}	/* D-TOR */

public:
	/** reset members to initialized state */
	virtual bool reset(double _sampleRate)
	{
		// --- if sample rate did not change
		if (sampleRate == _sampleRate)
		{
			// --- just flush buffer and return
			delayBuffer_L.flushBuffer();
			delayBuffer_R.flushBuffer();
			return true;
		}

		// --- create new buffer, will store sample rate and length(mSec)
		createDelayBuffers(_sampleRate, bufferLength_mSec);

		return true;
	}

	/** process MONO audio delay */
	/**
	\param xn input
	\return the processed sample
	*/
	virtual double processAudioSample(double xn)
	{
		// --- read delay
		double yn = delayBuffer_L.readBuffer(delayInSamples_L);

		// --- create input for delay buffer
		double dn = xn + (parameters.feedback_Pct / 100.0) * yn;

		// --- write to delay buffer
		delayBuffer_L.writeBuffer(dn);

		// --- form mixture out = dry*xn + wet*yn
		double output = dryMix * xn + wetMix * yn;

		return output;
	}

	/** return true: this object can also process frames */
	virtual bool canProcessAudioFrame() { return true; }

	/** process STEREO audio delay in frames */
	virtual bool processAudioFrame(const float* inputFrame,		/* ptr to one frame of data: pInputFrame[0] = left, pInputFrame[1] = right, etc...*/
		float* outputFrame,
		uint32_t inputChannels,
		uint32_t outputChannels)
	{
		// --- make sure we have input and outputs
		if (inputChannels == 0 || outputChannels == 0)
			return false;

		// --- make sure we support this delay algorithm
		if (parameters.algorithm != delayAlgorithm::kNormal &&
			parameters.algorithm != delayAlgorithm::kPingPong)
			return false;

		// --- if only one output channel, revert to mono operation
		if (outputChannels == 1)
		{
			// --- process left channel only
			outputFrame[0] = processAudioSample(inputFrame[0]);
			return true;
		}

		// --- if we get here we know we have 2 output channels
		//
		// --- pick up inputs
		//
		// --- LEFT channel
		double xnL = inputFrame[0];

		// --- RIGHT channel (duplicate left input if mono-in)
		double xnR = inputChannels > 1 ? inputFrame[1] : xnL;

		// --- read delay LEFT
		double ynL = delayBuffer_L.readBuffer(delayInSamples_L);

		// --- read delay RIGHT
		double ynR = delayBuffer_R.readBuffer(delayInSamples_R);

		/////// read delay 2
		double yn2 = delayBuffer_L.readBuffer(delayInSamples2);

		/////// read delay 3
		double yn3 = delayBuffer_L.readBuffer(delayInSamples3);
		/////// read delay 4
		double yn4 = delayBuffer_L.readBuffer(delayInSamples4);

		// --- create input for delay buffer with LEFT channel info
		double dnL = xnL + (parameters.feedback_Pct / 100.0) * ynL;

		// --- create input for delay buffer with RIGHT channel info
		double dnR = xnR + (parameters.feedback_Pct / 100.0) * ynR;

		// --- decode
		if (parameters.algorithm == delayAlgorithm::kNormal)
		{
			// --- write to LEFT delay buffer with LEFT channel info
			delayBuffer_L.writeBuffer(dnL);

			// --- write to RIGHT delay buffer with RIGHT channel info
			delayBuffer_R.writeBuffer(dnR);
		}
		else if (parameters.algorithm == delayAlgorithm::kPingPong)
		{
			// --- write to LEFT delay buffer with RIGHT channel info
			delayBuffer_L.writeBuffer(dnR);

			// --- write to RIGHT delay buffer with LEFT channel info
			delayBuffer_R.writeBuffer(dnL);
		}
		
		//delayBuffer_L.writeBuffer(xnL);

		// --- write to RIGHT delay buffer with RIGHT channel info
		//delayBuffer_R.writeBuffer(xnR);


		// --- form mixture out = dry*xn + wet*yn
		double outputL = dryMix * xnL + wetMix *(tap1Mix * ynL + tap2Mix * yn2 + tap3Mix * yn3 + tap4Mix * yn4);

		// --- form mixture out = dry*xn + wet*yn
		double outputR = dryMix * xnR + wetMix *(tap1Mix * ynR + tap2Mix * yn2 + tap3Mix * yn3 + tap4Mix * yn4);

		// --- form mixture out = dry*xn + wet*yn
		//double outputL = wetMix * ynL + tap2Mix * yn2 + tap3Mix * yn3 + tap4Mix * yn4;

		// --- form mixture out = dry*xn + wet*yn
		//double outputR =  wetMix * ynR + tap2Mix * yn2 + tap3Mix * yn3 + tap4Mix * yn4;
		// --- set left channel
		outputFrame[0] = outputL;

		// --- set right channel
		outputFrame[1] = outputR;

		return true;
	}

	/** get parameters: note use of custom structure for passing param data */
	/**
	\return AudioDelayParameters custom data structure
	*/
	connorStereoDelayParameters getParameters() { return parameters; }

	/** set parameters: note use of custom structure for passing param data */
	/**
	\param AudioDelayParameters custom data structure
	*/
	void setParameters(connorStereoDelayParameters _parameters)
	{
		// --- check mix in dB for calc
		//if (_parameters.dryLevel_dB != parameters.dryLevel_dB)
			dryMix = pow(10.0, _parameters.dryLevel_dB / 20.0);
		if (_parameters.wetLevel_dB != parameters.wetLevel_dB)
			wetMix = pow(10.0, _parameters.wetLevel_dB / 20.0);

		if (_parameters.tap1Level_db != parameters.tap1Level_db)
			tap1Mix = pow(10.0, _parameters.tap1Level_db / 20.0);
		if (_parameters.tap2Level_db != parameters.tap2Level_db)
			tap2Mix = pow(10.0, _parameters.tap2Level_db / 20.0);
		if (_parameters.tap3Level_db != parameters.tap3Level_db)
			tap3Mix = pow(10.0, _parameters.tap3Level_db / 20.0);
		if (_parameters.tap4Level_db != parameters.tap4Level_db)
			tap4Mix = pow(10.0, _parameters.tap4Level_db / 20.0);
		// --- save; rest of updates are cheap on CPU
		parameters = _parameters;

		double newDelayInSamples2 = parameters.delay2_mSec*(samplesPerMSec);
		delayInSamples2 = newDelayInSamples2;
		double newDelayInSamples3 = parameters.delay3_mSec*(samplesPerMSec);
		delayInSamples3 = newDelayInSamples3;
		double newDelayInSamples4 = parameters.delay4_mSec*(samplesPerMSec);
		delayInSamples4 = newDelayInSamples4;

		// --- check update type first:
		if (parameters.updateType == delayUpdateType::kLeftAndRight)
		{
			// --- set left and right delay times
			// --- calculate total delay time in samples + fraction
			double newDelayInSamples_L = parameters.leftDelay_mSec*(samplesPerMSec);
			double newDelayInSamples_R = parameters.rightDelay_mSec*(samplesPerMSec);

			// --- new delay time with fraction
			delayInSamples_L = newDelayInSamples_L;
			delayInSamples_R = newDelayInSamples_R;
		}
		else if (parameters.updateType == delayUpdateType::kLeftPlusRatio)
		{
			// --- get and validate ratio
			double delayRatio = parameters.delayRatio_Pct / 100.0;
			boundValue(delayRatio, 0.0, 1.0);

			// --- calculate total delay time in samples + fraction
			double newDelayInSamples = parameters.leftDelay_mSec*(samplesPerMSec);

			// --- new delay time with fraction
			delayInSamples_L = newDelayInSamples;
			delayInSamples_R = delayInSamples_L * delayRatio;
		}
	}

	/** creation function */
	void createDelayBuffers(double _sampleRate, double _bufferLength_mSec)
	{
		// --- store for math
		bufferLength_mSec = _bufferLength_mSec;
		sampleRate = _sampleRate;
		samplesPerMSec = sampleRate / 1000.0;

		// --- total buffer length including fractional part
		bufferLength = (unsigned int)(bufferLength_mSec*(samplesPerMSec)) + 1; // +1 for fractional part

																			   // --- create new buffer
		delayBuffer_L.createCircularBuffer(bufferLength);
		delayBuffer_R.createCircularBuffer(bufferLength);
	}

private:
	connorStereoDelayParameters parameters; ///< object parameters

	double sampleRate = 0.0;		///< current sample rate
	double samplesPerMSec = 0.0;	///< samples per millisecond, for easy access calculation
	double delayInSamples_L = 0.0;	///< double includes fractional part
	double delayInSamples_R = 0.0;	///< double includes fractional part
	double delayInSamples2 = 0.0;	///< double includes fractional part
	double delayInSamples3 = 0.0;	///< double includes fractional part
	double delayInSamples4 = 0.0;	///< double includes fractional part
	double bufferLength_mSec = 0.0;	///< buffer length in mSec
	unsigned int bufferLength = 0;	///< buffer length in samples
	double wetMix = 0.707; ///< wet output default = -3dB
	double dryMix = 0.707; ///< dry output default = -3dB
	double tap1Mix = 0.707; ///< dry output default = -3dB
	double tap2Mix = 0.0; ///< dry output default = -3dB
	double tap3Mix = 0.0;
	double tap4Mix = 0.0;

	// --- delay buffer of doubles
	CircularBuffer<double> delayBuffer_L;	///< LEFT delay buffer of doubles
	CircularBuffer<double> delayBuffer_R;	///< RIGHT delay buffer of doubles
};
