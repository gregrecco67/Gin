/*
 ==============================================================================

 This file is part of the GIN library.
 Copyright (c) 2019 - Roland Rabien.

 ==============================================================================
 */

#pragma once


struct Matrix {
	inline friend Matrix operator*(const Matrix& m, const float s) { // scalar multiplication
		return { m.a * s, m.b * s, m.c * s, m.d * s };
	}
	float a, b, c, d;
};

struct StereoMatrix {
	Matrix left, right;
	inline friend StereoMatrix operator*(const StereoMatrix& m, const float s) { // scalar multiplication
			return { m.left * s, m.right * s };
		}
};

struct StereoPosition {
	float xL{ 0.f }, yL{ 0.f }, xR{ 0.f }, yR{ 0.f };
	
	inline friend StereoPosition operator*(const StereoPosition& p, const StereoMatrix& m) { // apply matrix to position
		return {.xL = m.left.a * p.xL + m.left.b * p.yL,
				.yL = m.left.c * p.xL + m.left.d * p.yL,
				.xR = m.right.a * p.xR + m.right.b * p.yR,
				.yR = m.right.c * p.xR + m.right.d * p.yR };
	}

	inline friend StereoPosition operator*(const StereoPosition& p, const float s) { // scalar multiplication
			return { p.xL * s, p.yL * s, p.xR * s, p.yR * s };
		}

	inline StereoPosition operator+(const StereoPosition otherPos) {
		return { this->xL + otherPos.xL, this->yL + otherPos.yL,
			this->xR + otherPos.xR, this->yR + otherPos.yR };
	}
};

//==============================================================================
/** Virtual Analog Stereo oscillator.
*/
class StereoOscillator
{
public:
    StereoOscillator (BandLimitedLookupTables& bllt_) : bllt (bllt_) {}

    struct Params
    {
        Wave wave = Wave::sawUp;
        float leftGain = 1.0;
        float rightGain = 1.0;
        float pw = 0.5;
        float fold = 0.0f;
        float asym = 0.0f;
    };

    /*
    	struct Params {
        Wavetype wave = Wavetype::sine;
		float tones{ 1.0 }, pan{ 0.f };
    */

    void setSampleRate (double sr)  { sampleRate = sr; }
    void noteOn (float p = -1);

    void process (float note, const Params& params, juce::AudioSampleBuffer& buffer)
    {
        buffer.clear();
        processAdding (note, params, buffer);
    }

    void processAdding (float note, const Params& params, juce::AudioSampleBuffer& buffer)
    {
        float freq = float (std::min (sampleRate / 2.0, 440.0 * std::pow (2.0, (note - 69.0) / 12.0)));
        float delta = 1.0f / (float ((1.0f / freq) * sampleRate));

        int samps = buffer.getNumSamples();
        auto l = buffer.getWritePointer (0);
        auto r = buffer.getWritePointer (1);

        for (int i = 0; i < samps; i++)
        {
            auto s = bllt.process (params.wave, note, phase, params.pw);
            postProcess (params, s);

            *l++ += s * params.leftGain;
            *r++ += s * params.rightGain;

            phase += delta;
            while (phase >= 1.0f)
                phase -= 1.0f;
        }
    }

    float qrtPhs(const float phase) {
        float p2 = phase + 0.25f;
        while (p2 >= 1.0f) {
            p2 -= 1.0f;
        }
        return p2;
    }

    // TODO handle tones / harmonics spread

    void renderPositions(float note, const Params& params, StereoPosition *positions, const int numSamples)) {
        float freq = float (std::min (sampleRate / 2.0, 440.0 * std::pow (2.0, (note - 69.0) / 12.0)));
        float delta = 1.0f / (float ((1.0f / freq) * sampleRate));
        
        for (int i = 0; i < numSamples; i++)
        {
            auto x = bllt.process (params.wave, note, phase, params.pw);
            auto y = bllt.process (params.wave, note, qrtPhase(phase), params.pw);
            postProcess (params, x);
            postProcess (params, y);

            positions[i].xL = x * params.leftGain;
            positions[i].yL = y * params.leftGain;
            positions[i].xR = x * params.rightGain;
            positions[i].yR = y * params.rightGain;
            
            phase += delta;
            while (phase >= 1.0f)
                phase -= 1.0f;
        }
    }

private:
    template<typename T>
    void postProcess (const Params& params, T& v)
    {
        if (params.asym > 0)
            v = math::lerp (v, math::pow4 (v - 1.0f) * -1.0f + 1.0f, math::pow2 (params.asym));

        if (params.fold > 0)
        {
            const auto fold = math::pow2 (math::pow2 (1.0f - params.fold)) * 1.5f;
            v = (v - ((math::max (v, fold) - fold) * T(2.0f)) - ((math::min (v, -fold) + fold) * T(2.0f)));
        }
    }

    BandLimitedLookupTables& bllt;
    double sampleRate = 44100.0;
    float phase = 0.0f;
};

struct VoicedOscillatorParams
{
    int voices = 1;
    float pan = 0.0f;
    float spread = 0.0f;
    float detune = 0.0f;
    float gain = 1.0f;
};

//==============================================================================
/** Stereo Oscillator with multiples voices, pan, spread, detune, etc
*/
template<typename O, typename P>
class VoicedStereoOscillator
{
public:
    VoicedStereoOscillator() = default;

    void setSampleRate (double sr)
    {
        for (auto o : oscillators)
            o->setSampleRate (sr);
    }

    void noteOn (float phase = -1)
    {
        for (auto o : oscillators)
            o->noteOn (phase);
    }

    void noteOn (float phases[])
    {
        for (auto idx = 0; auto o : oscillators)
            o->noteOn (phases[idx++]);
    }

    void process (float note, const P& params, juce::AudioSampleBuffer& buffer)
    {
        buffer.clear();
        processAdding (note, params, buffer);
    }

    void processAdding (float note, const P& params, juce::AudioSampleBuffer& buffer)
    {
        typename O::Params p;
        params.init (p);

        if (params.voices == 1)
        {
            p.leftGain  = params.gain * (1.0f - params.pan);
            p.rightGain = params.gain * (1.0f + params.pan);

            oscillators[0]->processAdding (note, p, buffer);
        }
        else
        {
            float baseNote  = note - params.detune / 2;
            float noteDelta = params.detune / (params.voices - 1);

            float basePan = params.pan - params.spread;
            float panDelta = (params.spread * 2) / (params.voices - 1);

            for (int i = 0; i < params.voices; i++)
            {
                float pan = juce::jlimit (-1.0f, 1.0f, basePan + panDelta * i);

                p.leftGain  = params.gain * (1.0f - pan) / float (std::sqrt (params.voices));
                p.rightGain = params.gain * (1.0f + pan) / float (std::sqrt (params.voices));

                oscillators[i]->processAdding (baseNote + noteDelta * i, p, buffer);
            }
        }
    }

protected:
    juce::OwnedArray<O> oscillators;
};

struct VoicedStereoOscillatorParams : public VoicedOscillatorParams
{
    Wave wave   = Wave::sawUp;
    float pw    = 0.5;
    float fold  = 0.0f;
    float asym  = 0.0f;

    inline void init (StereoOscillator::Params& p) const
    {
        p.wave  = wave;
        p.pw    = pw;
        p.asym  = asym;
        p.fold  = fold;
    }
};

//==============================================================================
/** Stereo Oscillator with multiple voices, pan, spread, detune, etc
 */
class BLLTVoicedStereoOscillator : public VoicedStereoOscillator<StereoOscillator, VoicedStereoOscillatorParams>
{
public:
    BLLTVoicedStereoOscillator (BandLimitedLookupTables& bllt, int maxVoices = 8)
    {
        for (int i = 0; i < maxVoices; i++)
            oscillators.add (new StereoOscillator (bllt));
    }
};
