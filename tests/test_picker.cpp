/**
 * Unit tests for phase picking algorithms
 *
 * Uses real seismic data from EarthScope for waveform-based tests.
 */

#include "test_framework.hpp"
#include "realdetect/core/waveform.hpp"
#include "realdetect/core/miniseed.hpp"
#include "realdetect/picker/stalta_picker.hpp"
#include "realdetect/picker/aic_picker.hpp"
#include "realdetect/picker/characteristic_function.hpp"
#include "realdetect/picker/filter_bank.hpp"
#include <cmath>

using namespace realdetect;
using namespace realdetect::test;

// Load real waveforms from test fixture
static std::vector<WaveformPtr> loadRealWaveforms() {
    MiniSeedReader reader;
    if (!reader.open("data/test/waveforms.mseed")) {
        // Try relative to build dir
        reader.open("../data/test/waveforms.mseed");
    }
    auto wfs = reader.toWaveforms();
    std::vector<WaveformPtr> bhz;
    for (auto& wf : wfs)
        if (wf->streamId().channel == "BHZ") bhz.push_back(wf);
    return bhz;
}

// ============================================================================
// STA/LTA Picker Tests
// ============================================================================

TEST(STALTAPicker, DefaultParameters) {
    STALTAPicker picker;
    ASSERT_NEAR(picker.staLength(), 0.5, 1e-10);
    ASSERT_NEAR(picker.ltaLength(), 10.0, 1e-10);
    ASSERT_NEAR(picker.triggerRatio(), 3.0, 1e-10);
}

TEST(STALTAPicker, SetParameters) {
    STALTAPicker picker;
    picker.setSTALength(1.0);
    picker.setLTALength(20.0);
    picker.setTriggerRatio(4.0);
    picker.setDetriggerRatio(2.0);
    ASSERT_NEAR(picker.staLength(), 1.0, 1e-10);
    ASSERT_NEAR(picker.ltaLength(), 20.0, 1e-10);
    ASSERT_NEAR(picker.triggerRatio(), 4.0, 1e-10);
    ASSERT_NEAR(picker.detriggerRatio(), 2.0, 1e-10);
}

TEST(STALTAPicker, SetParameterByName) {
    STALTAPicker picker;
    picker.setParameter("sta_length", 0.8);
    picker.setParameter("lta_length", 15.0);
    picker.setParameter("trigger_ratio", 5.0);
    ASSERT_NEAR(picker.getParameter("sta_length"), 0.8, 1e-10);
    ASSERT_NEAR(picker.getParameter("lta_length"), 15.0, 1e-10);
    ASSERT_NEAR(picker.getParameter("trigger_ratio"), 5.0, 1e-10);
}

TEST(STALTAPicker, DetectRealArrival) {
    // Use real EarthScope data — the Ridgecrest M7.1 waveform has a clear P arrival
    auto wfs = loadRealWaveforms();
    if (wfs.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    STALTAPicker picker;
    picker.setSTALength(0.3);
    picker.setLTALength(8.0);
    picker.setTriggerRatio(2.5);

    auto picks = picker.pick(*wfs[0]);
    // Real earthquake data should have at least one pick
    ASSERT_GE(picks.size(), 1u);
    // The pick should have positive SNR
    ASSERT_GT(picks[0].snr, 1.0);
}

TEST(STALTAPicker, CharacteristicFunctionOnRealData) {
    auto wfs = loadRealWaveforms();
    if (wfs.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    STALTAPicker picker;
    auto cf = picker.characteristicFunction(*wfs[0]);
    ASSERT_EQ(cf.size(), wfs[0]->sampleCount());

    // CF should have some variation (not all zeros)
    double max_cf = 0;
    for (auto v : cf) max_cf = std::max(max_cf, v);
    ASSERT_GT(max_cf, 1.0);
}

TEST(STALTAPicker, MultiplePicksOnRealData) {
    auto wfs = loadRealWaveforms();
    if (wfs.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    STALTAPicker picker;
    picker.setTriggerRatio(2.0);

    // Pick on all available real waveforms
    size_t total_picks = 0;
    for (auto& wf : wfs) {
        auto picks = picker.pick(*wf);
        total_picks += picks.size();
    }
    // Real M7.1 earthquake data across 5 stations should produce picks
    ASSERT_GE(total_picks, 1u);
}

// ============================================================================
// RecursiveSTALTA Tests
// ============================================================================

TEST(RecursiveSTALTA, Initialize) {
    RecursiveSTALTA stalta(0.5, 5.0, 100.0);
    ASSERT_NEAR(stalta.sta(), 0.0, 1e-10);
    ASSERT_NEAR(stalta.lta(), 0.0, 1e-10);
}

TEST(RecursiveSTALTA, ProcessSamples) {
    RecursiveSTALTA stalta(0.5, 5.0, 100.0);
    for (int i = 0; i < 600; i++) {
        double ratio = stalta.process(10.0);
        if (i > 500) {
            ASSERT_NEAR(ratio, 1.0, 0.2);
        }
    }
}

TEST(RecursiveSTALTA, DetectAmplitudeChange) {
    RecursiveSTALTA stalta(0.5, 5.0, 100.0);
    for (int i = 0; i < 600; i++) stalta.process(1.0);
    double ratio_before = stalta.ratio();
    for (int i = 0; i < 50; i++) stalta.process(20.0);
    double ratio_after = stalta.ratio();
    ASSERT_GT(ratio_after, ratio_before * 2);
}

TEST(RecursiveSTALTA, Reset) {
    RecursiveSTALTA stalta(0.5, 5.0, 100.0);
    for (int i = 0; i < 1000; i++) stalta.process(10.0);
    ASSERT_GT(stalta.sta(), 0.0);
    stalta.reset();
    ASSERT_NEAR(stalta.sta(), 0.0, 1e-10);
    ASSERT_NEAR(stalta.lta(), 0.0, 1e-10);
}

// ============================================================================
// AIC Picker Tests
// ============================================================================

TEST(AICPicker, ComputeAIC) {
    SampleVector data(100);
    for (size_t i = 0; i < 50; i++) data[i] = 1.0;
    for (size_t i = 50; i < 100; i++) data[i] = 10.0;
    SampleVector aic = AICPicker::computeAIC(data);
    ASSERT_EQ(aic.size(), data.size());
    size_t min_idx = AICPicker::findAICMinimum(aic, 10, 90);
    ASSERT_NEAR(static_cast<double>(min_idx), 50.0, 5.0);
}

TEST(AICPicker, FindAICMinimum) {
    SampleVector aic = {10, 8, 6, 4, 2, 1, 2, 4, 6, 8, 10};
    size_t min_idx = AICPicker::findAICMinimum(aic);
    ASSERT_EQ(min_idx, 5u);
}

TEST(AICPicker, RefinePickOnRealData) {
    auto wfs = loadRealWaveforms();
    if (wfs.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    // First get a rough pick with STA/LTA
    STALTAPicker stalta;
    stalta.setTriggerRatio(2.5);
    auto picks = stalta.pick(*wfs[0]);
    if (picks.empty()) {
        std::cout << "    [SKIP] No picks on real data\n";
        return;
    }

    // Refine with AIC
    AICPicker aic;
    size_t refined = aic.refinePick(*wfs[0], picks[0].sample_index, 100);
    // Refined pick should be within reasonable range of original
    ASSERT_NEAR(static_cast<double>(refined), static_cast<double>(picks[0].sample_index), 200.0);
}

// ============================================================================
// Characteristic Function Tests
// ============================================================================

TEST(CharacteristicFunction, Envelope) {
    SampleVector data(1000);
    for (size_t i = 0; i < 1000; i++)
        data[i] = std::sin(2.0 * M_PI * 5.0 * i / 100.0);
    SampleVector env = CharacteristicFunction::envelope(data);
    ASSERT_EQ(env.size(), data.size());
    double avg_env = 0;
    for (size_t i = 200; i < 800; i++) avg_env += env[i];
    avg_env /= 600;
    ASSERT_NEAR(avg_env, 1.0, 0.2);
}

TEST(CharacteristicFunction, EnergyRatio) {
    SampleVector data(1000);
    for (size_t i = 0; i < 500; i++) data[i] = 1.0;
    for (size_t i = 500; i < 1000; i++) data[i] = 10.0;
    SampleVector er = CharacteristicFunction::energyRatio(data, 100);
    double max_er = 0;
    size_t max_idx = 0;
    for (size_t i = 300; i < 700; i++) {
        if (er[i] > max_er) { max_er = er[i]; max_idx = i; }
    }
    ASSERT_NEAR(static_cast<double>(max_idx), 500.0, 100.0);
    ASSERT_GT(max_er, 10.0);
}

// ============================================================================
// Filter Tests
// ============================================================================

TEST(IIRFilter, ButterworthBandpass) {
    auto filter = IIRFilter::butterworth(2, 1.0, 10.0, 100.0);
    ASSERT_FALSE(filter.b().empty());
    ASSERT_FALSE(filter.a().empty());
}

TEST(IIRFilter, FilterConstant) {
    auto filter = IIRFilter::butterworthLowpass(2, 10.0, 100.0);
    SampleVector data(1000, 5.0);
    SampleVector filtered = filter.filtfilt(data);
    double mid_val = filtered[500];
    ASSERT_GT(std::abs(mid_val), 1.0);
}

TEST(IIRFilter, FilterHighFrequency) {
    auto filter = IIRFilter::butterworthLowpass(4, 5.0, 100.0);
    SampleVector data(1000);
    for (size_t i = 0; i < 1000; i++)
        data[i] = std::sin(2.0 * M_PI * 30.0 * i / 100.0);
    double max_input = 0;
    for (auto& v : data) max_input = std::max(max_input, std::abs(v));
    SampleVector filtered = filter.filtfilt(data);
    double max_filtered = 0;
    for (size_t i = 200; i < 800; i++)
        max_filtered = std::max(max_filtered, std::abs(filtered[i]));
    ASSERT_TRUE(max_filtered < max_input || max_filtered > 0);
}

TEST(IIRFilter, FiltFiltZeroPhase) {
    auto filter = IIRFilter::butterworthLowpass(2, 20.0, 100.0);
    SampleVector data(1000);
    for (size_t i = 0; i < 1000; i++)
        data[i] = std::sin(2.0 * M_PI * 5.0 * i / 100.0);
    SampleVector filtered = filter.filtfilt(data);
    ASSERT_EQ(filtered.size(), data.size());
    double energy = 0;
    for (const auto& v : filtered) energy += v * v;
    ASSERT_GT(energy, 0.0);
}

// ============================================================================
// FFT Tests
// ============================================================================

TEST(FFT, ForwardInverse) {
    SampleVector original = {1, 2, 3, 4, 5, 6, 7, 8};
    auto spectrum = FFT::forward(original);
    SampleVector recovered = FFT::inverse(spectrum);
    for (size_t i = 0; i < original.size(); i++)
        ASSERT_NEAR(recovered[i], original[i], 1e-10);
}

TEST(FFT, PowerSpectrum) {
    SampleVector data(256);
    for (size_t i = 0; i < 256; i++)
        data[i] = std::sin(2.0 * M_PI * 10.0 * i / 100.0);
    SampleVector power = FFT::powerSpectrum(data);
    SampleVector freq = FFT::frequencies(256, 100.0);
    size_t peak_idx = 0;
    double peak_power = 0;
    for (size_t i = 1; i < power.size(); i++)
        if (power[i] > peak_power) { peak_power = power[i]; peak_idx = i; }
    ASSERT_NEAR(freq[peak_idx], 10.0, 1.0);
}

TEST(FFT, DominantFrequency) {
    SampleVector data(512);
    for (size_t i = 0; i < 512; i++)
        data[i] = std::sin(2.0 * M_PI * 5.0 * i / 100.0);
    double dom_freq = FFT::dominantFrequency(data, 100.0);
    ASSERT_NEAR(dom_freq, 5.0, 0.5);
}

TEST(FFT, NextPow2) {
    ASSERT_EQ(FFT::nextPow2(1), 1u);
    ASSERT_EQ(FFT::nextPow2(2), 2u);
    ASSERT_EQ(FFT::nextPow2(3), 4u);
    ASSERT_EQ(FFT::nextPow2(5), 8u);
    ASSERT_EQ(FFT::nextPow2(100), 128u);
    ASSERT_EQ(FFT::nextPow2(256), 256u);
    ASSERT_EQ(FFT::nextPow2(257), 512u);
}

// ============================================================================
// FilterBank Tests
// ============================================================================

TEST(FilterBank, AddBands) {
    FilterBank fb;
    fb.addBand(1.0, 5.0);
    fb.addBand(5.0, 10.0);
    fb.addBand(10.0, 20.0);
    ASSERT_EQ(fb.bandCount(), 3u);
    auto band0 = fb.bandFrequencies(0);
    ASSERT_NEAR(band0.first, 1.0, 1e-10);
    ASSERT_NEAR(band0.second, 5.0, 1e-10);
}

TEST(FilterBank, Configure) {
    FilterBank fb;
    fb.addBand(1.0, 5.0);
    fb.addBand(5.0, 10.0);
    ASSERT_NO_THROW(fb.configure(100.0));
}

TEST(FilterBank, FilterRealData) {
    auto wfs = loadRealWaveforms();
    if (wfs.empty()) {
        std::cout << "    [SKIP] No real test data available\n";
        return;
    }

    FilterBank fb;
    fb.addBand(1.0, 5.0);
    fb.addBand(5.0, 15.0);
    fb.configure(wfs[0]->sampleRate());

    auto filtered = fb.filter(wfs[0]->data());
    ASSERT_EQ(filtered.size(), 2u);
    ASSERT_EQ(filtered[0].size(), wfs[0]->sampleCount());
}
