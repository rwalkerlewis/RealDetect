#include "realdetect/picker/polarization_picker.hpp"
#include "realdetect/picker/aic_picker.hpp"
#include <cmath>
#include <algorithm>
#include <numeric>

namespace realdetect {

PolarizationPicker::PolarizationPicker()
    : sta_length_(0.3)
    , lta_length_(8.0)
    , trigger_ratio_(1.5)
    , detrigger_ratio_(0.5)
    , pol_window_(0.5)
    , rect_threshold_(0.4)
    , filter_low_(1.0)
    , filter_high_(15.0)
    , use_filter_(true)
{
}

// ── Parameter accessors ─────────────────────────────────────────

void PolarizationPicker::setParameter(const std::string& name, double value) {
    if      (name == "sta_length")      sta_length_      = value;
    else if (name == "lta_length")      lta_length_      = value;
    else if (name == "trigger_ratio")   trigger_ratio_   = value;
    else if (name == "detrigger_ratio") detrigger_ratio_ = value;
    else if (name == "pol_window")      pol_window_      = value;
    else if (name == "rect_threshold")  rect_threshold_  = value;
    else if (name == "filter_low")      filter_low_      = value;
    else if (name == "filter_high")     filter_high_     = value;
    else if (name == "use_filter")      use_filter_      = (value != 0);
}

double PolarizationPicker::getParameter(const std::string& name) const {
    if      (name == "sta_length")      return sta_length_;
    else if (name == "lta_length")      return lta_length_;
    else if (name == "trigger_ratio")   return trigger_ratio_;
    else if (name == "detrigger_ratio") return detrigger_ratio_;
    else if (name == "pol_window")      return pol_window_;
    else if (name == "rect_threshold")  return rect_threshold_;
    else if (name == "filter_low")      return filter_low_;
    else if (name == "filter_high")     return filter_high_;
    else if (name == "use_filter")      return use_filter_ ? 1.0 : 0.0;
    return 0;
}

// ── Bandpass filter (capped at 80 % Nyquist) ────────────────────

SampleVector PolarizationPicker::filterData(const SampleVector& data, double sr) const {
    if (!use_filter_ || filter_low_ <= 0 || filter_high_ <= 0 || data.size() < 100)
        return data;

    double nyquist = sr / 2.0;
    double hi = std::min(filter_high_, nyquist * 0.80);  // cap at 80% Nyquist
    if (hi <= filter_low_) return data;                   // nothing to filter

    auto bp = IIRFilter::butterworth(2, filter_low_, hi, sr);
    return bp.filtfilt(data);
}

// ── STA/LTA on energy ───────────────────────────────────────────

SampleVector PolarizationPicker::computeSTALTA(const SampleVector& data, double sr) const {
    size_t n = data.size();
    size_t sta_n = std::max(size_t(1), static_cast<size_t>(sta_length_ * sr));
    size_t lta_n = std::max(sta_n + 1, static_cast<size_t>(lta_length_ * sr));

    SampleVector energy(n);
    for (size_t i = 0; i < n; i++) energy[i] = data[i] * data[i];

    // Cumulative sum for fast windowed means
    SampleVector cum(n + 1, 0);
    for (size_t i = 0; i < n; i++) cum[i + 1] = cum[i] + energy[i];

    SampleVector cf(n, 0);
    for (size_t i = lta_n; i < n; i++) {
        double sta = (cum[i + 1] - cum[i + 1 - sta_n]) / sta_n;
        double lta = (cum[i + 1] - cum[i + 1 - lta_n]) / lta_n;
        cf[i] = (lta > 1e-20) ? sta / lta : 0;
    }
    return cf;
}

// ── Rectilinearity from 3C covariance eigenvalues ───────────────

SampleVector PolarizationPicker::computeRectilinearity(
    const SampleVector& z, const SampleVector& n, const SampleVector& e, double sr) const
{
    size_t len = std::min({z.size(), n.size(), e.size()});
    size_t win = std::max(size_t(3), static_cast<size_t>(pol_window_ * sr));
    SampleVector rect(len, 0);
    if (len < win) return rect;

    for (size_t i = win; i < len; i++) {
        // Build 3×3 covariance matrix in the window [i-win, i)
        double mean_z = 0, mean_n = 0, mean_e = 0;
        for (size_t j = i - win; j < i; j++) {
            mean_z += z[j]; mean_n += n[j]; mean_e += e[j];
        }
        mean_z /= win; mean_n /= win; mean_e /= win;

        double c00 = 0, c01 = 0, c02 = 0;
        double c11 = 0, c12 = 0, c22 = 0;
        for (size_t j = i - win; j < i; j++) {
            double dz = z[j] - mean_z;
            double dn = n[j] - mean_n;
            double de = e[j] - mean_e;
            c00 += dz * dz; c01 += dz * dn; c02 += dz * de;
            c11 += dn * dn; c12 += dn * de;
            c22 += de * de;
        }

        // Invariants of the 3×3 symmetric matrix
        double I1 = c00 + c11 + c22;                                     // trace
        double I2 = c00*c11 + c11*c22 + c00*c22 - c01*c01 - c12*c12 - c02*c02;

        // Rectilinearity ≈ 1 − 2√I2 / I1   (exact when λ1 >> λ2,λ3)
        if (I1 > 1e-20) {
            double r = 1.0 - 2.0 * std::sqrt(std::max(0.0, I2)) / I1;
            rect[i] = std::max(0.0, std::min(1.0, r));
        }
    }
    return rect;
}

// ── AIC-based pick refinement ───────────────────────────────────

size_t PolarizationPicker::refinePick(const SampleVector& data, size_t approx, double /*sr*/) const {
    // Use AIC in a window around the approximate pick
    size_t half = 50;
    size_t lo = (approx > half) ? approx - half : 0;
    size_t hi = std::min(approx + half, data.size());
    if (hi - lo < 10) return approx;

    SampleVector window(data.begin() + lo, data.begin() + hi);
    SampleVector aic = AICPicker::computeAIC(window);
    size_t min_idx = AICPicker::findAICMinimum(aic);
    return lo + min_idx;
}

// ── Fallback: 1-component pick (delegates to plain STA/LTA) ────

std::vector<PickResult> PolarizationPicker::pick(const Waveform& waveform) {
    // Without horizontals, just run STA/LTA on the vertical
    if (waveform.sampleCount() == 0) return {};

    SampleVector data = waveform.data();
    double mean = std::accumulate(data.begin(), data.end(), 0.0) / data.size();
    for (auto& s : data) s -= mean;
    data = filterData(data, waveform.sampleRate());

    SampleVector cf = computeSTALTA(data, waveform.sampleRate());

    std::vector<PickResult> picks;
    bool triggered = false;
    size_t trig_idx = 0;
    double max_cf = 0;

    for (size_t i = 0; i < cf.size(); i++) {
        if (!triggered) {
            if (cf[i] >= trigger_ratio_) {
                triggered = true; trig_idx = i; max_cf = cf[i];
            }
        } else {
            if (cf[i] > max_cf) max_cf = cf[i];
            if (cf[i] < detrigger_ratio_) {
                size_t idx = refinePick(data, trig_idx, waveform.sampleRate());
                PickResult p;
                p.time = waveform.timeAt(idx);
                p.phase_type = PhaseType::P;
                p.sample_index = idx;
                p.snr = max_cf;
                p.confidence = std::min(1.0, max_cf / 10.0);
                // Amplitude
                size_t end = std::min(idx + static_cast<size_t>(waveform.sampleRate()), data.size());
                double amp = 0;
                for (size_t j = idx; j < end; j++) amp = std::max(amp, std::abs(data[j]));
                p.amplitude = amp;
                picks.push_back(p);
                triggered = false; max_cf = 0;
            }
        }
    }
    return picks;
}

// ── 3-component pick — the coupled STA/LTA × Rectilinearity ────

std::vector<PickResult> PolarizationPicker::pick3C(
    const Waveform& z_wf, const Waveform& n_wf, const Waveform& e_wf)
{
    if (z_wf.sampleCount() == 0) return {};
    double sr = z_wf.sampleRate();

    // 1. Demean
    auto demean = [](SampleVector v) {
        double m = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        for (auto& s : v) s -= m;
        return v;
    };
    SampleVector zd = demean(z_wf.data());
    SampleVector nd = demean(n_wf.data());
    SampleVector ed = demean(e_wf.data());

    // 2. Bandpass filter all components (capped at 80% Nyquist)
    zd = filterData(zd, sr);
    nd = filterData(nd, sr);
    ed = filterData(ed, sr);

    // 3. Compute STA/LTA on Z
    SampleVector stalta = computeSTALTA(zd, sr);

    // 4. Compute rectilinearity on 3C
    SampleVector rect = computeRectilinearity(zd, nd, ed, sr);

    // 5. Coupled CF = STA/LTA × rectilinearity
    size_t n = std::min({stalta.size(), rect.size()});
    SampleVector coupled(n, 0);
    for (size_t i = 0; i < n; i++) {
        coupled[i] = stalta[i] * rect[i];
    }

    // 6. Trigger on coupled CF
    std::vector<PickResult> picks;
    bool triggered = false;
    size_t trig_idx = 0;
    double max_coupled = 0;
    double max_rect = 0;

    for (size_t i = 0; i < n; i++) {
        if (!triggered) {
            if (coupled[i] >= trigger_ratio_ && rect[i] >= rect_threshold_) {
                triggered = true;
                trig_idx = i;
                max_coupled = coupled[i];
                max_rect = rect[i];
            }
        } else {
            if (coupled[i] > max_coupled) {
                max_coupled = coupled[i];
                max_rect = rect[i];
            }
            if (coupled[i] < detrigger_ratio_ || rect[i] < rect_threshold_ * 0.5) {
                // 7. Refine with AIC on Z
                size_t idx = refinePick(zd, trig_idx, sr);

                PickResult p;
                p.time = z_wf.timeAt(idx);
                p.phase_type = PhaseType::P;
                p.sample_index = idx;
                p.snr = max_coupled;
                p.confidence = max_rect;   // rectilinearity as confidence
                // Amplitude on Z
                size_t end = std::min(idx + static_cast<size_t>(sr), zd.size());
                double amp = 0;
                for (size_t j = idx; j < end; j++) amp = std::max(amp, std::abs(zd[j]));
                p.amplitude = amp;

                picks.push_back(p);
                triggered = false;
                max_coupled = 0;
                max_rect = 0;
            }
        }
    }
    return picks;
}

// ── Diagnostic: return the coupled CF ───────────────────────────

SampleVector PolarizationPicker::coupledCF(
    const Waveform& z_wf, const Waveform& n_wf, const Waveform& e_wf) const
{
    double sr = z_wf.sampleRate();
    auto demean = [](SampleVector v) {
        double m = std::accumulate(v.begin(), v.end(), 0.0) / v.size();
        for (auto& s : v) s -= m;
        return v;
    };
    SampleVector zd = filterData(demean(z_wf.data()), sr);
    SampleVector nd = filterData(demean(n_wf.data()), sr);
    SampleVector ed = filterData(demean(e_wf.data()), sr);

    SampleVector stalta = computeSTALTA(zd, sr);
    SampleVector rect   = computeRectilinearity(zd, nd, ed, sr);

    size_t n = std::min({stalta.size(), rect.size()});
    SampleVector coupled(n, 0);
    for (size_t i = 0; i < n; i++) coupled[i] = stalta[i] * rect[i];
    return coupled;
}

} // namespace realdetect
