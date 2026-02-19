#pragma once

#include "picker.hpp"
#include "filter_bank.hpp"

namespace realdetect {

/**
 * PolarizationPicker - 3-component STA/LTA × Rectilinearity coupled picker
 *
 * Computes STA/LTA on the vertical component and rectilinearity from the
 * 3×3 covariance matrix of the 3C data.  The coupled characteristic function
 * (STA/LTA × rectilinearity) suppresses noise triggers that lack coherent
 * particle motion, yielding more reliable P-phase picks.
 *
 * Falls back to plain STA/LTA when only a single component is available.
 */
class PolarizationPicker : public BasePicker {
public:
    PolarizationPicker();

    std::vector<PickResult> pick(const Waveform& waveform) override;
    std::string name() const override { return "Polarization"; }

    void setParameter(const std::string& name, double value) override;
    double getParameter(const std::string& name) const override;

    /**
     * 3-component pick using coupled STA/LTA × rectilinearity.
     * @param z_wf  Vertical component
     * @param n_wf  North component
     * @param e_wf  East component
     */
    std::vector<PickResult> pick3C(
        const Waveform& z_wf, const Waveform& n_wf, const Waveform& e_wf);

    /** Return the coupled characteristic function (diagnostic). */
    SampleVector coupledCF(
        const Waveform& z_wf, const Waveform& n_wf, const Waveform& e_wf) const;

private:
    double sta_length_;
    double lta_length_;
    double trigger_ratio_;
    double detrigger_ratio_;
    double pol_window_;
    double rect_threshold_;
    double filter_low_;
    double filter_high_;
    bool   use_filter_;

    SampleVector filterData(const SampleVector& data, double sr) const;
    SampleVector computeSTALTA(const SampleVector& data, double sr) const;
    SampleVector computeRectilinearity(
        const SampleVector& z, const SampleVector& n,
        const SampleVector& e, double sr) const;
    size_t refinePick(const SampleVector& data, size_t approx, double sr) const;
};

} // namespace realdetect
