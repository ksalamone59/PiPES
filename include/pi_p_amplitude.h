#ifndef PI_P_AMPLITUDE_H
#define PI_P_AMPLITUDE_H

#include "physics_helpers.h"
#include "partial_wave.h"
#include "PhaseShiftLoader.h"

#include <array>
#include <complex>
#include <limits>
#include <memory>

using cdouble = std::complex<double>;
inline static constexpr cdouble I(0.0, 1.0);

// Tromborg Table I values: {q/mu, Sigma_0+, Sigma_1-, Sigma_1+}
// Coulomb phase shifts as function of q/mu
inline static constexpr std::array<std::array<double,4>, 10> coulomb_phase_shifts = {{
    {0.5,  3.08e-4,  1.32e-2,  1.43e-2},
    {1.0,  5.66e-4,  6.85e-3,  8.80e-3},
    {1.5,  7.38e-4,  4.38e-3,  6.93e-3},
    {2.0,  8.23e-4,  2.90e-3,  5.76e-3},
    {2.5,  8.41e-4,  1.90e-3,  4.87e-3},
    {3.0,  8.20e-4,  1.22e-3,  4.17e-3},
    {3.5,  7.79e-4,  7.63e-4,  3.60e-3},
    {4.0,  7.30e-4,  4.51e-4,  3.14e-3},
    {5.0,  6.33e-4,  9.45e-5,  2.46e-3},
    {6.0,  5.50e-4, -7.33e-5,  2.00e-3}
}};

class PiPAmplitude
{
    private:
        partial_wave phase_shifts;
        std::array<double,3> coulomb_phases{0.,0.,0.};
        std::shared_ptr<const PhaseShiftLoader> phase_shift_loader;
        double gamma_coulomb{0.}, cm_momentum{0.}, cm_momentum_fm{0.};
        double q_over_mu{0.}, alpha_kinematic{0.0};
        double momentum_lab{0.}; // MeV 
        double s{0.}; // Mandelstam s, MeV^2
        double gamma_relativsitc{0.}; 
        double beta_cm_pion{0.0}, beta_cm_proton{0.0}, beta_lab_cm_boost{0.0};
        double energy_lab{0.}, energy_cm{0.}, energy_cm_proton{0.}, W{0.};
        double phiC_integral_cache{std::numeric_limits<double>::quiet_NaN()};
        double phiC_integral_cache_k{-1.0};
        charge charge_polarity;
        bool verbose{false};
        int L_MAX{1};

        /**
         * @brief Returns the Coulomb amplitude 
         * @param theta_cm Center of mass scattering angle, in degrees
         * @returns The Coulomb amplitude
         */
        inline cdouble coulomb_no_spin_flip_amplitude(const double theta_cm) 
        {
            // Tromborg
            const double theta_cm_rad = physics_helpers::deg2rad(theta_cm);
            const double t = physics_helpers::t(cm_momentum, theta_cm_rad);
            const double M = physics_helpers::m_proton;
            cdouble ret = ((2 * cm_momentum * gamma_coulomb / t) + physics_helpers::alpha_em / (2. * W) * 
                ((M + W)/(energy_cm_proton + M))) * 
                physics_helpers::F1P(t);
            ret += (W - M + t / (4 * (energy_cm_proton + M))) * 
                    physics_helpers::alpha_em * physics_helpers::F2P(t)/(2.*W*physics_helpers::m_proton);
            ret *= physics_helpers::F_Pi(t) * std::exp(I * phiC(theta_cm));
            return ret * physics_helpers::hbarc; // MeV to fm
        }
        /**
         * @brief Returns the Coulomb spin flip amplitude 
         * @param theta_cm Center of mass scattering angle, in degrees
         * @returns The Coulomb spin flip amplitude
         */
        inline cdouble coulomb_spin_flip_amplitude(const double theta_cm) const
        {
            // Tromborg
            const double theta_cm_rad = physics_helpers::deg2rad(theta_cm);
            const double t = physics_helpers::t(cm_momentum, theta_cm_rad);
            const double M = physics_helpers::m_proton;
            cdouble ret = ((W + M) / (energy_cm_proton + M)) * physics_helpers::F1P(t) + 
                (W + t/(4.*(energy_cm_proton + M))) * physics_helpers::F2P(t)/physics_helpers::m_proton;
            ret *= physics_helpers::alpha_em * physics_helpers::F_Pi(t) / (2*W*std::tan(0.5*theta_cm_rad));
            return ret * physics_helpers::hbarc; // MeV to fm
        }
        /**
         * @brief Returns the integral of the Coulomb phase shift
         * @returns The integral of the Coulomb phase shift
         */
        inline double phiC_integral() 
        {
            // Cache: this is momentum-dependent but angle-independent
            if (cm_momentum == phiC_integral_cache_k) return phiC_integral_cache;
            phiC_integral_cache_k = cm_momentum;
            auto integrand = [](double t_) -> double {
                if(std::fabs(t_) < 1.e-6)
                {
                    const double h = 1e-4;
                    const double plus  = (1.0 - physics_helpers::F_Pi(h)  * physics_helpers::F1P(h))  / h;
                    const double minus = (1.0 - physics_helpers::F_Pi(-h) * physics_helpers::F1P(-h)) / (-h);
                    return 0.5 * (plus + minus); // central-difference estimate of the removable limit
                }
                return (1.0 - physics_helpers::F_Pi(t_) * physics_helpers::F1P(t_)) / t_;
            };
            phiC_integral_cache = physics_helpers::adaptive_simpson(-4.0 * cm_momentum * cm_momentum, 0.0, integrand, 1e-8);
            return phiC_integral_cache;
        }
        /**
         * @brief Returns the Coulomb phase shift for Coulomb amplitudes
         * @param theta_cm Center of mass scattering angle, in degrees 
         * @returns The Coulomb phase shift, in radians
         */
        inline double phiC(const double theta_cm) 
        {
            const double theta_cm_rad = physics_helpers::deg2rad(theta_cm);
            const double s2 = std::pow(std::sin(theta_cm_rad / 2.0), 2.0);
            return -gamma_coulomb * std::log(s2) + gamma_coulomb * phiC_integral();
        }
        /**
         * @brief Returns the cached Coulomb phase shift for hadronic amplitudes
         * @param l The angular momentum
         * @param plus_minus A boolean indicating which state we are in (+1/2 or -1/2)
         * @returns The cached Coulomb phase shift
         */
        inline double get_coulomb_phase(const int l, const bool plus_minus) const
        {
            if(l == 0) return coulomb_phases[0];
            return plus_minus ? coulomb_phases[2] : coulomb_phases[1];
        }
        /**
         * @brief Returns T matrix element divided by momentum (in fm)
         * @param delta The phase shift
         * @param eta_elas The elasticity parameter (default: 1.0)
         * @returns The T matrix element divided by momentum (in fm)
         */
        inline cdouble T(const double delta, const double eta_elas = 1.0) const
        {
            const cdouble S = eta_elas * std::exp(2.0 * I * delta);
            return ((S - 1.0) / (2.0 * I * cm_momentum_fm)); 
        }
    public: 
        PiPAmplitude() = delete;
        ~PiPAmplitude() = default;
        /**
         * @brief Constructor for the pion-proton amplitude class
         * @param phase_shift_loader_ A shared pointer to a PhaseShiftLoader object
         * @param momentum_lab_ Lab momentum
         * @param ch Charge of the pion
         * @param verbose_ Whether to print verbose output
         * @param L_MAX_ Maximum partial wave to calculate. Default: p-wave (L_MAX = 1)
         */
        PiPAmplitude(const std::shared_ptr<const PhaseShiftLoader> &phase_shift_loader_, const double momentum_lab_, const charge ch, bool verbose_ = false, const int L_MAX_ = 1) 
            : phase_shift_loader(std::move(phase_shift_loader_)), momentum_lab(momentum_lab_), charge_polarity(ch), verbose(verbose_), L_MAX(L_MAX_)
        {
            int charge_sign = static_cast<int>(charge_polarity);
            update_parameters_for_new_momentum(momentum_lab);
            if(verbose)
            {
                std::cout << "Partial Wave Analyis initialized with the following kinematic parameters:" << std::endl;
                std::cout << "Charge of pion: " << charge_sign << std::endl;
                std::cout << "Lab momentum: " << momentum_lab << "MeV\n";
                std::cout << "CM momentum: " << cm_momentum << "MeV\n";
                std::cout << "q/mu = " << cm_momentum / physics_helpers::m_pion << std::endl;
                std::cout << "Energy in lab: " << energy_lab << "MeV\n";
                std::cout << "Energy in CM: " << energy_cm << "MeV\n";
                std::cout << "s = " << s << " MeV^2. W = " << W << " MeV\n";
                std::cout << "Gamma relativistic: " << gamma_relativsitc << std::endl;
                std::cout << "Gamma Coulomb = " << gamma_coulomb << std::endl;
            }
        }
        /**
         * @brief Returns the CM boost factor for the gamma parameter
         * @returns The gamma parameter for the CM boost
         */
        inline double get_gamma_cm_boost() const {return gamma_relativsitc;}
        /**
         * @brief Returns the CM boost factor for the beta parameter
         * @returns The beta parameter for the CM boost
         */
        inline double get_beta_cm_boost() const {return beta_lab_cm_boost;} 
        /**
         * @brief Returns the kinematic alpha parameter
         * @returns The alpha parameter for the kinematic calculations
         */
        inline double get_alpha_kinematic() const {return alpha_kinematic;}
        /**
         * @brief Returns the momentum in the CM frame
         * @returns The momentum in the CM frame (in MeV)
         */
        inline double get_momentum_cm() const noexcept {return cm_momentum;}
        /**
         * @brief Calculates dsigma/dOmega in CM frame 
         * @param theta_cm: CM Scattering Angle, in degrees
         * @returns Differential cross-section in CM frame at the given theta_cm (in fm^2/sr)
         */ 
        inline double dsigma_domega_cm(const double theta_cm) 
        {   
            cdouble spin_no_flip_amplitude(0.0, 0.0), spin_flip_amplitude(0.0, 0.0);
            const double theta_rad = physics_helpers::deg2rad(theta_cm);
            const double costh = std::cos(theta_rad);
            cdouble a0_plus(0.0, 0.0);
            static const double one_third = 1.0/3.0, two_thirds = 2.0/3.0;
            if(charge_polarity == charge::plus)
            {
                a0_plus = T(phase_shifts[static_cast<std::size_t>(wave::S31)]);
            }
            else if(charge_polarity == charge::minus)
            {
                cdouble T_32 = T(phase_shifts[static_cast<std::size_t>(wave::S31)]);
                cdouble T_12 = T(phase_shifts[static_cast<std::size_t>(wave::S11)]);
                a0_plus = one_third * T_32 + two_thirds * T_12;
            }
            a0_plus *= std::exp(2.0 * I * get_coulomb_phase(0, charge_polarity == charge::plus));
            spin_no_flip_amplitude = a0_plus;
            for(int l=1;l<=L_MAX;l++)
            {
                cdouble a_plus(0.0,0.0), a_minus(0.0,0.0);
                std::size_t base = 2 + 4 * (l - 1);
                std::size_t idx_12_minus = base;
                std::size_t idx_12_plus = base + 1;
                std::size_t idx_32_minus = base + 2;
                std::size_t idx_32_plus = base + 3;
                if(charge_polarity == charge::plus)
                {
                    a_plus = T(phase_shifts[idx_32_plus]);
                    a_minus = T(phase_shifts[idx_32_minus]);
                }
                else if(charge_polarity == charge::minus)
                {
                    cdouble T_32_minus = T(phase_shifts[idx_32_minus]);
                    cdouble T_32_plus = T(phase_shifts[idx_32_plus]);
                    cdouble T_12_minus = T(phase_shifts[idx_12_minus]);
                    cdouble T_12_plus = T(phase_shifts[idx_12_plus]);
                    a_minus = one_third * T_32_minus + two_thirds * T_12_minus;
                    a_plus = one_third * T_32_plus + two_thirds * T_12_plus;
                }
                cdouble phase_c_plus = std::exp(2.0 * I * get_coulomb_phase(l, true));
                cdouble phase_c_minus = std::exp(2.0 * I * get_coulomb_phase(l, false));
                a_plus *= phase_c_plus;
                a_minus *= phase_c_minus;
                spin_no_flip_amplitude += ((l + 1.0) * a_plus + static_cast<double>(l) * a_minus) * std::legendre(l, costh);
                spin_flip_amplitude += (a_plus - a_minus) * physics_helpers::associated_legendre(l, theta_rad);
            }
            cdouble spin_no_flip_coulomb_amplitude = coulomb_no_spin_flip_amplitude(theta_cm);
            cdouble spin_flip_coulomb_amplitude = coulomb_spin_flip_amplitude(theta_cm);
            double dcs = std::norm(spin_no_flip_amplitude + spin_no_flip_coulomb_amplitude) + std::norm(spin_flip_amplitude + spin_flip_coulomb_amplitude);
            return dcs;
        }
        /**
         * @brief: Update PWA parameters for a given lab momentum 
         * @param momentum_lab The lab momentum in MeV
         */
        void update_parameters_for_new_momentum(const double momentum_lab)
        {
            this->momentum_lab = momentum_lab;
            phase_shifts = phase_shift_loader->set_s_p_phase_shifts(momentum_lab, verbose);
            energy_lab = physics_helpers::pion_lab_energy(momentum_lab);
            s = physics_helpers::s(energy_lab);
            cm_momentum = physics_helpers::p_lab_to_cm(s);
            W = std::sqrt(s);
            q_over_mu = cm_momentum / physics_helpers::m_pion;
            coulomb_phases = physics_helpers::interpolate(q_over_mu, coulomb_phase_shifts);
            cm_momentum_fm = cm_momentum / physics_helpers::hbarc;
            energy_cm_proton = physics_helpers::E_cm_p(cm_momentum);
            energy_cm = physics_helpers::E_cm_pi(cm_momentum);
            beta_cm_proton = physics_helpers::beta_proton_cm(cm_momentum);
            beta_cm_pion = physics_helpers::beta_pion_cm(cm_momentum);
            beta_lab_cm_boost = physics_helpers::beta_boost(momentum_lab, energy_lab);
            alpha_kinematic = physics_helpers::alpha_kin(beta_lab_cm_boost, beta_cm_pion);
            gamma_relativsitc = physics_helpers::gamma(beta_lab_cm_boost);
            gamma_coulomb = physics_helpers::alpha_em * (s - physics_helpers::m_proton_squared - physics_helpers::m_pion_squared) / (2.0 * cm_momentum * W);
        }
        inline const partial_wave &get_phase_shifts() const {return phase_shifts;}
};

#endif 