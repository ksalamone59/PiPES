    #ifndef PI_P_AMPLITUDE_H
    #define PI_P_AMPLITUDE_H

    #include "physics_helpers.h"
    #include "partial_wave.h"
    #include "PhaseShiftLoader.h"

    #include <array>
    #include <complex>
    #include <limits>
    #include <memory>
    #include <optional>
    #include <vector>

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
    // Tromborg Table II: pi+p, {q/mu, Delta(S1/2), Delta(P1/2), Delta(P3/2)}, degrees
    inline static constexpr std::array<std::array<double,4>, 18> delta_pip_deg = {{
        {0.5,  0.10, 0.01, -0.07}, 
        {0.8,  0.09, 0.02, -0.23}, 
        {1.0,  0.10, 0.04, -0.47},
        {1.1,  0.10, 0.04, -0.65}, 
        {1.2,  0.10, 0.05, -0.90}, 
        {1.3,  0.11, 0.06, -1.19},
        {1.35, 0.11, 0.07, -1.32}, 
        {1.4,  0.11, 0.07, -1.42}, 
        {1.45, 0.12, 0.07, -1.46},
        {1.5,  0.12, 0.08, -1.42}, 
        {1.6,  0.12, 0.09, -1.13}, 
        {1.7,  0.13, 0.09, -0.72},
        {1.8,  0.13, 0.10, -0.36}, 
        {1.9,  0.13, 0.11, -0.10}, 
        {2.0,  0.13, 0.12,  0.07},
        {2.2,  0.13, 0.13,  0.25}, 
        {2.5,  0.14, 0.15,  0.32}, 
        {3.0,  0.14, 0.20,  0.27}
    }};

    // Tromborg Table III: pi-p, {q/mu, Delta1S1/2, Delta3S1/2, Delta13S1/3, Delta1P3/2, Delta3P3/2, Delta13P3/2}; degrees
    // P1/2 is negligible
    inline static constexpr std::array<std::array<double,7>, 18> delta_pim_deg = {{
        {0.5, -0.31, 0.52, 0.06,  0.02, -0.74,  0.11}, 
        {0.8, -0.21, 0.35, 0.06,  0.05, -1.48,  0.15},
        {1.0, -0.15, 0.30, 0.07,  0.08, -2.13,  0.15}, 
        {1.1, -0.13, 0.28, 0.07,  0.10, -2.51,  0.14},
        {1.2, -0.10, 0.27, 0.08,  0.11, -2.87,  0.12}, 
        {1.3, -0.07, 0.26, 0.09,  0.13, -3.09,  0.09},
        {1.35,-0.06, 0.25, 0.09,  0.14, -3.03,  0.07}, 
        {1.4, -0.05, 0.25, 0.10,  0.14, -2.82,  0.05},
        {1.45,-0.03, 0.25, 0.10,  0.15, -2.44,  0.01}, 
        {1.5, -0.02, 0.24, 0.11,  0.16, -1.93, -0.02},
        {1.6,  0.00, 0.24, 0.12,  0.17, -0.73, -0.08}, 
        {1.7,  0.01, 0.23, 0.14,  0.19,  0.26, -0.13},
        {1.8,  0.03, 0.22, 0.15,  0.20,  0.87, -0.16}, 
        {1.9,  0.04, 0.22, 0.16,  0.22,  1.15, -0.17},
        {2.0,  0.05, 0.21, 0.17,  0.23,  1.23, -0.18}, 
        {2.2,  0.09, 0.20, 0.20,  0.27,  1.14, -0.17},
        {2.5,  0.14, 0.18, 0.24,  0.36,  0.83, -0.16}, 
        {3.0,  0.20, 0.13, 0.35,  0.47,  0.37, -0.13}
    }};
    // Tromborg Table IV: pi-p inelasticity corrections.
    // {q/mu, etabar1(S1/2), etabar3(S1/2), eta13(S1/2), etabar1(P3/2), etabar3(P3/2), eta13(P3/2)}
    inline static constexpr std::array<std::array<double,7>, 13> eta_pim_e4 = {{
        {0.8, 17e-4,  6e-4, 11e-4, 1e-4,   4e-4,  2e-4},
        {1.0, 22e-4,  8e-4, 14e-4, 2e-4,  14e-4,  4e-4},
        {1.2, 26e-4,  9e-4, 16e-4, 3e-4,  49e-4, 10e-4},
        {1.3, 27e-4,  9e-4, 17e-4, 3e-4,  69e-4, 14e-4},
        {1.4, 28e-4, 10e-4, 18e-4, 4e-4,  89e-4, 16e-4},
        {1.5, 29e-4, 10e-4, 18e-4, 4e-4, 109e-4, 19e-4},
        {1.6, 30e-4, 10e-4, 19e-4, 5e-4, 113e-4, 20e-4},
        {1.7, 31e-4, 11e-4, 19e-4, 5e-4, 105e-4, 20e-4},
        {1.8, 31e-4, 11e-4, 20e-4, 6e-4,  91e-4, 20e-4},
        {2.0, 32e-4, 12e-4, 21e-4, 7e-4,  61e-4, 18e-4},
        {2.2, 32e-4, 12e-4, 21e-4, 7e-4,  39e-4, 16e-4},
        {2.5, 31e-4, 14e-4, 22e-4, 9e-4,  20e-4, 13e-4},
        {3.0, 26e-4, 21e-4, 25e-4, 9e-4,   6e-4,  8e-4}
    }};

    class PiPAmplitude
    {
        private:
            partial_wave phase_shifts, hadronic_phase_shifts;
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
            double Delta13_S{0.0}, eta13_S{0.0};       
            double Delta13_P32{0.0}, eta13_P32{0.0};   
            double eta_S11{1.0}, eta_S31{1.0}, eta_P13{1.0}, eta_P33{1.0}; 
            charge charge_polarity;
            std::vector<std::array<std::optional<double>,2>> coulomb_phase_cache;
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
            inline cdouble coulomb_spin_flip_amplitude(const double theta_cm)
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
             * @brief Computes the isospin-mixing amplitude f_13 (Tromborg Eq. 10c)
             * @param delta_13 Mixing phase-shift correction Delta_13 (radians)
             * @param eta_13 Mixing inelasticity correction eta_13 (dimensionless, already x1e-4 applied)
             * @param delta_1 Isospin-1/2 nuclear phase shift delta^1 (Eq. 11a), radians
             * @param delta_3 Isospin-3/2 nuclear phase shift delta^3 (Eq. 11b), radians
             * @returns f_13
             */
            inline cdouble f13(const int l, const double delta_13, const double eta_13, const double delta_1, const double delta_3)
            {
                if(charge_polarity == charge::plus || l > 1) return cdouble(0.0, 0.0);
                return (2./3.) * std::numbers::sqrt2 * (eta_13 + I * delta_13) * std::exp(I * (delta_1 + delta_3));
            }
            /**
             * @brief Calculates the Coulomb phase shift integral \Sigma_{l\pm} (Eq. 12)
             * @param l The orbital angular momentum quantum number (l >= 0)
             * @param plus_minus true for l+ (J = l + 1/2), false for l- (J = l - 1/2)
             * @returns The Coulomb phase shift
             */
            inline double coulomb_phase_integral(const int l, const bool plus_minus)
            {
                if (l > L_MAX) throw std::out_of_range("Requested Coulomb phase above L_MAX");
                auto &cached_value = coulomb_phase_cache[l][plus_minus];
                if(cached_value) return *cached_value;
                std::cout << "Computing Coulomb phase integral for l = " << l << ", plus_minus = " << plus_minus << std::endl;
                using namespace physics_helpers;
                int sgn = plus_minus ? 1 : -1;
                if (l == 0 && !plus_minus) 
                {
                    throw std::domain_error("l- state (J = -1/2) does not exist for l = 0.");
                }
                const double J = l + sgn * 0.5;
                const int l_prime = l + sgn;
                auto first_integrand = [this, l](double z) -> double 
                {
                    z = std::clamp(z, -1.0, 1.0);
                    double t_z = -2.0 * cm_momentum * cm_momentum * (1.0 - z);
                    double term1 = 0.0;
                    if (std::abs(1.0 - z) < 1e-10) 
                    {
                        term1 = (l == 0) ? 0.0 : (gamma_coulomb / (2.0 * cm_momentum)) * l * (l + 1) * F1P(0.0);
                    } 
                    else
                    {
                        term1 = (2.0 * cm_momentum * gamma_coulomb / t_z) * F1P(t_z) * (legendre(l, z) - 1.0);
                    }
                    double term2 = (alpha_em / (2.0 * W)) * F2P(t_z) * legendre(l, z);
                    double ret = F_Pi(t_z) * (term1 - term2);
                    return ret;
                };
                auto second_integrand = [this, l, l_prime](double z) -> double 
                {
                    z = std::clamp(z, -1.0, 1.0);
                    double t_z = -2.0 * cm_momentum * cm_momentum * (1.0 - z);
                    double f1_term = (W + m_proton) / (energy_cm_proton + m_proton) * F1P(t_z);
                    double f2_term = (W / m_proton + t_z / (4.0 * m_proton * (energy_cm_proton + m_proton))) * F2P(t_z);
                    double P_prime_sum = legendre_derivative(l, z) + legendre_derivative(l_prime, z);
                    return F_Pi(t_z) * (f1_term + f2_term) * P_prime_sum;
                };
                double first_integral = 0.5 * cm_momentum * adaptive_simpson(-1.0, 1.0, first_integrand, 1e-13, 100);
                double second_integral = sgn * ((alpha_em * cm_momentum) / (4.0 * W * (J + 0.5))) 
                                        * adaptive_simpson(-1.0, 1.0, second_integrand, 1e-13, 100);
                cached_value = first_integral + second_integral;
                if(charge_polarity == charge::minus) *cached_value *= -1.0;
                return cached_value.value();
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
                coulomb_phase_cache.resize(L_MAX + 1);
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
                    std::cout << "Coulomb Phase Shifts: " << std::endl;
                    std::cout << "  l=0, +: " << coulomb_phase_cache[0][1].value() << " deg" << std::endl;
                    std::cout << "  l=1, -: " << coulomb_phase_cache[1][0].value() << " deg" << std::endl;
                    std::cout << "  l=1, +: " << coulomb_phase_cache[1][1].value() << " deg" << std::endl;
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
                const double sinth = std::sin(theta_rad);
                cdouble a0_plus(0.0, 0.0);
                static const double one_third = 1.0/3.0, two_thirds = 2.0/3.0;
                if(charge_polarity == charge::plus)
                {
                    a0_plus = T(phase_shifts[static_cast<std::size_t>(wave::S31)]);
                }
                else if(charge_polarity == charge::minus)
                {
                    cdouble T_32 = T(phase_shifts[static_cast<std::size_t>(wave::S31)], eta_S31);
                    cdouble T_12 = T(phase_shifts[static_cast<std::size_t>(wave::S11)], eta_S11);
                    cdouble f_13_S = f13(0, Delta13_S, eta13_S, phase_shifts[static_cast<std::size_t>(wave::S11)], phase_shifts[static_cast<std::size_t>(wave::S31)]);
                    a0_plus = one_third * T_32 + two_thirds * T_12 - two_thirds * std::numbers::sqrt2 * f_13_S / (2. * I * cm_momentum_fm);
                }
                a0_plus *= std::exp(2.0 * I * coulomb_phase_integral(0, true));
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
                        cdouble T_32_plus = T(phase_shifts[idx_32_plus], eta_P33);
                        cdouble T_12_minus = T(phase_shifts[idx_12_minus]);
                        cdouble T_12_plus = T(phase_shifts[idx_12_plus], eta_P13);
                        a_minus = one_third * T_32_minus + two_thirds * T_12_minus; // No mixing here
                        cdouble mixing_term = f13(l, Delta13_P32, eta13_P32, phase_shifts[idx_12_plus], phase_shifts[idx_32_plus]);
                        a_plus = one_third * T_32_plus + two_thirds * T_12_plus - two_thirds * std::numbers::sqrt2 * mixing_term / (2. * I * cm_momentum_fm);
                    }
                    cdouble phase_c_plus = std::exp(2.0 * I * coulomb_phase_integral(l, true));
                    cdouble phase_c_minus = std::exp(2.0 * I * coulomb_phase_integral(l, false));
                    a_plus *= phase_c_plus;
                    a_minus *= phase_c_minus;
                    spin_no_flip_amplitude += ((l + 1.0) * a_plus + static_cast<double>(l) * a_minus) * std::legendre(l, costh);
                    spin_flip_amplitude += (a_plus - a_minus) * physics_helpers::legendre_derivative(l, costh) * sinth;
                }
                cdouble spin_no_flip_coulomb_amplitude = coulomb_no_spin_flip_amplitude(theta_cm);
                cdouble spin_flip_coulomb_amplitude = coulomb_spin_flip_amplitude(theta_cm);
                // Couomb amplitudes are defined for pi+p, Tromborg: this gives pi-p
                if(charge_polarity == charge::minus)
                {
                    spin_no_flip_coulomb_amplitude = -std::conj(spin_no_flip_coulomb_amplitude);
                    spin_flip_coulomb_amplitude = -spin_flip_coulomb_amplitude;
                }
                double dcs = std::norm(spin_no_flip_amplitude + spin_no_flip_coulomb_amplitude) + std::norm(spin_flip_amplitude + spin_flip_coulomb_amplitude);
                return dcs;
            }
            /**
             * @brief: Update PWA parameters for a given lab momentum 
             * @param momentum_lab The lab momentum in MeV
             */
            inline void update_parameters_for_new_momentum(const double momentum_lab)
            {
                this->momentum_lab = momentum_lab;
                // SAID gives hadronic phase shifts only
                hadronic_phase_shifts = phase_shift_loader->set_s_p_phase_shifts(momentum_lab, L_MAX, verbose);
                energy_lab = physics_helpers::pion_lab_energy(momentum_lab);
                s = physics_helpers::s(energy_lab);
                cm_momentum = physics_helpers::p_lab_to_cm(s);
                W = std::sqrt(s);
                q_over_mu = cm_momentum / physics_helpers::m_pion;
                cm_momentum_fm = cm_momentum / physics_helpers::hbarc;
                energy_cm_proton = physics_helpers::E_cm_p(cm_momentum);
                energy_cm = physics_helpers::E_cm_pi(cm_momentum);
                beta_cm_proton = physics_helpers::beta_proton_cm(cm_momentum);
                beta_cm_pion = physics_helpers::beta_pion_cm(cm_momentum);
                beta_lab_cm_boost = physics_helpers::beta_boost(momentum_lab, energy_lab);
                alpha_kinematic = physics_helpers::alpha_kin(beta_lab_cm_boost, beta_cm_pion);
                gamma_relativsitc = physics_helpers::gamma(beta_lab_cm_boost);
                gamma_coulomb = physics_helpers::alpha_em * (s - physics_helpers::m_proton_squared - physics_helpers::m_pion_squared) / (2.0 * cm_momentum * W);
                if(charge_polarity == charge::minus) 
                {
                    auto delta = physics_helpers::interpolate(q_over_mu, delta_pim_deg);
                    Delta13_S = physics_helpers::deg2rad(delta[2]);
                    Delta13_P32 = physics_helpers::deg2rad(delta[5]);
                    physics_helpers::col_array<7UL> etas;
                    if(q_over_mu >= eta_pim_e4.front().front()) etas = physics_helpers::interpolate(q_over_mu, eta_pim_e4);
                    else etas.fill(0.0);
                    eta13_S = etas[2];
                    eta13_P32 = etas[5];
                    eta_S11 = 1. - etas[0];
                    eta_S31 = 1. - etas[1];
                    eta_P13 = 1. - etas[3];
                    eta_P33 = 1. - etas[4];
                    std::cout << "eta_S11 = " << eta_S11 << ", eta_S31 = " << eta_S31 << ", eta_P13 = " << eta_P13 << ", eta_P33 = " << eta_P33 << std::endl;
                    std::cout << "Delta13_S = " << Delta13_S << ", eta13_S = " << eta13_S << ", Delta13_P32 = " << Delta13_P32 << ", eta13_P32 = " << eta13_P32 << std::endl;
                    phase_shifts[static_cast<std::size_t>(wave::S11)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::S11)] - (2./3.) * physics_helpers::deg2rad(delta[0]);
                    phase_shifts[static_cast<std::size_t>(wave::S31)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::S31)] - (1./3.) * physics_helpers::deg2rad(delta[1]);
                    if(L_MAX >= 1)
                    {
                        phase_shifts[static_cast<std::size_t>(wave::P33)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::P33)] - (1./3.) * physics_helpers::deg2rad(delta[4]);
                        phase_shifts[static_cast<std::size_t>(wave::P11)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::P11)];
                        phase_shifts[static_cast<std::size_t>(wave::P31)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::P31)];
                        phase_shifts[static_cast<std::size_t>(wave::P13)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::P13)] - (2./3.) * physics_helpers::deg2rad(delta[3]);
                        if(L_MAX >= 2)
                        {
                            phase_shifts[static_cast<std::size_t>(wave::D33)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::D33)];
                            phase_shifts[static_cast<std::size_t>(wave::D35)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::D35)];
                            phase_shifts[static_cast<std::size_t>(wave::D13)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::D13)];
                            phase_shifts[static_cast<std::size_t>(wave::D15)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::D15)];
                        }
                    }
                }
                else 
                {
                    auto delta = physics_helpers::interpolate(q_over_mu, delta_pip_deg);
                    phase_shifts[static_cast<std::size_t>(wave::S11)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::S11)];
                    phase_shifts[static_cast<std::size_t>(wave::S31)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::S31)] + physics_helpers::deg2rad(delta[0]);
                    if(L_MAX >= 1) 
                    {
                        phase_shifts[static_cast<std::size_t>(wave::P31)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::P31)] + physics_helpers::deg2rad(delta[1]);
                        phase_shifts[static_cast<std::size_t>(wave::P33)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::P33)] + physics_helpers::deg2rad(delta[2]);
                        phase_shifts[static_cast<std::size_t>(wave::P11)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::P11)];
                        phase_shifts[static_cast<std::size_t>(wave::P13)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::P13)];
                        if(L_MAX >= 2) 
                        {
                            phase_shifts[static_cast<std::size_t>(wave::D33)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::D33)];
                            phase_shifts[static_cast<std::size_t>(wave::D35)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::D35)];
                            phase_shifts[static_cast<std::size_t>(wave::D13)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::D13)];
                            phase_shifts[static_cast<std::size_t>(wave::D15)] = hadronic_phase_shifts[static_cast<std::size_t>(wave::D15)];
                        }
                    }
                }
                // Fill the cache 
                coulomb_phase_cache.assign(L_MAX + 1, {});
                auto sigmas = physics_helpers::interpolate(cm_momentum/physics_helpers::m_pion, coulomb_phase_shifts);
                coulomb_phase_cache[0][1] = sigmas[0];
                coulomb_phase_cache[1][0] = sigmas[1];
                coulomb_phase_cache[1][1] = sigmas[2];
            }
            inline const partial_wave &get_phase_shifts() const {return phase_shifts;}
    };

    #endif 