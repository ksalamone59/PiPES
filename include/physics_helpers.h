/**
    * @brief Custom made physics helper functions, including boost kinematics, parametrizations of Form Factors. Written in natural units unless otherwise specified
    * @author Kyle Salamone 
    * @date July 2026
*/

#ifndef PHYSICS_HELPERS_H
#define PHYSICS_HELPERS_H

#include <iostream>
#include <cmath>
#include <numbers>
#include <array>
#include <algorithm>
#include <stdexcept>

enum class charge{minus = -1, plus = 1};

struct threeVector
{
    double x{0.}, y{0.}, z{0.};
    /**
     * @brief Returns the magnitude squared of the three-vector
     * @returns The magnitude squared of the three-vector
     */
    double mag2() const noexcept
    {
        return x*x + y*y + z*z;
    }
    /**
     * @brief Normalizes the three-vector
     * @returns The normalized three vector
     */
    threeVector normalize() const 
    {
        double mag = std::sqrt(mag2());
        if(mag == 0.0)
        {
            throw std::runtime_error("Cannot normalize zero vector");
        }
        return {x/mag, y/mag, z/mag};
    }
    /**
     * @brief Dot product operator of two three vectors
     * @returns The dot product of the two three vectors
     */
    double operator*(const threeVector &other) const noexcept
    {
        return x * other.x + y * other.y + z * other.z;
    }
};

struct fourVector
{
    double E{0.}, px{0.}, py{0.}, pz{0.};
    fourVector() = default;
    /**
     * @brief Constructs a four-vector with the given components
     */
    fourVector(double E_, double px_, double py_, double pz_) : E(E_), px(px_), py(py_), pz(pz_) {}
    double p() const noexcept 
    {
        return std::hypot(px, py, pz);
    }
};

namespace physics_helpers
{
    constexpr double m_pion = 139.57; // MeV/c^2
    constexpr double m_pion_squared = m_pion * m_pion;
    constexpr double m_proton = 938.27; // MeV/c^2
    constexpr double m_proton_squared = m_proton * m_proton;
    constexpr double magnetic_moment_proton =  2.792847344; // Nuclear magnetons
    constexpr double reduced_mass = (m_pion * m_proton) / (m_pion + m_proton); // MeV/c^2
    constexpr double hbarc = 197.327; // MeV*fm
    constexpr double alpha_em = 1.0 / 137.036;
    const double e_charge = std::sqrt(4*std::numbers::pi*alpha_em);

    inline constexpr double deg2rad(double x) {return x * std::numbers::pi / 180.;}
    inline constexpr double rad2deg(double x) {return x * 180. / std::numbers::pi;}
    /**
     * @brief Converts lab frame momentum to CM frame momentum
     * @param s: Mandelstam s variable
     * @returns The CM frame momentum
     */
    inline double p_lab_to_cm(const double s)
    {
        return std::sqrt((s - std::pow(m_proton + m_pion, 2)) * (s - std::pow(m_proton - m_pion, 2)) ) / (2.0 * std::sqrt(s));
    }
    /**
     * @brief Calculates the lab frame energy of the pion
     * @param momentum_lab: Lab frame momentum
     * @returns The lab frame energy of the pion
     */
    inline double pion_lab_energy(const double momentum_lab)
    {
        return std::sqrt(momentum_lab*momentum_lab + m_pion*m_pion);
    }
    /**
     * @brief Calculates the Mandelstam s variable
     * @param energy_lab: Lab frame energy
     * @returns The Mandelstam s variable
     */
    inline double s(const double energy_lab)
    {
        return m_proton*m_proton + m_pion*m_pion + 2.0*m_proton*energy_lab;
    }
    /**
     * @brief Calculates the Mandelstam t variable
     * @param k_cm: CM momentum
     * @param theta_cm: Scattering angle in CM frame, in radians
     * @returns The Mandelstam t variable
     */
    inline double t(const double k_cm, const double theta_cm)
    {
        return -2.0 * k_cm * k_cm * (1.0 - std::cos(theta_cm));
    }
    /**
     * @brief Calculates the Lorentz factor
     * @param beta: Velocity divided by the speed of light
     * @returns The Lorentz factor
     */
    inline double gamma(const double beta)
    {
        return 1.0 / std::sqrt(1.0 - beta*beta);
    }
    /**
     * @brief Calculates the velocity of the CM frame in the lab frame
     * @param momentum_lab: Lab frame momentum
     * @param energy_lab: Lab frame energy
     * @returns The velocity of the CM frame in the lab frame
     */
    inline double beta_boost(const double momentum_lab, const double energy_lab) // Speed of CM frame in lab
    {
        return momentum_lab / (energy_lab + m_proton);
    }
    /**
     * @brief Calculates the energy of the pion in the CM frame
     * @param k_cm: CM momentum
     * @returns The energy of the pion in the CM frame
     */
    inline double E_cm_pi(const double k_cm)
    {
        return std::sqrt(k_cm*k_cm + m_pion*m_pion);
    }
    /**
     * @brief Calculates the energy of the proton in the CM frame
     * @param k_cm: CM momentum
     * @returns The energy of the proton in the CM frame
     */
    inline double E_cm_p(const double k_cm)
    {
        return std::sqrt(k_cm*k_cm + m_proton*m_proton);
    }
    /**
     * @brief Calculates the velocity of the pion in the CM frame
     * @param momentum_cm: CM momentum
     * @returns The velocity of the pion in the CM frame
     */
    inline double beta_pion_cm(const double momentum_cm)
    {
        return momentum_cm / std::sqrt(momentum_cm*momentum_cm + m_pion*m_pion);
    }
    /**
     * @brief Calculates the velocity of the proton in the CM frame
     * @param momentum_cm: CM momentum
     * @returns The velocity of the proton in the CM frame
     */
    inline double beta_proton_cm(const double momentum_cm)
    {
        return momentum_cm / std::sqrt(momentum_cm*momentum_cm + m_proton*m_proton);
    }
    /**
     * @brief Calculates the kinematic parameter alpha
     * @param beta_cm: Velocity of the CM frame divided by the speed of light
     * @param beta_pion: Velocity of the pion divided by the speed of light
     * @returns The kinematic parameter alpha
     */
    inline double alpha_kin(const double beta_cm, const double beta_pion)
    {
        return beta_cm / beta_pion;
    }
    /**
     * @brief Converts lab frame scattering angle to CM frame angle
     * @param lab_theta: Scattering angle in lab frame, in degrees
     * @param gamma_boost: Relativistic factor
     * @param alpha_kinematic: Kinematic parameter
     * @returns Scattering angle in CM frame, in degrees
     */
    inline double theta_lab_to_cm(double lab_theta, const double gamma_boost, const double alpha_kinematic)
    {
        const double lab_rad = deg2rad(lab_theta);
        const double sin_lab = std::sin(lab_rad);
        const double cos_lab = std::cos(lab_rad);
        double theta_cm = lab_rad * 1.2; // initial guess
        for(int i = 0; i < 50; i++)
        {
            const double sin_cm = std::sin(theta_cm);
            const double cos_cm = std::cos(theta_cm);
            const double f  =  sin_cm * cos_lab - sin_lab * gamma_boost * (cos_cm + alpha_kinematic);
            const double df =  cos_cm * cos_lab + sin_lab * gamma_boost * sin_cm;
            const double step = f / df;
            theta_cm -= step;
            if(std::fabs(step) < 1e-12) break;
        }
        return rad2deg(theta_cm);
    }
    /**
     * @brief Converts CM frame scattering angle to lab frame angle
     * @param cm_theta: Scattering angle in CM frame, in degrees
     * @param gamma_boost: Relativistic factor
     * @param alpha_kinematic: Kinematic parameter
     * @returns Scattering angle in lab frame, in degrees
     */
    inline double theta_cm_to_lab(const double cm_theta, const double gamma_boost, const double alpha_kinematic)
    {
        const double cm_rad = deg2rad(cm_theta);
        const double numerator = std::sin(cm_rad);
        const double denominator = gamma_boost * (std::cos(cm_rad) + alpha_kinematic);
        double theta_lab = std::atan2(numerator, denominator);
        if(theta_lab < 0) theta_lab += std::numbers::pi;
        return rad2deg(theta_lab);
    }
    /**
     * @brief Calculates the Jacobian for the transformation from lab to CM frame
     * @param lab_theta: Scattering angle in lab frame, in degrees
     * @param gamma_boost: Relativistic factor
     * @param alpha_kinematic: Kinematic parameter
     * @returns The Jacobian dOmega_cm/dOmega_lab
     */
    inline double jacobian(const double lab_theta, const double gamma_boost, const double alpha_kinematic) 
    {
        const double eps = 1e-6;
        double t1 = theta_lab_to_cm(lab_theta + eps, gamma_boost, alpha_kinematic);
        double t0 = theta_lab_to_cm(lab_theta - eps, gamma_boost, alpha_kinematic);

        double dcm_dlab = (deg2rad(t1) - deg2rad(t0)) / (2 * eps * std::numbers::pi / 180.0);
        double th_cm = deg2rad(theta_lab_to_cm(lab_theta, gamma_boost, alpha_kinematic));
        double lab = deg2rad(lab_theta);

        return (std::sin(th_cm) / std::sin(lab)) * dcm_dlab;
    }
    /**
     * @brief Calculates the associated Legendre polynomial
     * @param l: Degree of the polynomial
     * @param theta: Angle in radians
     * @returns The value of the associated Legendre polynomial
     */
    inline double associated_legendre(const int l, const double theta)
    {
        const double cos_theta = std::cos(theta);
        const double sin_theta = std::sin(theta);
        double ret = 0.0;
        switch(l)
        {
            case 0:
                return 0.;
            case 1:
                ret = 1.;
                break;
            case 2:
                ret = 3. * cos_theta;
                break;
            case 3:
                ret = 0.5 * (15. * cos_theta * cos_theta - 3.);
                break;
            default:
                throw std::invalid_argument("Invalid associated Legendre index " + std::to_string(l));
        }
        return -1. * sin_theta * ret;
    }
    /**
     * @brief Proton Dirac Form Factor
     * @param t_ Mandelstam t (internally converted to t/m_pi^2)
     * @returns The value of the Proton Dirac Form Factor
     */
    inline double F1P(double t_)
    {
        t_ /= m_pion_squared;
        return (0.47737 + 0.04497/std::pow(1.-t_/18.2, 2.))/(1.-t_/27.523)
            + 0.147/(1.-t_/91.) - 0.16934/(1.-t_/130.) + 0.92/(1.-t_/31.468)
            - 0.352/(1.-t_/70.) - 0.068/(1.-t_/140.);
    }
    /**
     * @brief Proton Pauli form factor 
     * @param t_ Mandelstam t (in units of m_pi^2 = 1)
     * @returns The value of the Proton Pauli form factor
     */
    inline double F2P(double t_)
    {
        t_ /= m_pion_squared;
        return (2.6676 + 0.48086/(1.-t_/13.758))/(1.-t_/30.943) - 2.2/(1.-t_/91.)
            + 0.90454/(1.-t_/130.) - 0.072/(1.-t_/31.468)
            - 0.008/(1.-t_/70.) + 0.02/(1.-t_/140.);
    }
    /**
     * @brief Pion Form Factor 
     * @param t_ Mandelstam t (in units of m_pi^2 = 1)
     * @returns The value of the Pion Form Factor
     */
    inline double F_Pi(double t_)
    {
        t_ /= m_pion_squared;
        return 1.0538/(1-(t_/31.093)) - 0.0538/(1-(t_/76.8)); 
    }
    /**
     * @brief Simpson's rule for numerical integration
     * @param lower Lower limit of integration
     * @param upper Upper limit of integration
     * @param f_lower Function value at lower limit
     * @param f_middle Function value at middle point
     * @param f_upper Function value at upper limit
     * @return Integrated value
     */
    inline double simpson_integrate(const double lower, const double upper, const double f_lower, const double f_middle, const double f_upper)
    {
        return ((upper - lower) / 6.0) * (f_lower + 4.0 * f_middle + f_upper);
    }
    /**
     * @brief Adaptive Simpson's rule for numerical integration
     * @param lower Lower limit of integration
     * @param upper Upper limit of integration
     * @param f Function to integrate
     * @param tol Tolerance for convergence
     * @param N_ITER Maximum number of iterations
     * @return Integrated value
     */
    template<typename Func>
    inline double adaptive_simpson(const double lower, const double upper, Func&& f, double tol, int N_ITER = 100)
    {
        const double mid = 0.5 * (lower + upper);
        const double f_low = f(lower), f_mid = f(mid), f_upper = f(upper);
        const double whole_integral = simpson_integrate(lower, upper, f_low, f_mid, f_upper);
        const double left_middle = 0.5 * (lower + mid);
        const double right_middle = 0.5 * (mid + upper);
        const double f_left_mid = f(left_middle), f_right_mid = f(right_middle);
        const double left_integral = simpson_integrate(lower, mid, f_low, f_left_mid, f_mid);
        const double right_integral = simpson_integrate(mid, upper, f_mid, f_right_mid, f_upper);
        if(std::fabs(left_integral + right_integral - whole_integral) < 15.*tol || --N_ITER <= 0) 
        {
            return left_integral + right_integral + (left_integral + right_integral - whole_integral) / 15.0;
        }
        return adaptive_simpson(lower, mid, f, tol/2., N_ITER) + adaptive_simpson(mid, upper, f, tol/2., N_ITER);
    }
    template<std::size_t N_COLS>
    using col_array = std::array<double, N_COLS - 1>;
    /**
     * @brief Interpolates values in a table, written specificall for tables with row values in q/mu 
     * @param x_extract: The value to interpolate at
     * @param table: The table of values to interpolate from
     * @returns The interpolated values
     */
    template<std::size_t N_ROWS, std::size_t N_COLS>
    inline col_array<N_COLS> interpolate(const double x_extract, const std::array<std::array<double,N_COLS>, N_ROWS> &table)
    {
        col_array<N_COLS> res{0.};
        try
        {
            static_assert(N_COLS > 1, "Table must have at least two columns");
            auto extract = [](const std::array<double, N_COLS> &row) -> col_array<N_COLS>
            {
                col_array<N_COLS> out;
                for(std::size_t i=0;i<N_COLS-1;i++)
                {
                    out[i] = row[i+1];
                }
                return out;
            };
            if(x_extract <= table.front()[0])
            {
                if(x_extract < table.front()[0])
                {
                    throw std::invalid_argument("x_extract is out of bounds: too small");
                }
                return extract(table.front());
            }
            if(x_extract >= table.back()[0])
            {
                if(x_extract > table.back()[0])
                {
                    throw std::invalid_argument("x_extract is out of bounds: too large");
                }
                return extract(table.back());
            }
            auto it = std::upper_bound(table.begin(), table.end(), x_extract, 
                [](double val, const std::array<double,N_COLS> &row)
                {
                    return val < row[0];
                });
            const std::size_t idx1 = static_cast<std::size_t>(std::distance(table.begin(), it));
            const std::size_t idx0 = idx1 - 1;

            const double x0 = table[idx0][0];
            const double x1 = table[idx1][0];
            double frac = (x_extract - x0) / (x1 - x0);

            for(std::size_t i=0;i<N_COLS-1;i++)
            {
                res[i] = table[idx0][i + 1] + frac * (table[idx1][i + 1] - table[idx0][i + 1]);
            }
        }
        catch(const std::exception &e)
        {
            std::cerr << "Error in interpolate: " << e.what() << std::endl;
        }
        return res;
    }
    /**
     * @brief Boosts a given four vector from CM to lab frame 
     * @param vec: The four vector to boost
     * @param beta: Three vector representing the boost
     * @returns The boosted four vector in the lab frame
     */
    fourVector boost_cm_to_lab(const fourVector &vec, const threeVector &beta)
    {
        double beta2 = beta.mag2();
        if(beta2 == 0.) return vec;
        double gamma = 1.0 / std::sqrt(1.0 - beta2);
        double bp = beta.x * vec.px + beta.y * vec.py + beta.z * vec.pz;
        double gamma2 = ((gamma - 1.0)/beta2)*bp + gamma*vec.E;
        return {
            gamma * (vec.E + bp),
            vec.px + beta.x * gamma2,
            vec.py + beta.y * gamma2,
            vec.pz + beta.z * gamma2
        };
    }
    /**
     * @brief Rotates a given four vector around a specified axis
     * @param fVec: The four vector to rotate, passsed as reference
     * @param axis: The axis to rotate around
     */
    void rotate_four_vector(fourVector &fVec, const threeVector &axis)
    {
        static const threeVector zAxis{0., 0., 1.};
        threeVector target = (axis.mag2() == 1) ? axis : axis.normalize();
        double cos_theta = zAxis * target;
        if(std::fabs(cos_theta - 1.0) < 1e-6) return;
        else if(std::fabs(cos_theta + 1.0) < 1.e-6)
        {
            fVec.px = -fVec.px;
            fVec.py = -fVec.py;
            fVec.pz = -fVec.pz;
            return;
        }
        threeVector rotation_axis = {-target.y, target.x, 0.};
        rotation_axis = rotation_axis.normalize();
        double sin_theta = std::sqrt(1.0 - cos_theta * cos_theta);
        double ux = rotation_axis.x;
        double uy = rotation_axis.y;    
        double uz = rotation_axis.z;
        double px = fVec.px;
        double py = fVec.py;
        double pz = fVec.pz;
        double dot = rotation_axis * threeVector{px, py, pz};
        fVec.px = px * cos_theta + (uy * pz - uz * py) * sin_theta + ux * dot * (1.0 - cos_theta);
        fVec.py = py * cos_theta + (uz * px - ux * pz) * sin_theta + uy * dot * (1.0 - cos_theta);
        fVec.pz = pz * cos_theta + (ux * py - uy * px) * sin_theta + uz * dot * (1.0 - cos_theta);
    }
    /**
     * @brief Returns spherical angles of a vector relative to an arbitrary axis.
     *
     * The polar angle theta is measured from axis.
     * phi=0 is defined by the +x direction projected perpendicular to the axis.
     * phi increases clockwise when looking along the axis direction.
     *
     * @param vec Momentum vector
     * @param axis Reference axis (beam direction)
     * @returns pair(theta,phi) in radians
     */
    inline std::pair<double,double> spherical_angles_relative_to_axis(const threeVector &vec, const threeVector &axis)
    {
        threeVector z = axis.normalize();
        threeVector x_ref{1.,0.,0.};
        if(std::fabs(z*x_ref) > 0.999)
        {
            x_ref = {0.,1.,0.};
        }
        x_ref = {x_ref.x - (x_ref*z)*z.x, x_ref.y - (x_ref*z)*z.y, x_ref.z - (x_ref*z)*z.z};
        x_ref = x_ref.normalize();
        threeVector y_ref{z.y*x_ref.z - z.z*x_ref.y, z.z*x_ref.x - z.x*x_ref.z, z.x*x_ref.y - z.y*x_ref.x };
        threeVector p = vec.normalize();
        double theta = std::acos(std::clamp(p*z,-1.0,1.0));
        double phi = std::atan2(p*y_ref, p*x_ref);
        if(phi < 0) phi += 2.0*std::numbers::pi;
        return {theta,phi};
    }
}

#endif 