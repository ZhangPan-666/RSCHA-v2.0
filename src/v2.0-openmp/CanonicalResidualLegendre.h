#pragma once
#include "RSCHACompatibility.h"

// Independent canonical-residual associated-Legendre kernel.
// No external reference evaluator, table, data file, or runtime call is used.

#include <algorithm>
#include <array>
#include <cmath>
#include <complex>
#include <limits>

namespace crli {

constexpr double pi = 3.141592653589793238462643383279502884;

enum class Region { evanescent, turning, oscillatory };
enum class Status {
    ok,
    nonfinite_input,
    outside_real_ferrers_domain,
    outside_current_airy_chart,
    airy_kernel_failure
};

struct SignedLogValue {
    double value = std::numeric_limits<double>::quiet_NaN();
    double log_abs = std::numeric_limits<double>::quiet_NaN();
    int sign = 0;
    bool value_representable = false;
};

struct Result {
    Status status = Status::nonfinite_input;
    SignedLogValue p_tilde;
    SignedLogValue q_tilde;
    SignedLogValue p_theta;
    bool p_theta_available = false;
    double theta_turning = std::numeric_limits<double>::quiet_NaN();
};

namespace detail {

constexpr long double airy_ai_0 =
    0.35502805388781723926006318600418317639797917419918L;
constexpr long double airy_aip_0 =
   -0.25881940379280679840518356018920396347909113835493L;
constexpr long double airy_bi_0 =
    0.61492662744600073515092236909361355359472818864860L;
constexpr long double airy_bip_0 =
    0.44828835735382635791482371039882839086622679921226L;

struct AiryScaled {
    // For x>0: Ai=ai_scaled*exp(-exponent), Bi=bi_scaled*exp(exponent).
    // For x<=0 exponent is zero and the fields are ordinary values.
    double ai_scaled = 0.0;
    double bi_scaled = 0.0;
    double aip_scaled = 0.0;
    double bip_scaled = 0.0;
    double exponent = 0.0;
};

inline double airy_u(int k) {
    long double product = 1.0L;
    for (int j = 2 * k + 1; j <= 6 * k - 1; j += 2) product *= j;
    long double factorial = 1.0L;
    for (int j = 2; j <= k; ++j) factorial *= j;
    long double power = 1.0L;
    for (int j = 0; j < k; ++j) power *= 216.0L;
    return static_cast<double>(product / (factorial * power));
}

struct SeriesPair { double value, derivative; };

inline SeriesPair airy_series(double x, long double a0, long double a1) {
    std::array<long double, 256> a{};
    a[0] = a0;
    a[1] = a1;
    long double x_power = 1.0L;
    long double sum = a0;
    long double derivative = a1;
    int consecutive_small = 0;
    for (int k = 1; k < static_cast<int>(a.size()); ++k) {
        if (k >= 3) a[k] = a[k - 3] / (static_cast<long double>(k) * (k - 1));
        x_power *= static_cast<long double>(x);
        const long double term = a[k] * x_power;
        sum += term;
        if (k >= 2 && x != 0.0) derivative += k * term / x;
        if (a[k] != 0.0L && k > 20
            && std::abs(term) < std::abs(sum) * 1.0e-20L) {
            ++consecutive_small;
        } else if (a[k] != 0.0L) {
            consecutive_small = 0;
        }
        if (consecutive_small >= 3) break;
    }
    return {static_cast<double>(sum), static_cast<double>(derivative)};
}

inline AiryScaled airy(double x,
                       double canonical_action = std::numeric_limits<double>::quiet_NaN()) {
    if (x > -6.0 && x < 5.0) {
        const SeriesPair ai = airy_series(x, airy_ai_0, airy_aip_0);
        const SeriesPair bi = airy_series(x, airy_bi_0, airy_bip_0);
        return {ai.value, bi.value, ai.derivative, bi.derivative, 0.0};
    }

    if (x >= 5.0) {
        const double xi = (2.0 / 3.0) * x * std::sqrt(x);
        const double prefactor = 1.0 / (std::sqrt(pi) * std::pow(x, 0.25));
        double sum_ai = 1.0;
        double sum_bi = 1.0;
        double derivative_ai = 0.0;
        double derivative_bi = 0.0;
        double last = std::numeric_limits<double>::infinity();
        for (int k = 1; k <= 15; ++k) {
            const double term = airy_u(k) / std::pow(xi, k);
            if (!(term < last)) break; // optimal truncation of the divergent series
            sum_ai += (k & 1) ? -term : term;
            sum_bi += term;
            derivative_ai += (k & 1) ? k * term / xi : -k * term / xi;
            derivative_bi += -k * term / xi;
            last = term;
        }
        const double ai_prime = 0.5 * prefactor / x
            * (-0.25 * sum_ai + 1.5 * xi * (derivative_ai - sum_ai));
        const double bi_prime = prefactor / x
            * (-0.25 * sum_bi + 1.5 * xi * (derivative_bi + sum_bi));
        const double exponent = std::isfinite(canonical_action)
            ? canonical_action : xi;
        return {0.5 * prefactor * sum_ai, prefactor * sum_bi,
                ai_prime, bi_prime, exponent};
    }

    const double X = -x;
    const double xi = (2.0 / 3.0) * X * std::sqrt(X);
    const double phase = (std::isfinite(canonical_action)
                         ? canonical_action : xi) - 0.25 * pi;
    const double prefactor = 1.0 / (std::sqrt(pi) * std::pow(X, 0.25));
    double even = 1.0;
    double odd = 0.0;
    double even_derivative = 0.0;
    double odd_derivative = 0.0;
    double last = std::numeric_limits<double>::infinity();
    for (int k = 1; k <= 15; ++k) {
        const double term = airy_u(k) / std::pow(xi, k);
        if (!(term < last)) break;
        if ((k & 1) == 0) {
            const double signed_term = ((k / 2) & 1) ? -term : term;
            even += signed_term;
            even_derivative -= k * signed_term / xi;
        } else {
            const double signed_term = (((k - 1) / 2) & 1) ? -term : term;
            odd += signed_term;
            odd_derivative -= k * signed_term / xi;
        }
        last = term;
    }
    const double cp = std::cos(phase), sp = std::sin(phase);
    const double fai = cp * even + sp * odd;
    const double fbi = -sp * even + cp * odd;
    const double fai_xi = -sp * even + cp * odd
                        + cp * even_derivative + sp * odd_derivative;
    const double fbi_xi = -cp * even - sp * odd
                        - sp * even_derivative + cp * odd_derivative;
    const double derivative_prefactor = prefactor / X;
    return {prefactor * fai, prefactor * fbi,
            derivative_prefactor * (0.25 * fai - 1.5 * xi * fai_xi),
            derivative_prefactor * (0.25 * fbi - 1.5 * xi * fbi_xi), 0.0};
}

inline long double oscillatory_action_extended(long double a, long double p) {
    const long double t = std::tan(p) / a;
    if (std::abs(t) < 0.08) {
        const long double t2 = t * t;
        long double power = t * t2;
        long double sum = 0.0L;
        // atan(alpha*t)-alpha*atan(t), with the cancelled linear term omitted.
        for (int k = 1; k <= 8; ++k) {
            const long double numerator = std::pow(a, 2 * k + 1) - a;
            const long double term = ((k & 1) ? -1.0L : 1.0L)
                              * numerator * power / (2 * k + 1);
            sum += term;
            power *= t2;
        }
        return sum;
    }
    return p-a*std::atan2(std::sin(p),a*std::cos(p));
}

inline double oscillatory_action(long double a, long double p) {
    return static_cast<double>(oscillatory_action_extended(a,p));
}

inline long double evanescent_action_extended(long double a, long double tt) {
    if (std::abs(tt) < 0.08L) {
        const long double t2 = tt * tt;
        long double power = tt * t2;
        long double sum = 0.0L;
        // alpha*atanh(t)-atanh(alpha*t), with the linear term omitted.
        for (int k = 1; k <= 12; ++k) {
            sum += (a - std::pow(a, 2 * k + 1))
                 * power / (2 * k + 1);
            power *= t2;
        }
        return sum;
    }
    return a*std::atanh(tt)-std::atanh(a*tt);
}

inline double evanescent_action(long double a, long double tt) {
    return static_cast<double>(evanescent_action_extended(a,tt));
}

struct CanonicalMap {
    Region region;
    double q;
    double zeta;
    double action;
    double zeta_over_q;
};

template <std::size_t N>
struct SymmetricGaussRule {
    static_assert(N % 2 == 0, "Gauss order must be even");
    std::array<double, N / 2> nodes{};
    std::array<double, N / 2> weights{};
};

template <std::size_t N>
inline const SymmetricGaussRule<N>& gauss_legendre_rule() {
    static const SymmetricGaussRule<N> rule = [] {
        SymmetricGaussRule<N> out;
        for (std::size_t i = 0; i < N / 2; ++i) {
            double x = std::cos(pi * (static_cast<double>(i) + 0.75)
                              / (static_cast<double>(N) + 0.5));
            double derivative = 0.0;
            for (int iteration = 0; iteration < 20; ++iteration) {
                double p0 = 1.0, p1 = x;
                for (std::size_t k = 2; k <= N; ++k) {
                    const double pk = ((2.0 * k - 1.0) * x * p1
                                     - (k - 1.0) * p0) / k;
                    p0 = p1;
                    p1 = pk;
                }
                derivative = N * (x * p1 - p0) / (x * x - 1.0);
                const double next = x - p1 / derivative;
                if (std::abs(next - x) <= 4.0 * std::numeric_limits<double>::epsilon()
                                             * std::abs(next)) {
                    x = next;
                    break;
                }
                x = next;
            }
            double p0 = 1.0, p1 = x;
            for (std::size_t k = 2; k <= N; ++k) {
                const double pk = ((2.0 * k - 1.0) * x * p1
                                 - (k - 1.0) * p0) / k;
                p0 = p1;
                p1 = pk;
            }
            derivative = N * (x * p1 - p0) / (x * x - 1.0);
            out.nodes[i] = x;
            out.weights[i] = 2.0 / ((1.0 - x * x) * derivative * derivative);
        }
        return out;
    }();
    return rule;
}

#pragma float_control(precise, on, push)
struct DoubleDouble { double hi=0.0,lo=0.0; };
inline DoubleDouble dd_normalize(double h,double l){double s=h+l;return {s,l-(s-h)};}
inline DoubleDouble operator+(DoubleDouble a,DoubleDouble b){
    double s=a.hi+b.hi,v=s-a.hi,e=(a.hi-(s-v))+(b.hi-v)+a.lo+b.lo;return dd_normalize(s,e);
}
inline DoubleDouble operator-(DoubleDouble a,DoubleDouble b){return a+DoubleDouble{-b.hi,-b.lo};}
inline DoubleDouble operator*(DoubleDouble a,DoubleDouble b){
    double p=a.hi*b.hi,e=std::fma(a.hi,b.hi,-p)+a.hi*b.lo+a.lo*b.hi;return dd_normalize(p,e);
}
inline DoubleDouble operator*(DoubleDouble a,double b){
    double p=a.hi*b,e=std::fma(a.hi,b,-p)+a.lo*b;return dd_normalize(p,e);
}
inline DoubleDouble operator/(DoubleDouble a,DoubleDouble b){
    double q1=a.hi/b.hi;DoubleDouble r=a-b*DoubleDouble{q1,0.0};double q2=(r.hi+r.lo)/b.hi;
    return dd_normalize(q1,q2);
}
inline DoubleDouble dd_sqrt_positive(DoubleDouble a){
    const double seed=std::sqrt(a.hi+a.lo);
    DoubleDouble value{seed,0.0};
    value=(value+a/value)*0.5;
    return value;
}
inline DoubleDouble dd_log_positive(DoubleDouble a){
    int exponent=0;
    std::frexp(a.hi,&exponent);
    DoubleDouble mantissa{std::ldexp(a.hi,-exponent),
                          std::ldexp(a.lo,-exponent)};
    if(mantissa.hi<0.70710678118654752440){
        mantissa=mantissa*2.0;
        --exponent;
    }
    const DoubleDouble one{1.0,0.0};
    const DoubleDouble z=(mantissa-one)/(mantissa+one);
    const DoubleDouble z2=z*z;
    DoubleDouble term=z,sum=z;
    for(int denominator=3;denominator<=63;denominator+=2){
        term=term*z2;
        sum=sum+term*(1.0/static_cast<double>(denominator));
    }
    constexpr DoubleDouble ln2{0.69314718055994528623,
                               2.3190468138462995584e-17};
    return sum*2.0+ln2*static_cast<double>(exponent);
}
#pragma float_control(pop)

inline CanonicalMap canonical_map(double alpha, double theta, double alpha_tail=0.0) {
    const long double a = static_cast<long double>(alpha);
    const long double th = static_cast<long double>(theta);
    const long double s_ld = std::sin(th);
    const long double c_ld = std::cos(th);
    const long double sigma_ld = std::sqrt((1.0L-a)*(1.0L+a));
    const DoubleDouble exact_alpha{alpha,alpha_tail};
    const DoubleDouble sine{static_cast<double>(s_ld),0.0};
    const DoubleDouble cosine{static_cast<double>(c_ld),0.0};
    const DoubleDouble sine_squared=alpha>0.9 ?
        DoubleDouble{1.0,0.0}-cosine*cosine : sine*sine;
    const DoubleDouble potential=(sine_squared-exact_alpha*exact_alpha)/sine_squared;
    const long double q_ld=potential.hi+potential.lo;
    const double q = static_cast<double>(q_ld);
    const double tolerance = 64.0 * std::numeric_limits<double>::epsilon();
    CanonicalMap out{Region::turning, q, 0.0, 0.0, 0.0};
    if (std::abs(q) <= tolerance) {
        out.zeta_over_q = static_cast<double>(std::pow(
            2.0L*sigma_ld/a,-2.0L/3.0L));
        return out;
    }
    if (q > 0.0) {
        out.region = Region::oscillatory;
        const long double local_t=s_ld*std::sqrt(q_ld)/(a*c_ld);
        const long double phi=std::atan(a*local_t);
        out.action = oscillatory_action(a, phi);
        out.zeta = static_cast<double>(std::pow(1.5L*std::max(
            static_cast<long double>(out.action),0.0L),2.0L/3.0L));
    } else {
        out.region = Region::evanescent;
        const long double y=s_ld*std::sqrt(-q_ld)/c_ld;
        const long double t = RSCHACompat::clamp(y/a,0.0L,
            1.0L-4.0L*std::numeric_limits<long double>::epsilon());
        out.action = evanescent_action(a, t);
        out.zeta = -static_cast<double>(std::pow(1.5L*std::max(
            static_cast<long double>(out.action),0.0L),2.0L/3.0L));
    }
    out.zeta_over_q = (std::abs(q) > 1.0e-8 && std::abs(out.zeta) > 1.0e-16)
        ? out.zeta / q : static_cast<double>(std::pow(
            2.0L*sigma_ld/a,-2.0L/3.0L));
    return out;
}

inline double airy_residual_r0(double alpha) {
    const double a2 = alpha * alpha;
    return -std::cbrt(2.0) * (4.0 * a2 - 1.0) * (a2 - 4.0)
         / (280.0 * std::pow(alpha, 4.0 / 3.0)
            * std::pow(1.0 - a2, 4.0 / 3.0));
}

inline double airy_residual_r1(double alpha) {
    const double a2 = alpha * alpha;
    return (a2 + 1.0) * (8.0 * a2 * a2 - 9.0 * a2 + 8.0)
         / (300.0 * a2 * (1.0 - a2) * (1.0 - a2));
}

inline double airy_residual_r2(double alpha) {
    const double a2 = alpha * alpha;
    return -std::pow(2.0, 2.0 / 3.0) * (a2 - 4.0) * (4.0 * a2 - 1.0)
         * (69.0 * a2 * a2 + 13.0 * a2 + 69.0)
         / (53900.0 * std::pow(alpha, 8.0 / 3.0)
            * std::pow(1.0 - a2, 8.0 / 3.0));
}

inline double airy_residual(double alpha, double theta, const CanonicalMap& m,
                            bool potential_absorbed=false) {
    // The explicit series prevents catastrophic cancellation in the analytic
    // Schwarzian expression at the simple turning point.
    if (std::abs(m.zeta) < 0.01) {
        double value=airy_residual_r0(alpha)+airy_residual_r1(alpha)*m.zeta
             +airy_residual_r2(alpha)*m.zeta*m.zeta;
        if(potential_absorbed){const double s=std::sin(theta);
            value-=0.25*m.zeta_over_q/(s*s);}
        return value;
    }
    const double s = std::sin(theta);
    const double c = std::cos(theta);
    const double a2 = alpha * alpha;
    const double q_prime = 2.0 * a2 * c / (s * s * s);
    const double q_second = -2.0 * a2 / (s * s)
                          - 6.0 * a2 * c * c / (s * s * s * s);
    const double zeta_prime = std::sqrt(m.q / m.zeta);
    const double zeta_second = 0.5 * zeta_prime
        * (q_prime / m.q - zeta_prime / m.zeta);
    const double log_h_prime = 0.25
        * (zeta_prime / m.zeta - q_prime / m.q);
    const double log_h_second = 0.25
        * (zeta_second / m.zeta
           - zeta_prime * zeta_prime / (m.zeta * m.zeta)
           - q_second / m.q + q_prime * q_prime / (m.q * m.q));
    const double h_second_over_h = log_h_second + log_h_prime * log_h_prime;
    const double lower_order_potential = potential_absorbed?0.0:0.25 / (s * s);
    return (lower_order_potential + h_second_over_h) * m.zeta / m.q;
}

inline double theta_from_zeta(double alpha, double target_zeta,
                              double theta_endpoint) {
    const double theta_turning = std::asin(alpha);
    if (target_zeta == 0.0) return theta_turning;
    double lower = std::min(theta_turning, theta_endpoint);
    double upper = std::max(theta_turning, theta_endpoint);
    const double target_t = std::sqrt(std::abs(target_zeta));
    const CanonicalMap endpoint = canonical_map(alpha, theta_endpoint);
    const double endpoint_t = std::sqrt(std::abs(endpoint.zeta));
    double theta = theta_turning
        + (theta_endpoint - theta_turning) * target_t / endpoint_t;
    for (int iteration = 0; iteration < 18; ++iteration) {
        const CanonicalMap m = canonical_map(alpha, theta);
        const double residual = m.zeta - target_zeta;
        if (std::abs(residual) <= 8.0 * std::numeric_limits<double>::epsilon()
                                      * (1.0 + std::abs(target_zeta))) return theta;
        if (residual < 0.0) lower = theta; else upper = theta;
        double next = 0.5 * (lower + upper);
        if (m.zeta != 0.0 && m.q != 0.0) {
            const double derivative = std::sqrt(m.q / m.zeta);
            const double newton = theta - residual / derivative;
            if (std::isfinite(newton) && newton > lower && newton < upper)
                next = newton;
        }
        theta = next;
    }
    return theta;
}

inline double airy_transport_b0(double alpha, double theta,
                                const CanonicalMap& target,
                                bool potential_absorbed=false) {
    if (!potential_absorbed&&std::abs(target.zeta) < 0.01) {
        return airy_residual_r0(alpha)
             + airy_residual_r1(alpha) * target.zeta / 3.0
             + airy_residual_r2(alpha) * target.zeta * target.zeta / 5.0;
    }
    // Integrate in t=sqrt(|zeta|), where the canonical residual is analytic.
    // This removes the endpoint loss seen when theta itself is the quadrature
    // variable and makes a moderate-order rule substantially more accurate.
    const auto& rule = gauss_legendre_rule<8>();
    const double T = std::sqrt(std::abs(target.zeta));
    const double zeta_sign = std::signbit(target.zeta) ? -1.0 : 1.0;
    double average = 0.0;
    for (std::size_t i = 0; i < rule.nodes.size(); ++i) {
        for (double side : {-1.0, 1.0}) {
            const double parameter = 0.5 * (1.0 + side * rule.nodes[i]);
            const double t = parameter * T;
            const double zeta = zeta_sign * t * t;
            const double point = theta_from_zeta(alpha, zeta, theta);
            const CanonicalMap m = canonical_map(alpha, point);
            average += 0.5 * rule.weights[i]
                     * airy_residual(alpha, point, m,potential_absorbed);
        }
    }
    return average;
}

inline double airy_transport_a1(double alpha, double theta,
                                bool potential_absorbed=false) {
    // A1 is two asymptotic orders smaller; a 32-point outer rule combined
    // with the high-accuracy B0 kernel is sufficient and avoids O(N^2) cost.
    const auto& rule = gauss_legendre_rule<32>();
    const double theta_turning = std::asin(alpha);
    const CanonicalMap target = canonical_map(alpha, theta);
    const double target_R = airy_residual(alpha, theta, target,potential_absorbed);
    const double target_B = airy_transport_b0(alpha, theta, target,potential_absorbed);
    const double B_prime = std::abs(target.zeta) < 0.01
        ? airy_residual_r1(alpha) / 3.0
          + 2.0 * airy_residual_r2(alpha) * target.zeta / 5.0
        : (target_R - target_B) / (2.0 * target.zeta);
    const double half = 0.5 * (theta - theta_turning);
    const double middle = 0.5 * (theta + theta_turning);
    double integral = 0.0;
    for (std::size_t i = 0; i < rule.nodes.size(); ++i) {
        for (double side : {-1.0, 1.0}) {
            const double point = middle + half * side * rule.nodes[i];
            const CanonicalMap m = canonical_map(alpha, point);
            const double zeta_prime = std::sqrt(m.q / m.zeta);
            const double R = airy_residual(alpha, point, m,potential_absorbed);
            const double B = airy_transport_b0(alpha, point, m,potential_absorbed);
            integral += rule.weights[i] * R * B * zeta_prime;
        }
    }
    // Integrated form of A1'=-(B0''+R B0)/2.  It avoids a numerically
    // destructive second derivative of B0 and includes A1(0)=-R'(0)/6.
    return -0.5 * B_prime - 0.5 * half * integral;
}

struct DirectOscillatoryPair {
    bool ok = false;
    double p = 0.0;
    double q = 0.0;
    double p_log_scale = 0.0;
    double q_log_scale = 0.0;
};


inline long double sin_half_pi_dd(DoubleDouble x){
    const double cycles=std::nearbyint(x.hi/4.0);
    const DoubleDouble r=x-DoubleDouble{4.0*cycles,0.0};
    constexpr double halfpi_hi=1.5707963267948965580;
    constexpr double halfpi_lo=6.1232339957367660359e-17;
    const double angle_hi=halfpi_hi*r.hi;
    const double angle_lo=std::fma(halfpi_hi,r.hi,-angle_hi)
        +halfpi_hi*r.lo+halfpi_lo*r.hi;
    return static_cast<long double>(std::sin(angle_hi)*std::cos(angle_lo)
        +std::cos(angle_hi)*std::sin(angle_lo));
}

// Stable log(Gamma(x)/Gamma(x+1/2)).  Computing the two lgamma values
// separately loses several digits when x is large.  The Bernoulli-polynomial
// quotient expansion contains only odd inverse powers; recurrence first moves
// small arguments into its rapidly convergent range.
inline double log_gamma_half_ratio(double x){
    if(!(x>0.0) || !std::isfinite(x))
        return std::numeric_limits<double>::quiet_NaN();
    double correction=0.0;
    while(x<16.0){
        correction+=std::log1p(0.5/x);
        x+=1.0;
    }
    const double z=1.0/x, z2=z*z;
    // c_k=(2-2^{-k}) B_{k+1}/(k(k+1)), k=1,3,... .
    const double series=z*(1.0/8.0+z2*(-1.0/192.0+z2*(1.0/640.0
        +z2*(-17.0/14336.0+z2*(31.0/18432.0
        +z2*(-691.0/180224.0))))));
    return correction-0.5*std::log(x)+series;
}

struct SignedLogGamma {
    long double log_abs = 0.0L;
    int sign = 1;
    bool pole = false;
};

inline SignedLogGamma signed_log_gamma(long double x) {
    if (x > 0.0L) return {std::lgamma(x), 1, false};
    const long double sine = std::sin(static_cast<long double>(pi) * x);
    if (std::abs(sine) < 64.0L * std::numeric_limits<long double>::epsilon())
        return {0.0L, 0, true};
    return {std::log(static_cast<long double>(pi))-std::log(std::abs(sine))
            -std::lgamma(1.0L-x), sine < 0.0L ? -1 : 1, false};
}

struct PiOverTwoData { long double value, derivative; bool ok; };

// Canonically normalized sqrt(sin(theta))*P_nu^{-order}(cos(theta)) and
// its theta derivative at theta=pi/2, from the exact Gamma values.
inline PiOverTwoData normalized_p_at_pi_over_two(long double nu,
                                                  long double order) {
    const long double u=nu+0.5L;
    const long double a=0.5L*(nu-order+1.0L);
    const long double b=0.5L*(nu+order+1.0L);
    if(!(a>0.0)||!(b>0.0))return {0.0L,0.0L,false};
    const double ra=log_gamma_half_ratio(static_cast<double>(a));
    const double rb=log_gamma_half_ratio(static_cast<double>(b));
    if(!std::isfinite(ra)||!std::isfinite(rb))return {0.0L,0.0L,false};
    // Reflection followed by duplication reduces the exactly normalized
    // center data to half-step Gamma ratios.  This form has no subtraction
    // of O(nu log nu) quantities, even at nu=10^6.
    const long double amplitude=std::sqrt(u/static_cast<long double>(pi));
    const DoubleDouble phase=DoubleDouble{static_cast<double>(nu),0.0}
        -DoubleDouble{static_cast<double>(order),0.0};
    const long double sa=sin_half_pi_dd(phase+DoubleDouble{1.0,0.0});
    const long double sd=sin_half_pi_dd(phase);
    long double value=amplitude*sa*std::exp(0.5*(ra+rb));
    long double derivative=-2.0L*amplitude*sd*std::exp(-0.5*(ra+rb));
    // The representable double c=0.5*pi lies halfpi_lo below mathematical
    // pi/2.  Shift the exact Gamma anchor onto that actual ODE grid point;
    // otherwise u=1e6 inherits an avoidable O(1e-10) phase error.
    constexpr long double center_delta=-6.1232339957367660359e-17L;
    const long double potential=u*u+0.25L-order*order;
    const long double old_value=value;
    value+=center_delta*derivative;
    derivative-=center_delta*potential*old_value;
    if (!std::isfinite(value) || !std::isfinite(derivative))
        return {0.0L,0.0L,false};
    return {value,derivative,true};
}

struct ScaledOdeValue {
    bool ok = false;
    double value = 0.0;
    double log_scale = 0.0;
    double derivative_over_u = 0.0;
};

struct LogScalar { int sign=0; double log_abs=-std::numeric_limits<double>::infinity(); bool ok=false; };
struct EndpointScalar { LogScalar value,derivative; bool ok=false; };

inline long double endpoint_gamma_ratio_log(double nu,double order){
    const double nearest=std::round(order);
    if(order>=0.0&&nearest<=128.0&&std::abs(order-nearest)<1.0e-12){
        const int m=static_cast<int>(nearest);
        long double sum=0.0L,compensation=0.0L;
        for(int j=-m+1;j<=m;++j){
            const long double argument=static_cast<long double>(nu)+j;
            if(!(argument>0.0L))return std::numeric_limits<long double>::quiet_NaN();
            const long double term=std::log(argument)-compensation;
            const long double updated=sum+term;
            compensation=(updated-sum)-term;
            sum=updated;
        }
        return sum;
    }
    return std::lgamma(nu+order+1.0L)-std::lgamma(nu-order+1.0L);
}

inline LogScalar endpoint_normalized_p(double nu,double order,double theta){
    const long double gamma_ratio=endpoint_gamma_ratio_log(nu,order);
    if(!std::isfinite(gamma_ratio))return {};
    const SignedLogGamma g0=signed_log_gamma(1.0L+order);
    if(g0.pole)return {};
    long double term=g0.sign, sum=term, sum_abs=std::abs(term);
    long double log_scale=-g0.log_abs;
    const long double z=std::sin(0.5L*theta)*std::sin(0.5L*theta);
    const bool transformed=order>=0.0;
    bool converged=false;
    for(int k=0;k<600;++k){
        const long double den=(k+order+1.0L)*(k+1.0L);
        if(den==0.0L)return {};
        const long double numerator=transformed
            ?(k+order-nu)*(k+order+nu+1.0L)
            :(k-nu)*(k+nu+1.0L);
        term*=(numerator/den)*z;
        sum+=term;
        sum_abs+=std::abs(term);
        const long double mag=std::max(std::abs(term),std::abs(sum));
        if(mag>1.0e100L){term/=mag;sum/=mag;sum_abs/=mag;log_scale+=std::log(mag);}
        if(std::abs(term)<16.0L*std::numeric_limits<long double>::epsilon()*std::abs(sum)){
            converged=true;break;
        }
    }
    if(!converged||sum==0.0L||!std::isfinite(sum)
       ||(transformed&&sum_abs>1.0e6L*std::abs(sum)))return {};
    const long double logn=0.5L*(std::log(nu+0.5L)+gamma_ratio);
    const long double logv=logn+0.5L*std::log(std::sin(theta))
        +order*(transformed?std::log(0.5L*std::sin(theta))
                           :std::log(std::tan(0.5L*theta)))
        +log_scale+std::log(std::abs(sum));
    return {sum<0?-1:1,static_cast<double>(logv),std::isfinite(logv)};
}

inline EndpointScalar endpoint_normalized_p_with_derivative(double nu,double order,double theta){
    const long double gamma_ratio=endpoint_gamma_ratio_log(nu,order);
    const SignedLogGamma g0=signed_log_gamma(1.0L+order);
    if(!std::isfinite(gamma_ratio)||g0.pole)return {};
    const long double z=std::sin(0.5L*theta)*std::sin(0.5L*theta);
    long double term=g0.sign,sum=term,dzsum=0.0L,sum_abs=std::abs(term),
        dzsum_abs=0.0L,log_scale=-g0.log_abs;
    const bool transformed=order>=0.0;
    bool converged=false;
    for(int k=0;k<600;++k){
        const long double den=(k+order+1.0L)*(k+1.0L);
        if(den==0.0L)return {};
        const long double numerator=transformed
            ?(k+order-nu)*(k+order+nu+1.0L)
            :(k-nu)*(k+nu+1.0L);
        term*=(numerator/den)*z;
        sum+=term; dzsum+=(k+1.0L)*term/z;
        sum_abs+=std::abs(term);dzsum_abs+=std::abs((k+1.0L)*term/z);
        const long double mag=std::max({std::abs(term),std::abs(sum),std::abs(dzsum)});
        if(mag>1e100L){term/=mag;sum/=mag;dzsum/=mag;sum_abs/=mag;
            dzsum_abs/=mag;log_scale+=std::log(mag);}
        if(std::abs(term)<16.0L*std::numeric_limits<long double>::epsilon()
                         *std::max(1.0L,std::abs(sum))){converged=true;break;}
    }
    if(!converged||sum==0.0L||!std::isfinite(sum)||!std::isfinite(dzsum)
       ||(transformed&&(sum_abs>1.0e6L*std::abs(sum)
          ||dzsum_abs>1.0e6L*std::max(1.0L,std::abs(dzsum)))))return {};
    const long double logn=0.5L*(std::log(nu+0.5L)+gamma_ratio);
    const long double logv=logn+0.5L*std::log(std::sin(theta))
        +order*(transformed?std::log(0.5L*std::sin(theta))
                           :std::log(std::tan(0.5L*theta)))
        +log_scale+std::log(std::abs(sum));
    const long double ld=0.5L/std::tan(theta)
        +order*(transformed?1.0L/std::tan(theta):1.0L/std::sin(theta))
        +(dzsum/sum)*0.5L*std::sin(theta);
    LogScalar v{sum<0?-1:1,static_cast<double>(logv),std::isfinite(logv)};
    LogScalar d;
    if(ld==0.0L)d={0,-std::numeric_limits<double>::infinity(),true};
    else d={v.sign*(ld<0?-1:1),static_cast<double>(logv+std::log(std::abs(ld))),true};
    return {v,d,v.ok&&d.ok};
}

inline LogScalar log_linear_combination(const LogScalar&a,double ca,
                                         const LogScalar&b,double cb){
    if(!a.ok||!b.ok)return {};
    const int sa=ca==0?0:(ca<0?-a.sign:a.sign), sb=cb==0?0:(cb<0?-b.sign:b.sign);
    const double la=ca==0?-std::numeric_limits<double>::infinity():a.log_abs+std::log(std::abs(ca));
    const double lb=cb==0?-std::numeric_limits<double>::infinity():b.log_abs+std::log(std::abs(cb));
    const double m=std::max(la,lb);
    const double scaled=sa*std::exp(la-m)+sb*std::exp(lb-m);
    if(scaled==0.0||!std::isfinite(scaled))return {};
    return {scaled<0?-1:1,m+std::log(std::abs(scaled)),true};
}

// Removable integer-order limit of the normalized Q connection formula.
// Symmetric sampling cancels odd order derivatives.  Deliberately wide nodes
// avoid the severe cancellation of the connection quotient next to an
// integer; interpolation in h^2 then takes the removable limit at h=0.
inline LogScalar endpoint_integer_q_limit(double nu,int m,double theta){
    constexpr std::array<double,5> steps{{0.12,0.10,0.08,0.06,0.04}};
    if(m<1||static_cast<double>(m)+steps[0]>nu)return {};
    auto noninteger_q=[&](double order)->LogScalar{
        const LogScalar plus=endpoint_normalized_p(nu,order,theta);
        const LogScalar minus=endpoint_normalized_p(nu,-order,theta);
        if(!plus.ok||!minus.ok)return {};
        const double s=std::sin(pi*order);
        if(s==0.0)return {};
        LogScalar q=log_linear_combination(minus,1.0,plus,-std::cos(pi*order));
        if(!q.ok)return {};
        q.sign*=s<0.0?-1:1;
        q.log_abs-=std::log(std::abs(s));
        return q;
    };
    const double endpoint_parameter=std::pow(
        (nu+0.5)*std::sin(0.5*theta),2)/(1.0+m);
    if(endpoint_parameter<1.0e-3){
        constexpr double h=1.0e-3;
        std::array<LogScalar,6> q{};
        for(int j=0;j<3;++j){
            const double step=std::ldexp(h,-j);
            q[2*j]=noninteger_q(m+step);
            q[2*j+1]=noninteger_q(m-step);
            if(!q[2*j].ok||!q[2*j+1].ok)return {};
        }
        double scale=-std::numeric_limits<double>::infinity();
        for(const auto&v:q)scale=std::max(scale,v.log_abs);
        std::array<double,3> even{};
        for(int j=0;j<3;++j)even[j]=0.5*(
            q[2*j].sign*std::exp(q[2*j].log_abs-scale)
           +q[2*j+1].sign*std::exp(q[2*j+1].log_abs-scale));
        const double extrapolated=(even[0]-20.0*even[1]+64.0*even[2])/45.0;
        if(extrapolated==0.0||!std::isfinite(extrapolated))return {};
        return {extrapolated<0.0?-1:1,
            scale+std::log(std::abs(extrapolated)),true};
    }
    std::array<LogScalar,10> q{};
    for(int j=0;j<5;++j){
        q[2*j]=noninteger_q(m+steps[j]);
        q[2*j+1]=noninteger_q(m-steps[j]);
        if(!q[2*j].ok||!q[2*j+1].ok)return {};
    }
    double scale=-std::numeric_limits<double>::infinity();
    for(const auto&v:q)scale=std::max(scale,v.log_abs);
    if(!std::isfinite(scale))return {};
    std::array<double,5> even{};
    for(int j=0;j<5;++j)even[j]=0.5*(
        q[2*j].sign*std::exp(q[2*j].log_abs-scale)
       +q[2*j+1].sign*std::exp(q[2*j+1].log_abs-scale));
    double extrapolated=0.0;
    for(int j=0;j<5;++j){
        const double tj=steps[j]*steps[j];
        double weight=1.0;
        for(int k=0;k<5;++k)if(k!=j){
            const double tk=steps[k]*steps[k];
            weight*=(-tk)/(tj-tk);
        }
        extrapolated+=weight*even[j];
    }
    if(extrapolated==0.0||!std::isfinite(extrapolated))return {};
    return {extrapolated<0.0?-1:1,
        scale+std::log(std::abs(extrapolated)),true};
}

inline LogScalar endpoint_integer_q(double nu,int m,double theta){
    if(m<1||m>nu)return {};
    constexpr int terms=18;
    std::array<long double,terms> inverse_sinc_squared{};
    inverse_sinc_squared[0]=1.0L;
    for(int n=1;n<terms;++n){
        long double convolution=0.0L;
        for(int k=1;k<=n;++k){
            const long double g=(k&1?-1.0L:1.0L)
                *std::pow(2.0L,2*k+1)/std::tgamma(2*k+3.0L);
            convolution+=g*inverse_sinc_squared[n-k];
        }
        inverse_sinc_squared[n]=-convolution;
    }
    const long double u=nu+0.5L,c=0.25L-static_cast<long double>(m)*m;
    std::array<long double,terms> a{}; a[0]=1.0L;
    for(int n=1;n<terms;++n){
        if(n>=m)break;
        long double rhs=u*u*a[n-1];
        for(int j=0;j<n;++j)rhs+=c*inverse_sinc_squared[j+1]*a[n-1-j];
        a[n]=-rhs/(4.0L*n*(n-m));
    }
    const long double x=theta*theta; long double sum=0.0L;
    for(int n=std::min(terms-1,m-1);n>=0;--n)sum=sum*x+a[n];
    if(!(sum>0.0L)||!std::isfinite(sum))return {};
    const long double logn=0.5L*(std::log(u)-endpoint_gamma_ratio_log(nu,m));
    const long double logq=m*std::log(2.0L)-std::log(static_cast<long double>(pi))
        +std::lgamma(static_cast<long double>(m))+logn
        +(0.5L-m)*std::log(theta)+std::log(sum);
    return {1,static_cast<double>(logq),std::isfinite(logq)};
}

inline double digamma_positive(double x){
    double correction=0.0;
    while(x<12.0){correction-=1.0/x;x+=1.0;}
    const double z=1.0/x,z2=z*z;
    return correction+std::log(x)-0.5*z-z2*(1.0/12.0-z2*(1.0/120.0
        -z2*(1.0/252.0-z2*(1.0/240.0-z2*(1.0/132.0)))));
}

// At order zero the endpoint connection formula is 0/0.  Differentiate its
// regular hypergeometric series analytically with respect to order:
// Q_tilde=-(2/pi) d(P_tilde)/dmu at mu=0.  This supplies the pole-transition
// layer without importing a Bessel evaluator.
inline DirectOscillatoryPair endpoint_zero_order_pair(double nu,double theta){
    const long double z=std::sin(0.5L*theta)*std::sin(0.5L*theta);
    long double term=1.0L,sum=1.0L,dsum=0.0L,harmonic=0.0L;
    bool converged=false;
    for(int k=0;k<300;++k){
        term*=((k-nu)*(k+nu+1.0L)/((k+1.0L)*(k+1.0L)))*z;
        harmonic+=1.0L/(k+1.0L);
        sum+=term;
        dsum-=harmonic*term;
        if(std::abs(term)<16.0L*std::numeric_limits<long double>::epsilon()
                          *std::max(1.0L,std::abs(sum))){converged=true;break;}
    }
    if(!converged||!std::isfinite(sum)||!std::isfinite(dsum))return {};
    const long double pref=std::sqrt((nu+0.5L)*std::sin(theta));
    const long double dlog=digamma_positive(nu+1.0)+std::log(std::tan(0.5L*theta))
        -digamma_positive(1.0);
    const long double p=pref*sum;
    const long double q=-2.0L/static_cast<long double>(pi)*pref*(dlog*sum+dsum);
    if(!std::isfinite(p)||!std::isfinite(q))return {};
    return {true,static_cast<double>(p),static_cast<double>(q),0.0,0.0};
}

inline DirectOscillatoryPair direct_endpoint_pair(double nu,double mu,double theta){
    const long double s=std::sin(static_cast<long double>(pi)*mu);
    const LogScalar plus=endpoint_normalized_p(nu,mu,theta);
    const double nearest=std::round(mu);
    if(std::abs(mu-nearest)<1.0e-12&&nearest>=1.0){
        LogScalar q=endpoint_integer_q_limit(nu,static_cast<int>(nearest),theta);
        if(!q.ok)q=endpoint_integer_q(nu,static_cast<int>(nearest),theta);
        if(!plus.ok||!q.ok)return {};
        return {true,static_cast<double>(plus.sign),static_cast<double>(q.sign),
                plus.log_abs,q.log_abs};
    }
    if(std::abs(s)<1.0e-3L)return {};
    const LogScalar minus=endpoint_normalized_p(nu,-mu,theta);
    if(!plus.ok||!minus.ok)return {};
    LogScalar q=log_linear_combination(minus,1.0,plus,-std::cos(pi*mu));
    if(!q.ok)return {};
    q.sign*=s<0?-1:1; q.log_abs-=std::log(std::abs(static_cast<double>(s)));
    return {true,static_cast<double>(plus.sign),static_cast<double>(q.sign),
            plus.log_abs,q.log_abs};
}

struct CenterQData { long double value, derivative; bool ok; };

// The connection quotient for Q has a removable singularity at every integer
// order and is badly conditioned nearby.  Reflection and duplication reduce
// it identically to the same half-step ratios used for P, with sine replaced
// by cosine.  Thus no numerical l'Hopital limit or large-lgamma subtraction
// is required.
inline CenterQData normalized_q_at_pi_over_two(long double nu,
                                                long double order) {
    const long double u=nu+0.5L;
    const long double a=0.5L*(nu-order+1.0L);
    const long double b=0.5L*(nu+order+1.0L);
    if(!(a>0.0)||!(b>0.0))return {0.0L,0.0L,false};
    const double ra=log_gamma_half_ratio(static_cast<double>(a));
    const double rb=log_gamma_half_ratio(static_cast<double>(b));
    if(!std::isfinite(ra)||!std::isfinite(rb))return {0.0L,0.0L,false};
    const long double amplitude=std::sqrt(u/static_cast<long double>(pi));
    const DoubleDouble phase=DoubleDouble{static_cast<double>(nu),0.0}
        -DoubleDouble{static_cast<double>(order),0.0};
    const long double ca=sin_half_pi_dd(phase+DoubleDouble{2.0,0.0});
    const long double sa=sin_half_pi_dd(phase+DoubleDouble{1.0,0.0});
    long double value=amplitude*ca*std::exp(0.5*(ra+rb));
    long double derivative=-2.0L*amplitude*sa*std::exp(-0.5*(ra+rb));
    constexpr long double center_delta=-6.1232339957367660359e-17L;
    const long double potential=u*u+0.25L-order*order;
    const long double old_value=value;
    value+=center_delta*derivative;
    derivative-=center_delta*potential*old_value;
    return {value,derivative,std::isfinite(value)&&std::isfinite(derivative)};
}

// Propagate one normal-form solution y''+[u^2+(1/4-mu^2)csc^2(theta)]y=0.
// Dynamic rescaling retains its logarithm in deep forbidden regions.
inline ScaledOdeValue propagate_scaled(double u, double mu, double theta,
                                       double target, double value,
                                       double derivative_over_u,
                                       double initial_log_scale,
                                       double rtol=2.0e-14,
                                       double atol=2.0e-16) {
    std::array<double,2> y{value,derivative_over_u};
    double log_scale = initial_log_scale;
    double h = std::copysign(std::min(0.02/u,
        std::max(1.0e-6,std::abs(target-theta)/20.0)),target-theta);
    auto rhs=[&](double x,const std::array<double,2>& s){
        const double sn=std::sin(x);
        const double potential=u*u+(0.25-mu*mu)/(sn*sn);
        return std::array<double,2>{u*s[1],-potential*s[0]/u};
    };
    for(int step=0;step<1000000;++step){
        if((target-theta)*h<=0.0 || std::abs(h)>std::abs(target-theta)) h=target-theta;
        if(h==0.0) return {true,y[0],log_scale,y[1]};
        const auto k1=rhs(theta,y); auto s=y;
        for(int i=0;i<2;++i)s[i]=y[i]+h*k1[i]/5.0;
        const auto k2=rhs(theta+h/5.0,s);
        for(int i=0;i<2;++i)s[i]=y[i]+h*(3*k1[i]/40.0+9*k2[i]/40.0);
        const auto k3=rhs(theta+3*h/10.0,s);
        for(int i=0;i<2;++i)s[i]=y[i]+h*(44*k1[i]/45.0-56*k2[i]/15.0+32*k3[i]/9.0);
        const auto k4=rhs(theta+4*h/5.0,s);
        for(int i=0;i<2;++i)s[i]=y[i]+h*(19372*k1[i]/6561.0-25360*k2[i]/2187.0+64448*k3[i]/6561.0-212*k4[i]/729.0);
        const auto k5=rhs(theta+8*h/9.0,s);
        for(int i=0;i<2;++i)s[i]=y[i]+h*(9017*k1[i]/3168.0-355*k2[i]/33.0+46732*k3[i]/5247.0+49*k4[i]/176.0-5103*k5[i]/18656.0);
        const auto k6=rhs(theta+h,s); std::array<double,2> fifth{};
        for(int i=0;i<2;++i)fifth[i]=y[i]+h*(35*k1[i]/384.0+500*k3[i]/1113.0+125*k4[i]/192.0-2187*k5[i]/6784.0+11*k6[i]/84.0);
        const auto k7=rhs(theta+h,fifth); double error=0.0;
        for(int i=0;i<2;++i){
            const double fourth=y[i]+h*(5179*k1[i]/57600.0+7571*k3[i]/16695.0+393*k4[i]/640.0-92097*k5[i]/339200.0+187*k6[i]/2100.0+k7[i]/40.0);
            const double scale=atol+rtol*std::max(std::abs(y[i]),std::abs(fifth[i]));
            error=std::max(error,std::abs(fifth[i]-fourth)/scale);
        }
        if(error<=1.0){
            theta+=h; y=fifth;
            const double magnitude=std::max(std::abs(y[0]),std::abs(y[1]));
            if(magnitude>1.0e100 || (magnitude>0.0 && magnitude<1.0e-100)){
                y[0]/=magnitude; y[1]/=magnitude; log_scale+=std::log(magnitude);
            }
        }
        const double factor=error==0.0?5.0:RSCHACompat::clamp(0.9*std::pow(error,-0.2),0.2,5.0);
        h*=factor;
        if(std::abs(h)<2.0e-15 || !std::isfinite(y[0]) || !std::isfinite(y[1])) return {};
    }
    return {};
}

inline DirectOscillatoryPair endpoint_propagated_pair(double nu,double mu,double target){
    const double u=nu+0.5;
    if(std::abs(std::sin(pi*mu))<1.0e-3)return {};
    const double argument=std::sqrt(4.0*(1.0+mu))/u;
    if(!(argument>0.0&&argument<1.0))return {};
    const double theta0=2.0*std::asin(argument);
    if(!(theta0<target))return {};
    const EndpointScalar plus=endpoint_normalized_p_with_derivative(nu,mu,theta0);
    const EndpointScalar minus=endpoint_normalized_p_with_derivative(nu,-mu,theta0);
    if(!plus.ok||!minus.ok)return {};
    const double c=-std::cos(pi*mu),s=std::sin(pi*mu);
    LogScalar q=log_linear_combination(minus.value,1.0,plus.value,c);
    LogScalar qd=log_linear_combination(minus.derivative,1.0,plus.derivative,c);
    if(!q.ok||!qd.ok)return {};
    q.sign*=s<0?-1:1;qd.sign*=s<0?-1:1;
    q.log_abs-=std::log(std::abs(s));qd.log_abs-=std::log(std::abs(s));
    auto advance=[&](const LogScalar&v,const LogScalar&d){
        if(v.sign==0)return ScaledOdeValue{};
        const double derivative_scaled=d.sign==0?0.0:
            static_cast<double>(d.sign)*std::exp(d.log_abs-v.log_abs)/u;
        return propagate_scaled(u,mu,theta0,target,static_cast<double>(v.sign),
                                derivative_scaled,v.log_abs);
    };
    const ScaledOdeValue p=advance(plus.value,plus.derivative);
    // The regular P solution is stable outward.  The dominant Q solution is
    // instead retained from its endpoint connection series at the target;
    // outward propagation of Q is exponentially ill-conditioned.
    const DirectOscillatoryPair target_series=direct_endpoint_pair(nu,mu,target);
    if(!p.ok||!target_series.ok)return {};
    return {true,p.value,target_series.q,p.log_scale,target_series.q_log_scale};
}

// In the forbidden region the recessive P solution is propagated outward
// from its exact regular-endpoint hypergeometric germ.  Q is propagated in
// the stable direction from pi/2.  Thus neither solution is obtained by an
// ill-conditioned extraction from the dominant solution.
inline DirectOscillatoryPair direct_evanescent_pair(double nu, double mu,
                                                     double target_theta) {
    const double u=nu+0.5;
    const CenterQData center_q=normalized_q_at_pi_over_two(nu,mu);
    if(!center_q.ok) return {};
    const long double q0=center_q.value, qd0=center_q.derivative;
    double qscale=std::max(std::abs(static_cast<double>(q0)),
                           std::abs(static_cast<double>(qd0/u)));
    if(!(qscale>0.0) || !std::isfinite(qscale)) return {};
    const ScaledOdeValue q=propagate_scaled(u,mu,0.5*pi,target_theta,
        static_cast<double>(q0)/qscale,static_cast<double>(qd0/u)/qscale,
        std::log(qscale));
    if(!q.ok) return {};

    const double theta0=std::min(target_theta,0.02/u);
    const long double z=std::sin(0.5L*theta0)*std::sin(0.5L*theta0);
    long double term=1.0L, hyper=1.0L, hyper_derivative=0.0L;
    for(int k=0;k<80;++k){
        term*=((k-nu)*(k+nu+1.0L)/((k+mu+1.0L)*(k+1.0L)))*z;
        hyper+=term;
        hyper_derivative+=(k+1.0L)*term/z;
        if(std::abs(term)<8.0L*std::numeric_limits<long double>::epsilon()*std::abs(hyper)) break;
        if(k==79 || !std::isfinite(hyper)) return {};
    }
    if(!(hyper>0.0L)) return {};
    const long double log_normalization=0.5L*(std::log(static_cast<long double>(u))
        +std::lgamma(nu+mu+1.0L)-std::lgamma(nu-mu+1.0L));
    const long double log_p0=log_normalization+0.5L*std::log(std::sin(theta0))
        +mu*std::log(std::tan(0.5L*theta0))-std::lgamma(1.0L+mu)+std::log(hyper);
    const long double logarithmic_derivative=0.5L/std::tan(theta0)+mu/std::sin(theta0)
        +(hyper_derivative/hyper)*0.5L*std::sin(theta0);
    const ScaledOdeValue p=propagate_scaled(u,mu,theta0,target_theta,1.0,
        static_cast<double>(logarithmic_derivative/u),static_cast<double>(log_p0));
    if(!p.ok) return {};
    return {true,p.value,q.value,p.log_scale,q.log_scale};
}

// Exact center initialization with independently scaled propagation in both
// theta directions.
inline DirectOscillatoryPair direct_center_scaled_pair(double nu,double mu,
                                                        double target_theta){
    const double u=nu+0.5;
    const PiOverTwoData p0=normalized_p_at_pi_over_two(nu,mu);
    const CenterQData q0=normalized_q_at_pi_over_two(nu,mu);
    if(!p0.ok||!q0.ok)return {};
    auto propagate=[&](long double value,long double derivative)->ScaledOdeValue{
        double scale=std::max(std::abs(static_cast<double>(value)),
                              std::abs(static_cast<double>(derivative/u)));
        if(!(scale>0.0)||!std::isfinite(scale))return {};
        return propagate_scaled(u,mu,0.5*pi,target_theta,
            static_cast<double>(value)/scale,
            static_cast<double>(derivative/u)/scale,std::log(scale));
    };
    const ScaledOdeValue p=propagate(p0.value,p0.derivative);
    const ScaledOdeValue q=propagate(q0.value,q0.derivative);
    if(!p.ok||!q.ok)return {};
    return {true,p.value,q.value,p.log_scale,q.log_scale};
}

// Pole-separated canonical residual chart.  For
//   y''+[u^2+lambda csc^2(theta)]y=0,
// write Z=exp(i*u*(theta-pi/2))*sum a_n/(2*i*u)^n.  With t=cot(theta),
// the transport polynomials are generated without numerical quadrature by
//   d a_{n+1}/dt=(1+t^2)a_n''+2t a_n'+lambda a_n.
// Their integration constants are zero, and exact center data fixes the two
// real linear combinations.  This is an independent residual expansion, not
// a call to a Bessel or reference Legendre evaluator.
inline DirectOscillatoryPair pole_residual_pair(double nu,double mu,double theta,
                                                 double* p_derivative=nullptr,
                                                 double* q_derivative=nullptr,
                                                 double* truncation_error=nullptr){
    constexpr int count=30, width=62;
    using Poly=std::array<long double,width>;
    std::array<Poly,count> a{};
    a[0][0]=1.0L;
    const long double lambda=0.25L-static_cast<long double>(mu)*mu;
    for(int n=0;n<count-1;++n){
        Poly rhs{};
        for(int k=1;k<width;++k)rhs[k]+=2.0L*k*a[n][k];
        for(int k=2;k<width;++k){
            rhs[k-2]+=k*(k-1)*a[n][k];
            rhs[k]+=k*(k-1)*a[n][k];
        }
        for(int k=0;k<width;++k)rhs[k]+=lambda*a[n][k];
        for(int k=0;k<width-1;++k)a[n+1][k+1]=rhs[k]/(k+1.0L);
    }
    struct Basis{std::complex<double> z,zd;double error;bool ok;};
    const double u=nu+0.5;
    auto basis=[&](double x)->Basis{
        const long double t=1.0L/std::tan(static_cast<long double>(x));
        std::complex<long double> sum{0.0L,0.0L},sum_theta{0.0L,0.0L};
        std::complex<long double> factor{1.0L,0.0L};
        const std::complex<long double> step{0.0L,-0.5L/u};
        long double last_value=0.0L,last_derivative=0.0L;
        for(int n=0;n<count;++n){
            long double p=0.0L,dp=0.0L;
            for(int k=width-1;k>=0;--k)p=p*t+a[n][k];
            for(int k=width-1;k>=1;--k)dp=dp*t+k*a[n][k];
            sum+=factor*p;
            sum_theta+=factor*(-(1.0L+t*t)*dp);
            if(n==count-1){
                last_value=std::abs(factor*p);
                last_derivative=std::abs(factor*
                    (std::complex<long double>(0.0L,u)*p-(1.0L+t*t)*dp))/u;
            }
            factor*=step;
        }
        // Form and reduce the carrier in extended precision.  At u=1e6 a
        // double phase product alone leaves O(1e-10) value errors near zeros.
        const long double pi_ld=std::acos(-1.0L);
        long double phase=0.0L;
        if(nu==std::floor(nu)){
            const auto degree=static_cast<unsigned long long>(nu);
            phase=static_cast<long double>(u)*static_cast<long double>(x)
                -static_cast<long double>(degree&3ULL)*(0.5L*pi_ld)
                -0.25L*pi_ld;
        }else{
            phase=static_cast<long double>(u)*
                (static_cast<long double>(x)-0.5L*pi_ld);
        }
        const long double reduced=std::remainder(phase,2.0L*pi_ld);
        const std::complex<long double> carrier{std::cos(reduced),std::sin(reduced)};
        const std::complex<long double> z=carrier*sum;
        const std::complex<long double> zd=carrier*
            (std::complex<long double>(0.0L,u)*sum+sum_theta);
        const long double value_scale=std::max(std::abs(sum),1.0e-300L);
        const long double derivative_scale=std::max(
            std::abs(std::complex<long double>(0.0L,1.0L)*sum
                     +sum_theta/static_cast<long double>(u)),
            1.0e-300L);
        const double error=static_cast<double>(
            std::max(last_value/value_scale,last_derivative/derivative_scale));
        Basis out{{static_cast<double>(z.real()),static_cast<double>(z.imag())},
                  {static_cast<double>(zd.real()),static_cast<double>(zd.imag())},
                  error,true};
        out.ok=std::isfinite(out.z.real())&&std::isfinite(out.z.imag())
            &&std::isfinite(out.zd.real())&&std::isfinite(out.zd.imag())
            &&std::isfinite(out.error);
        return out;
    };
    const Basis center=basis(0.5*pi), target=basis(theta);
    const PiOverTwoData p0=normalized_p_at_pi_over_two(nu,mu);
    const CenterQData q0=normalized_q_at_pi_over_two(nu,mu);
    if(!center.ok||!target.ok||!p0.ok||!q0.ok)return {};
    const double det=center.z.real()*center.zd.imag()-center.z.imag()*center.zd.real();
    if(!(std::abs(det)>0.0)||!std::isfinite(det))return {};
    auto transfer=[&](long double y,long double yd){
        const double A=(static_cast<double>(y)*center.zd.imag()
                       -static_cast<double>(yd)*center.z.imag())/det;
        const double B=(center.z.real()*static_cast<double>(yd)
                       -center.zd.real()*static_cast<double>(y))/det;
        return std::array<double,2>{A*target.z.real()+B*target.z.imag(),
                                    A*target.zd.real()+B*target.zd.imag()};
    };
    const auto p=transfer(p0.value,p0.derivative);
    const auto q=transfer(q0.value,q0.derivative);
    if(!std::isfinite(p[0])||!std::isfinite(q[0])||!std::isfinite(p[1])
       ||!std::isfinite(q[1]))return {};
    if(p_derivative)*p_derivative=p[1];
    if(q_derivative)*q_derivative=q[1];
    if(truncation_error)*truncation_error=target.error;
    return {true,p[0],q[0],0.0,0.0};
}

inline DirectOscillatoryPair direct_oscillatory_pair(double nu, double mu,
                                                       double target_theta) {
    const PiOverTwoData plus = normalized_p_at_pi_over_two(nu, mu);
    const CenterQData center_q=normalized_q_at_pi_over_two(nu,mu);
    if (!plus.ok || !center_q.ok) return {};
    std::array<double,4> y{
        static_cast<double>(plus.value),
        static_cast<double>(plus.derivative/(nu+0.5)),
        static_cast<double>(center_q.value),
        static_cast<double>(center_q.derivative/(nu+0.5))};
    for (double value : y) if (!std::isfinite(value)) return {};

    const double u = nu + 0.5;
    auto rhs = [&](double theta, const std::array<double,4>& state) {
        const double s = std::sin(theta);
        const double potential = u*u + (0.25-mu*mu)/(s*s);
        return std::array<double,4>{u*state[1],-potential*state[0]/u,
                                    u*state[3],-potential*state[2]/u};
    };
    double theta = 0.5*pi;
    double h = -std::min(0.02/u, std::max(1.0e-5,(theta-target_theta)/20.0));
    constexpr double relative_tolerance = 2.0e-12;
    constexpr double absolute_tolerance = 2.0e-14;
    for (int step = 0; step < 1000000 && theta > target_theta; ++step) {
        if (theta+h < target_theta) h = target_theta-theta;
        const auto k1=rhs(theta,y);
        auto stage=y;
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(k1[i]/5.0);
        const auto k2=rhs(theta+h/5.0,stage);
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(3*k1[i]/40.0+9*k2[i]/40.0);
        const auto k3=rhs(theta+3*h/10.0,stage);
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(44*k1[i]/45.0-56*k2[i]/15.0+32*k3[i]/9.0);
        const auto k4=rhs(theta+4*h/5.0,stage);
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(19372*k1[i]/6561.0-25360*k2[i]/2187.0+64448*k3[i]/6561.0-212*k4[i]/729.0);
        const auto k5=rhs(theta+8*h/9.0,stage);
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(9017*k1[i]/3168.0-355*k2[i]/33.0+46732*k3[i]/5247.0+49*k4[i]/176.0-5103*k5[i]/18656.0);
        const auto k6=rhs(theta+h,stage);
        std::array<double,4> fifth{};
        for(int i=0;i<4;++i)fifth[i]=y[i]+h*(35*k1[i]/384.0+500*k3[i]/1113.0+125*k4[i]/192.0-2187*k5[i]/6784.0+11*k6[i]/84.0);
        const auto k7=rhs(theta+h,fifth);
        std::array<double,4> fourth{};
        double error=0.0;
        for(int i=0;i<4;++i){
            fourth[i]=y[i]+h*(5179*k1[i]/57600.0+7571*k3[i]/16695.0+393*k4[i]/640.0-92097*k5[i]/339200.0+187*k6[i]/2100.0+k7[i]/40.0);
            const double scale=absolute_tolerance+relative_tolerance*std::max(std::abs(y[i]),std::abs(fifth[i]));
            error=std::max(error,std::abs(fifth[i]-fourth[i])/scale);
        }
        if(error<=1.0){theta+=h;y=fifth;}
        const double factor=error==0.0?5.0:RSCHACompat::clamp(0.9*std::pow(error,-0.2),0.2,5.0);
        h*=factor;
        if(std::abs(h)<2.0e-15) return {};
        if(theta<=target_theta){
            return {true,y[0],y[2],0.0,0.0};
        }
    }
    return {};
}

struct AiryBasisValue { bool ok=false; double a=0.0,b=0.0; };

inline double airy_transport_b0_fast(double alpha, double theta,
                                     const CanonicalMap &target) {
    if (std::abs(target.zeta) < 0.01)
        return airy_residual_r0(alpha)
             + airy_residual_r1(alpha) * target.zeta / 3.0
             + airy_residual_r2(alpha) * target.zeta * target.zeta / 5.0;
    const auto &rule = gauss_legendre_rule<4>();
    const double sine = std::sin(theta);
    const double cosine = std::cos(theta);
    const double sigma = std::sqrt((1.0 - alpha) * (1.0 + alpha));
    const bool oscillatory = target.region == Region::oscillatory;
    const double endpoint = oscillatory
        ? std::atan2(std::sqrt(std::max(0.0, sine * sine - alpha * alpha)), cosine)
        : std::sqrt(std::max(0.0, alpha * alpha - sine * sine))
            / (alpha * cosine);
    double average = 0.0;
    for (std::size_t i = 0; i < rule.nodes.size(); ++i) {
        for (double side : {-1.0, 1.0}) {
            const double parameter = 0.5 * (1.0 + side * rule.nodes[i]);
            const double wanted = parameter * parameter * parameter
                                * target.action;
            double lo = 0.0, hi = endpoint;
            double coordinate = parameter * endpoint;
            for (int iteration = 0; iteration < 32; ++iteration) {
                const double action = oscillatory
                    ? oscillatory_action(alpha, coordinate)
                    : evanescent_action(alpha, coordinate);
                if (std::abs(action - wanted) <=
                    64.0 * std::numeric_limits<double>::epsilon() *
                        (1.0 + std::abs(wanted)))
                    break;
                if (action < wanted) lo = coordinate; else hi = coordinate;
                double derivative = 0.0;
                if (oscillatory) {
                    const double sp = std::sin(coordinate), cp = std::cos(coordinate);
                    derivative = sp * sp * (1.0 - alpha * alpha)
                               / (alpha * alpha * cp * cp + sp * sp);
                } else {
                    const double t2 = coordinate * coordinate;
                    derivative = alpha * t2 * (1.0 - alpha * alpha)
                               / ((1.0 - t2) * (1.0 - alpha * alpha * t2));
                }
                const double candidate = coordinate - (action - wanted) / derivative;
                coordinate = std::isfinite(candidate) && candidate > lo && candidate < hi
                    ? candidate : 0.5 * (lo + hi);
            }
            const double point = oscillatory
                ? std::acos(sigma * std::cos(coordinate))
                : std::acos(sigma / std::sqrt(
                      1.0 - alpha * alpha * coordinate * coordinate));
            const CanonicalMap map = canonical_map(alpha, point);
            average += 0.5 * rule.weights[i]
                     * airy_residual(alpha, point, map);
        }
    }
    return average;
}

inline double airy_transport_a1_fast(double alpha, double theta) {
    const auto &rule = gauss_legendre_rule<8>();
    const double turning = std::asin(alpha);
    const CanonicalMap target = canonical_map(alpha, theta);
    const double target_r = airy_residual(alpha, theta, target);
    const double target_b = airy_transport_b0_fast(alpha, theta, target);
    const double b_prime = std::abs(target.zeta) < 0.01
        ? airy_residual_r1(alpha) / 3.0
          + 2.0 * airy_residual_r2(alpha) * target.zeta / 5.0
        : (target_r - target_b) / (2.0 * target.zeta);
    const double half = 0.5 * (theta - turning);
    const double middle = 0.5 * (theta + turning);
    double integral = 0.0;
    for (std::size_t i = 0; i < rule.nodes.size(); ++i)
        for (double side : {-1.0, 1.0}) {
            const double point = middle + half * side * rule.nodes[i];
            const CanonicalMap map = canonical_map(alpha, point);
            const double zeta_prime = std::sqrt(map.q / map.zeta);
            const double residual = airy_residual(alpha, point, map);
            const double b0 = airy_transport_b0_fast(alpha, point, map);
            integral += rule.weights[i] * residual * b0 * zeta_prime;
        }
    return -0.5 * b_prime - 0.5 * half * integral;
}

inline AiryBasisValue airy_basis_value_fast(double nu, double mu,
                                             double theta) {
    const double u = nu + 0.5;
    const double alpha = mu / u;
    const CanonicalMap map = canonical_map(alpha, theta,std::fma(-alpha,u,mu)/u);
    if (!(map.zeta_over_q > 0.0) || !std::isfinite(map.zeta_over_q)) return {};
    const double u23 = std::pow(u, 2.0 / 3.0);
    const AiryScaled ab = airy(-u23 * map.zeta, u * map.action);
    const double derivative_scale = std::pow(u, -4.0 / 3.0)
                                  * airy_transport_b0_fast(alpha, theta, map);
    const double relative_a1 = u < 4096.0
        ? airy_transport_a1_fast(alpha, theta) / (u * u) : 0.0;
    const double value_scale = 1.0 + (std::abs(relative_a1) <= 0.01
                                   ? relative_a1 : 0.0);
    const double ca = value_scale * ab.ai_scaled
                    - derivative_scale * ab.aip_scaled;
    const double cb = value_scale * ab.bi_scaled
                    - derivative_scale * ab.bip_scaled;
    if (ab.exponent > 30.0) return {};
    const double pref = std::sqrt(2.0) * std::pow(u, 1.0 / 6.0)
                      * std::pow(map.zeta_over_q, 0.25);
    const double a = pref * ca * std::exp(-ab.exponent);
    const double b = pref * cb * std::exp(ab.exponent);
    return {std::isfinite(a) && std::isfinite(b), a, b};
}

inline AiryBasisValue airy_basis_value(double nu,double mu,double theta,
                                       double effective_alpha=-1.0,
                                       bool potential_absorbed=false){
    const double u=nu+0.5;
    const double alpha=effective_alpha>=0.0?effective_alpha:mu/u;
    const CanonicalMap map=canonical_map(alpha,theta,effective_alpha<0.0 ? std::fma(-alpha,u,mu)/u : 0.0);
    if(!(map.zeta_over_q>0.0)||!std::isfinite(map.zeta_over_q))return {};
    const double u23=std::pow(u,2.0/3.0),arg=-u23*map.zeta;
    const AiryScaled ab=airy(arg,u*map.action);
    const double b0=airy_transport_b0(alpha,theta,map,potential_absorbed);
    const double a1=airy_transport_a1(alpha,theta,potential_absorbed),rel=a1/(u*u);
    const double value_scale=1.0+(std::abs(rel)<=0.01?rel:0.0);
    const double derivative_scale=std::pow(u,-4.0/3.0)*b0;
    const double ca=value_scale*ab.ai_scaled-derivative_scale*ab.aip_scaled;
    const double cb=value_scale*ab.bi_scaled-derivative_scale*ab.bip_scaled;
    if(ab.exponent>30.0)return {};
    const double pref=std::sqrt(2.0)*std::pow(u,1.0/6.0)*std::pow(map.zeta_over_q,0.25);
    const double a=pref*ca*std::exp(-ab.exponent),b=pref*cb*std::exp(ab.exponent);
    return {std::isfinite(a)&&std::isfinite(b),a,b};
}

inline DirectOscillatoryPair empirically_calibrated_airy_pair(double nu,double mu,double theta,
                                                             double* p_derivative=nullptr,
                                                             double* q_derivative=nullptr){
    const double u=nu+0.5,alpha=mu/u;
    const bool high_carrier=nu>10000.0;
    // A tenth-order stencil on the local Airy scale avoids both the
    // carrier-scale cancellation and the truncation of a wide sixth-order one.
    const double h=high_carrier ? 0.08*std::cbrt(mu/(2.0*std::sqrt((1.0-alpha)*(1.0+alpha))))/u : 0.02/u;
    const int radius=high_carrier?5:3;
    std::array<AiryBasisValue,11> v;
    for(int j=-radius;j<=radius;++j){v[j+5]=airy_basis_value(nu,mu,theta+j*h);if(!v[j+5].ok)return {};}
    auto derivative=[&](bool second){
        auto f=[&](int j){return second?v[j+5].b:v[j+5].a;};
        if(high_carrier) return ((5.0/6.0)*(f(1)-f(-1))-(5.0/21.0)*(f(2)-f(-2))
            +(5.0/84.0)*(f(3)-f(-3))-(5.0/504.0)*(f(4)-f(-4))
            +(1.0/1260.0)*(f(5)-f(-5)))/h;
        return (-f(-3)+9*f(-2)-45*f(-1)+45*f(1)-9*f(2)+f(3))/(60*h);
    };
    const double ad=derivative(false),bd=derivative(true);
    double wa=v[5].a*bd-ad*v[5].b;
    const auto center_map=canonical_map(alpha,theta,std::fma(-alpha,u,mu)/u);
    if(high_carrier && std::abs(center_map.zeta)<0.01 && center_map.zeta!=0.0) {
        // The analytic Airy Wronskian removes differencing noise from the
        // value calibration. Differentiate the same polynomial transport
        // coefficients used in this core; no sampled flux enters this scale.
        const double zp=std::sqrt(center_map.q/center_map.zeta);
        const double u23=std::pow(u,2.0/3.0),arg=-u23*center_map.zeta;
        const double r0=airy_residual_r0(alpha),r1=airy_residual_r1(alpha),r2=airy_residual_r2(alpha);
        const double zeta=center_map.zeta;
        const double residual=r0+r1*zeta+r2*zeta*zeta;
        const double b0=r0+r1*zeta/3.0+r2*zeta*zeta/5.0;
        const double a1=airy_transport_a1(alpha,theta);
        const bool applied=std::abs(a1/(u*u))<=0.01;
        const double a=1.0+(applied?a1/(u*u):0.0);
        const double b=std::pow(u,-4.0/3.0)*b0;
        const double at=applied ? -0.5*(2.0*r2/5.0+residual*b0)*zp/(u*u) : 0.0;
        const double bt=std::pow(u,-4.0/3.0)*(r1/3.0+2.0*r2*zeta/5.0)*zp;
        wa=(-2.0*u/pi)*(a*a-arg*b*b+(b*at-a*bt)/(-u23*zp));
    }
    const PiOverTwoData p0=normalized_p_at_pi_over_two(nu,mu);
    const CenterQData q0=normalized_q_at_pi_over_two(nu,mu);
    if(!p0.ok||!q0.ok||wa==0.0||!std::isfinite(wa))return {};
    const double we=(static_cast<double>(p0.value)*static_cast<double>(q0.derivative)
                    -static_cast<double>(p0.derivative)*static_cast<double>(q0.value));
    if(we/wa<=0.0||!std::isfinite(we))return {};
    // Empirical value calibration retained for compatibility with the
    // independently tested regular solution.  Scaling both columns by
    // we/wa does NOT itself impose the exact Wronskian.
    const double scale=we/wa;
    const double p=scale*v[5].a,q=scale*v[5].b;
    if(!std::isfinite(p)||!std::isfinite(q))return {};
    if(p_derivative)*p_derivative=scale*ad;
    if(q_derivative)*q_derivative=scale*bd;
    return {true,p,q,0.0,0.0};
}

// Preserve the calibrated values.  Degree recurrence supplies derivatives
// independently of the noisy finite-difference flux used by the value
// calibration; a minimum-norm derivative correction enforces the known
// Wronskian at this evaluation point.  This invariant is a consistency
// condition, not a bound on the Airy approximation's function error.
inline DirectOscillatoryPair wronskian_normalized_airy_pair(double nu,double mu,double theta,
                                                             double* p_derivative=nullptr,
                                                             double* q_derivative=nullptr) {
    double pd=0.0,qd=0.0;
    const auto value=empirically_calibrated_airy_pair(nu,mu,theta,&pd,&qd);
    if(!value.ok || (!p_derivative && !q_derivative)) return value;
    if(nu<=10000.0) {
        const auto previous=empirically_calibrated_airy_pair(nu-1.0,mu,theta);
        if(!previous.ok || !(nu>mu)) return {};
        const double u=nu+0.5, sine=std::sin(theta), cosine=std::cos(theta);
        const double coupling=std::sqrt((nu-mu)*(nu+mu)*u/(u-1.0));
        pd=((nu+0.5)*cosine*value.p-coupling*previous.p)/sine;
        qd=((nu+0.5)*cosine*value.q-coupling*previous.q)/sine;
    }
    const auto p0=normalized_p_at_pi_over_two(nu,mu);
    const auto q0=normalized_q_at_pi_over_two(nu,mu);
    if(!p0.ok || !q0.ok) return {};
    const double expected=static_cast<double>(p0.value*q0.derivative-p0.derivative*q0.value);
    const double denominator=value.p*value.p+value.q*value.q;
    if(!(denominator>0.0) || !std::isfinite(denominator)) return {};
    const double correction=(expected-std::fma(value.p,qd,-pd*value.q))/denominator;
    pd-=value.q*correction;
    qd+=value.p*correction;
    if(!std::isfinite(pd) || !std::isfinite(qd)) return {};
    if(p_derivative) *p_derivative=pd;
    if(q_derivative) *q_derivative=qd;
    return value;
}

// Continue the independently normalized turning-point solution with the
// original Legendre ODE.  The integration starts at a fixed O(1) Airy action,
// so its cost is proportional to the local phase distance from the turning
// point, not to nu.  This removes the phase drift of a truncated outer Airy
// series precisely where a component is close to a zero.
inline DirectOscillatoryPair locally_continued_airy_pair(double nu,double mu,
                                                          double target_theta,
                                                          double* p_derivative=nullptr){
    const double u=nu+0.5,alpha=mu/u;
    const CanonicalMap target=canonical_map(alpha,target_theta,std::fma(-alpha,u,mu)/u);
    if(target.region!=Region::oscillatory)return {};
    const double seed_scaled_action=nu>10000.0 ? 1.0 : 0.002;
    const double target_scaled_action=u*target.action;
    if(!(target_scaled_action>seed_scaled_action))return {};
    const double seed_action=seed_scaled_action/u;
    const double seed_zeta=std::pow(1.5*seed_action,2.0/3.0);
    const double seed_theta=theta_from_zeta(alpha,seed_zeta,target_theta);
    double pd=0.0,qd=0.0;
    DirectOscillatoryPair seed=
        wronskian_normalized_airy_pair(nu,mu,seed_theta,&pd,&qd);
    if(!seed.ok||!std::isfinite(pd)||!std::isfinite(qd))return {};

    // Replace the approximate recessive Airy datum by the convergent regular
    // endpoint solution.  Q remains supplied by the turning chart, where its
    // endpoint connection formula would be ill-conditioned.  Finally impose
    // the exact global Wronskian on this hybrid pair.
    const double endpoint_argument=2.0/u;
    if(mu<=128.0 && endpoint_argument>0.0&&endpoint_argument<1.0){
        const double endpoint_theta=2.0*std::asin(endpoint_argument);
        if(endpoint_theta<seed_theta){
            const EndpointScalar ep=endpoint_normalized_p_with_derivative(
                nu,mu,endpoint_theta);
            if(ep.ok&&ep.value.sign!=0){
                const double ds=ep.derivative.sign==0?0.0:
                    static_cast<double>(ep.derivative.sign)
                    *std::exp(ep.derivative.log_abs-ep.value.log_abs)/u;
                const ScaledOdeValue exact_p=propagate_scaled(u,mu,endpoint_theta,
                    seed_theta,static_cast<double>(ep.value.sign),ds,ep.value.log_abs);
                if(exact_p.ok&&std::abs(exact_p.log_scale)<700.0){
                    const double factor=std::exp(exact_p.log_scale);
                    const double p=exact_p.value*factor;
                    const double pder=u*exact_p.derivative_over_u*factor;
                    const PiOverTwoData p0=normalized_p_at_pi_over_two(nu,mu);
                    const CenterQData q0=normalized_q_at_pi_over_two(nu,mu);
                    const double exact_w=(p0.ok&&q0.ok)?
                        static_cast<double>(p0.value*q0.derivative-p0.derivative*q0.value):0.0;
                    const double hybrid_w=p*qd-pder*seed.q;
                    if(exact_w/hybrid_w>0.0&&std::isfinite(exact_w/hybrid_w)){
                        const double qscale=exact_w/hybrid_w;
                        seed.p=p;pd=pder;seed.q*=qscale;qd*=qscale;
                    }
                }
            }
        }
    }
    std::array<double,4> y{seed.p,pd/u,seed.q,qd/u};
    auto rhs=[&](double theta,const std::array<double,4>& state){
        const double s=std::sin(theta);
        const double potential=u*u+(0.25-mu*mu)/(s*s);
        return std::array<double,4>{u*state[1],-potential*state[0]/u,
                                    u*state[3],-potential*state[2]/u};
    };
    double theta=seed_theta;
    double h=std::min(0.02/u,std::max(1.0e-14,(target_theta-theta)/20.0));
    constexpr double rtol=2.0e-14,atol=2.0e-16;
    for(int step=0;step<200000&&theta<target_theta;++step){
        if(theta+h>target_theta)h=target_theta-theta;
        const auto k1=rhs(theta,y);auto stage=y;
        for(int i=0;i<4;++i)stage[i]=y[i]+h*k1[i]/5.0;
        const auto k2=rhs(theta+h/5.0,stage);
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(3*k1[i]/40.0+9*k2[i]/40.0);
        const auto k3=rhs(theta+3*h/10.0,stage);
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(44*k1[i]/45.0-56*k2[i]/15.0+32*k3[i]/9.0);
        const auto k4=rhs(theta+4*h/5.0,stage);
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(19372*k1[i]/6561.0-25360*k2[i]/2187.0+64448*k3[i]/6561.0-212*k4[i]/729.0);
        const auto k5=rhs(theta+8*h/9.0,stage);
        for(int i=0;i<4;++i)stage[i]=y[i]+h*(9017*k1[i]/3168.0-355*k2[i]/33.0+46732*k3[i]/5247.0+49*k4[i]/176.0-5103*k5[i]/18656.0);
        const auto k6=rhs(theta+h,stage);
        std::array<double,4> fifth{};
        for(int i=0;i<4;++i)fifth[i]=y[i]+h*(35*k1[i]/384.0+500*k3[i]/1113.0+125*k4[i]/192.0-2187*k5[i]/6784.0+11*k6[i]/84.0);
        const auto k7=rhs(theta+h,fifth);double error=0.0;
        for(int i=0;i<4;++i){
            const double fourth=y[i]+h*(5179*k1[i]/57600.0+7571*k3[i]/16695.0+393*k4[i]/640.0-92097*k5[i]/339200.0+187*k6[i]/2100.0+k7[i]/40.0);
            const double scale=atol+rtol*std::max(std::abs(y[i]),std::abs(fifth[i]));
            error=std::max(error,std::abs(fifth[i]-fourth)/scale);
        }
        if(error<=1.0){theta+=h;y=fifth;}
        h*=error==0.0?5.0:RSCHACompat::clamp(0.9*std::pow(error,-0.2),0.2,5.0);
        if(std::abs(h)<2.0e-16)return {};
        if(theta>=target_theta){
            if(!std::isfinite(y[0])||!std::isfinite(y[2]))return {};
            if(p_derivative) *p_derivative=u*y[1];
            return {true,y[0],y[2],0.0,0.0};
        }
    }
    return {};
}

// In the moderate-order overlap, initialize at the rigorously decreasing
// pole-residual series and propagate backward to the Airy core.  This fixes
// the one remaining connection constant without importing any external
// evaluator.  The order cap bounds the O(mu^2) phase distance and runtime.
inline DirectOscillatoryPair pole_continued_pair(double nu,double mu,
                                                  double target_theta,
                                                  double match_parameter=0.55){
    if(!(mu>=10.0&&mu<=100.0))return {};
    const double u=nu+0.5;
    const double s=(mu*mu+0.25)/(2.0*u*match_parameter);
    if(!(s>0.0&&s<1.0))return {};
    const double match=std::asin(s);
    if(!(match>target_theta&&match<0.5*pi))return {};
    double pd=0.0,qd=0.0;
    const DirectOscillatoryPair seed=pole_residual_pair(nu,mu,match,&pd,&qd);
    if(!seed.ok||!std::isfinite(pd)||!std::isfinite(qd))return {};
    const ScaledOdeValue p=propagate_scaled(u,mu,match,target_theta,
        seed.p,pd/u,0.0,5.0e-14,5.0e-16);
    const ScaledOdeValue q=propagate_scaled(u,mu,match,target_theta,
        seed.q,qd/u,0.0,5.0e-14,5.0e-16);
    if(!p.ok||!q.ok)return {};
    double qvalue=q.value,qlog=q.log_scale;
    const PiOverTwoData p0=normalized_p_at_pi_over_two(nu,mu);
    const CenterQData q0=normalized_q_at_pi_over_two(nu,mu);
    const double numerical_w=u*(p.value*q.derivative_over_u
        -p.derivative_over_u*q.value)*std::exp(p.log_scale+q.log_scale);
    const double exact_w=(p0.ok&&q0.ok)?static_cast<double>(
        p0.value*q0.derivative-p0.derivative*q0.value):0.0;
    const double flux_scale=exact_w/numerical_w;
    if(flux_scale>0.0&&std::isfinite(flux_scale))qlog+=std::log(flux_scale);
    return {true,p.value,qvalue,p.log_scale,qlog};
}

// Exact nonoscillatory envelope/phase equations.  In x=u(theta-pi/2), write
// P+iQ=exp(a+i(phi0+sigma*x+delta)).  The Ermakov system is
// a'=b, b'=w^2-k-b^2, delta'=w-sigma, w=C*exp(-2a).
// It contains no carrier oscillations, so its step count follows variation of
// k rather than the degree.  This chart is used only where both components
// remain well conditioned; recessive forbidden-region P continues to come
// from the endpoint chart.
inline DirectOscillatoryPair phase_amplitude_pair(double nu,double mu,double target){
    const double u=nu+0.5,c=0.5*pi;
    const PiOverTwoData p0=normalized_p_at_pi_over_two(nu,mu);
    const CenterQData q0=normalized_q_at_pi_over_two(nu,mu);
    if(!p0.ok||!q0.ok)return {};
    const double p=static_cast<double>(p0.value),q=static_cast<double>(q0.value);
    const double pd=static_cast<double>(p0.derivative),qd=static_cast<double>(q0.derivative);
    const double lambda=0.25-mu*mu,k0=1.0+lambda/(u*u);
    if(!(k0>0.0))return {};
    // Select the nonoscillatory Kummer branch at the symmetry point.  Its
    // logarithmic-amplitude derivative is exactly zero there.  The displayed
    // k'' correction is the first center WKB correction to the phase speed.
    const double w0=std::sqrt(k0-lambda/(2.0*std::pow(u,4)*k0));
    const double a0=-0.5*std::log(w0),rho0=std::exp(a0);
    constexpr double cx=1.0,sigma=1.0;
    std::array<double,3> y{a0,0.0,0.0};
    double x=0.0,xt=u*(target-c),h=std::copysign(std::min(100.0,std::abs(xt)/20.0),xt);
    auto rhs=[&](double xx,const std::array<double,3>&v){
        const double th=c+xx/u,s=std::sin(th);
        const double k=1.0+(0.25-mu*mu)/(u*u*s*s);
        const double w=cx*std::exp(-2.0*v[0]);
        return std::array<double,3>{v[1],w*w-k-v[1]*v[1],w-sigma};
    };
    constexpr double rtol=2e-12,atol=2e-12;
    for(int step=0;step<500000;++step){
        if((xt-x)*h<=0.0||std::abs(h)>std::abs(xt-x))h=xt-x;
        if(h==0.0)break;
        const auto k1=rhs(x,y);auto s=y;
        for(int i=0;i<3;i++)s[i]=y[i]+h*k1[i]/5;const auto k2=rhs(x+h/5,s);
        for(int i=0;i<3;i++)s[i]=y[i]+h*(3*k1[i]/40+9*k2[i]/40);const auto k3=rhs(x+3*h/10,s);
        for(int i=0;i<3;i++)s[i]=y[i]+h*(44*k1[i]/45-56*k2[i]/15+32*k3[i]/9);const auto k4=rhs(x+4*h/5,s);
        for(int i=0;i<3;i++)s[i]=y[i]+h*(19372*k1[i]/6561-25360*k2[i]/2187+64448*k3[i]/6561-212*k4[i]/729);const auto k5=rhs(x+8*h/9,s);
        for(int i=0;i<3;i++)s[i]=y[i]+h*(9017*k1[i]/3168-355*k2[i]/33+46732*k3[i]/5247+49*k4[i]/176-5103*k5[i]/18656);const auto k6=rhs(x+h,s);
        std::array<double,3> fifth{};for(int i=0;i<3;i++)fifth[i]=y[i]+h*(35*k1[i]/384+500*k3[i]/1113+125*k4[i]/192-2187*k5[i]/6784+11*k6[i]/84);
        const auto k7=rhs(x+h,fifth);double error=0;
        for(int i=0;i<3;i++){double fourth=y[i]+h*(5179*k1[i]/57600+7571*k3[i]/16695+393*k4[i]/640-92097*k5[i]/339200+187*k6[i]/2100+k7[i]/40);
            double scale=atol+rtol*std::max(std::abs(y[i]),std::abs(fifth[i]));error=std::max(error,std::abs(fifth[i]-fourth)/scale);}
        if(error<=1){x+=h;y=fifth;}
        h*=error==0?5.0:RSCHACompat::clamp(0.9*std::pow(error,-0.2),0.2,5.0);
        if(std::abs(h)<1e-12||!std::isfinite(y[0])||!std::isfinite(y[1])||!std::isfinite(y[2]))return {};
    }
    if(std::abs(x-xt)>1e-8*std::max(1.0,std::abs(xt)))return {};
    const double carrier=sigma*xt;
    const double hi=carrier,lo=std::fma(sigma,xt,-hi);
    constexpr double twopi_hi=6.283185307179586232,twopi_lo=2.449293598294706414e-16;
    const double ncycle=std::nearbyint(hi/twopi_hi),prod=ncycle*twopi_hi;
    const double prodlo=std::fma(ncycle,twopi_hi,-prod)+ncycle*twopi_lo;
    const double phase=(hi-prod)+(lo-prodlo)+y[2];
    const double amp=std::exp(y[0]);
    const double cs=std::cos(phase),sn=std::sin(phase);
    const double Ap=p/rho0,Bp=pd/(u*w0*rho0);
    const double Aq=q/rho0,Bq=qd/(u*w0*rho0);
    const double pv=amp*(Ap*cs+Bp*sn),qv=amp*(Aq*cs+Bq*sn);
    if(!std::isfinite(pv)||!std::isfinite(qv))return {};
    return {true,pv,qv,0.0,0.0};
}

inline SignedLogValue make_signed_log(double scaled, double log_scale) {
    SignedLogValue out;
    if (!std::isfinite(scaled) || !std::isfinite(log_scale)) return out;
    if (scaled == 0.0) {
        out.value = 0.0;
        out.log_abs = -std::numeric_limits<double>::infinity();
        out.sign = 0;
        out.value_representable = true;
        return out;
    }
    out.sign = std::signbit(scaled) ? -1 : 1;
    out.log_abs = std::log(std::abs(scaled)) + log_scale;
    constexpr double log_max = 709.782712893383973096;
    constexpr double log_min_subnormal = -744.44007192138121809;
    if (out.log_abs <= log_max && out.log_abs >= log_min_subnormal) {
        out.value = static_cast<double>(out.sign) * std::exp(out.log_abs);
        out.value_representable = true;
    } else {
        out.value = out.log_abs > log_max
            ? static_cast<double>(out.sign) * std::numeric_limits<double>::infinity()
            : static_cast<double>(out.sign) * 0.0;
        out.value_representable = false;
    }
    return out;
}

}

// Canonically normalized global Airy chart for P_tilde and Q_tilde.
// Current separation conditions avoid the pole-turning and turning-endpoint
// coalescence regimes.  They are explicit mathematical chart conditions, not
// a fallback to another evaluator.
inline Result evaluate(double nu, double mu, double theta, bool with_derivative=false) {
    Result out;
    if (!std::isfinite(nu) || !std::isfinite(mu) || !std::isfinite(theta)) {
        out.status = Status::nonfinite_input;
        return out;
    }
    if (nu < 0.0 || nu>1.0e8 || mu < 0.0 || mu > nu || theta <= 0.0 || theta >= pi) {
        out.status = Status::outside_real_ferrers_domain;
        return out;
    }

    const double u = nu + 0.5;
    const double alpha = mu / u;

    const double exact_turning_ratio=mu>0.5?std::sqrt(mu*mu-0.25)/u:0.0;
    out.theta_turning = exact_turning_ratio<1.0?std::asin(exact_turning_ratio):0.5*pi;

    // Direct chart for the complete low/mid-degree domain.
    // It precedes all Airy guards and therefore also covers mu=0, nu<10 and
    // the reflected half interval.  Only stable-direction endpoint P is used
    // in the left forbidden region; all other cases use scaled center data.
    if(nu<=1000.0){
        const bool left_evan=mu>0.5&&exact_turning_ratio<1.0&&theta<out.theta_turning;
        const detail::DirectOscillatoryPair direct=left_evan
            ?detail::direct_evanescent_pair(nu,mu,theta)
            :detail::direct_center_scaled_pair(nu,mu,theta);
        if(direct.ok){

            out.p_tilde=detail::make_signed_log(direct.p,direct.p_log_scale);
            out.q_tilde=detail::make_signed_log(direct.q,direct.q_log_scale);

            out.status=Status::ok;
            return out;
        }
    }

    // Exact Gamma anchor with only a short O(1/u) propagation.  This removes
    // avoidable Airy phase and normalization error around theta=pi/2.
    if(std::abs(theta-0.5*pi)*u<=0.25){
        const detail::DirectOscillatoryPair center=
            detail::direct_center_scaled_pair(nu,mu,theta);
        if(center.ok){

            out.p_tilde=detail::make_signed_log(center.p,center.p_log_scale);
            out.q_tilde=detail::make_signed_log(center.q,center.q_log_scale);

            out.status=Status::ok;
            return out;
        }
        out.status=Status::outside_current_airy_chart;
        return out;
    }

    // In the Weber coalescence layer the useful interval around pi/2 is
    // narrow in theta even at very high degree.  A budgeted exact-center
    // propagation is faster and more accurate than applying a nonuniform
    // Airy expansion there.  The 3000 carrier-unit cap bounds the work.
    if((theta>=out.theta_turning||u*(out.theta_turning-theta)<=5.0)
       &&u*(1.0-alpha)<=500.0
       &&u*std::abs(theta-0.5*pi)<=7000.0){
        const detail::DirectOscillatoryPair center=
            detail::direct_center_scaled_pair(nu,mu,theta);
        if(center.ok){

            out.p_tilde=detail::make_signed_log(center.p,center.p_log_scale);
            out.q_tilde=detail::make_signed_log(center.q,center.q_log_scale);

            out.status=Status::ok;
            return out;
        }
    }

    // Regular endpoint series.  The conservative u*theta bound makes the
    // hypergeometric terms rapidly decreasing even at degree one million.
    const double endpoint_series_parameter=
        std::pow(u*std::sin(0.5*theta),2)/(1.0+mu);
    if(mu==0.0 && endpoint_series_parameter<=8.0){
        const detail::DirectOscillatoryPair endpoint=
            detail::endpoint_zero_order_pair(nu,theta);
        if(endpoint.ok){

            out.p_tilde=detail::make_signed_log(endpoint.p,0.0);
            out.q_tilde=detail::make_signed_log(endpoint.q,0.0);

            out.status=Status::ok;
            return out;
        }
    }
    if(endpoint_series_parameter<=8.0){
        const detail::DirectOscillatoryPair endpoint=
            detail::direct_endpoint_pair(nu,mu,theta);
        if(endpoint.ok){

            out.p_tilde=detail::make_signed_log(endpoint.p,endpoint.p_log_scale);
            out.q_tilde=detail::make_signed_log(endpoint.q,endpoint.q_log_scale);

            out.status=Status::ok;
            return out;
        }
        // Near an integer order the endpoint connection needs its logarithmic
        // Frobenius companion.  If that chart declines the point, allow the
        // separated residual/Airy charts below to try it.
    }
    const double extended_endpoint_limit=std::min(20.0,2.0+0.08*mu);
    if(endpoint_series_parameter>8.0
       &&endpoint_series_parameter<=extended_endpoint_limit){
        const detail::DirectOscillatoryPair extended=
            detail::endpoint_propagated_pair(nu,mu,theta);
        if(extended.ok){

            out.p_tilde=detail::make_signed_log(extended.p,extended.p_log_scale);
            out.q_tilde=detail::make_signed_log(extended.q,extended.q_log_scale);

            out.status=Status::ok;
            return out;
        }
    }

    if(theta>0.5*pi){
        out.status=Status::outside_current_airy_chart;
        return out;
    }

    // Small-order/high-degree pole chart.  The expansion parameter controls
    // the canonical residual itself; points outside it remain assigned to the
    // turning-point charts below.
    const double pole_parameter=(mu*mu+0.25)/(2.0*u*std::sin(theta));
    if(theta>out.theta_turning && u*std::sin(theta)>=16.0
       && pole_parameter<=(mu<=33.0?8.0:2.6)){
        const detail::DirectOscillatoryPair pole=
            detail::pole_residual_pair(nu,mu,theta);
        if(pole.ok){

            out.p_tilde=detail::make_signed_log(pole.p,0.0);
            out.q_tilde=detail::make_signed_log(pole.q,0.0);

            out.status=Status::ok;
            return out;
        }
    }

    // An Airy chart requires both separations to be resolved.
    if (nu < 10.0 || alpha >= 1.0 || mu < 10.0
        || u * (1.0 - alpha) < 10.0) {
        out.status = Status::outside_current_airy_chart;
        return out;
    }

    const detail::CanonicalMap map = detail::canonical_map(alpha, theta,std::fma(-alpha,u,mu)/u);

    if (!(map.zeta_over_q > 0.0) || !std::isfinite(map.zeta_over_q)) {
        out.status = Status::airy_kernel_failure;
        return out;
    }

    if(u*std::abs(map.action)<=20.0){
        double pd=0.0;
        const detail::DirectOscillatoryPair normalized=
            detail::wronskian_normalized_airy_pair(nu,mu,theta,with_derivative?&pd:nullptr);
        if(normalized.ok){
            out.p_tilde=detail::make_signed_log(normalized.p,0.0);
            if(with_derivative) {out.p_theta=detail::make_signed_log(pd,0.0);out.p_theta_available=true;}
            out.q_tilde=detail::make_signed_log(normalized.q,0.0);

            out.status=Status::ok;
            return out;
        }
    }

    // Outside the Airy core, retain exact phase by a short continuation.  The
    // cap bounds runtime; farther out the outer expansion is used.
    if(map.region==Region::oscillatory&&u*map.action<=1000.0){
        const detail::DirectOscillatoryPair pole_continued=
            detail::pole_continued_pair(nu,mu,theta);
        if(pole_continued.ok){
            // P is recessive at the endpoint and is more accurately supplied
            // by the regular endpoint/Airy hybrid; Q is the component whose
            // Airy connection constant benefits from pole continuation.
            const detail::DirectOscillatoryPair endpoint_hybrid=
                detail::locally_continued_airy_pair(nu,mu,theta);
            const detail::DirectOscillatoryPair& psource=
                endpoint_hybrid.ok&&std::abs(pole_continued.p)<0.1
                    ?endpoint_hybrid:pole_continued;
            out.p_tilde=detail::make_signed_log(psource.p,psource.p_log_scale);
            out.q_tilde=detail::make_signed_log(pole_continued.q,pole_continued.q_log_scale);

            out.status=Status::ok;
            return out;
        }
        if(u*map.action<=200.0){
            double pd=0.0;
            const detail::DirectOscillatoryPair continued=
                detail::locally_continued_airy_pair(nu,mu,theta,with_derivative?&pd:nullptr);
            if(continued.ok){
                out.p_tilde=detail::make_signed_log(continued.p,0.0);
                if(with_derivative) {out.p_theta=detail::make_signed_log(pd,0.0);out.p_theta_available=true;}
                out.q_tilde=detail::make_signed_log(continued.q,0.0);

                out.status=Status::ok;
                return out;
            }
        }
    }

    // The direct chart is universal through degree 1000.  A narrow extension
    // covers the small-alpha boundary layer where the Airy transport
    // coefficients converge unusually slowly, without imposing O(nu) work
    // throughout the full high-degree domain.
    if (nu <= 1000.0 || (nu <= 3000.0 && alpha <= 0.03)) {
        const detail::DirectOscillatoryPair direct = map.region==Region::evanescent
            ? detail::direct_evanescent_pair(nu,mu,theta)
            : detail::direct_oscillatory_pair(nu,mu,theta);
        if (direct.ok) {
            out.p_tilde = detail::make_signed_log(direct.p, direct.p_log_scale);
            out.q_tilde = detail::make_signed_log(direct.q, direct.q_log_scale);

            out.status = Status::ok;
            return out;
        }
    }

    const double u_two_thirds = std::pow(u, 2.0 / 3.0);
    const double airy_argument = -u_two_thirds * map.zeta;
    const detail::AiryScaled ab = detail::airy(airy_argument, u * map.action);
    const double transport_b0 = detail::airy_transport_b0(alpha, theta, map);
    const double transport_a1 = detail::airy_transport_a1(alpha, theta);
    const double derivative_scale = std::pow(u, -4.0 / 3.0) * transport_b0;
    // Asymptotic-series truncation rule: do not append A1 once that term has
    // ceased to be perturbative in an unresolved coalescence layer.
    const double a1_relative = transport_a1 / (u * u);
    const bool apply_a1 = std::abs(a1_relative) <= 0.01;
    const double value_scale = 1.0
        + (apply_a1 ? a1_relative : 0.0);
    const double corrected_ai = value_scale * ab.ai_scaled
                              - derivative_scale * ab.aip_scaled;
    const double corrected_bi = value_scale * ab.bi_scaled
                              - derivative_scale * ab.bip_scaled;
    const double prefactor = std::sqrt(2.0) * std::pow(u, 1.0 / 6.0)
                           * std::pow(map.zeta_over_q, 0.25);

    out.p_tilde = detail::make_signed_log(prefactor * corrected_ai, -ab.exponent);
    out.q_tilde = detail::make_signed_log(prefactor * corrected_bi,  ab.exponent);

    if (!std::isfinite(out.p_tilde.log_abs) || !std::isfinite(out.q_tilde.log_abs)) {
        out.status = Status::airy_kernel_failure;
        return out;
    }
    out.status = Status::ok;
    return out;
}

}
