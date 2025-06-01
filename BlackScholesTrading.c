#include <stdio.h>
#include <math.h>


// Function
double ln(double x);
double sqr(double x);
double cdf(double x);

double BlackScholesCall(float S, float K, float sigma, float rf, float t);


int main() {
    const double S = 100.0;     // Current Underlying Price
    const double K = 100.0;     // Strike Price
    const double sigma = 0.3;   // Volatility (σ)
    const double rf = 0.03;     // Risk-Free Interest Rate
    const double T = 1.0;       // Time to Maturity

    double CallOptionPrice = BlackScholesCall(S, K, sigma, rf, T);
    
    
    printf("Black-Scholes Call Option Price: %.5f\n", CallOptionPrice);
}


// Math Function
double ln(double x) {
    return log(x);
}

double sqr(double x) {
    return x * x;
}

double cdf(double x) {
    return 0.5 * (1 + erf(x * M_SQRT1_2));
}


// Black-Scholes Call
// C = Call Option Price
// S = Current Underlying Price
// K = Strike Price
// sigma = Volatility (σ)
// rf = Risk-Free Interest Rate
// t = Time to Maturity
double BlackScholesCall(float S, float K, float sigma, float rf, float t) {
    double d1 = (ln(S / K) + ((rf + (sqr(sigma) / 2)) * t)) / (sigma * sqrt(t));
    double d2 = d1 - (sigma * sqrt(t));

    double C = S * cdf(d1) - exp(-rf * t) * K * cdf(d2);

    return C;
}