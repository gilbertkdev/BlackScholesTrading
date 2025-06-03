#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <time.h>

#define DAYS 252

// Function
double ln(double x);
double sqr(double x);
double cdf(double x);

double BlackScholesCall(float S, float K, float sigma, float rf, float t);
double BrownianMotion(double S, double mu, double sigma, double dt, double path[]);
double BoxMuller();


int main() {
    const double S = 100.0;         // Current Underlying Price
    const double K = 100.0;         // Strike Price
    const double sigma = 0.3;       // Volatility (σ)
    const double rf = 0.03;         // Risk-Free Interest Rate
    const double T = 1.0;           // Time to Maturity

    const double mu = 0.05;         // Expected return
    const double dt = 1.0 / DAYS;   // Time Step (1 Trading Day)
    double path[DAYS];

    srand(time(NULL));

    double CallOptionPrice = BlackScholesCall(S, K, sigma, rf, T);
    BrownianMotion(S, mu, sigma, dt, path);


    FILE *fout = fopen("BrownianMotion.txt", "w");

    for (int i = 0; i < DAYS; i++) {
        fprintf(fout, "%d %.5f\n", i, path[i]);
    }

    fclose(fout);
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

// Brownian Motion
// Z = Standard Normal Random Variable
double BrownianMotion(double S, double mu, double sigma, double dt, double path[]) {
    path[0] = S;

    for (int i = 1; i < DAYS; i++) {
        double Z = BoxMuller();
        double dW = sqrt(dt) * Z;

        S *= exp((mu - 0.5 * sigma * sigma) * dt + sigma * dW);

        path[i] = S;
    }

    return S;
}

// Box-Muller
double BoxMuller() {
    double U1 = (double)rand() / RAND_MAX;
    double U2 = (double)rand() / RAND_MAX;

    return sqrt(-2.0 * log(U1)) * cos(2.0 * M_PI * U2);
}