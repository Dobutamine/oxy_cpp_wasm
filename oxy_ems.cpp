/*

Compile to WASM with Emscripten for use in web applications using:

emcc oxy_ems.cpp -o oxy.js -O3 -s ENVIRONMENT="web,worker" --std=c++17  --bind 

-o specifies the outfile name. The extension .js is important. Emscripten generates a javascript glue file and a wasm file.
-O3 = maximal optimization (O0-O3, Oz = minimal size, Os=optimized and minimal size)
-s ENVIRONMENT="web,worker” is also optional as the default is ENVIRONMENT = "web,webview,worker,node”. It omits code for web view and node environments and makes the glue code size smaller.
—std=c++17 tells the c++ version
—bind is important as it provides code for the bindings of c++ datatypes with javascript datatypes
-s WASM=1 -s MODULARIZE=1 -s EXPORT_ES6=1  are optional as the default for all s options = 1

*/

#include <cmath>
#include <iostream>
#include <stdexcept>

// include emscripten dependencies
#include <emscripten.h>
#include <emscripten/bind.h>

struct O2Result {
    double PO2;   // Partial pressure of O2 in mmHg
    double SaO2;  // O2 saturation fraction (0-1)
    int counter; // Iteration count
};

/**
 * Calculate arterial PO2 and O2 saturation from blood gas parameters in SI units.
 * Uses Newton-Raphson for fast convergence.
 *
 * @param pH          Blood pH
 * @param pCO2        Arterial PCO2 in mmHg
 * @param temperature Blood temperature in °C
 * @param dpg         2,3-DPG concentration in mmol/L
 * @param TO2         Total O2 content in mmol O2/L
 * @param Hb          Hemoglobin concentration in mmol Hb/L
 * @param tol         Convergence tolerance (default 1e-6)
 * @param max_iter    Maximum iterations for Newton-Raphson (default 50)
 * @return            Struct containing PO2 (mmHg) SaO2 (fraction) and iteration count
 */


O2Result calculateO2Parameters(double pH,
                                double pCO2,
                                double temperature,
                                double dpg,
                                double TO2,
                                double Hb,
                                double prev_PO2 = -1.0,
                                double P50_0 = 26.8,
                                double tol = 1e-8,
                                int max_iter = 50) {

    // physical constants
    const double n = 2.7;             // Hill coefficient
    const double alpha = 1.38e-5;     // O2 solubility in mmol/L/mmHg

    // Compute shifts from standard conditions
    double dpH = pH - 7.40;             //Bohr effect: ↓pH → right shift → ↑P₅₀)
    double dpCO2 = pCO2 - 40.0;         // Haldane effect: ↑pCO2 → right shift → ↑P₅₀
    double dT = temperature - 37.0;     // ↑T → right shift → ↑P₅₀
    double dDPG = dpg - 5.0;            // ↑DPG → right shift → ↑P₅₀

    // Adjust P50 on log10 scale
    double log10_p50 = std::log10(P50_0)
                      - 0.48 * dpH
                      + 0.014 * dpCO2
                      + 0.024 * dT
                      + 0.051 * dDPG;
    double P50 = std::pow(10.0, log10_p50);
    double P50_n = std::pow(P50, n);

    // Initial PO2 guess is the P50
    
    double po2 = P50;
    // if the previous PO2 is not zero then use it as an educated guess
    if (prev_PO2 > 0.0) {
        po2 = prev_PO2;
    } 
    double po2 = P50;
    double sat = 0.0;

    int counter = 0;
    // Newton-Raphson iteration
    for (int i = 0; i < max_iter; ++i) {
        counter++;
        double po2_n = std::pow(po2, n);
        double denom = po2_n + P50_n;
        sat = po2_n / denom;
        double f = Hb * sat + alpha * po2 - TO2;
        if (std::abs(f) < tol) break;
        double df_sat = n * P50_n * std::pow(po2, n - 1) / (denom * denom);
        double df = Hb * df_sat + alpha;
        po2 -= f / df;
        if (po2 < tol) {
            po2 = tol;
            break;
        }
    }

    return {po2, sat, counter};
}

// set the emscripten bindings
EMSCRIPTEN_BINDINGS(my_module) {
    emscripten::function("calculateO2Parameters", &calculateO2Parameters);
    emscripten::value_object<O2Result>("O2Result")
        .field("PO2", &O2Result::PO2)
        .field("SaO2", &O2Result::SaO2)
        .field("counter", &O2Result::counter);
}