#pragma once

#include <mitsuba/mitsuba.h>
#include <drjit/math.h>

NAMESPACE_BEGIN(mitsuba)
NAMESPACE_BEGIN(acoustic)

// -----------------------------------------------------------------------
//! @{ \name Speed of sound calculation (adapts to the type of input parameters)
// -----------------------------------------------------------------------

/**
 * \brief Calculation methods and automatic method selector for the speed of sound
 *
 * This Function calculates the speed of sound in air 
 *
 * \param temperature
 *      The temperature in degree Celsius.
 * \param relative_humidity
 *     Relative humidity in the range of 0 to 1.
 * \param atmospheric_pressure
 *    Atmospheric pressure in Pascal
 * \param saturation_vapor_pressure
 *    Saturation vapor pressure in Pascal.
 * \param co2_ppm
 *   CO2 concentration in parts per million (ppm).
 * \param method
 *  The method to use for the calculation. Possible values are:
 *   - "auto" (default): automatically selects the method based on the types of input parameters.
 *   - "simple" The calculation follows ISO 9613-1 (Formula A.5).
 *   - "ideal_gas": calculates based on chapter 6.3 in V. E. Ostashev and D. K. Wilson, Acoustics in Moving Inhomogeneous Media, 2nd ed. London: CRC Press, 2015. doi: 10.1201/b18922.
 *   - "cramer": calculates using Cramers method described in O. Cramer, “The variation of the specific heat ratio and the speed of sound in air with temperature, pressure, humidity, and CO2 concentration,” The Journal of the Acoustical Society of America, vol. 93, no. 5, pp. 2510-2516, May 1993, doi: 10.1121/1.405827.
 * 
 * \return
 *      The speed of sound in meters per second
 */
template <typename Value>
Value speed_of_sound(Value temperature, 
                    Value relative_humidity, 
                    Value atmospheric_pressure, 
                    Value saturation_vapor_pressure, 
                    Value co2_ppm, 
                    std::string method = "auto") {
    
    // input validation - at least temperature must be provided
    if (temperature == null) {
        throw std::invalid_argument("Temperature must be provided.");
    }

    std::string selected_method = method;

    // selection logic
    if (selected_method == "auto") {
        if (relative_humidity == null) {
            selected_method = "simple";
        } else if (co2_ppm != null) {
            selected_method = "cramer";
        } else {
            selected_method = "ideal_gas";
        }
    }

    if (selected_method == "simple") {
        if (temperature < -20.0f || temperature > 50.0f) {
            throw std::invalid_argument("Temperature out of range for simple method (-20°C to 50°C).");
        }
        return 343.2f * dr::sqrt((temperature + 273.15f) / 293.15f);
    } else if (selected_method == "ideal_gas") {
        // Ideal gas calculation based on Ostashev and Wilson
        float R = 8.314f; // J/(mol*K)
        float gamma_a = 1.400f;
        float gamma_w = 1.330f;
        float mu_a = 28.97*1e-3f; //kg/mol
        float mu_w = 18.02*1e-3f; //kg/mol
        float R_a = R / mu_a;
        float p;

        if (saturation_vapor_pressure == null) {
            float e_s = 6.1094 * dr::exp((17.625 * temperature) / (temperature + 243.04));
            p = 100 * e_s;
        }
        else {
            p = saturation_vapor_pressure;
        }
        
        float e = p * relative_humidity;
        float alpha = alpha = mu_a / mu_w;
        float delta = (1 - (1/gamma_a)) / (1 - (1/gamma_w));
        float nu = (gamma_a - 1) / (gamma_w - 1);
        float C = (e/P) / (alpha * (1 - e/P));

        return dr::sqrt((gamma_a * R_a * temperature_kelvin * (1 + (alpha * (1 + delta - nu) - 1) * C)));

        return;
    } else if (selected_method == "cramer") {
        // Cramer's method
        return;
    } else {
        throw std::invalid_argument("Invalid method specified for speed of sound calculation.");
    }

}

NAMESPACE_END(acoustic)
NAMESPACE_END(mitsuba)