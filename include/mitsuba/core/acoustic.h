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
    }

    else if (relative_humidity < 0.0f || relative_humidity > 1.0f) {
        throw std::invalid_argument("Relative humidity must be in the range of 0 to 1.");
    }
    else if (atmospheric_pressure < 0.0f) {
        throw std::invalid_argument("Atmospheric pressure must be non-negative.");
    }

    else if (selected_method == "ideal_gas") {

        // Ideal gas calculation based on Ostashev and Wilson
        float R = 8.314f; // J/(mol*K)
        float gamma_a = 1.400f;
        float gamma_w = 1.330f;
        float mu_a = 28.97*1e-3f; //kg/mol
        float mu_w = 18.02*1e-3f; //kg/mol
        float R_a = R / mu_a;
        float p;

        if (saturation_vapor_pressure == null) {
            p = saturation_vapor_pressure_magnus(temperature); // compute saturation vapor pressure using Magnus formula if not provided
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
    } 
    else if (selected_method == "cramer") {
        // Cramer's specific bounds for temperature and pressure
        if (temperature < 0.0f || temperature > 30.0f) {
            throw std::invalid_argument("Temperature out of range for Cramer's method (0°C to 30°C).");
        }
        else if (atmospheric_pressure < 75000.0f || atmospheric_pressure > 102000.0f) {
            throw std::invalid_argument("Atmospheric pressure out of range for Cramer's method (75,000 Pa to 102,000 Pa).");
        }
        else if (co2_ppm < 0.0f || co2_ppm > 10000.0f) {
            throw std::invalid_argument("CO2 concentration out of range for Cramer's method (0 ppm to 10,000 ppm).");
        }
        if (atmospheric_pressure == null) {
            atmospheric_pressure = 101325.0;
        }

        float x_c = co2_ppm * 1e-6f; // Convert ppm to mole fraction
        float T = temperature + 273.15f; // Convert to Kelvin
        float p = atmospheric_pressure; // in Pa

        float f = 1.00062+3.14e-8f * p + 5.6e-7f*temperature*temperature;
        float p_sv = dr::exp(1.2811805e-5f*T*T-1.9509874e-2f*T+34.04926034-6.3536311e3f/T);

        float x_w = relative_humidity * p_sv / p; // Mole fraction of water vapor

        if (x_w < 0.0f || x_w > 0.06f) {
            throw std::invalid_argument("Calculated mole fraction of water vapor is out of range (0 to 0.06). Check input parameters.");
        }

        float a0 = 331.5024;
        float a1 = 0.603055;
        float a2 = -0.000528;
        float a3 = 51.471935;
        float a4 = 0.1495874;
        float a5 = -0.000782;
        float a6 = -1.82e-7;
        float a7 = 3.73e-8;
        float a8 = -2.93e-10;
        float a9 = -85.20931;
        float a10 = -0.228525;
        float a11 = 5.91e-5;
        float a12 = -2.835149;
        float a13 = -2.15e-13;
        float a14 = 29.179762;
        float a15 = 0.000486;

        return  a0 + a1*temperature + a2*temperature*temperature 
                + (a3 + a4*temperature + a5*temperature*temperature) * x_w 
                + (a6 + a7*temperature + a8*temperature*temperature) * p 
                + (a9 + a10*temperature + a11*temperature*temperature) * x_c
                + (a12*x_w*x_w + a13*p*p + a14*x_c*x_c + a15*x_c*p*x_w);

    } else {
        throw std::invalid_argument("Invalid method specified for speed of sound calculation.");
    }
}

NAMESPACE_END(acoustic)
NAMESPACE_END(mitsuba)

float saturation_vapor_pressure_magnus (float temperature) {
    float e_s = 6.1094 * dr::exp((17.625 * temperature) / (temperature + 243.04));
    return 100 * e_s; // Convert from hPa to Pa
}