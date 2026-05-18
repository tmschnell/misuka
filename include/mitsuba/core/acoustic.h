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
        // ISO 9613-1 (Formula A.5)
        return 331.3f * dr::sqrt(1 + temperature / 273.15f);
    } else if (selected_method == "ideal_gas") {
        // Ideal gas calculation based on Ostashev and Wilson
        Value T = temperature + 273.15f; // Convert to Kelvin
        Value R = 287.05f; // Specific gas constant for dry air in J/(kg*K)
        Value gamma = 1.4f; // Adiabatic index for dry air
        return dr::sqrt(gamma * R * T);
    } else if (selected_method == "cramer") {
        // Cramer's method
        Value T = temperature + 273.15f; // Convert to Kelvin
        Value P = atmospheric_pressure; // Pressure in Pascal
        Value H = relative_humidity; // Relative humidity (0 to 1)
        Value S = saturation_vapor_pressure; // Saturation vapor pressure in Pascal
        Value C = co2_ppm; // CO2 concentration in ppm

        Value a = 331.3f; // Speed of sound at 0°C in m/s
        Value b = 0.6f; // Temperature coefficient in m/s per °C
        Value c = 0.0124f; // Humidity coefficient in m/s per %RH
        Value d = -0.0001f; // Pressure coefficient in m/s per Pa
        Value e = -0.00001f; // CO2 coefficient in m/s per ppm

        return a + b * temperature + c * H * 100 + d * P + e * C;
    } else {
        throw std::invalid_argument("Invalid method specified for speed of sound calculation.");
    }

}

NAMESPACE_END(acoustic)
NAMESPACE_END(mitsuba)