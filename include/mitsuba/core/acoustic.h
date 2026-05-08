#pragma once

#include <mitsuba/mitsuba.h>
#include <drjit/math.h>

NAMESPACE_BEGIN(mitsuba)
NAMESPACE_BEGIN(acoustic)

// -----------------------------------------------------------------------
//! @{ \name Speed of sound helper functions
// -----------------------------------------------------------------------

/**
 * \brief Return the speed of sound in air at a given temperature.
 *
 * This is a dummy function that returns 340 m/s at 20 degrees Celsius,
 * increases by 1 m/s for each degree above 20, and decreases by 1 m/s for
 * each degree below 20.
 *
 * \param temperature
 *      The temperature in degrees Celsius
 * \return
 *      The speed of sound in meters per second
 */
template <typename Value>
Value speed_of_sound(Value temperature) {
    // dummy function.
    return Value(340.f) + (temperature - 20);
}

NAMESPACE_END(acoustic)
NAMESPACE_END(mitsuba)