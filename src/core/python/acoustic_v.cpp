#include <mitsuba/core/acoustic.h>
#include <mitsuba/python/python.h>

MI_PY_EXPORT(acoustic) {
    MI_PY_IMPORT_TYPES()

    m.def("speed_of_sound",
          acoustic::speed_of_sound<Float>,
          "temperature"_a,D(acoustic, speed_of_sound));
}