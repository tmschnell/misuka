import mitsuba as mi


def test01_speed_of_sound_simple(variants_all_acoustic):
    assert hasattr(mi, "acoustic")
    assert mi.acoustic.speed_of_sound(20.0) == 340
    assert mi.acoustic.speed_of_sound(21.0) == 341