def test_plumerise_briggs():
    from ..utils import plumerise_briggs
    dz = plumerise_briggs(
        1, 1, 288.15, pres_a=101325., temp_a=288.15
    )
    assert dz == 0
