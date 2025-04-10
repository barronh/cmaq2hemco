def test_open_gdemis():
    from ..mp2022 import open_gdemis
    import pandas as pd
    date = pd.to_datetime('2022-07-01')
    sct = 'merged_nobeis_norwc'
    f = open_gdemis(date, sct, cache=False)
    assert len(f.sizes) == 3
    assert 'time' in f.sizes
    assert 'ROW' in f.sizes
    assert 'COL' in f.sizes


def test_open_ptemis():
    from ..mp2022 import open_ptemis
    import pandas as pd
    date = pd.to_datetime('2022-07-01')
    sct = 'ptegu'
    f = open_ptemis(date, sct, cache=False)
    assert len(f.sizes) == 2
    assert 'time' in f.sizes
    assert 'stack' in f.sizes
