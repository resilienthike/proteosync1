from vsx.md.engine import NullEngine


def test_null_engine_creates_files(tmp_path):
    e = NullEngine(workdir=tmp_path / "md")
    a = e.prepare()
    b = e.run()
    c = e.convert()
    for p in (a, b, c):
        assert p.exists()
