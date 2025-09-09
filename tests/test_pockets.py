from vsx.feats.pocket_def import init_pocket_file, load_pocket_residues, validate_pocket


def test_init_and_validate(tmp_path, monkeypatch):
    monkeypatch.setattr("vsx.feats.pocket_def.POCKETS_DIR", tmp_path)
    p = init_pocket_file("GLP1R", "semaglutide")
    assert p.exists()
    assert validate_pocket("GLP1R", "semaglutide").startswith("EMPTY")
    p.write_text("# header\nA\t123\nA\t127\n")
    res = load_pocket_residues("GLP1R", "semaglutide")
    assert res == [("A", 123), ("A", 127)]
    assert validate_pocket("GLP1R", "semaglutide").startswith("OK")
