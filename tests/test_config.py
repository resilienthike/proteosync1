from vsx.utils.config import load_config
def test_load_config():
    cfg = load_config()
    assert "project" in cfg and "paths" in cfg
