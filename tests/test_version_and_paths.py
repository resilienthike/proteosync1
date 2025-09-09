import vsx
from vsx.utils.paths import REPO_ROOT, ARTIFACTS_DIR, RAW_DIR, STD_DIR


def test_version_string():
    assert isinstance(vsx.__version__, str)
    assert vsx.__version__.count(".") >= 1


def test_paths_exist():
    assert (REPO_ROOT / "pyproject.toml").exists()
    for d in (ARTIFACTS_DIR, RAW_DIR, STD_DIR):
        assert d.exists() and d.is_dir()
