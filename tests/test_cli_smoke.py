import subprocess
import sys


def test_cli_help_runs():
    res = subprocess.run(
        [sys.executable, "-m", "svphaser", "--help"], capture_output=True, text=True
    )
    assert res.returncode == 0
    assert "SvPhaser" in res.stdout or "Usage" in res.stdout


def test_cli_version_runs():
    res = subprocess.run(
        [sys.executable, "-m", "svphaser", "--version"], capture_output=True, text=True
    )
    assert res.returncode == 0
    assert res.stdout.strip()  # prints something (tag version in wheels, "0+unknown" in raw dev)
