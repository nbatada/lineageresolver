from __future__ import annotations

import subprocess
import sys


def test_cli_module_help_smoke() -> None:
    completed = subprocess.run(
        [sys.executable, "-m", "lineageresolver.cli", "--help"],
        check=False,
        capture_output=True,
        text=True,
    )
    assert completed.returncode == 0
    assert "lineageresolver" in completed.stdout
