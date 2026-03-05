"""Task configuration validation and loading utilities."""

from __future__ import annotations

import json
from pathlib import Path
from typing import Any


def load_task_config(task_config: Any) -> dict[str, Any]:
    """Load task config from dict, JSON path, or simple YAML path."""
    if isinstance(task_config, dict):
        return task_config
    if isinstance(task_config, Path):
        return _load_task_config_path(task_config)
    if isinstance(task_config, str):
        return _load_task_config_path(Path(task_config))
    raise ValueError("task_config must be a dict or a path to YAML/JSON.")


def _load_task_config_path(path: Path) -> dict[str, Any]:
    if not path.exists():
        raise FileNotFoundError(f"task_config file not found: {path}")

    suffix = path.suffix.lower()
    if suffix == ".json":
        with open(path, "r", encoding="utf-8") as handle:
            loaded = json.load(handle)
    elif suffix in {".yaml", ".yml"}:
        with open(path, "r", encoding="utf-8") as handle:
            loaded = _parse_simple_yaml(handle.read())
    else:
        raise ValueError(f"Unsupported task_config file extension: {path.suffix}")

    if not isinstance(loaded, dict):
        raise ValueError("task_config must be a mapping/object.")
    return loaded


def _parse_simple_yaml(text: str) -> dict[str, Any]:
    """Parse a minimal YAML subset needed for project config files."""
    root: dict[str, Any] = {}
    stack: list[tuple[int, dict[str, Any]]] = [(-1, root)]

    for raw_line in text.splitlines():
        line = raw_line.rstrip()
        if not line or line.lstrip().startswith("#"):
            continue

        indent = len(line) - len(line.lstrip(" "))
        stripped = line.strip()
        if ":" not in stripped:
            raise ValueError(f"Invalid YAML line: {raw_line}")
        key, value = stripped.split(":", 1)
        key = key.strip()
        value = value.strip()

        while stack and indent <= stack[-1][0]:
            stack.pop()
        if not stack:
            raise ValueError("Invalid YAML indentation.")
        parent = stack[-1][1]

        if value == "":
            child: dict[str, Any] = {}
            parent[key] = child
            stack.append((indent, child))
        else:
            parent[key] = _parse_scalar(value)

    return root


def _parse_scalar(value: str) -> Any:
    lower = value.lower()
    if lower == "true":
        return True
    if lower == "false":
        return False

    if value.startswith("[") and value.endswith("]"):
        inner = value[1:-1].strip()
        if not inner:
            return []
        return [_parse_scalar(token.strip()) for token in _split_csv_like(inner)]

    if value.startswith("{") and value.endswith("}"):
        inner = value[1:-1].strip()
        if not inner:
            return {}
        result: dict[str, Any] = {}
        for token in _split_csv_like(inner):
            if ":" not in token:
                raise ValueError(f"Invalid inline mapping token: {token}")
            key, raw_val = token.split(":", 1)
            result[key.strip()] = _parse_scalar(raw_val.strip())
        return result

    if value.startswith(("\"", "'")) and value.endswith(("\"", "'")) and len(value) >= 2:
        return value[1:-1]

    try:
        if "." in value or "e" in lower:
            return float(value)
        return int(value)
    except ValueError:
        return value


def _split_csv_like(value: str) -> list[str]:
    """Split comma-separated inline values, respecting [] and {} nesting."""
    tokens: list[str] = []
    current: list[str] = []
    depth = 0
    for char in value:
        if char in "[{":
            depth += 1
        elif char in "]}":
            depth -= 1
        if char == "," and depth == 0:
            tokens.append("".join(current).strip())
            current = []
            continue
        current.append(char)
    if current:
        tokens.append("".join(current).strip())
    return tokens
