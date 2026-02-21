# iDiffIR

iDiffIR identifies differential intron retention from RNA-seq data.

## Python support

Python 2 is **formally dropped**. This project supports Python:

- 3.10
- 3.11
- 3.12

## Installation and development with `uv`

### One-time setup

```bash
uv sync --group dev
```

This creates/uses a local virtual environment and installs project + development dependencies from `uv.lock`.

### Run tests

```bash
uv run pytest
```

### Run project scripts/tools

Use `uv run` for project commands, for example:

```bash
uv run python scripts/idiffir.py --help
```

Use `uv tool` for standalone CLI tools that are not project dependencies, for example:

```bash
uv tool run ruff --version
```

## Contributing

Contribution guidelines are in `CONTRIBUTING.md`.
Community expectations are in `CODE_OF_CONDUCT.md`.

## Packaging

The project uses PEP 621 metadata in `pyproject.toml`.

Runtime dependencies are declared in `project.dependencies`, with optional groups in:

- `[dependency-groups].dev`

## Documentation

- Project docs source lives under `doc/`.
- Read the Docs hosting is tracked in backlog issue `#47`.

## License

This project is distributed under the GNU GPL; see `LICENSE`.

## Release Notes

### 0.3.2

- Hardened legacy CLI scripts to avoid import-time argument parsing and side effects.
- Added smoke coverage for all script `--help` entrypoints and no-side-effect checks.
- Added strict compile gate (`-W error::SyntaxWarning`) in tests and CI.
- Cleaned remaining invalid escape sequences that blocked strict compile checks.
- Aligned runtime package version with `pyproject.toml`.
