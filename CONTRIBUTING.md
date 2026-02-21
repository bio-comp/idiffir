# Contributing to iDiffIR

## Workflow

1. Open a GitHub issue for each discrete change.
2. Create a dedicated feature branch from `master`.
3. Open a PR referencing the issue.
4. Run the required local checks before requesting review.

## Community Conduct

Follow `CODE_OF_CONDUCT.md` in all project interactions.

## Local Setup

```bash
uv sync --group dev
```

## Required Verification

```bash
uv run pytest -q
uv run python -W error::SyntaxWarning -m compileall -f -q iDiffIR scripts tests
```

## Scope and Compatibility

- Keep CLI flags and output contracts stable unless the issue explicitly allows breaking changes.
- Prefer small, reviewable PRs with focused scope.
- For migration work, prioritize behavior parity and add regression tests for each bug or edge case fixed.

## Style and Refactoring

- Follow repository conventions from `pyproject.toml` and local `AGENTS.md`.
- Use `snake_case` in touched legacy code when safe.
- Add concise docstrings for touched functions when behavior is non-obvious.
