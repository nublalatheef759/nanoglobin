#!/usr/bin/env python3
"""One-shot branch bootstrap; removed from the generated commit."""
from __future__ import annotations

import base64
import io
from pathlib import Path
import shutil
import subprocess
import tarfile


def run(*args: str) -> None:
    subprocess.run(args, check=True)


def main() -> None:
    chunks = sorted(Path(".bootstrap").glob("payload-*.b64"))
    if not chunks:
        raise RuntimeError("bootstrap payload chunks are missing")
    encoded = "".join(path.read_text(encoding="ascii").strip() for path in chunks)
    data = base64.b64decode(encoded, validate=True)

    with tarfile.open(fileobj=io.BytesIO(data), mode="r:gz") as archive:
        for member in archive.getmembers():
            if not member.isfile():
                continue
            path = Path(member.name)
            if path.is_absolute() or ".." in path.parts:
                raise RuntimeError(f"unsafe archive path: {member.name}")
            source = archive.extractfile(member)
            if source is None:
                raise RuntimeError(f"missing archive member: {member.name}")
            path.parent.mkdir(parents=True, exist_ok=True)
            path.write_bytes(source.read())
            path.chmod(member.mode)

    shutil.rmtree(".bootstrap")
    Path("bootstrap_foundation.py").unlink(missing_ok=True)
    Path(".github/workflows/bootstrap-foundation.yml").unlink(missing_ok=True)

    run("git", "config", "user.name", "github-actions[bot]")
    run("git", "config", "user.email", "41898282+github-actions[bot]@users.noreply.github.com")
    run("git", "add", "-A")
    if subprocess.run(["git", "diff", "--cached", "--quiet"]).returncode == 0:
        return
    run("git", "commit", "-m", "Add cohort and assay compiler foundations")
    run("git", "push", "origin", "HEAD")


if __name__ == "__main__":
    main()
