#!/usr/bin/env python3
"""Build a self-contained manuscript-figure source package and labeled panels."""
import argparse
import csv
import hashlib
import io
import json
import zipfile
from pathlib import Path
from validate_publication import validate

CODE_FILES = [".gitattributes", "README.md", "Figures.Rmd", "Figures.md", "Figures.html",
              "Rfunctions/SELEXseqFunctions.R", "Rfunctions/SELEXseqFunctions_response.R",
              "Rfunctions/motifDNAScoring.R", "Rfunctions/scatter_forward_rc.R",
              "Rfunctions/validate_projection_coordinates.R", "Rfunctions/manuscript_panels.csv",
              "scripts/render_figures.R", "scripts/validate_publication.py", "scripts/package_publication.py"]


def package(root):
    root = root.resolve()
    validation = validate(root)
    files = set(CODE_FILES)
    for directory in ["CnDmelHox/PSAMs", "CnHoxMonomer/PSAMs"]:
        files.update(p.relative_to(root).as_posix() for p in (root / directory).rglob("*PSAM*.txt"))
    for directory in ["outputs/publication_panels", "outputs/manuscript_source_data"]:
        files.update(p.relative_to(root).as_posix() for p in (root / directory).iterdir() if p.is_file())
    assert all((root / name).is_file() for name in files)
    assert not any("forward_RC" in name or "validation/" in name or "comparison/" in name for name in files)
    checksums = "".join(hashlib.sha256((root / name).read_bytes()).hexdigest() + "  " + name + "\n" for name in sorted(files))
    destination = root / "dist"
    destination.mkdir(exist_ok=True)
    source_zip = destination / "SELEX_manuscript_source_data.zip"
    with zipfile.ZipFile(source_zip, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        for name in sorted(files):
            archive.write(root / name, "CnidarianHox/" + name)
        archive.writestr("CnidarianHox/SHA256SUMS.txt", checksums)
    with zipfile.ZipFile(source_zip) as archive:
        assert archive.testzip() is None
        for name in files:
            assert archive.read("CnidarianHox/" + name) == (root / name).read_bytes()
    with (root / "Rfunctions/manuscript_panels.csv").open(encoding="utf-8-sig", newline="") as handle:
        panels = list(csv.DictReader(handle))
    panel_zip = destination / "SELEX_manuscript_panels.zip"
    with zipfile.ZipFile(panel_zip, "w", zipfile.ZIP_DEFLATED, compresslevel=6) as archive:
        table = io.StringIO(newline="")
        writer = csv.writer(table)
        writer.writerow(["Panel", "Description", "PNG", "PDF"])
        for panel in panels:
            basename = panel["Panel"].replace(" ", "_")
            targets = {ext: f"{ext}/{basename}.{ext.lower()}" for ext in ["PNG", "PDF"]}
            for ext in targets:
                archive.write(root / panel[ext], targets[ext])
            writer.writerow([panel["Panel"], panel["Description"], targets["PNG"], targets["PDF"]])
        archive.writestr("panel_manifest.csv", table.getvalue())
    with zipfile.ZipFile(panel_zip) as archive:
        assert archive.testzip() is None and len(archive.namelist()) == 121
    summary = {"validation": validation, "source_files": len(files), "archives": [
        {"path": path.relative_to(root).as_posix(), "bytes": path.stat().st_size,
         "sha256": hashlib.sha256(path.read_bytes()).hexdigest()} for path in [source_zip, panel_zip]]}
    (destination / "publication_package.json").write_text(json.dumps(summary, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(summary))


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    package(parser.parse_args().root)
