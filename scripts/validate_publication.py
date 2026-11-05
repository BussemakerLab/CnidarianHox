#!/usr/bin/env python3
"""Validate manuscript figure data and exports with the Python standard library."""
import argparse
import base64
import csv
import hashlib
import json
import math
import re
import struct
from collections import defaultdict
from html.parser import HTMLParser
from pathlib import Path


def read_csv(path):
    with path.open(encoding="utf-8-sig", newline="") as handle:
        return list(csv.DictReader(handle))


def sha256(path):
    return hashlib.sha256(path.read_bytes()).hexdigest()


def close(actual, expected, tolerance=1e-9):
    assert math.isfinite(actual) and math.isfinite(expected)
    assert math.isclose(actual, expected, rel_tol=tolerance, abs_tol=tolerance), (actual, expected)


def reverse_complement(sequence):
    return sequence.translate(str.maketrans("ACGTN", "TGCAN"))[::-1]


def score_windows(sequence, columns):
    sequence = "NNNN" + sequence + "NNNN"
    return math.fsum(math.exp(math.fsum(
        min(column.values()) if base == "N" else column[base]
        for base, column in zip(sequence[start:start + len(columns)], columns)
    )) for start in range(len(sequence) - len(columns) + 1))


class ReportImages(HTMLParser):
    def __init__(self):
        super().__init__()
        self.images = []

    def handle_starttag(self, tag, attrs):
        if tag == "img":
            self.images.append(dict(attrs).get("src", ""))


def validate(root):
    root = root.resolve()
    source = root / "outputs/manuscript_source_data"
    panels = read_csv(source / "panel_manifest.csv")
    assert panels == read_csv(root / "Rfunctions/manuscript_panels.csv")
    assert len(panels) == 60 and sum(p["Kind"] == "scatter" for p in panels) == 28
    assert len({p["Panel"] for p in panels}) == 60
    for extension in ["PNG", "PDF"]:
        expected = {p[extension] for p in panels}
        actual = {p.relative_to(root).as_posix() for p in
                  (root / "outputs/publication_panels").glob("*." + extension.lower())}
        assert actual == expected, ("Unexpected exports", sorted(actual ^ expected))
    for panel in panels:
        png = (root / panel["PNG"]).read_bytes()
        assert png[:8] == b"\x89PNG\r\n\x1a\n"
        width, height = struct.unpack(">II", png[16:24])
        assert width > 500 and height > 500
        pdf = (root / panel["PDF"]).read_bytes()
        assert pdf.startswith(b"%PDF-") and b"%%EOF" in pdf[-100:]
    document = ReportImages()
    document.feed((root / "Figures.html").read_text(encoding="utf-8"))
    assert len(document.images) == 60
    embedded_hashes = []
    for src in document.images:
        assert src.startswith("data:image/png;base64,"), src[:80]
        embedded_hashes.append(hashlib.sha256(base64.b64decode(src.split(",", 1)[1])).hexdigest())
    assert embedded_hashes == [sha256(root / p["PNG"]) for p in panels]
    markdown = (root / "Figures.md").read_text(encoding="utf-8")
    assert re.findall(r'<img\s+src="([^"]+)"', markdown) == [p["PNG"] for p in panels]
    rmd = (root / "Figures.Rmd").read_text(encoding="utf-8")
    for label in ["Figure S3b", "Figure S4b", "Figure S5b", "Figure S11c"]:
        assert f"**{label}**" in rmd
    assert "Figure S11d" not in rmd and "supp11-monomer-tree" not in rmd
    assert "figures_forward_rc/" not in rmd and "Rfunctions/scatter_forward_rc.R" in rmd
    for entry in read_csv(source / "input_files.csv"):
        path = (root / entry["Path"]).resolve()
        assert root in path.parents and sha256(path) == entry["SHA256"]
    sequences = read_csv(source / "scatter_sequences.csv")
    assert len(sequences) == 2000 and len({r["Sequence"] for r in sequences}) == 2000
    assert all(len(r["Sequence"]) == 12 and set(r["Sequence"]) <= set("ACGT") for r in sequences)
    assert all(reverse_complement(r["Sequence"]) == r["Reverse_complement"] for r in sequences)
    assert [r["Sequence"] for r in sequences[:5]] == [
        "AGTTGATAAATG", "ACATGATAAATT", "TGATTTATTAAG", "AGTGATCGATCA", "CTTCTGATCGAT"]
    coefficients = defaultdict(dict)
    for row in read_csv(source / "scatter_model_coefficients.csv"):
        key = (row["Model"], int(row["Model_width"]))
        coefficients[key].setdefault(int(row["Position"]), {})[row["Base"]] = float(row["Scoring_coefficient"])
    for key, positions in coefficients.items():
        assert sorted(positions) == list(range(1, key[1] + 1))
        for values in positions.values():
            assert set(values) == set("ACGT")
            close(max(values.values()), 0)
    affinities = {}
    by_model = defaultdict(list)
    for row in read_csv(source / "scatter_model_affinities.csv"):
        key = (row["Model"], int(row["Model_width"]), row["Kmer"])
        assert key not in affinities
        affinities[key] = row
        by_model[key[:2]].append(row)
        close(float(row["Affinity"]), float(row["ForwardAffinity"]) + float(row["ReverseComplementAffinity"]))
    assert len(by_model) == 16 and len(affinities) == 32000
    max_scoring_error = 0.0
    for model, rows in by_model.items():
        assert len(rows) == 2000
        columns = [coefficients[model][i] for i in range(1, model[1] + 1)]
        rc_columns = [{base: col[reverse_complement(base)] for base in "ACGT"} for col in reversed(columns)]
        for row in rows[:50] + rows[50::53]:
            sequence = row["Kmer"]
            forward = score_windows(sequence, columns)
            reverse = score_windows(sequence, rc_columns)
            close(reverse, score_windows(reverse_complement(sequence), columns))
            for actual, expected in [(forward, float(row["ForwardAffinity"])),
                                     (reverse, float(row["ReverseComplementAffinity"])),
                                     (forward + reverse, float(row["Affinity"]))]:
                close(actual, expected)
                max_scoring_error = max(max_scoring_error, abs(actual - expected))
    scatter_panels = {row["Panel"]: row for row in read_csv(source / "scatter_panels.csv")}
    assert set(scatter_panels) == {p["Panel"] for p in panels if p["Kind"] == "scatter"}
    coordinates = read_csv(source / "scatter_coordinates.csv")
    groups = defaultdict(list)
    panel_sequences = defaultdict(set)
    for row in coordinates:
        spec = scatter_panels[row["Panel"]]
        sequence = row["Kmer"]
        assert row["probe"] in sequence
        x_key = (spec["X_model"], int(spec["X_model_width"]), sequence)
        y_key = (spec["Y_model"], int(spec["Y_model_width"]), sequence)
        close(float(row["set1Affinities"]), float(affinities[x_key]["Affinity"]))
        close(float(row["set2Affinities"]), float(affinities[y_key]["Affinity"]))
        groups[(row["Panel"], row["probe"])].append((float(row["set1Affinities"]), float(row["set2Affinities"])))
        panel_sequences[row["Panel"]].add(sequence)
    assert len(groups) == 196 and all(len(s) == 2000 for s in panel_sequences.values())
    fits = read_csv(source / "scatter_regression_fits.csv")
    assert len(fits) == 196
    for row in fits:
        points = groups[(row["Panel"], row["probe"])]
        assert len(points) == int(row["n"])
        x_mean = math.fsum(x for x, _ in points) / len(points)
        y_mean = math.fsum(y for _, y in points) / len(points)
        xx = math.fsum((x - x_mean) ** 2 for x, _ in points)
        yy = math.fsum((y - y_mean) ** 2 for _, y in points)
        xy = math.fsum((x - x_mean) * (y - y_mean) for x, y in points)
        close(xy / xx, float(row["slope"]))
        close(y_mean - xy / xx * x_mean, float(row["intercept"]))
        close(xy * xy / (xx * yy), float(row["r_squared"]))
    basis = read_csv(source / "projection_error_basis.csv")
    assert len(basis) == 36
    for a in ["u1", "u2"]:
        for b in ["u1", "u2"]:
            close(math.fsum(float(r[a]) * float(r[b]) for r in basis), float(a == b))
    differences = read_csv(source / "projection_error_differences.csv")
    pairs = {r["Pair"]: r for r in read_csv(source / "replicate_pair_distances.csv")}
    assert len(differences) == len(pairs) == 4
    distances = []
    for row in differences:
        projected = [math.fsum(float(row[b["Coordinate"]]) * float(b[axis]) for b in basis) for axis in ["u1", "u2"]]
        distance = math.hypot(*projected)
        close(distance, float(pairs[row["Pair"]]["Projected_Euclidean_distance"]))
        distances.append(distance)
    for row in read_csv(source / "projection_error_bars.csv"):
        assert int(row["Input_coordinates"]) == 36 and int(row["Complete_pairs"]) == 4
        close(float(row["Half_width"]), math.fsum(distances) / 4)
        close(float(row["Total_length"]), 2 * math.fsum(distances) / 4)
    return {"passed": True, "publication_panels": 60, "scatter_panels": 28,
            "unique_sequences": 2000, "scored_models": 16,
            "model_sequence_scores": len(affinities), "scatter_source_rows": len(coordinates),
            "regression_fits": len(fits), "maximum_independent_scoring_error": max_scoring_error,
            "html_images_match_exports": True, "markdown_images_match_exports": True,
            "projection_error_input_coordinates": 36, "projection_error_pairs": 4}


if __name__ == "__main__":
    parser = argparse.ArgumentParser(description=__doc__)
    parser.add_argument("--root", type=Path, default=Path(__file__).resolve().parents[1])
    args = parser.parse_args()
    result = validate(args.root)
    destination = args.root / "outputs/manuscript_panel_validation/publication_validation.json"
    destination.parent.mkdir(parents=True, exist_ok=True)
    destination.write_text(json.dumps(result, indent=2) + "\n", encoding="utf-8")
    print(json.dumps(result))
