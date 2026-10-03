"""NEO-CRISPR Core Engine -- local FastAPI server.

Serves a fully offline web UI (no CDN; see ``static/``) and a small JSON API
that drives the bioinformatics engine in :mod:`engine`. The heavy genome scan
runs in a worker thread so the server stays responsive.
"""

import os
import threading
import webbrowser

import uvicorn
from fastapi import FastAPI
from fastapi.concurrency import run_in_threadpool
from fastapi.responses import FileResponse
from fastapi.staticfiles import StaticFiles

import engine

app = FastAPI(title="NEO-CRISPR Core Engine")

BASE_DIR = os.path.dirname(os.path.abspath(__file__))
DATA_DIR = os.path.join(BASE_DIR, "data")
STATIC_DIR = os.path.join(BASE_DIR, "static")
MAX_CANDIDATES = 15

app.mount("/static", StaticFiles(directory=STATIC_DIR), name="static")


# -----------------------------------------------------------------------------
# Species discovery
# -----------------------------------------------------------------------------
def scan_genomes():
    """Find plug-and-play species directories under ``data/``.

    A directory qualifies if it contains at least one FASTA (``.fa``/``.fasta``)
    file. A ``.gtf`` annotation is optional; if present, off-targets are
    labelled with the gene they fall in.
    """
    species = []
    if not os.path.exists(DATA_DIR):
        os.makedirs(DATA_DIR)
    for item in sorted(os.listdir(DATA_DIR)):
        path = os.path.join(DATA_DIR, item)
        if not os.path.isdir(path):
            continue
        fa = [f for f in os.listdir(path) if f.endswith((".fa", ".fasta"))]
        gtf = [f for f in os.listdir(path) if f.endswith(".gtf")]
        if not fa:
            continue
        species.append({
            "id": item,
            "name": item.replace("_", " "),
            "fasta": os.path.join(path, fa[0]),
            "gtf": os.path.join(path, gtf[0]) if gtf else None,
            "annotated": bool(gtf),
        })
    return species


# -----------------------------------------------------------------------------
# Design (runs off the event loop)
# -----------------------------------------------------------------------------
def _design(seq: str, species: dict) -> dict:
    cleaned = engine.clean_sequence(seq)
    sites = engine.find_guide_sites(cleaned)
    if not sites:
        return {
            "candidates": [], "input_sequence": cleaned,
            "device": engine.device_name(), "annotated": species["annotated"],
        }

    spacers = [s.spacer for s in sites]
    off_lists = engine.scan_off_targets(spacers, species["fasta"])

    gtf_index = engine.load_gtf(species["gtf"]) if species["gtf"] else None

    candidates = []
    for site, offs in zip(sites, off_lists):
        for ot in offs:
            if gtf_index is not None:
                ot.gene = engine.annotate(gtf_index, ot.chrom, ot.pos)
        candidates.append({
            "target": site.target,
            "spacer": site.spacer,
            "pam": site.pam,
            "strand": site.strand,
            "start": site.start,
            "end": site.end,
            "gc": f"{engine.gc_percent(site.spacer)}%",
            "score_raw": engine.heuristic_score(site.spacer),
            "score": f"{engine.heuristic_score(site.spacer)}",
            "off_targets": [
                {
                    "loc": ot.loc, "seq": ot.seq, "mm": ot.mm,
                    "seed_mm": ot.seed_mm, "pam": ot.pam,
                    "gene": ot.gene,
                }
                for ot in offs
            ],
        })

    # Rank: fewest off-targets first, then highest heuristic score.
    candidates.sort(key=lambda c: (len(c["off_targets"]), -c["score_raw"]))
    if candidates:
        candidates[0]["is_best"] = True

    return {
        "candidates": candidates[:MAX_CANDIDATES],
        "input_sequence": cleaned,
        "device": engine.device_name(),
        "annotated": species["annotated"],
    }


# -----------------------------------------------------------------------------
# API endpoints
# -----------------------------------------------------------------------------
@app.get("/api/species")
def get_species():
    return [
        {"id": s["id"], "name": s["name"], "annotated": s["annotated"]}
        for s in scan_genomes()
    ]


@app.post("/api/design")
async def design_guides(data: dict):
    seq = data.get("sequence", "")
    species_id = data.get("species_id", "")
    species = next((s for s in scan_genomes() if s["id"] == species_id), None)
    if not species:
        return {"error": "Target species configuration missing."}
    if not engine.clean_sequence(seq):
        return {"error": "No sequence provided."}
    # Offload the CPU/GPU-heavy scan so the event loop is not blocked.
    return await run_in_threadpool(_design, seq, species)


@app.get("/")
def index():
    return FileResponse(os.path.join(STATIC_DIR, "index.html"))


def open_browser():
    webbrowser.open("http://127.0.0.1:8000")


if __name__ == "__main__":
    threading.Timer(1.5, open_browser).start()
    uvicorn.run(app, host="127.0.0.1", port=8000)
