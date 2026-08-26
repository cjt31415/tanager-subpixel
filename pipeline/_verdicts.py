"""Exp 35 stage helpers (copied from experiment 34, unchanged): the verdict ledger, and a log file next to the outputs.

Five stages each decide part of the experiment, so `metrics/verdicts.csv` is written
incrementally rather than in one place (experiment 32 did the same thing by hand in
`draw_footprints.py`). :func:`upsert_verdict` makes that an operation instead of a
copy-pasted concat, so re-running one stage replaces exactly its own rows and leaves the
others alone — a stale verdict from a previous run is the one thing this file must never
be able to keep.
"""
from __future__ import annotations

import logging
from pathlib import Path

import pandas as pd


def configure_logging(log_level: str = "info") -> None:
    """Root logger in drift's line format, without importing drift.

    These stages are republished standalone (see ``just exp35-build-study``), where there
    is no drift to import, so the experiment uses the same plain logger in both trees —
    one code path. UTC timestamps, unlike drift's local-time formatter.
    """
    logging.basicConfig(level=getattr(logging, log_level.upper(), logging.INFO),
                        format="%(asctime)s - %(levelname)s - %(message)s",
                        datefmt="%Y-%m-%d %H:%M:%S", force=True)


#: Column order of `metrics/verdicts.csv`, matching experiments 31 and 32 exactly so the
#: three tables can be concatenated when the arc is written up.
VERDICT_COLUMNS = ["hypothesis", "statistic", "value", "threshold", "verdict", "note"]


def verdict_row(hypothesis: str, statistic: str, value: str, threshold: str,
                verdict: str, note: str = "") -> dict[str, str]:
    """One verdict, with every field spelled out. ``verdict`` is the human-readable call.

    ``hypothesis`` must start with the tag the stage owns (``H1``, ``G2``, …); that prefix
    is the upsert key.
    """
    return {"hypothesis": hypothesis, "statistic": statistic, "value": value,
            "threshold": threshold, "verdict": verdict, "note": note}


def upsert_verdict(metrics_dir: Path, rows: list[dict[str, str]]) -> pd.DataFrame:
    """Replace every existing row whose tag matches one of ``rows``, then append.

    The tag is the leading whitespace-delimited token of ``hypothesis`` (``"H1 the cells
    are organised"`` → ``H1``), so a stage owns its tags and cannot disturb another's.
    """
    path = metrics_dir / "verdicts.csv"
    incoming = pd.DataFrame(rows, columns=VERDICT_COLUMNS)
    tags = {str(h).split(maxsplit=1)[0] for h in incoming["hypothesis"]}
    if path.exists():
        existing = pd.read_csv(path).fillna("")
        keep = ~existing["hypothesis"].astype(str).str.split().str[0].isin(tags)
        incoming = pd.concat([existing[keep], incoming], ignore_index=True)
    incoming.to_csv(path, index=False)
    return incoming


def configure_stage_logging(log_path: Path, log_level: str = "info") -> None:
    """Console logging as everywhere else, plus a fresh file next to the outputs.

    The file is truncated per run on purpose: it documents the run that produced the
    outputs sitting beside it, and an appended file would blur two runs into one story.
    """
    configure_logging(log_level)
    log_path.parent.mkdir(parents=True, exist_ok=True)
    handler = logging.FileHandler(log_path, mode="w")
    handler.setFormatter(logging.Formatter("%(asctime)s - %(levelname)s - %(message)s"))
    logging.getLogger().addHandler(handler)
