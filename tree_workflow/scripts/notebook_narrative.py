"""Helpers that turn numbers into the short findings printed at the top and bottom of the notebooks.

The static explanation in a notebook is hand-written; the *findings* are generated here so they follow
each run's numbers. Every judgement comes from an explicit threshold (a :class:`Rule`), listed in the
notebook's ``THRESHOLDS`` cell, so changing what counts as "a problem" means editing one number.

    OK    nothing to do
    LOOK  worth a look, probably fine
    BAD   a problem to fix before using the result
"""
from dataclasses import dataclass
from typing import Iterable, List, Optional

OK, LOOK, BAD = "OK", "LOOK", "BAD"
MARK = {OK: "\u2705", LOOK: "\u26a0\ufe0f", BAD: "\u274c"}   # shown as emoji in the rendered notebook
_ORDER = {OK: 0, LOOK: 1, BAD: 2}


@dataclass(frozen=True)
class Rule:
    """Thresholds for a measure. ``higher_is_worse`` False means small values are the problem."""
    look: float
    bad: float
    higher_is_worse: bool = True

    def status(self, value: Optional[float]) -> str:
        if value is None or value != value:   # None / NaN
            return LOOK
        if self.higher_is_worse:
            return BAD if value >= self.bad else LOOK if value >= self.look else OK
        return BAD if value <= self.bad else LOOK if value <= self.look else OK


def worst(statuses: Iterable[str]) -> str:
    return max(statuses, key=_ORDER.__getitem__, default=OK)


def pct(x: Optional[float], digits: int = 0) -> str:
    return "n/a" if x is None or x != x else f"{100 * x:.{digits}f}%"


def count(n: int, singular: str, plural: Optional[str] = None) -> str:
    return f"{n:,} {singular if n == 1 else (plural or singular + 's')}"


def across(rule, values, fmt=lambda v: f"{v:.2f}"):
    """(status, text) for one measure that has a value per variant.

    The status is the worst over the variants; the text is the range, naming the variants at both ends
    (just the value when they all print the same).
    """
    status = worst(rule.status(v) for v in values.values())
    lo, hi = min(values, key=values.get), max(values, key=values.get)
    if len(values) == 1 or fmt(values[lo]) == fmt(values[hi]):
        return status, fmt(values[lo])
    return status, f"{fmt(values[lo])} (`{lo}`) to {fmt(values[hi])} (`{hi}`)"


class Findings:
    """Collects findings (status, sentence, optional action) and renders them as markdown."""

    def __init__(self):
        self.items: List[dict] = []

    def add(self, status: str, text: str, evidence: str = "", action: str = "") -> str:
        self.items.append({"status": status, "text": text, "evidence": evidence, "action": action})
        return status

    def say(self, status: str, text: str, evidence: str = "", action: str = "") -> str:
        """Record a finding and show it as a one-line takeaway right where it is called (under a figure)."""
        self.add(status, text, evidence, action)
        show(f"{MARK[status]} **Takeaway.** {text}")
        return status

    def overall(self) -> str:
        return worst(i["status"] for i in self.items)

    def bullets(self) -> str:
        return "\n".join(f"- {MARK[i['status']]} {i['text']}" for i in self.items)

    def verdict(self) -> str:
        """One sentence summing up all findings."""
        counts = {k: sum(i["status"] == k for i in self.items) for k in (BAD, LOOK, OK)}
        head = {OK: "No problems found.", LOOK: "Nothing is broken, a few things are worth a look.",
                BAD: "There are problems to fix before relying on this result."}[self.overall()]
        return (f"**{MARK[self.overall()]} {head}** {len(self.items)} checks: "
                f"{counts[OK]} fine, {counts[LOOK]} worth a look, {counts[BAD]} problems.")

    def table(self):
        """Findings as a list of dicts (for writing to a TSV)."""
        return [dict(i) for i in self.items]

    def actions(self) -> str:
        """Markdown table of the findings that need attention (not OK), worst first."""
        rows = sorted((i for i in self.items if i["status"] != OK), key=lambda i: -_ORDER[i["status"]])
        if not rows:
            return "Nothing needs attention."
        lines = ["| | finding | evidence | suggested action |", "|---|---|---|---|"]
        for i in rows:
            lines.append(f"| {MARK[i['status']]} | {_cell(i['text'])} | {_cell(i['evidence'])} | {_cell(i['action'])} |")
        return "\n".join(lines)


def _cell(text: str) -> str:
    return (text or "").replace("|", "/").replace("\n", " ")


def status_table(rows: List[dict], columns: List[str]) -> str:
    """Markdown table; a ``status`` key in a row is shown as the first column."""
    head = ["", *columns] if rows and "status" in rows[0] else columns
    lines = ["| " + " | ".join(head) + " |", "|" + "---|" * len(head)]
    for row in rows:
        cells = ([MARK[row["status"]]] if "status" in row else []) + [_cell(str(row.get(c, ""))) for c in columns]
        lines.append("| " + " | ".join(cells) + " |")
    return "\n".join(lines)


def show(markdown: str) -> None:
    """Display generated markdown in a notebook (falls back to print outside Jupyter)."""
    try:
        from IPython.display import Markdown, display
        display(Markdown(markdown))
    except ImportError:  # pragma: no cover
        print(markdown)
