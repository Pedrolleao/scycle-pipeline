"""
cycle_ui.py — terminal presentation for the cycle pipelines (mcycle / ncycle / scycle),
in the Vee Lab website palette (vee-lab.eu, same scheme as archaeahq-update).

Palette (CSS tokens of vee-lab.eu/assets/site.css):
  surface #0e0e0e · primary #c084fc · secondary #54c7fc · tertiary #f3ffcd · lime #c5fe00
  on-surface-variant #adaaaa · outline #767575 · outline-variant #484847

The module uses `rich` when it is importable and degrades to plain text otherwise, so
the launcher can call it before the conda env (which ships rich) is active. `--plain`,
the NO_COLOR convention and a non-terminal stdout all give the plain rendering too.

Presentation only: nothing here touches the pipeline's inputs, databases or calls.
"""

from __future__ import annotations

import os
import sys
import time
from contextlib import contextmanager
from typing import Dict, Iterable, Optional, Sequence

try:
    from rich import box
    from rich.console import Console
    from rich.panel import Panel
    from rich.progress import (BarColumn, MofNCompleteColumn, Progress, SpinnerColumn,
                               TextColumn, TimeElapsedColumn, TimeRemainingColumn)
    from rich.table import Table
    from rich.text import Text
    from rich.theme import Theme
    HAVE_RICH = True
except ImportError:  # outside the conda env, before the bootstrap
    HAVE_RICH = False

PALETTE = {
    "primary": "#c084fc",     # purple — titles, finished bars
    "secondary": "#54c7fc",   # cyan — stage rules, table headers, spinners
    "tertiary": "#f3ffcd",    # cream — accents, hints
    "lime": "#c5fe00",        # OK marks
    "muted": "#adaaaa",       # on-surface-variant
    "faint": "#767575",       # outline
    "outline": "#484847",     # outline-variant — borders, rules
    "danger": "#ff6b6b",
}

CREDIT_LINE = "Created at Vee Lab - Radboud University"

_console = None
_plain = True
_tool = "cycle"


# ── Setup ─────────────────────────────────────────────────────────────────────

def init(tool: str, plain: bool = False):
    """Choose the rendering. Plain when asked, when rich is missing, when NO_COLOR is
    set, or when stdout is not a terminal (logs, CI)."""
    global _console, _plain, _tool
    _tool = tool
    want_plain = plain or not HAVE_RICH or bool(os.environ.get("NO_COLOR")) or not sys.stdout.isatty()
    _plain = want_plain
    if HAVE_RICH:
        theme = Theme({
            "primary": PALETTE["primary"],
            "secondary": PALETTE["secondary"],
            "tertiary": PALETTE["tertiary"],
            "lime": PALETTE["lime"],
            "muted": PALETTE["muted"],
            "faint": PALETTE["faint"],
            "ok": f"bold {PALETTE['lime']}",
            "warn": f"bold {PALETTE['tertiary']}",
            "err": f"bold {PALETTE['danger']}",
            "stage": f"bold {PALETTE['secondary']}",
            "title": f"bold {PALETTE['primary']}",
            "rule.line": PALETTE["outline"],
            "progress.percentage": PALETTE["secondary"],
            "progress.elapsed": PALETTE["faint"],
            "progress.remaining": PALETTE["faint"],
            "bar.complete": PALETTE["primary"],
            "bar.finished": PALETTE["lime"],
            "bar.pulse": PALETTE["secondary"],
        })
        _console = Console(theme=theme, highlight=False, no_color=want_plain,
                           emoji=not want_plain)
    return _console


def console():
    if _console is None:
        init(_tool, plain=True)
    return _console


def is_plain() -> bool:
    return _plain


def is_live() -> bool:
    """True when live widgets (progress bars, spinners) can be drawn."""
    return HAVE_RICH and not _plain and _console is not None and _console.is_terminal


def _out(s: str) -> None:
    print(s, flush=True)


# ── Messages ──────────────────────────────────────────────────────────────────

def banner(subcommand: str, version: str = "", subtitle: str = "") -> None:
    if _plain:
        head = f"{_tool}-pipeline" + (f"  {version}" if version else "") + (f"  ·  {subcommand}" if subcommand else "")
        _out(head)
        if subtitle:
            _out(subtitle)
        _out(CREDIT_LINE)
        _out("")
        return
    c = console()
    title = Text()
    title.append(_tool[0], style=f"bold {PALETTE['primary']}")
    title.append(_tool[1:], style=f"bold {PALETTE['secondary']}")
    title.append("  pipeline", style=f"bold {PALETTE['tertiary']}")
    tail = "   " + "  ·  ".join(x for x in (version, subcommand) if x)
    if tail.strip():
        title.append(tail, style="muted")
    body = Text(subtitle, style="faint") if subtitle else Text("")
    credit = Text(CREDIT_LINE, style="faint")
    parts = [title, "\n", body, "\n", credit] if subtitle else [title, "\n", credit]
    c.print(Panel(Text.assemble(*parts), border_style=PALETTE["outline"], box=box.ROUNDED, padding=(0, 2)))


def info(msg: str) -> None:
    if _plain:
        _out(f"  . {msg}")
    else:
        console().print(f"[faint]·[/] {msg}")


def ok(msg: str) -> None:
    if _plain:
        _out(f"  [OK] {msg}")
    else:
        console().print(f"[ok]✓[/] {msg}")


def warn(msg: str) -> None:
    if _plain:
        _out(f"  [WARN] {msg}")
    else:
        console().print(f"[warn]![/] {msg}")


def err(msg: str) -> None:
    if _plain:
        _out(f"  [ERROR] {msg}")
    else:
        console().print(f"[err]✗[/] {msg}")


def kv(key: str, value, key_width: int = 18) -> None:
    if _plain:
        _out(f"  {key:<{key_width}} {value}")
    else:
        console().print(f"  [muted]{key:<{key_width}}[/] {value}")


def ask(question: str) -> str:
    """Free-text prompt on the terminal."""
    if _plain:
        return input(f"? {question} ").strip()
    return console().input(f"[secondary]?[/] {question} ").strip()


def choose(question: str, choices: Dict[str, str], default: str) -> str:
    """Prompt for one of `choices` (key -> label); the first letter of the answer picks.
    Returns the chosen key; an empty answer or end of input gives `default`."""
    hint = " / ".join(f"[{k[0]}]{k[1:]}" for k in choices)
    try:
        ans = ask(f"{question} {hint}").lower()
    except (EOFError, KeyboardInterrupt):
        _out("")
        return default
    for k in choices:
        if ans.startswith(k[0]):
            return k
    return default


@contextmanager
def stage(index: int, total: int, title: str):
    """Stage header, body, elapsed time."""
    t0 = time.time()
    if _plain:
        _out("")
        _out(f"[{index}/{total}] {title}")
    else:
        c = console()
        c.print()
        c.rule(Text(f"[{index}/{total}]  {title}", style="stage"), align="left")
    try:
        yield
    finally:
        dt = time.time() - t0
        if _plain:
            _out(f"  finished in {fmt_duration(dt)}")
        else:
            console().print(f"[faint]  finished in {fmt_duration(dt)}[/]")


def fmt_duration(seconds: float) -> str:
    if seconds < 60:
        return f"{seconds:.1f} s"
    m, s = divmod(int(seconds), 60)
    if m < 60:
        return f"{m} min {s:02d} s"
    h, m = divmod(m, 60)
    return f"{h} h {m:02d} min"


def fmt_size(n_bytes: int) -> str:
    for unit in ("B", "kB", "MB", "GB"):
        if n_bytes < 1000 or unit == "GB":
            return f"{n_bytes:.0f} {unit}" if unit == "B" else f"{n_bytes:.1f} {unit}"
        n_bytes /= 1000
    return f"{n_bytes:.1f} GB"


# ── Tables ────────────────────────────────────────────────────────────────────

def table(title: str, columns: Sequence[str], rows: Iterable[Sequence],
          justify: Optional[Dict[str, str]] = None, row_styles: Optional[Sequence[str]] = None) -> None:
    rows = [[str(x) for x in r] for r in rows]
    justify = justify or {}
    if _plain:
        widths = [max([len(c)] + [len(r[i]) for r in rows]) for i, c in enumerate(columns)]
        def fmt(cells):
            return "  " + "  ".join(
                (c.ljust(w) if justify.get(columns[i], "right" if i else "left") == "left" else c.rjust(w))
                for i, (c, w) in enumerate(zip(cells, widths)))
        _out("")
        _out(f"  {title}")
        _out(fmt(columns))
        _out("  " + "  ".join("-" * w for w in widths))
        for r in rows:
            _out(fmt(r))
        return
    t = Table(title=title, title_style="title", box=box.SIMPLE_HEAD, header_style="secondary",
              border_style=PALETTE["outline"], pad_edge=False, show_edge=False, title_justify="left")
    for i, col in enumerate(columns):
        default = "left" if i == 0 else "right"
        t.add_column(col, justify=justify.get(col, default), no_wrap=(i == 0))
    for i, r in enumerate(rows):
        style = row_styles[i] if row_styles and i < len(row_styles) else None
        t.add_row(*r, style=style)
    console().print(t)


# ── Progress ──────────────────────────────────────────────────────────────────

class _PlainProgress:
    """Stand-in for rich.progress.Progress in plain mode: one line per update of the
    completed count, no redraws."""

    def __init__(self):
        self._tasks = {}
        self._n = 0
        self.console = _PlainConsole()

    def __enter__(self):
        return self

    def __exit__(self, *a):
        return False

    def add_task(self, description, total=None, **kw):
        self._n += 1
        self._tasks[self._n] = {"desc": description, "total": total, "done": 0}
        return self._n

    def update(self, task_id, advance=0, completed=None, total=None, description=None, **kw):
        t = self._tasks[task_id]
        before = t["done"]
        if total is not None:
            t["total"] = total
        if completed is not None:
            t["done"] = completed
        t["done"] += advance
        if description is not None:
            t["desc"] = description
        if t["done"] != before:
            tot = f"/{t['total']}" if t["total"] else ""
            _out(f"  {t['done']}{tot}  {t['desc']}")

    def advance(self, task_id, n=1):
        self.update(task_id, advance=n)


class _PlainConsole:
    def print(self, *a, **kw):
        import re
        text = " ".join(str(x) for x in a)
        _out(re.sub(r"\[/?[a-z#0-9 ]*\]", "", text))


def progress(transient: bool = False):
    if not is_live():
        return _PlainProgress()
    return Progress(
        SpinnerColumn(style=PALETTE["secondary"]),
        TextColumn("[muted]{task.description}"),
        BarColumn(bar_width=28),
        MofNCompleteColumn(),
        TimeElapsedColumn(),
        TimeRemainingColumn(),
        console=console(),
        transient=transient,
    )


@contextmanager
def spinner(message: str):
    if not is_live():
        _out(f"  {message} ...")
        yield
        return
    with console().status(f"[muted]{message}", spinner="dots", spinner_style=PALETTE["secondary"]):
        yield


# ── Panels ────────────────────────────────────────────────────────────────────

def error_panel(title: str, detail: str, tail: str = "", hint: str = "") -> None:
    if _plain:
        _out("")
        _out(f"[ERROR] {title}: {detail}")
        if tail:
            _out(tail.strip()[-2500:])
        if hint:
            _out(hint)
        return
    body = Text()
    body.append(detail + "\n", style="bold")
    if tail:
        body.append("\n" + tail.strip()[-2500:] + "\n", style="faint")
    if hint:
        body.append("\n" + hint, style="tertiary")
    console().print(Panel(body, title=f"[err]{title}[/]", title_align="left", border_style=PALETTE["danger"],
                          box=box.ROUNDED, expand=False, padding=(0, 2)))


def results_panel(title: str, headline: str, files: Sequence[Sequence[str]], hint: str = "") -> None:
    """Closing panel: a headline, then (label, path) rows, then an optional hint. The
    border hugs the content; a path longer than the terminal folds under its column."""
    if _plain:
        _out("")
        _out(f"{title}: {headline}")
        for label, path in files:
            _out(f"  {label:<18} {path}")
        if hint:
            _out(f"  {hint}")
        return
    from rich.console import Group
    parts = [Text(headline, style=f"bold {PALETTE['lime']}")]
    if files:
        grid = Table.grid(padding=(0, 2))
        grid.add_column(style="muted", no_wrap=True)
        grid.add_column(style="tertiary", overflow="fold")
        for label, path in files:
            grid.add_row(label, str(path))
        parts += [Text(""), Text("Output files", style="secondary"), grid]
    if hint:
        parts += [Text(""), Text(hint, style="faint")]
    console().print(Panel(Group(*parts), title=f"[title]{title}[/]", title_align="left",
                          border_style=PALETTE["primary"], box=box.ROUNDED, expand=False, padding=(0, 2)))
