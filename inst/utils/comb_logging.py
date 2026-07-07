#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Uniform logging for the Python pangen steps (comb_06_align, comb_07_mafft,
comb_09_mafft_combine), mirroring the R checkpoint-logging convention defined
in inst/utils/chunk_logging.R.

For every step that receives --path.log this writes:

    <path.log>/script.log   -- the main, human-readable log (always written)

Content: a run header (key parameters), progress, an itemised list of failed
loci (timeouts / bad alignments / missing inputs) and a final summary with
counts and elapsed time. Per-OK-locus lines are intentionally NOT written so
the log stays bounded on inputs with millions of loci; failures always are.

Log levels match the R side (0..3):
    * the file is written whenever --path.log is given, regardless of level;
    * console (stderr) echo of ordinary "main" messages is gated by level>=2,
      matching ll.main in chunk_logging.R; per-item ("loop") messages by >=3.
Callers may force echo (echo=True) for summaries that should always be visible.

No third-party dependencies; safe to import when scripts run standalone
(no --path.log -> logging silently disabled, behaviour unchanged).
"""

import os
import sys
import time

# Mirror ll.main / ll.loop from chunk_logging.R
LL_MAIN = 2
LL_LOOP = 3


def _coerce_level(log_level):
    if log_level is None:
        return 0
    try:
        return int(log_level)
    except (TypeError, ValueError):
        try:
            return int(float(log_level))
        except (TypeError, ValueError):
            return 0


class Logger:
    """Minimal main-process logger writing <path_log>/script.log."""

    def __init__(self, path_log=None, log_level=None, script=None):
        self.path_log = path_log or None
        self.level = _coerce_level(log_level)
        self.echo_main = self.level >= LL_MAIN
        self.echo_loop = self.level >= LL_LOOP
        self.file = None
        self._t0 = time.time()

        if self.path_log:
            try:
                os.makedirs(self.path_log, exist_ok=True)
                self.file = os.path.join(self.path_log, "script.log")
                # create / truncate, matching chunk_logging.R (fresh script.log)
                with open(self.file, "w", encoding="utf-8"):
                    pass
            except OSError as e:
                print("[log-warn] cannot open log dir %s: %s" % (self.path_log, e),
                      file=sys.stderr)
                self.file = None

        if script is not None:
            self.log("=== %s ===" % script, echo=False)

    def _emit(self, msg, echo):
        if self.file is not None:
            try:
                with open(self.file, "a", encoding="utf-8") as fh:
                    fh.write(msg + "\n")
            except OSError:
                pass
        if echo:
            print(msg, file=sys.stderr, flush=True)

    def log(self, msg, echo=None):
        """Main-level message. File always; echo gated by level>=2 unless forced."""
        self._emit(msg, self.echo_main if echo is None else bool(echo))

    def loop(self, msg, echo=None):
        """Per-item message. File always; echo gated by level>=3 unless forced."""
        self._emit(msg, self.echo_loop if echo is None else bool(echo))

    def elapsed(self):
        return time.time() - self._t0


def add_log_args(parser):
    """Attach the standard --path.log / --log.level options (dotted + dashed)."""
    parser.add_argument("--path.log", "--path-log", dest="path_log", default=None,
                        help="Directory for script.log (uniform with the R steps).")
    parser.add_argument("--log.level", "--log-level", dest="log_level", default=None,
                        help="Log level 0..3 (console echo gating).")
    return parser
