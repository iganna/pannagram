# Checkpoint-logging pattern (canonical)

Reference spec for parallel R scripts in `inst/pangen/`, `inst/analys/`, etc.
Goal: **logging doubles as the checkpoint system**, with a **minimal number of files**.
Apply this pattern identically in every parallel `%dopar%` script.

Reference implementation: `inst/pangen/comb_01_one_ref.R`.
Shared helpers live in `inst/utils/chunk_logging.R` (sourced by every R script).

---

## Core principles (invariants)

1. **One log file per worker**, not per item → file count = `num.cores`, not thousands.
   File name: `core_<workerID>.log`.
2. **The log file *is* the checkpoint ledger.** A processed item ends with a line
   `Done: <id>`. Progress is rebuilt by scanning these lines.
3. **`markDone` goes AFTER the real output write** — the marker means "fully done".
4. **Writes are idempotent** (`try(<delete>) -> write`), because a crash can happen
   between the output write and `markDone`.
5. **Never wipe output on resume**: create output only if missing.
6. **Checkpoint decisions are log-driven only** (approach A). Accepts one edge:
   "logs survived but output deleted by hand" → silent gap. Fine as long as `clean`
   removes logs and output together.
7. **Resume is robust to a different core count**: `getDoneSet()` globs *all*
   `core_*.log`, so a rerun with 10 cores sees markers written by a 30-core run.

---

## Shared helpers (already in `chunk_logging.R`)

```r
loop.done.prefix <- 'Done:'

# Per-worker log file (bounded number of files). Stable worker id (.worker.id),
# PID fallback for the sequential path.
initLoopLog <- function(path.log, worker.id = NULL) {
  if (is.null(path.log)) return(NULL)
  if (is.null(worker.id)) {
    worker.id <- if (exists('.worker.id', envir = .GlobalEnv))
      get('.worker.id', envir = .GlobalEnv) else Sys.getpid()
  }
  f <- paste0(path.log, 'core_', worker.id, '.log')
  if (!file.exists(f)) invisible(file.create(f))
  return(f)
}

# Set of completed ids from ALL core_*.log (independent of core count).
getDoneSet <- function(path.log, prefix = loop.done.prefix) {
  if (is.null(path.log)) return(character(0))
  files <- list.files(path = path.log, pattern = '^core_.*\\.log$', full.names = TRUE)
  if (length(files) == 0) return(character(0))
  pat <- paste0('^\\s*', prefix, '\\s*')
  done <- unlist(lapply(files, function(f) {
    lines <- grep(pat, readLines(f, warn = FALSE), value = TRUE)
    trimws(sub(pat, '', lines))
  }), use.names = FALSE)
  unique(done)
}

# Write the checkpoint marker "Done: <id>".
markDone <- function(id, file, echo = FALSE) {
  pokaz(paste0(loop.done.prefix, ' ', id), file = file, echo = echo)
}
```

---

## Per-script integration (copy this shape)

`ITEM` = the outer loop variable (e.g. `s.comb`, `f.blast`, `f.maj`, `acc`).
`SUB`  = the inner loop variable if nested checkpoints are wanted (e.g. `acc`).

### 1. loop.function — accept `done.set`, skip done items, open worker log, end with marker

```r
loop.function <- function(ITEM, done.set = character(0), echo.loop = T){

  # Outer checkpoint: skip a fully completed item
  item.id <- ITEM
  if(item.id %in% done.set){
    return(NULL)
  }

  file.log.loop = initLoopLog(path.log)          # one file per worker
  pokaz('...', file=file.log.loop, echo=echo.loop)

  # ... work ...

  markDone(item.id, file=file.log.loop, echo=echo.loop)   # AFTER the work
  return(NULL)
}
```

### 2. Nested checkpoint (optional) — inner loop + idempotent output

```r
  # create output only if missing (never wipe on resume)
  if(!file.exists(file.out)){
    h5createFile(file.out); h5createGroup(file.out, gr.accs.e)
  }

  for(SUB in subitems){
    sub.id <- paste0(item.id, '_', SUB)
    if(sub.id %in% done.set){ next }               # inner checkpoint

    # ... compute ...

    suppressMessages({
      try(h5delete(file.out, dsname(SUB)), silent = TRUE)   # idempotent
      h5write(value, file.out, dsname(SUB))
    })
    markDone(sub.id, file=file.log.loop, echo=echo.loop)    # AFTER the write
  }
```

Resulting log content:
```
Combination: 1 1
Accession acc1 ...
Done: 1_1_acc1
Done: 1_1
```

### 3. Driver — build done.set, assign stable worker ids, pass done.set in

```r
done.set <- getDoneSet(path.log)

if(num.cores == 1){
  assign('.worker.id', 1, envir = .GlobalEnv)          # single stable core_1.log
  for(ITEM in items){
    loop.function(ITEM, done.set = done.set, echo.loop = echo.loop)
  }
} else {
  myCluster <- makeCluster(num.cores, type = "PSOCK")
  registerDoParallel(myCluster)

  # stable worker id 1..N -> bounded, reused file names
  parallel::clusterApply(myCluster, seq_len(num.cores),
                         function(i) assign('.worker.id', i, envir = .GlobalEnv))

  foreach(ITEM = items,
          .packages = c('rhdf5','crayon'),
          .export   = c('done.set')) %dopar% {
    loop.function(ITEM, done.set = done.set, echo.loop = echo.loop)
  }
  stopCluster(myCluster)
}
```

---

## Checklist when converting a file

- [ ] `loop.function` takes `done.set = character(0)`.
- [ ] outer id defined and `if(id %in% done.set) return(NULL)` at the top.
- [ ] `file.log.loop = initLoopLog(path.log)` (drop old `paste0(path.log,'loop_...')`
      and per-file `checkDone`).
- [ ] every output write is idempotent (`try(<delete>) -> write`).
- [ ] `markDone(id, ...)` AFTER the write; nested ids as `paste0(outer,'_',sub)`.
- [ ] output created only if missing (no wipe on resume).
- [ ] driver: `done.set <- getDoneSet(path.log)`, `.worker.id` assigned (seq + cluster),
      `done.set` passed as arg and `.export`ed.

## Notes / gotchas

- `pokaz(..., file=NULL)` silently writes nothing → when `path.log` is NULL there are
  no checkpoints (logging off). Same behavior as before.
- Item ids must be unique strings. For numeric loop vars, `paste`-join them
  (`paste0(query.chr,'_',base.chr)`), matching whatever id the driver iterates on.
- Items skipped for reasons other than completion (e.g. missing input, `next` before
  any write) get NO `Done:` marker — correct, since nothing was produced.
- `h5delete` requires rhdf5 >= ~2.30 (env pins 2.50.0). For non-h5 outputs use the
  matching idempotent removal (`file.remove`, overwrite, etc.).
