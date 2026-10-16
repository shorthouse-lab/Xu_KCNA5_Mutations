## ---------------------------------------------------------------------------
## Timestamped, immediately-flushed progress messages. Plain message()/cat()
## calls get block-buffered when Rscript's stdout is redirected to a log
## file, so nothing appears until the buffer fills or the process exits --
## this forces every line out immediately so `tail -f` on the log is useful.
## ---------------------------------------------------------------------------
log_progress <- function(...) {
  cat(sprintf("[%s] %s\n", format(Sys.time(), "%H:%M:%S"), paste0(...)))
  flush(stdout())
}
