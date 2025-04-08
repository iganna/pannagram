# Common options

## Help

In order to display help information about Pannagrams' usage and available options, use `-h` for a short description and `-help` for an extended description.

## Cores

The flag `-cores` specifies the number of CPU cores to use for parallel processing.

## Restarting Interrupted Processes

If the run is interrupted or stops unexpectedly, restarting **the same command again** will resume execution from where it was interrupted.

> 🧨 **Risky Tips:** You can manually specify the steps to execute.  
> But this is **not recommended** unless you are sure you need to do so.
> - Start at step `N`: use `-s N`.
> - Run only step `N`: use `-s N -one_step`.
> - Run from step `N` to `M`: use `-s N -e M`.
>
> If you want to clean up the output files from the previous run, add the `-clean` flag.

## Logging

Everything is logged into the log files; however, only logs with the specified level are shown in the command line.  
There are four levels of logging for the command line output:

- `-log 0`: nothing is logged.
- `-log 1`: only step progress is logged (default).
- `-log 2`: detailed information about all steps is logged.

If Pannagram encounters any issues, to help us understand the problem, please provide the `${PATH_OUT}logs/` folder.
