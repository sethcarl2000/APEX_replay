# APEX_replay

This code is meant to take raw APEX CODA data, and output `.root` files which contain VDC tracks. 

## Replay Workflow

> [!NOTE] 
> All the following terminal commands are executed in the **bash** shell (by default as of writing, ifarm uses **tsch**; to use bash instead, just type the command `bash`). All commands are executed on ifarm from the following directory: `/work/halla/apex/disk1/sethhall/apex_replay`. These commands _should_ still work from another directory, as long as but I have not tested this. 

First, all code is meant to be executed from the following directory: 

`/work/halla/apex/disk1/sethhall/apex_replay`

This should work for another directory; but some parameters must be modified. See `set_apex_replay_env.sh` for lines which must be changed for install in another directory 

### Compiling Source Code
When you clone this repository (on ifarm), you can clone it:
```bash
git clone https://github.com/sethcarl2000/APEX_replay.git
```
And check to make sure you have a compatible version of the [analyzer](https://github.com/JeffersonLab/analyzer). The September 9 2026 replay pass used version `Podd 1.7.15 git@Release-1715-0-ga6e08d60 23 Jul 2025`. Check which version you have pre-installed on ifarm, use: 
```bash
analyzer --version
```
I have only checked the code in its current state with `Podd 1.7.15`. Earlier versions have a known bug with shower-module decoding. I have not tested later versions. 

To compile the source code: 
```bash
cmake -B build -S . 
cmake --build build -j 4
```
No `cmake --install` step is necessary. 


### Setting APEX env 
You must run the command: 
```bash
source set_apex_replay_env.sh
```
to set necessary env variables. See this script for more details 
(FIXME: have cmake auto-generate this script, and remove it from version control). 

### Caching APEX Rawfiles (CODA files)
First, raw files must be cached. This is done using the command: 
```bash
./scripts/cache_rawFiles.sh [first-run] [last-run]
```
Where the integers `[first-run]` and `[last-run]` are the desired run range, the usual choice for APEX production runs is run `3875` to `5025`. For a full list of APEX runs, see [this list](https://hallaweb.jlab.org/experiment/APEX/runlist/runlist.php). 

I get an email when the filesystem is done copying files to tape; to have it sent to another email see the `cache_rawFiles.sh` script. 

For more details on where rawfiles are placed, and how the caching process works, see [cache_rawFiles.sh docs](docs/cache_rawFiles.sh)


### Executing the Replay

To execute the replay, use the following command: 
```bash
./scripts/run_full_replay.sh [first-run] [last-run]
```

This script automatically creates and submits the **slurm task array**. Slurm is the batch-job software that ifarm uses. Specifically, a **task array** is when many jobs are submitted as a single request, each job being indexed by a unique **task-id**.  

Again the typical choice of run range is `3875` to `5025`. It is usually a good idea to submit a smaller subrange, perhaps `3875` to `4050`, to make sure there are no problems with the slurm jobs. 

Upon completon, the _replayed_ output files will be stored under the following directory: 

`/volatile/halla/apex/full_replay/production/replay.[run].rawfile-[rawfile].segment-[segment].root`

Where: 
1. `[run]` is the run number
2. `[rawfile]` is the 0-indexed rawfile number. Each run may have 1 or many rawfiles
3. `[segment]` is the 0-indexed segment number. Each rawfile starts with segment `0`, and increments each time a rawfile needs to be replayed in several chunks. For details on how these chunks are configured, see [`execute_specific_run.C`](execute_specific_run.C). 

The slurm-job output (`stdout`) is saved to: 

`/volatile/halla/apex/full_replay/slurm/full-replay/apex_replay_[first-run]_[last-run].[task-id].[task-array-id].out`

And error printing (`stderr`) is saved to: 

`/volatile/halla/apex/full_replay/slurm/full-replay/apex_replay_[first-run]_[last-run].[task-id].[task-array-id].err`

Where: 

1. `[first-run]` and `[last-run]` are same as before
2. `[task-id]` is the 0-indexed slurm task. See [`scripts/run_full_replay.sh`](scripts/run_full_replay.sh) for details. 
3. `[task-array-id]` is the unique slurm job-id assigned to the entire array. Each slurm task has its own slurm job id, but this identifies the entire array (so you can see which slurms tasks were part of the same task array). 

### Farm Resource Requirements

For the full replay pass executed on 9 September 2026 (most recent, as of writing), the entire pass in the range ~ `3875-5025` was split by `run_full_replay.sh` into **~700 distinct slurm tasks**. To see the specific resource requests, see [`scripts/slurm-replay-task`](scripts/slurm-replay-task#L3-L12). 

In summary, each task requests: 
1. 4 CPUS
2. 5 GB of RAM
3. 10 hours of runtime (even though most use nowhere near this amount)
4. 2500 MB of `/scratch` disk space. See [`slurm-replay-task`](scripts/slurm-replay-task) and [`execute_specific_run.C`](execute_specific_run.C) to see how this is- used.  

For a rough estimate, the 9 Sep 2026 pass used **~4200 CPU-hours** total, which implies a mean walltime of 1-2 hours (but this many not be close to the _median_ walltime). 

### Zipped Code Payload

To ensure that the _exact_ code used for a particular replay pass is saved, each slurm task uses a copy of a compressed archive of the code, generated by `run_full_replay.sh`. These can be found in the auto-generated directory `slurm_payloads`. 

For a full copy of all code / parameters used for the 9 Sep 2026 replay can be found on ifarm under: 

```
slurm_payloads/payload__2026-Sep-09_22.09.27__run-range_3875-4050.tar.gz
slurm_payloads/payload__2026-Sep-09_23.29.16__run-range_4051-5020.tar.gz
```
Here, you can see the code is stored as a `.tar.gz` compressed archive, identified with the run-range and timestamp of the replay passes' submission. 