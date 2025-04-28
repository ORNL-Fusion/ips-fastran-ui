# Run environment

```bash
source /fusion/projects/codes/ips/local/xfastran/env.setup.bash
```

# Run command

```
xfastran.py [options]
```

# Input

- EFIT directory 
- GAProfile directory

```
test samples:

EFIT directory: /fusion/projects/codes/ips/samples/153648/fit00
GAProfiles directory: /fusion/projects/codes/ips/samples/153648/fit02
```

# Output 

- `instate_<shot>`: initial condtion or profiles for a time slice to be analized
- `innubeam_<shot>`: NUBEAM input
- `intoray_<shot>`: TORAY input
- `t<shot>.nc`: time trace of ne, Te, Ti, rotation, zeff, plasma boundary, PNB, PEC, ... for time dependent analysis

# Example - snapshot analysis

```
xfastran_d3d.py \
--snap \
--shot=153648 \
--efitdir='/fusion/projects/codes/ips/samples/153648/fit00' \
--profdir='/fusion/projects/codes/ips/samples/153648/fit02' \
--time=4250 \
--dtavg=500 \
--rdir='SIMULATION' \
--nb \
--dtbeam_avg=100 \
--ec \
--dtech_avg=100
```

- Profiles will be ~averaged over t = tmin ~ tmax
	- collect profiles between t = tmin and tmax
	- sum of profiles / (number of profiles)

# Example - time dependent analysis

```
xfastran_d3d.py \
--timetrace \
--shot=153648 \
--efitdir='/fusion/projects/codes/ips/samples/153648/fit00' \
--profdir='/fusion/projects/codes/ips/samples/153648/fit02' \
--tmin=4000 \
--tmax=4500 \
--rdir='SIMULATION' \
--nb \
--dtbeam=20 \
--ec \
--dtech=20
```

# Options

`--shot`: shot number, integer

`--efitdir`: EFIT GEQDSK/AEQDSK directory

`--profdir`: GAProfiles directory

`--tmin`: start time [msec, integer]

`--tmax`: end time [msec, integer]

`--dtavg`: profile average time [msec, integer]

`--nb'`: write NB power timetrace

`--dtbeam'`: time interval for beam power output [msec, integer]

`--dtbeam_avg`: time interval for beam power average [msec, integer]

`--ec'`: write ECH power timetrace

`--dtech'`: time interval for ECH power output [msec, integer]

`--dtech_avg`: time interval for ECH power average [msec, integer]
