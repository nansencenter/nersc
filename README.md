The following briefly goes through the steps of installing FABM-ECOSMO to HYCOM model.

# BGC MODEL SETUP

## Compile FABM coupler and its models

### Clone the models

Create a BGC top folder.

```bash
mkdir -p ${HOME}/FABM
```

You will need the following clones to be able to run ECOSMO with HYCOM:

```bash
cd ${HOME}/FABM/
git clone https://github.com/fabm-model/fabm.git
git clone https://github.com/pmlmodelling/ersem.git
git clone https://github.com/nansencenter/nersc.git
git clone -b v6.0 --recurse-submodules https://github.com/gotm-model/code.git gotm
```

### Installation

```bash
rm -rf ${HOME}/FABM/build
mkdir ${HOME}/FABM/build && cd ${HOME}/FABM/build/
```

For coupling to HYCOM:

```bash
cmake ${HOME}/FABM/fabm \
    -DFABM_HOST=hycom \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DFABM_INSTITUTES="ersem;nersc;gotm" \
    -DFABM_NERSC_BASE=${HOME}/FABM/nersc \
    -DFABM_ERSEM_BASE=${HOME}/FABM/ersem
make install
```

For coupling to GOTM (1-d model):

```bash
cmake ${HOME}/FABM/gotm \
    -DFABM_BASE=${HOME}/FABM/fabm \
    -DCMAKE_Fortran_COMPILER=ifort \
    -DFABM_INSTITUTES="ersem;nersc;gotm" \
    -DFABM_NERSC_BASE=${HOME}/FABM/nersc \
    -DFABM_ERSEM_BASE=${HOME}/FABM/ersem
make install
```

## Experiment folder setup to run BGC model

### EXPT.src file

Add the following lines to EXPT.src file:

```
export COMPILE_BIOMODEL="yes"
```

**IMPORTANT** If you are running the operational model or activating sea-ice algae, add the following line as well:

```
export IA_DRIFT="yes"
```

### blkdat.input file

Make sure the following entries in blkdat.input files are as follows:

```
1      'ntracr' = number of tracers (0=none,negative (i.e. -1) to initialize from climatology)
1      'trcrlx' = activate lat. bound. tracer nudging  (0=F,1=T)
```

### hycom_fabm.nml file

If you do not have hycom_fabm.nml file, create it. Make sure hycom_fabm.nml looks like the following, **especially** for the operational model. For other use cases, for example if you are not running with nesting, remove *nested_variables* line. If you are not running with ice-algae, set *do_icealgae = .false.*. Other lines can be set to *.false.* if you are debugging the code.

```
&hycom_fabm
  do_vertical_movement = .true.
  do_interior_sources = .true.
  do_bottom_sources = .true.
  do_surface_sources = .true.
  do_icealgae = .true.
  nested_variables = 'ECO_no3','ECO_pho','ECO_sil'
/
```

### fabm.yaml file

Copy the *fabm.yaml* file into experiment folder.
If you are running the **operational model**:

```bash
cp $HOME/FABM/nersc/ecosmo/fabm.yaml.operational ./fabm.yaml
```

If you are **NOT** running the operational model:

```bash
cp $HOME/FABM/nersc/ecosmo/fabm.yaml ./fabm.yaml
```

## Compile HYCOM-FABM-ECOSMO
```bash
$HOME/NERSC-HYCOM-CICE/bin/compile_model.sh ifort -u
```
## Boundary conditions and forcing files (needs revision)

You are expected to have a working copy of relaxation and river forcing files. If this is not your first time, just copy the relaxation and river old experiment number folder as you current experiment number folder. Otherwise, execute:

```bash
./create_ref_case.sh
```

### Nesting instructions
