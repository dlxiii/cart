## readme

The initial values of concentration C(t, x) and age concentration β(t, x) were both set to 0 in the whole domain. At the releasing point (x=xr), the concentration C(t, x) was always set to 1, that is, the water particle was released continuously at x=xr . The age concentration at x=x r was set to 0, resulting in a zero age of water particle at x=xr  (Bolin and Rodhe, 1973; Takeoka, 1984). At the two ends (x=0 and x=L), the concentration C(t, x) and age concentrationβ(t, x) both were set to 0, indicating that the water particle could not re-enter the model domain.

## code modification

### mod_dye.F

* add variables of dye_age

~~~
REAL(SP), ALLOCATABLE, TARGET :: DYE_AGE(:,:)
REAL(SP), ALLOCATABLE         :: DYEF_AGE(:,:)
~~~

If defined "MDATA", also ""DYE_S_AGE(:,:) and "DYE_SF_AGE(:,:)" should be allocated.

* allocate variables of dye_age

~~~
ALLOCATE(DYE_AGE(0:MT,KB))   ;DYE_AGE  = ZERO
ALLOCATE(DYEF_AGE(0:MT,KB))  ;DYEF_AGE = ZERO
~~~

* add several subroutines

~~~
SUBROUTINE ADV_DYE_AGE
SUBROUTINE VDIF_DYE_AGE(F_AGE)
SUBROUTINE BCOND_DYE_AGE
~~~

### mod_main.F

* add variables for partial difference

~~~
REAL(SP), ALLOCATABLE,TARGET :: PFPXB_AGE(:)   
REAL(SP), ALLOCATABLE,TARGET :: PFPYB_AGE(:) 
~~~

* allocate variables for partial difference

~~~
ALLOCATE(PFPXB_AGE(MT))  ;PFPXB_AGE = ZERO
ALLOCATE(PFPYB_AGE(MT))  ;PFPYB_AGE = ZERO
~~~

### mod_ncdio.F

* allocate variables for output

~~~
allocate(DYE_AGE(MGL,KB),stat=status)
IF (STATUS /=0 ) CALL FATAL_ERROR("COULD NOT ALLOCATE MEMORY ON IO PROC FOR OUTPUT DATA:EL")
DYE_AGE = 0.0_SP
~~~

* print out values of variable

~~~
! Dye_AGE
VAR  => NC_MAKE_AVAR(name='DYE_AGE',&
        & values=DYE_AGE, DIM1= DIM_node, DIM2= DIM_siglay, DIM3=DIM_time)

ATT  => NC_MAKE_ATT(name='long_name',values='DYE_AGE') 
VAR  => ADD(VAR,ATT)

ATT  => NC_MAKE_ATT(name='units',values='meters') 
VAR  => ADD(VAR,ATT)

ATT  => NC_MAKE_ATT(name='positive',values='up') 
VAR  => ADD(VAR,ATT)

ATT  => NC_MAKE_ATT(name='standard_name',values='DYE_AGE') 
VAR  => ADD(VAR,ATT)

# if defined(UCF)
    ATT  => NC_MAKE_ATT(name='grid',values='grid4') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values="lat lon") 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='grid_location',values='node') 
    VAR  => ADD(VAR,ATT)
# else
    ATT  => NC_MAKE_ATT(name='grid',values='SigmaLayer_Mesh') 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='coordinates',values=CoordVar) 
    VAR  => ADD(VAR,ATT)

    ATT  => NC_MAKE_ATT(name='type',values='data') 
    VAR  => ADD(VAR,ATT)
# endif

NCF  => ADD(NCF,VAR)
~~~

### internal_step.F

* call subroutines for dye age

~~~
CALL ADV_DYE_AGE
CALL VDIF_DYE(DYEF_AGE)
CALL BCOND_DYE_AGE
DYE_AGE = DYEF_AGE 