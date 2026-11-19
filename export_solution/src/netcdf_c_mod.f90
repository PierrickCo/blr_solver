!===============================================================================
! Interface minimale vers l'API C de NetCDF (libnetcdf), par ISO_C_BINDING.
! N'utilise PAS l'interface Fortran de NetCDF (netcdf.mod, libnetcdff) : elle
! doit etre compilee avec le meme compilateur que le programme et manque
! souvent ; la bibliotheque C suffit ici.
!
! Constantes : valeurs de netcdf.h (stables depuis NetCDF 4.0).
! Indices et dimensions en convention C : la derniere dimension varie le
! plus vite, et les indices commencent a 0.
!===============================================================================
MODULE NETCDF_C_MOD
  USE ISO_C_BINDING
  IMPLICIT NONE
  PRIVATE
  PUBLIC :: NC_CHECK, NC_CREATE_FILE, NC_DEF_DIM_F, NC_DEF_VAR_F, NC_PUT_ATT_INT_F, &
            NC_PUT_ATT_TEXT_F, NC_ENDDEF_F, NC_PUT_VARA_FLOAT_F, NC_PUT_VARA_INT_F, NC_CLOSE_F
  PUBLIC :: NC_INT, NC_FLOAT, NC_GLOBAL

  INTEGER(C_INT), PARAMETER :: NC_NOERR = 0, NC_CLOBBER = 0, NC_NETCDF4 = 4096
  INTEGER(C_INT), PARAMETER :: NC_INT = 4, NC_FLOAT = 5, NC_CHAR = 2, NC_GLOBAL = -1

  INTERFACE
    INTEGER(C_INT) FUNCTION nc_create(path, cmode, ncidp) BIND(C, NAME='nc_create')
      IMPORT :: C_INT, C_CHAR
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: path(*)
      INTEGER(C_INT), VALUE :: cmode
      INTEGER(C_INT), INTENT(OUT) :: ncidp
    END FUNCTION nc_create
    INTEGER(C_INT) FUNCTION nc_def_dim(ncid, name, len, idp) BIND(C, NAME='nc_def_dim')
      IMPORT :: C_INT, C_CHAR, C_SIZE_T
      INTEGER(C_INT), VALUE :: ncid
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: name(*)
      INTEGER(C_SIZE_T), VALUE :: len
      INTEGER(C_INT), INTENT(OUT) :: idp
    END FUNCTION nc_def_dim
    INTEGER(C_INT) FUNCTION nc_def_var(ncid, name, xtype, ndims, dimidsp, varidp) BIND(C, NAME='nc_def_var')
      IMPORT :: C_INT, C_CHAR
      INTEGER(C_INT), VALUE :: ncid, xtype, ndims
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: name(*)
      INTEGER(C_INT), INTENT(IN) :: dimidsp(*)
      INTEGER(C_INT), INTENT(OUT) :: varidp
    END FUNCTION nc_def_var
    INTEGER(C_INT) FUNCTION nc_put_att_int(ncid, varid, name, xtype, len, op) BIND(C, NAME='nc_put_att_int')
      IMPORT :: C_INT, C_CHAR, C_SIZE_T
      INTEGER(C_INT), VALUE :: ncid, varid, xtype
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: name(*)
      INTEGER(C_SIZE_T), VALUE :: len
      INTEGER(C_INT), INTENT(IN) :: op(*)
    END FUNCTION nc_put_att_int
    INTEGER(C_INT) FUNCTION nc_put_att_text(ncid, varid, name, len, op) BIND(C, NAME='nc_put_att_text')
      IMPORT :: C_INT, C_CHAR, C_SIZE_T
      INTEGER(C_INT), VALUE :: ncid, varid
      CHARACTER(KIND=C_CHAR), INTENT(IN) :: name(*), op(*)
      INTEGER(C_SIZE_T), VALUE :: len
    END FUNCTION nc_put_att_text
    INTEGER(C_INT) FUNCTION nc_enddef(ncid) BIND(C, NAME='nc_enddef')
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE :: ncid
    END FUNCTION nc_enddef
    INTEGER(C_INT) FUNCTION nc_put_vara_float(ncid, varid, startp, countp, op) BIND(C, NAME='nc_put_vara_float')
      IMPORT :: C_INT, C_SIZE_T, C_FLOAT
      INTEGER(C_INT), VALUE :: ncid, varid
      INTEGER(C_SIZE_T), INTENT(IN) :: startp(*), countp(*)
      REAL(C_FLOAT), INTENT(IN) :: op(*)
    END FUNCTION nc_put_vara_float
    INTEGER(C_INT) FUNCTION nc_put_vara_int(ncid, varid, startp, countp, op) BIND(C, NAME='nc_put_vara_int')
      IMPORT :: C_INT, C_SIZE_T
      INTEGER(C_INT), VALUE :: ncid, varid
      INTEGER(C_SIZE_T), INTENT(IN) :: startp(*), countp(*)
      INTEGER(C_INT), INTENT(IN) :: op(*)
    END FUNCTION nc_put_vara_int
    INTEGER(C_INT) FUNCTION nc_close(ncid) BIND(C, NAME='nc_close')
      IMPORT :: C_INT
      INTEGER(C_INT), VALUE :: ncid
    END FUNCTION nc_close
    TYPE(C_PTR) FUNCTION nc_strerror(ncerr) BIND(C, NAME='nc_strerror')
      IMPORT :: C_INT, C_PTR
      INTEGER(C_INT), VALUE :: ncerr
    END FUNCTION nc_strerror
    INTEGER(C_SIZE_T) FUNCTION c_strlen(s) BIND(C, NAME='strlen')
      IMPORT :: C_SIZE_T, C_PTR
      TYPE(C_PTR), VALUE :: s
    END FUNCTION c_strlen
  END INTERFACE

CONTAINS

  ! .TRUE. si STATUS = NC_NOERR ; sinon affiche le message NetCDF et QUOI.
  LOGICAL FUNCTION NC_CHECK(STATUS, QUOI)
    INTEGER(C_INT), INTENT(IN) :: STATUS
    CHARACTER(LEN=*), INTENT(IN) :: QUOI
    CHARACTER(KIND=C_CHAR), POINTER :: MSG(:)
    TYPE(C_PTR) :: P
    INTEGER :: L
    NC_CHECK = (STATUS .EQ. NC_NOERR)
    IF (NC_CHECK) RETURN
    P = nc_strerror(STATUS)
    L = INT(c_strlen(P))
    CALL C_F_POINTER(P, MSG, [L])
    WRITE(0,'(A,A,A,*(A))') ' ERREUR NetCDF (', QUOI, ') : ', MSG(1:L)
  END FUNCTION NC_CHECK

  FUNCTION CSTR(S) RESULT(C)
    CHARACTER(LEN=*), INTENT(IN) :: S
    CHARACTER(KIND=C_CHAR, LEN=LEN_TRIM(S)+1) :: C
    C = TRIM(S)//C_NULL_CHAR
  END FUNCTION CSTR

  ! Fichier NetCDF-4 (format HDF5 : pas de limite de taille de variable).
  INTEGER(C_INT) FUNCTION NC_CREATE_FILE(PATH, NCID)
    CHARACTER(LEN=*), INTENT(IN) :: PATH
    INTEGER(C_INT), INTENT(OUT) :: NCID
    NC_CREATE_FILE = nc_create(CSTR(PATH), IOR(NC_CLOBBER, NC_NETCDF4), NCID)
  END FUNCTION NC_CREATE_FILE

  INTEGER(C_INT) FUNCTION NC_DEF_DIM_F(NCID, NAME, LEN, DIMID)
    INTEGER(C_INT), INTENT(IN) :: NCID
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    INTEGER, INTENT(IN) :: LEN
    INTEGER(C_INT), INTENT(OUT) :: DIMID
    NC_DEF_DIM_F = nc_def_dim(NCID, CSTR(NAME), INT(LEN, C_SIZE_T), DIMID)
  END FUNCTION NC_DEF_DIM_F

  ! DIMIDS en convention C (derniere dimension = la plus rapide).
  INTEGER(C_INT) FUNCTION NC_DEF_VAR_F(NCID, NAME, XTYPE, DIMIDS, VARID)
    INTEGER(C_INT), INTENT(IN) :: NCID, XTYPE, DIMIDS(:)
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    INTEGER(C_INT), INTENT(OUT) :: VARID
    NC_DEF_VAR_F = nc_def_var(NCID, CSTR(NAME), XTYPE, INT(SIZE(DIMIDS), C_INT), DIMIDS, VARID)
  END FUNCTION NC_DEF_VAR_F

  INTEGER(C_INT) FUNCTION NC_PUT_ATT_INT_F(NCID, VARID, NAME, VAL)
    INTEGER(C_INT), INTENT(IN) :: NCID, VARID
    CHARACTER(LEN=*), INTENT(IN) :: NAME
    INTEGER, INTENT(IN) :: VAL
    NC_PUT_ATT_INT_F = nc_put_att_int(NCID, VARID, CSTR(NAME), NC_INT, 1_C_SIZE_T, [INT(VAL, C_INT)])
  END FUNCTION NC_PUT_ATT_INT_F

  INTEGER(C_INT) FUNCTION NC_PUT_ATT_TEXT_F(NCID, VARID, NAME, TXT)
    INTEGER(C_INT), INTENT(IN) :: NCID, VARID
    CHARACTER(LEN=*), INTENT(IN) :: NAME, TXT
    NC_PUT_ATT_TEXT_F = nc_put_att_text(NCID, VARID, CSTR(NAME), INT(LEN_TRIM(TXT), C_SIZE_T), CSTR(TXT))
  END FUNCTION NC_PUT_ATT_TEXT_F

  INTEGER(C_INT) FUNCTION NC_ENDDEF_F(NCID)
    INTEGER(C_INT), INTENT(IN) :: NCID
    NC_ENDDEF_F = nc_enddef(NCID)
  END FUNCTION NC_ENDDEF_F

  ! START (0-based) et COUNT en convention C.
  INTEGER(C_INT) FUNCTION NC_PUT_VARA_FLOAT_F(NCID, VARID, START, COUNT, DATA)
    INTEGER(C_INT), INTENT(IN) :: NCID, VARID
    INTEGER, INTENT(IN) :: START(:), COUNT(:)
    REAL(C_FLOAT), INTENT(IN) :: DATA(*)
    NC_PUT_VARA_FLOAT_F = nc_put_vara_float(NCID, VARID, INT(START, C_SIZE_T), INT(COUNT, C_SIZE_T), DATA)
  END FUNCTION NC_PUT_VARA_FLOAT_F

  INTEGER(C_INT) FUNCTION NC_PUT_VARA_INT_F(NCID, VARID, START, COUNT, DATA)
    INTEGER(C_INT), INTENT(IN) :: NCID, VARID
    INTEGER, INTENT(IN) :: START(:), COUNT(:)
    INTEGER(C_INT), INTENT(IN) :: DATA(*)
    NC_PUT_VARA_INT_F = nc_put_vara_int(NCID, VARID, INT(START, C_SIZE_T), INT(COUNT, C_SIZE_T), DATA)
  END FUNCTION NC_PUT_VARA_INT_F

  INTEGER(C_INT) FUNCTION NC_CLOSE_F(NCID)
    INTEGER(C_INT), INTENT(IN) :: NCID
    NC_CLOSE_F = nc_close(NCID)
  END FUNCTION NC_CLOSE_F

END MODULE NETCDF_C_MOD
