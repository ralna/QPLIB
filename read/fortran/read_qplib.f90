! THIS VERSION: 31/01/2017 AT 11:50 GMT
! components also from GALAHAD SMT and GALAHAD QPT

!-*-*-*-*-*-*-*-*-*- G A L A H A D _ R P D   M O D U L E -*-*-*-*-*-*-*-*-

!  Copyright reserved, Gould/Orban/Toint, for GALAHAD productions
!  Principal author: Nick Gould

!  History -
!   originally released with GALAHAD Version 2.0. January 22nd 2006
!   extraced from GALAHAD 3.0, 14th December 2016

!  For full documentation, see
!   http://galahad.rl.ac.uk/galahad-www/specs.html

   MODULE GALAHAD_RPD_double

!     -------------------------------------------------
!     |                                               |
!     | Read data for the linear program (LP)         |
!     |                                               |
!     |    minimize           g(T) x + f              |
!     |    subject to     c_l <= A x <= c_u           |
!     |                   x_l <=  x  <= x_u           |
!     |                                               |
!     | the linear program with quadratic             |
!     | constraints (QCP)                             |
!     |                                               |
!     |    minimize       g(T) x + f                  |
!     |    subject to c_l <= A x +                    |
!     |            1/2 vec( x . H_c . x ) <= c_u      |
!     |                   x_l <=  x  <= x_u           |
!     |                                               |
!     | the bound-constrained quadratic program (BQP) |
!     |                                               |
!     |    minimize     1/2 x(T) H x + g(T) x + f     |
!     |    subject to     x_l <=  x  <= x_u           |
!     |                                               |
!     | the quadratic program (QP)                    |
!     |                                               |
!     |    minimize     1/2 x(T) H x + g(T) x + f     |
!     |    subject to     c_l <= A x <= c_u           |
!     |                   x_l <=  x  <= x_u           |
!     |                                               |
!     | or the quadratic program with quadratic       |
!     | constraints (QCQP)                            |
!     |                                               |
!     |    minimize     1/2 x(T) H x + g(T) x + f     |
!     |    subject to c_l <= A x +                    |
!     |            1/2 vec( x . H_c . x ) <= c_u      |
!     |                   x_l <=  x  <= x_u           |
!     |                                               |
!     | where vec( x . H_c . x ) is the vector        |
!     | whose ith component is  x(T) (H_c)_i x for    |
!     | the i-th constraint, from a problem-data file |
!     |                                               |
!     -------------------------------------------------

      IMPLICIT NONE

      PRIVATE
      PUBLIC :: RPD_read_problem_data

!--------------------
!   P r e c i s i o n
!--------------------

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )

!----------------------
!   P a r a m e t e r s
!----------------------

      REAL ( KIND = wp ), PARAMETER :: zero = 0.0_wp
      INTEGER, PARAMETER :: input_line_length = 256
      INTEGER, PARAMETER :: qp = 1
      INTEGER, PARAMETER :: qcqp = 2
      INTEGER, PARAMETER :: bqp = 3
      INTEGER, PARAMETER :: lp = 4
      INTEGER, PARAMETER :: qcp = 5

!-------------------------------------------------
!  D e r i v e d   t y p e   d e f i n i t i o n s
!-------------------------------------------------

!  ----------------------------
!  a sparse matrix derived type
!  ----------------------------

      TYPE, PUBLIC :: SMT_type
        INTEGER :: m                                 !  row dimension
        INTEGER :: n                                 !  column dimension
        INTEGER :: ne                                !  number of entries
        CHARACTER, ALLOCATABLE, DIMENSION( : ) :: id !  matrix id (not used)
        CHARACTER, ALLOCATABLE, DIMENSION( : ) :: type  !  matrix type
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: row !  row indices of entries
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: col  ! column indices of entries
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: ptr  ! array of other indices
                                                     ! (e.g. constraint numbers)
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: val ! values
      END TYPE

!  ------------------------------------
!  a quadratic programming derived type
!  ------------------------------------

      TYPE, PUBLIC :: QPT_problem_type
        INTEGER :: n   ! number of variables
        INTEGER :: m   ! number of constraints
        REAL ( KIND = wp ) :: f = 0.0_wp    ! constant term
        REAL ( KIND = wp ) :: infinity = ( 10.0_wp ) ** 20 ! bound infinity

!  allocatable arrays

        CHARACTER, ALLOCATABLE, DIMENSION( : ) :: name   !  problem name
        CHARACTER ( len = 10 ), ALLOCATABLE, DIMENSION( : ) :: X_names
        CHARACTER ( len = 10 ), ALLOCATABLE, DIMENSION( : ) :: C_names
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: G   ! linear terms
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X_l ! variable
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X_u ! bounds
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: C_l ! constraint
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: C_u ! bounds
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: X  ! primal variables
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y  ! Lagrange mults
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Z  ! dual variables
        TYPE ( SMT_type ) :: A    !  constraint matrix
        TYPE ( SMT_type ) :: H    !  objective Hessian
        TYPE ( SMT_type ) :: H_c  !  constraint Hessians

!  values and arrays not used in read_qplib

        INTEGER :: x_free, x_l_start, x_l_end, x_u_start, x_u_end
        INTEGER :: h_diag_end_free, h_diag_end_nonneg, h_diag_end_nonpos
        INTEGER :: h_diag_end_lower, h_diag_end_range, h_diag_end_upper
        INTEGER :: h_diag_end_fixed, c_equality, c_l_end, c_u_start
        INTEGER :: c_u_end, Hessian_kind, target_kind, gradient_kind
        REAL ( KIND = wp ) :: df, q, theta_max, theta, rho_g, rho_b
        LOGICAL :: new_problem_structure
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: C_status
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: X_status
        INTEGER, ALLOCATABLE, DIMENSION( : ) :: X_type
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DX_l, DX_u
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DC_l, DC_u
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Z_l, Z_u
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DZ_l, DZ_u
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: Y_l, Y_u
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DY_l, DY_u
        REAL ( KIND = wp ), ALLOCATABLE, DIMENSION( : ) :: DG, C, WEIGHT

!  arrays not used here

      END TYPE

!  --------------------
!  returned information
!  --------------------

      TYPE, PUBLIC :: RPD_inform_type
        INTEGER :: status   ! 0 = OK, -2 = allocation failure, -3 = end of file,
                            ! -4 = other read error, -5 = unrecognised type
        INTEGER :: alloc_status  ! status from last allocation attempt
        INTEGER :: io_status     ! status from last read attempt
        INTEGER :: line          ! number of last line read
        CHARACTER ( LEN = 10 ) :: bad_alloc = REPEAT( ' ', 10 ) ! last array
                                                        ! allocation attempt
      END TYPE

   CONTAINS

!-*-*-   R P D _ R E A D _ P R O B L E M _ D A T A   S U B R O U T I N E   -*-*-

      SUBROUTINE RPD_read_problem_data( input, prob, inform )
      INTEGER, INTENT( IN ) :: input
      TYPE ( QPT_problem_type ), INTENT( INOUT ) :: prob
      TYPE ( RPD_inform_type ), INTENT( OUT ) :: inform

!  Read the problem-data file from unit input into the derived type prob
!  (see above for components of inform, and GALAHAD_qpt for those of prob)

!  ****************************************************************************

!  For the linear program (LP)

!    minimize           g(T) x + f

!    subject to     c_l <= A x <= c_u
!                   x_l <= x <= x_u

!  the linear program with quadratic constraints (QCP)

!    minimize           g(T) x + f

!    subject to     c_l <= A x + 1/2 vec( x . H_c . x ) <= c_u
!                   x_l <= x <= x_u

!  the bound-constrained quadratic program (BQP)

!    minimize     1/2 x(T) H x + g(T) x + f

!    subject to     x_l <= x <= x_u

!  the quadratic program (QP)

!    minimize     1/2 x(T) H x + g(T) x + f

!    subject to     c_l <= A x <= c_u
!                   x_l <= x <= x_u

!  or the quadratic program with quadratic constraints (QCQP)

!    minimize     1/2 x(T) H x + g(T) x + f

!    subject to     c_l <= A x + 1/2 vec( x . H_c . x ) <= c_u
!                   x_l <= x <= x_u

!  where vec( x . H_c . x ) is the vector whose ith component is  x(T) H_c x
!  for the i-th constraint. Variables may be continuous or integer

!  ****************************************************************************

!  The data should be input in a file on unit 5. The data is in free format
!  (blanks separate values), but must occur in the order given here (depending
!  on the precise form of problem under consideration, certain data is not
!  required and should not be provided, see below). Any blank lines, or lines
!  starting with any of the  characters "!", "%" or "#" are ignored. Each term
!  in "quotes" denotes a required value. Any strings beyond those required on
!  a given lines will be regarded as comments and ignored.

!  "problem name"
!  "problem type"
!  "problem sense" i.e. one of the words minimize or maximize (case irrelevant)
!  "number variables, n"
!  "number general linear constraints, m"                                   [1]
!  "number of nonzeros in upper triangle of H"                              [2]
!  "row" "column" "value" for each entry of H (if any), one triple on each line
!  "default value for entries in g"
!  "number of non-default entries in g"
!  "index" "value" for each non-default term in g (if any), one pair per line
!  "value of f"
!  "number of nonzeros in upper triangles of H_c"                         [1,3]
!  "constraint" "row" "column" "value" for each entry of H_c (if any),
!    one quadruple on each line
!  "number of nonzeros in A"                                                [1]
!  "row" "column" "value" for each entry of A (if any), one triple on each line
!  "value for infinity" for bounds - any bound exceeding this in absolute value
!     is infinite
!  "default value for entries in c_l"                                       [1]
!  "number of non-default entries in c_l"                                   [1]
!  "index" "value" for each non-default term in c_l (if any), one pair per line
!  "default value for entries in c_u"                                       [1]
!  "number of non-default entries in c_u"                                   [1]
!  "index" "value" for each non-default term in c_u (if any), one pair per line
!  "default value for entries in x_l"
!  "number of non-default entries in x_l"
!  "index" "value" for each non-default term in x_l (if any), one pair per line
!  "default value for entries in x_u"
!  "number of non-default entries in x_u"
!  "index" "value" for each non-default term in x_u (if any), one pair per line
!  "default variable type"  (0 for a continuous variable, 1 for an integer one
!     and 2 for a binary one)                                               [4]
!  "number of non-default variables"                                        [4]
!  "index" "value" for each non-default variable type (if any), one pair/line
!  "default value for starting value for variables x"
!  "number of non-default starting entries in x"
!  "index" "value" for each non-default term in x (if any), one pair per line
!  "default value for starting value for Lagrange multipliers y for constraints"
!                                                                           [1]
!  "number of non-default starting entries in y"                            [1]
!  "index" "value" for each non-default term in y (if any), one pair per line
!  "default value for starting value for dual variables z for simple bounds"
!  "number of non-default starting entries in z"
!  "index" "value" for each non-default term in z (if any), one pair per line
!  "number of non-default names of variables" - default for variable i is "xi"
!  "index" "name" for each non-default name for variable x_i with index i
!    (if any)
!  "number of non-default names of constraints" - default for constraint i is
!    "ci"
!  "index" "name" for each non-default name for constraint with index i (if any)

!  The "problem type" is a string of three characters.

!  The first character indicates the type of objective function used.
!  It must be one of the following:

!   L  a linear objective function
!   D  a convex quadratic objective function whose Hessian is a diagonal matrix
!   C  a convex quadratic objective function
!   Q  a quadratic objective function whose Hessian may be indefinite

!  The second character indicates the types of variables that are present.
!  It must be one of the following:

!   C  all the variables are continuous
!   B  all the variables are binary (0-1)
!   M  the variables are a mix of continuous and binary
!   I  all the variables are integer
!   G  the variables are a mix of continuous, binary and integer

!  The third character indicates the type of the (most extreme)
!  constraint function used; other constraints may be of a lesser type.
!  It must be one of the following:

!   N  there are no constraints
!   B  some of the variables lie between lower and upper bounds (box constraint)
!   L  the constraint functions are linear
!   D  the constraint functions are convex quadratics with diagonal Hessians
!   C  the constraint functions are convex quadratics
!   Q  the constraint functions are quadratics whose Hessians may be indefinite

!  Thus for continuous problems, we would have

!    LCL            a linear program
!    LCC or LCQ     a linear program with quadratic constraints
!    CCB or QCB     a bound-constrained quadratic program
!    CCL or QCL     a quadratic program
!    CCC or CCQ or  a quadratic program with quadratic constraints
!    QCC or QCQ

!  For integer problems, the second character would be I rather than C,
!  and for mixed integer problems, the second character would by M or G.

!  [1] for bound-constrained QPs, these sections are omitted
!  [2] for linear program with quadratic constraints, this section is omitted
!  [3] for problems without quadratic constraints, this section is omitted
!  [4] for purely-continuous, purely-integer & binary problems, this section
!      is omitted.

!  *****************************************************************************

!  Local variables

     INTEGER :: i, ic, j, k, A_ne, H_ne, H_c_ne, nnzx_0, nnzy_0, nnzz_0
     INTEGER :: nnzg, nnzc_l, nnzc_u, nnzx_l, nnzx_u, smt_stat, ip, i_default
     INTEGER :: problem_type
     REAL ( KIND = wp ) :: rv, default
     LOGICAL :: objmax, binary
     CHARACTER ( LEN = 10 ) :: pname, cv
     CHARACTER ( LEN = 24 ) :: p_type
     CHARACTER ( LEN = input_line_length ) :: input_line, blank_line

     inform%line = 0
     inform%alloc_status = 0

!    DO i = 1, input_line_length
!      blank_line( i : i ) = ' '
!    END DO
     blank_line = REPEAT( ' ', input_line_length )

!  Determine the problem name

!    pname = '          '
     pname = REPEAT( ' ', 10 )
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) pname
     ALLOCATE( prob%name( 10 ) )
     prob%name = TRANSFER( pname, prob%name )

!  Determine the problem type

     p_type = REPEAT( ' ', 24 )
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) p_type
     p_type = ADJUSTL( p_type )
     CALL STRING_lower_word( p_type( 1 : 3 ) )

     IF ( p_type( 2 : 2 ) == 'm' .OR. p_type( 2 : 2 ) == 'g' ) THEN
       ip = 2
     ELSE IF ( p_type( 2 : 2 ) == 'i' .OR. p_type( 2 : 2 ) == 'b' ) THEN
       ip = 1
     ELSE
       ip = 0
     END IF
     binary =  p_type( 2 : 2 ) == 'b'

     IF ( p_type( 3 : 3 ) == 'q' .OR. p_type( 3 : 3 ) == 'c' .OR.              &
          p_type( 3 : 3 ) == 'd' ) THEN
       IF ( p_type( 1 : 1 ) == 'l' ) THEN
         problem_type = qcp
       ELSE
         problem_type = qcqp
       END IF
     ELSE IF ( p_type( 3 : 3 ) == 'l' ) THEN
       IF ( p_type( 1 : 1 ) == 'l' ) THEN
         problem_type = lp
       ELSE
         problem_type = qp
       END IF
     ELSE IF ( p_type( 3 : 3 ) == 'b' ) THEN
       problem_type = bqp
     ELSE
       GO TO 950
     END IF

!  Determine if the problem is a minimization or maximization one

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO

     CALL STRING_lower_word( input_line( 1 : 8 ) )
     IF ( input_line( 1 : 8 ) == 'maximize' ) THEN
       objmax = .TRUE.
     ELSE
       objmax = .FALSE.
     END IF

!  Determine the number of variables and constraints

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) prob%n
      IF ( problem_type /= bqp ) THEN
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) prob%m
     ELSE
       prob%m = 0
     END IF

!  Allocate suitable arrays

     ALLOCATE( prob%X( prob%n ), prob%X_l( prob%n ), prob%X_u( prob%n ),       &
               prob%G( prob%n ), prob%Z( prob%n ), STAT = inform%alloc_status )
     IF ( inform%alloc_status /= 0 ) THEN
       inform%status = - 2 ; inform%bad_alloc = 'X' ; RETURN
     END IF

     ALLOCATE( prob%C_l( prob%m ), prob%C_u( prob%m ), prob%Y( prob%m ),       &
               prob%C( prob%m ), STAT = inform%alloc_status )
     IF ( inform%alloc_status /= 0 ) THEN
       inform%status = - 2 ; inform%bad_alloc = 'Y' ; RETURN
     END IF

!  Fill component H

     IF ( problem_type == qp .OR. problem_type == bqp .OR.                     &
           problem_type == qcqp ) THEN
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO

       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) prob%H%ne
       ALLOCATE( prob%H%row( prob%H%ne + prob%n ),                             &
                 prob%H%col( prob%H%ne + prob%n ),                             &
                 prob%H%val( prob%H%ne + prob%n ), STAT = inform%alloc_status )
       IF ( inform%alloc_status /= 0 ) THEN
         inform%status = - 2 ; inform%bad_alloc = 'H' ; RETURN
       END IF

       H_ne = 0
       DO k = 1, prob%H%ne
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, j, rv
         IF ( rv == zero ) CYCLE
         IF ( objmax ) rv = - rv
         H_ne = H_ne + 1 ; prob%H%val( H_ne ) = rv
         IF ( i >= j ) THEN
           prob%H%row( H_ne ) = i
           prob%H%col( H_ne ) = j
         ELSE
           prob%H%row( H_ne ) = j
           prob%H%col( H_ne ) = i
         END IF
       END DO
     ELSE
       H_ne = 0
     END IF
     prob%H%ne = H_ne ; prob%H%m = prob%n ; prob%H%n = prob%n
     IF ( ALLOCATED( prob%H%type ) ) DEALLOCATE( prob%H%type )
     CALL SMT_put( prob%H%type, 'COORDINATE', smt_stat )

!  Fill component g

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     IF ( objmax ) default = - default
     prob%G = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzg
     DO k = 1, nnzg
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       IF ( objmax ) rv = - rv
       prob%G( i ) = rv
     END DO

!  Fill component f

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) prob%f
     IF ( objmax ) prob%f = - prob%f

!  Fill component H_c

     IF ( problem_type == qcqp .OR. problem_type == qcp ) THEN
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO

       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) prob%H_c%ne
       ALLOCATE( prob%H_c%row( prob%H_c%ne ), prob%H_c%col( prob%H_c%ne ),     &
                 prob%H_c%ptr( prob%H_c%ne ), prob%H_c%val( prob%H_c%ne ),     &
                 STAT = inform%alloc_status )
       IF ( inform%alloc_status /= 0 ) THEN
         inform%status = - 2 ; inform%bad_alloc = 'H_c' ; RETURN
       END IF

       H_c_ne = 0
       DO k = 1, prob%H_c%ne
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) ic, i, j, rv
         IF ( rv == zero ) CYCLE
         H_c_ne = H_c_ne + 1 ; prob%H_c%val( H_c_ne ) = rv
         prob%H_c%ptr( H_c_ne ) = ic
         IF ( i >= j ) THEN
           prob%H_c%row( H_c_ne ) = i
           prob%H_c%col( H_c_ne ) = j
         ELSE
           prob%H_c%row( H_c_ne ) = j
           prob%H_c%col( H_c_ne ) = i
         END IF
       END DO
     ELSE
       H_c_ne = 0
     END IF
     prob%H_c%ne = H_c_ne ; prob%H_c%m = prob%n ; prob%H_c%n = prob%n
     IF ( ALLOCATED( prob%H_c%type ) ) DEALLOCATE( prob%H_c%type )
     CALL SMT_put( prob%H_c%type, 'COORDINATE', smt_stat )

!  Fill component A

     IF ( problem_type /= bqp ) THEN
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) prob%A%ne
       ALLOCATE( prob%A%row( prob%A%ne ), prob%A%col( prob%A%ne ),             &
                 prob%A%val( prob%A%ne ), STAT = inform%alloc_status )
       IF ( inform%alloc_status /= 0 ) THEN
         inform%status = - 2 ; inform%bad_alloc = 'A' ; RETURN
       END IF

       A_ne = 0
       DO k = 1, prob%A%ne
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, j, rv
         IF ( rv == zero ) CYCLE
         A_ne = A_ne + 1 ; prob%A%val( A_ne ) = rv
         prob%A%row( A_ne ) = i ; prob%A%col( A_ne ) = j
       END DO
     ELSE
       A_ne = 0
       ALLOCATE( prob%A%row( A_ne ), prob%A%col( A_ne ),                       &
                 prob%A%val( A_ne ), STAT = inform%alloc_status )
       IF ( inform%alloc_status /= 0 ) THEN
         inform%status = - 2 ; inform%bad_alloc = 'A' ; RETURN
       END IF
     END IF
     prob%A%ne = A_ne ; prob%A%m = prob%m ; prob%A%n = prob%n
     IF ( ALLOCATED( prob%A%type ) ) DEALLOCATE( prob%A%type )
     CALL SMT_put( prob%A%type, 'COORDINATE', smt_stat )

!  Fill component infinity

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) prob%infinity

!  Fill component c_l

     IF ( problem_type /= bqp ) THEN
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) default
       prob%C_l = default
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) nnzc_l
       DO k = 1, nnzc_l
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, rv
         prob%C_l( i ) = rv
       END DO

!  Fill component c_u

       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) default
       prob%C_u = default
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) nnzc_u
       DO k = 1, nnzc_u
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, rv
         prob%C_u( i ) = rv
       END DO
     END IF

!  for binary problems, set the lower and upper bounds

     IF ( binary ) THEN
       prob%X_l( : prob%n ) = 0.0_wp
       prob%X_u( : prob%n ) = 1.0_wp

!  otherwise, fill component x_l

     ELSE
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) default
       prob%X_l = default
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) nnzx_l
       DO k = 1, nnzx_l
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, rv
         prob%X_l( i ) = rv
       END DO

!  Fill component x_u

       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) default
       prob%X_u = default
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) nnzx_u
       DO k = 1, nnzx_u
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, rv
         prob%X_u( i ) = rv
       END DO
     END IF

!  Fill component x_type

     ALLOCATE( prob%X_type( prob%n ), STAT = inform%alloc_status )
     IF ( inform%alloc_status /= 0 ) THEN
       inform%status = - 2 ; inform%bad_alloc = 'X_type' ; RETURN
     END IF
     IF ( ip == 2 ) THEN
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i_default
       prob%X_type = i_default
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) nnzx_0
       DO k = 1, nnzx_0
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, j
         prob%X_type( i ) = j
       END DO
     ELSE IF ( ip == 1 ) THEN
       prob%X_type = 1
     ELSE
       prob%X_type = 0
     END IF

!  Fill component x

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%X = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzx_0
     DO k = 1, nnzx_0
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%X( i ) = rv
     END DO

!  Fill component y

     IF ( problem_type /= bqp ) THEN
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) default
       prob%Y = default
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) nnzy_0
       DO k = 1, nnzy_0
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, rv
         prob%Y( i ) = rv
       END DO
     END IF

!  Fill component z

     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) default
     prob%Z = default
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzz_0
     DO k = 1, nnzz_0
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, rv
       prob%Z( i ) = rv
     END DO

!  Fill component x_names

     ALLOCATE( prob%X_names( prob%n ), prob%C_names( prob%m ),                 &
               STAT = inform%alloc_status )
     IF ( inform%alloc_status /= 0 ) THEN
       inform%status = - 2 ; inform%bad_alloc = 'X_names' ; RETURN
     END IF

     DO i = 1, prob%n
       prob%X_names( i ) = 'x' // REPEAT( ' ', 9 )
       WRITE( prob%X_names( i )( 2 : 10 ), "( I0 )" ) i
     END DO
     DO
       inform%line = inform%line + 1
       input_line = blank_line
       READ( input, "( A )", END = 930, ERR = 940 ) input_line
       IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
     END DO
     READ( input_line, *, IOSTAT = inform%io_status,                           &
           END = 930, ERR = 940 ) nnzx_0
     DO k = 1, nnzx_0
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) i, cv
       prob%X_names( i ) = cv
     END DO

!  Fill component c_names

     IF ( problem_type /= bqp ) THEN
       DO i = 1, prob%m
         prob%C_names( i ) = 'c' // REPEAT( ' ', 9 )
         WRITE( prob%C_names( i )( 2 : 10 ), "( I0 )" ) i
       END DO
       DO
         inform%line = inform%line + 1
         input_line = blank_line
         READ( input, "( A )", END = 930, ERR = 940 ) input_line
         IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
       END DO
       READ( input_line, *, IOSTAT = inform%io_status,                         &
             END = 930, ERR = 940 ) nnzy_0
       DO k = 1, nnzy_0
         DO
           inform%line = inform%line + 1
           input_line = blank_line
           READ( input, "( A )", END = 930, ERR = 940 ) input_line
           IF ( .NOT. RPD_ignore_string( input_line ) ) EXIT
         END DO
         READ( input_line, *, IOSTAT = inform%io_status,                       &
               END = 930, ERR = 940 ) i, cv
         prob%C_names( i ) = cv
       END DO
     END IF

!  - successful execution

     inform%status = 0
     RETURN

!  Error returns

!  - end of file encountered

 930 CONTINUE
     inform%status = 3
     RETURN

!  - other read error encountered

 940 CONTINUE
     inform%status = 4
     RETURN

!  - problem type unrecognised

 950 CONTINUE
     inform%status = 5
     RETURN

!  End of RPD_read_problem_data

     END SUBROUTINE RPD_read_problem_data

!-*-*-*-*-*-   R P D _ I G N O R E _ S T R I N G   F U N C T I O N   -*-*-*-*-*-

     FUNCTION RPD_ignore_string( input_line )
     LOGICAL :: RPD_ignore_string

!  Ignore a string if it is (a) blank or (b) starts with "!", "%" or "#"

     CHARACTER ( LEN = input_line_length ), INTENT( IN ) :: input_line

!  Local variables

     INTEGER :: i, length_string

     length_string = LEN_TRIM( input_line )
     IF ( length_string <= 0 ) THEN
       RPD_ignore_string = .TRUE.
       RETURN
     END IF

     DO i = 1, length_string
       IF ( input_line( i : i ) == ' ' ) CYCLE
       IF ( input_line( i : i ) == '!' .OR. input_line( i : i ) == '#' .OR.    &
            input_line( i : i ) == '%' .OR. input_line( i : i ) == '|' ) THEN
         RPD_ignore_string = .TRUE.
         RETURN
       END IF
       EXIT
     END DO
     RPD_ignore_string = .FALSE.

     RETURN

!  End of RPD_ignore_string

     END FUNCTION RPD_ignore_string

!-*-*-*-*-*-*-*-*-*-  S M T _ P U T   S U B R O U T I N E  -*-*-*-*-*-*-*-*-*-

     SUBROUTINE SMT_put( array, string,stat )
     CHARACTER, ALLOCATABLE :: array( : )
     CHARACTER( * ), INTENT( in ) ::  string
     INTEGER, INTENT( OUT ) ::  stat

!  put string into an allocatble array

     INTEGER :: i, l

     l = LEN_TRIM( string )
     IF ( ALLOCATED( array ) ) THEN
        DEALLOCATE( array, STAT = stat )
        IF ( stat /=0 ) RETURN
     END IF
     ALLOCATE( array( l ), STAT = stat )
     IF ( stat /= 0 ) RETURN
     DO i = 1, l
       array( i ) = string( i : i )
     END DO

!  end of subroutine SMT_put

     END SUBROUTINE SMT_put

!-*-*-*-*-*-*-*-   S T R I N G _ l o w e r   S U B R O U T I N E  -*-*-*-*-*-*-

     SUBROUTINE STRING_lower( string )

!  Convert a character variable from upper to lower case

!--------------------------------
!   D u m m y   A r g u m e n t
!--------------------------------

     CHARACTER, INTENT( INOUT ) :: string

!  Local variables

     INTEGER :: letter
     CHARACTER, DIMENSION( 26 ) :: LOWER, UPPER

     DATA LOWER / 'a', 'b', 'c', 'd', 'e', 'f', 'g', 'h', 'i', 'j',            &
                  'k', 'l', 'm', 'n', 'o', 'p', 'q', 'r', 's', 't',            &
                  'u', 'v', 'w', 'x', 'y', 'z' /
     DATA UPPER / 'A', 'B', 'C', 'D', 'E', 'F', 'G', 'H', 'I', 'J',            &
                  'K', 'L', 'M', 'N', 'O', 'P', 'Q', 'R', 'S', 'T',            &
                  'U', 'V', 'W', 'X', 'Y', 'Z' /

!  See if the current letter is upper case. If so replace it by its
!  lower case counterpart

     DO letter = 1, 26
       IF ( string == UPPER( letter ) ) THEN
         string = LOWER( letter )
         EXIT
       END IF
     END DO

     RETURN

!  End of subroutine STRING_lower

     END SUBROUTINE STRING_lower

!-*-*-*-*-   S T R I N G _ l o w e r  _ w o r d   S U B R O U T I N E  -*-*-*-*-

     SUBROUTINE STRING_lower_word( word )

!  Convert a word of character strings from upper to lower case

!--------------------------------
!   D u m m y   A r g u m e n t
!--------------------------------

     CHARACTER (LEN = * ), INTENT( INOUT ) :: word

!  Local variables

     INTEGER :: i

!  Change the word, letter by letter, to lower case

     DO i = 1,  LEN_TRIM( word )
       CALL STRING_lower( word( i : i ) )
     END DO

     RETURN

!  End of subroutine STRING_lower_word

     END SUBROUTINE STRING_lower_word

!  End of module RPD

   END MODULE GALAHAD_RPD_double
