!  qplib2lp - convert a qplib file to ILOG lp format
!  The qplib file should be read from standard input (unit 5)
!  and the lp format will be written on standard output (unit 6)
!  NB. The qplib file is assumed to be from a continuous problem
!  Nick Gould, 10th Oct 2014, revised 2nd Feb 2017

  PROGRAM QPLIB2LP

  INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
  INTEGER, PARAMETER :: dp = KIND( 1.0D+0 )
  REAL( KIND = dp ), PARAMETER :: teneps_d = 10.0_dp * EPSILON( 1.0_dp )

  CHARACTER ( len = 80 ) :: default, line, name, index, line_80, out_80
  CHARACTER ( len = 80 ) :: blank_80 = REPEAT( ' ', 80 )
  CHARACTER ( len = 40 ) :: default_val = REPEAT( ' ', 40 )
  CHARACTER ( len = 6 ) :: type = REPEAT( ' ', 6 )
  CHARACTER ( len = 8 ) :: minimize = REPEAT( ' ', 6 )
  INTEGER :: i, i_current, l, n, m, nnz_h, nnz_ch, nnz_a, nnz_v
  INTEGER :: j, j_up, b_type, status
  REAL( KIND = wp ) :: infinity, value, lower, upper
  LOGICAL :: hessian = .TRUE., constraints = .TRUE.
  LOGICAL :: c_hessian = .FALSE.
  LOGICAL :: name_set, no_bounds_needed, h_set, linear, quadratic

  CHARACTER ( len = 40 ) :: val, f
  INTEGER, ALLOCATABLE, DIMENSION( : ) :: A_start
  INTEGER, ALLOCATABLE, DIMENSION( : ) :: CH_start
  CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: G
  CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: X_l
  CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: X_u
  CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: C_l
  CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: C_u
  CHARACTER ( len = 16 ), ALLOCATABLE, DIMENSION( : ) :: A_row
  CHARACTER ( len = 16 ), ALLOCATABLE, DIMENSION( : ) :: A_col
  CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: A_val
  CHARACTER ( len = 16 ), ALLOCATABLE, DIMENSION( : ) :: H_row
  CHARACTER ( len = 16 ), ALLOCATABLE, DIMENSION( : ) :: H_col
  CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: H_val
  CHARACTER ( len = 16 ), ALLOCATABLE, DIMENSION( : ) :: CH_ind
  CHARACTER ( len = 16 ), ALLOCATABLE, DIMENSION( : ) :: CH_row
  CHARACTER ( len = 16 ), ALLOCATABLE, DIMENSION( : ) :: CH_col
  CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: CH_val
! CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: X_name
! CHARACTER ( len = 40 ), ALLOCATABLE, DIMENSION( : ) :: C_name

!  =================================
!  read the data from the QPLIB file
!  =================================

!  read the name

  CALL FIND_DATA_LINE( name )
  CALL FIND_DATA_LINE( line )

!  find and assess the type

  READ( line, "( A6 )" ) type
  SELECT CASE ( TRIM( type ) )
  CASE( 'QCB', 'CCB', 'DCB' )
    constraints = .FALSE.
  CASE( 'LCL' )
    hessian = .FALSE.
  CASE( 'LCQ', 'LCC', 'LCD' )
    hessian = .FALSE.
    c_hessian = .TRUE.
  CASE( 'QCQ', 'QCC', 'QCD', 'CCQ', 'CCC', 'CCD', 'DCQ', 'DCC', 'DCD' )
    c_hessian = .TRUE.
  END SELECT

!  minimize or maximize?

  CALL FIND_DATA_LINE( line )
  READ( line, "( A8 )" ) minimize

!  find the number of variables

  CALL FIND_DATA_LINE( line )
  READ( line, * ) n
  ALLOCATE( G( n ), X_l( n ), X_u( n ), STAT = status )

!  find the number of constraints, if any

  IF ( constraints ) THEN
    CALL FIND_DATA_LINE( line )
    READ( line, * ) m
    ALLOCATE( C_l( m ), C_u( m ), A_start( m + 1 ), STAT = status )
    IF ( c_hessian ) ALLOCATE( CH_start( m + 1 ), STAT = status )
  END IF

!  read the Hessian, if any

  IF ( hessian ) THEN
    CALL FIND_DATA_LINE( line )
    READ( line, * ) nnz_h
    ALLOCATE( H_row( nnz_h ), H_col( nnz_h ), H_val( nnz_h ), STAT = status )
    DO l = 1, nnz_h
      CALL FIND_DATA_LINE( line )
      READ( line, * ) H_row( l ), H_col( l ), H_val( l )
    END DO
  END IF

!  read g

    CALL FIND_DATA_LINE( default )
    READ( default, * ) val
    G = val
    CALL FIND_DATA_LINE( line )
    READ( line, * ) nnz_v
    DO l = 1, nnz_v
      CALL FIND_DATA_LINE( line )
      val = default_val
      READ( line, * ) i, val
      G( i ) = val
    END DO

!  read f

    CALL FIND_DATA_LINE( line )
    READ( line, * ) f

!  read the constraint Hessians, if any

  IF ( c_hessian ) THEN
    i_current = 0
    CALL FIND_DATA_LINE( line )
    READ( line, * ) nnz_ch
    ALLOCATE( CH_ind( nnz_ch ), CH_row( nnz_ch ), CH_col( nnz_ch ),            &
              CH_val( nnz_ch ), STAT = status )
    DO l = 1, nnz_ch
      CALL FIND_DATA_LINE( line )
      READ( line, * ) CH_ind( l ), CH_row( l ), CH_col( l ), CH_val( l )

!  double the stored value as the idiotic format does not have a 1/2 for
!  constraint Hessians ... duh!

      READ( CH_val( l ), * )  value
      value = 2.0_wp * value
      CH_val( l ) = TRIM_VALUE( value ) // REPEAT( ' ', 16 )

!  store the starting address for each Hessian

      READ( CH_ind( l ), * ) i
      IF ( i /= i_current ) THEN
        CH_start( i_current + 1 : i ) = l
        i_current = i
      END IF
    END DO
    CH_start( i_current + 1 : m + 1 ) = nnz_ch + 1
  END IF

!  read the constraint Jacobian, if any

  IF ( constraints ) THEN
    i_current = 0
    CALL FIND_DATA_LINE( line )
    READ( line, * ) nnz_a
    ALLOCATE( A_row( nnz_a ), A_col( nnz_a ), A_val( nnz_a ), STAT = status )
    DO l = 1, nnz_a
      CALL FIND_DATA_LINE( line )
      READ( line, * ) A_row( l ), A_col( l ), A_val( l )

!  store the starting address for each row

      READ( A_row( l ), * ) i
      IF ( i /= i_current ) THEN
        A_start( i_current + 1 : i ) = l
        i_current = i
      END IF
    END DO
    A_start( i_current + 1 : m + 1 ) = nnz_a + 1
!write(6,*) ' A_start ', A_start
  END IF
!stop
!  read the value of infinity

  CALL FIND_DATA_LINE( line )
  READ( line, * ) infinity

!  read the bounds on the constraints

  IF ( constraints ) THEN

!  lower bounds

    CALL FIND_DATA_LINE( default )
    READ( default, * ) val
    C_l = val
    CALL FIND_DATA_LINE( line )
    READ( line, * ) nnz_v
    DO l = 1, nnz_v
      CALL FIND_DATA_LINE( line )
      val = default_val
      READ( line, * ) i, val
      C_l( i ) = val
    END DO

!  upper bounds

    CALL FIND_DATA_LINE( default )
    READ( default, * ) val
    C_u = val
    CALL FIND_DATA_LINE( line )
    READ( line, * ) nnz_v
    DO l = 1, nnz_v
      CALL FIND_DATA_LINE( line )
      val = default_val
      READ( line, * ) i, val
      C_u( i ) = val
    END DO
  END IF

!  read the bounds on the variables

!  lower bounds

  CALL FIND_DATA_LINE( default )
  READ( default, * ) val
  READ( val, * ) lower
  X_l = val
  CALL FIND_DATA_LINE( line )
  READ( line, * ) nnz_v
  DO l = 1, nnz_v
    CALL FIND_DATA_LINE( line )
    val = default_val
    READ( line, * ) i, val
    X_l( i ) = val
  END DO

  no_bounds_needed = lower == 0.0_wp .AND. nnz_v == 0

!  upper bounds

  CALL FIND_DATA_LINE( default )
  READ( default, * ) val
  READ( val, * ) upper
  X_u = val
  CALL FIND_DATA_LINE( line )
  READ( line, * ) nnz_v
  DO l = 1, nnz_v
    CALL FIND_DATA_LINE( line )
    val = default_val
    READ( line, * ) i, val
    X_u( i ) = val
  END DO

  no_bounds_needed = no_bounds_needed .AND.                                    &
   ( upper > infinity .AND. nnz_v == 0 )

!  read the starting values for x (not used)

  CALL FIND_DATA_LINE( default )
  CALL FIND_DATA_LINE( line )
  READ( line, * ) nnz_v
  DO l = 1, nnz_v
    CALL FIND_DATA_LINE( line )
  END DO

!  read the starting values for y (not used)

  IF ( constraints ) THEN
    CALL FIND_DATA_LINE( default )
    CALL FIND_DATA_LINE( line )
    READ( line, * ) nnz_v
    DO l = 1, nnz_v
      CALL FIND_DATA_LINE( line )
    END DO
  END IF

!  read the starting values for z (not used)

  CALL FIND_DATA_LINE( default )
  CALL FIND_DATA_LINE( line )
  READ( line, * ) nnz_v
  DO l = 1, nnz_v
    CALL FIND_DATA_LINE( line )
  END DO

!  ================================
!  write the data in LP file format
!  ================================

  WRITE( 6, "( '\ ', A )" ) TRIM( name )
  WRITE( 6, "( '\ problem type: ', A )" ) TRIM( type )
! WRITE( 6, "( '' )" )
  WRITE( 6, "( A )" ) minimize

!  write the linear term of the objective function

  out_80 = blank_80
  name_set = .FALSE.
  DO i = 1, n
    READ( G( i ), * ) value
    WRITE( index, * ) i
    index = ADJUSTL( index )
    IF ( value /= 0.0_wp ) THEN
      IF ( name_set ) THEN
        IF ( value > 0.0_wp ) THEN
          IF ( value == 1.0_wp ) THEN
            line_80 = blank_80
            WRITE( line_80, "( '+ x', A )" ) TRIM( index )
            CALL LINE_manipulate( out_80, line_80 )
          ELSE
            line_80 = blank_80
            WRITE( line_80, "( '+ ', A, ' ', 'x', A )" )                       &
              TRIM( G( i ) ), TRIM( index )
            CALL LINE_manipulate( out_80, line_80 )
          END IF
        ELSE
          IF ( value == - 1.0_wp ) THEN
            line_80 = blank_80
            WRITE( line_80, "( '- x', A )" ) TRIM( index )
            CALL LINE_manipulate( out_80, line_80 )
          ELSE
            line_80 = blank_80
            WRITE( line_80, "( '- ', A, ' ', 'x', A )" )                       &
              TRIM( G( i )( 2 : ) ), TRIM( index )
            CALL LINE_manipulate( out_80, line_80 )
          END IF
        END IF
      ELSE
        IF ( value == 1.0_wp ) THEN
          line_80 = blank_80
          WRITE( line_80, "( 'obj: x', A )" ) TRIM( index )
          CALL LINE_manipulate( out_80, line_80 )
        ELSE IF ( value == - 1.0_wp ) THEN
          line_80 = blank_80
          WRITE( line_80, "( 'obj: - x', A )" ) TRIM( index )
          CALL LINE_manipulate( out_80, line_80 )
        ELSE
          IF ( value > 0.0_wp ) THEN
            line_80 = blank_80
            WRITE( line_80, "( 'obj: ', A, ' ', 'x', A )" )                   &
              TRIM( G( i ) ), TRIM( index )
            CALL LINE_manipulate( out_80, line_80 )
          ELSE
            line_80 = blank_80
            WRITE( line_80, "( 'obj: - ', A, ' ', 'x', A )" )                 &
              TRIM( G( i )( 2 : ) ), TRIM( index )
            CALL LINE_manipulate( out_80, line_80 )
          END IF
        END IF
        name_set = .TRUE.
      END IF
    END IF
  END DO

!  write the quadratic term of the objective function

  IF ( hessian .AND. nnz_h > 0 ) THEN
    h_set = .FALSE.
    DO l = 1, nnz_h
      READ( H_val( l ), * ) value
      IF ( value /= 0.0_wp ) THEN
        IF ( h_set ) THEN
          IF ( value > 0.0_wp ) THEN
            IF ( value == 1.0_wp ) THEN
              IF ( TRIM( ADJUSTL( H_row( l ) ) ) ==                            &
                   TRIM( ADJUSTL( H_col( l ) ) ) ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '+ x', A, '^2' )" )                         &
                  TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE
                line_80 = blank_80
                WRITE( line_80, "( '+ x', A, ' * x', A )" )                    &
                  TRIM( ADJUSTL( H_row( l ) ) ), TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              END IF
            ELSE
              IF ( TRIM( ADJUSTL( H_row( l ) ) ) ==                            &
                   TRIM( ADJUSTL( H_col( l ) ) ) ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '+ ', A, ' ', 'x', A, '^2' )" )             &
                  TRIM( H_val( l ) ), TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE
                line_80 = blank_80
                WRITE( line_80, "( '+ ', A, ' ', 'x',                          &
               &                   A, ' * x', A )" )    &
                  TRIM( H_val( l ) ), TRIM( ADJUSTL( H_row( l ) ) ),           &
                  TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              END IF
            END IF
          ELSE
            IF ( value == - 1.0_wp ) THEN
              IF ( TRIM( ADJUSTL( H_row( l ) ) ) ==                            &
                   TRIM( ADJUSTL( H_col( l ) ) ) ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '- x', A, '^2' )" )                         &
                  TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE
                line_80 = blank_80
                WRITE( line_80, "( '- x', A, ' * x', A )" )                    &
                  TRIM( ADJUSTL( H_row( l ) ) ), TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              END IF
            ELSE
              IF ( TRIM( ADJUSTL( H_row( l ) ) ) ==                            &
                   TRIM( ADJUSTL( H_col( l ) ) ) ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '- ', A, ' ', 'x', A, '^2' )" )             &
                  TRIM( H_val( l )( 2 : ) ), TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE
                line_80 = blank_80
                WRITE( line_80, "( '- ', A, ' ', 'x',                          &
               &                   A, ' * x', A )" )    &
                  TRIM( H_val( l )( 2 : ) ), TRIM( ADJUSTL( H_row( l ) ) ),    &
                  TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              END IF
            END IF
          END IF
        ELSE
          IF ( name_set ) THEN
            IF ( TRIM( ADJUSTL( H_row( l ) ) ) ==                              &
                 TRIM( ADJUSTL( H_col( l ) ) ) ) THEN
              IF ( value == 1.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '+ [ x', A, '^2' )" )                       &
                  TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE IF ( value == - 1.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '+ [ - x', A, '^2' )" )                     &
                  TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE IF ( value > 0.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '+ [ ', A, ' ', 'x', A, '^2' )" )           &
                  TRIM( H_val( l ) ), TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE
                line_80 = blank_80
                WRITE( line_80, "( '+ [ - ', A, ' ', 'x', A, '^2' )" )         &
                  TRIM( H_val( l )( 2 : ) ), TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              END IF
            ELSE
              IF ( value == 1.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '+ [  x', A, ' * x', A )" )                 &
                  TRIM( ADJUSTL( H_row( l ) ) ), TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE IF ( value == - 1.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '+ [ - x', A, ' * x', A )" )                &
                  TRIM( ADJUSTL( H_row( l ) ) ), TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE IF ( value > 0.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( '+ [ ', A, ' ', 'x', A, ' * x', A )" )      &
                  TRIM( H_val( l ) ), TRIM( ADJUSTL( H_row( l ) ) ),           &
                  TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE
                line_80 = blank_80
                WRITE( line_80, "( '+ [ - ', A, ' ', 'x',                      &
               &                    A, ' * x', A )" )    &
                  TRIM( H_val( l )( 2 : ) ), TRIM( ADJUSTL( H_row( l ) ) ),    &
                  TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              END IF
            END IF
          ELSE
            IF ( TRIM( ADJUSTL( H_row( l ) ) ) ==                              &
                 TRIM( ADJUSTL( H_col( l ) ) ) ) THEN
              IF ( value == 1.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( 'obj: [  x', A, '^2' )" )                   &
                  TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE IF ( value == - 1.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( 'obj: [ - x', A, '^2' )" )                  &
                  TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE IF ( value > 0.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( 'obj: [ ', A, ' ', 'x', A, '^2' )" )        &
                  TRIM( H_val( l ) ), TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE
                line_80 = blank_80
                WRITE( line_80, "( 'obj: [ - ', A, ' ', 'x', A, '^2' )" )      &
                  TRIM( H_val( l )( 2 : ) ), TRIM( ADJUSTL( H_row( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              END IF
            ELSE
              IF ( value == 1.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( 'obj: [  x', A, ' * x', A )" )              &
                  TRIM( ADJUSTL( H_row( l ) ) ), TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE IF ( value == - 1.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( 'obj: [ - x', A, ' * x', A )" )             &
                  TRIM( ADJUSTL( H_row( l ) ) ), TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE IF ( value > 0.0_wp ) THEN
                line_80 = blank_80
                WRITE( line_80, "( 'obj: [ ', A, ' ', 'x', A, ' * x', A )" )   &
                  TRIM( H_val( l ) ), TRIM( ADJUSTL( H_row( l ) ) ),           &
                  TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              ELSE
                line_80 = blank_80
                WRITE( line_80, "( 'obj: [ - ', A, ' ', 'x', A, ' * x', A )" ) &
                  TRIM( H_val( l )( 2 : ) ), TRIM( ADJUSTL( H_row( l ) ) ),    &
                  TRIM( ADJUSTL( H_col( l ) ) )
                CALL LINE_manipulate( out_80, line_80 )
              END IF
            END IF
            name_set = .TRUE.
          END IF
          h_set = .TRUE.
        END IF
      END IF
    END DO
    IF ( h_set ) THEN
      line_80 = blank_80
      WRITE( line_80, "( '] / 2' )" )
      CALL LINE_manipulate( out_80, line_80 )
    END IF
  END IF
  CALL LINE_flush( out_80 )

!  write the linear constraints, if any

  WRITE( 6, "( 'Subject to' )" )
  IF ( constraints ) THEN
    DO i = 1, m
      linear = A_start( i + 1 ) > A_start( i )
      IF ( c_hessian ) THEN
        quadratic = CH_start( i + 1 ) > CH_start( i )
      ELSE
        quadratic = .FALSE.
      END IF
      READ( C_l( i ), * ) lower
      READ( C_u( i ), * ) upper
      j_up = 1
      IF ( lower == upper ) THEN
        b_type = 0
      ELSE IF ( lower > - infinity .AND. upper < infinity ) THEN
        b_type = 3
        j_up = 2
      ELSE IF ( lower > - infinity ) THEN
        b_type = 1
      ELSE
        b_type = 2
      END IF
      DO j = 1, j_up
        name_set = .FALSE. ; h_set = .FALSE.
        IF ( linear ) THEN
          DO l = A_start( i ), A_start( i + 1 ) - 1
            READ( A_val( l ), * ) value
            IF ( value /= 0.0_wp ) THEN
              IF ( name_set ) THEN
                IF ( value > 0.0_wp ) THEN
                  IF ( value == 1.0_wp ) THEN
                    line_80 = blank_80
                    WRITE( line_80, "( '+ x', A )" )                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  ELSE
                    line_80 = blank_80
                    WRITE( line_80, "( '+ ', A, ' x', A )" )                   &
                      TRIM( ADJUSTL( A_val( l ) ) ),                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  END IF
                ELSE
                  IF ( value == - 1.0_wp ) THEN
                    line_80 = blank_80
                    WRITE( line_80, "( '- x', A )" )                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  ELSE
                    line_80 = blank_80
                    WRITE( line_80, "( '- ', A, ' x', A )" )                   &
                      TRIM( ADJUSTL( A_val( l )( 2 : ) ) ),                    &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  END IF
                END IF
              ELSE
                IF ( j == 1 ) THEN
                  IF ( value == 1.0_wp ) THEN
                    line_80 = blank_80
                    WRITE( line_80, "( 'c', A, ': x', A )" )                   &
                      TRIM( ADJUSTL( A_row( l ) ) ),                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  ELSE IF ( value == - 1.0_wp ) THEN
                     line_80 = blank_80
                    WRITE( line_80, "( 'c', A, ': - x', A )" )                 &
                      TRIM( ADJUSTL( A_row( l ) ) ),                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  ELSE IF ( value > 0.0_wp ) THEN
                    line_80 = blank_80
                    WRITE( line_80, "( 'c', A, ': ', A, ' x', A )" )           &
                      TRIM( ADJUSTL( A_row( l ) ) ),                           &
                      TRIM( ADJUSTL( A_val( l ) ) ),                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  ELSE
                    line_80 = blank_80
                    WRITE( line_80, "( 'c', A, ': - ', A, ' x', A )" )         &
                      TRIM( ADJUSTL( A_row( l ) ) ),                           &
                      TRIM( ADJUSTL( A_val( l )( 2 : ) ) ),                    &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  END IF
                ELSE
                  IF ( value == 1.0_wp ) THEN
                    line_80 = blank_80
                    WRITE( line_80, "( 'c', A, '_range : x', A )" )            &
                      TRIM( ADJUSTL( A_row( l ) ) ),                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  ELSE IF ( value == - 1.0_wp ) THEN
                    line_80 = blank_80
                    WRITE( line_80, "( 'c', A, '_range : - x', A )" )          &
                      TRIM( ADJUSTL( A_row( l ) ) ),                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  ELSE IF ( value > 0.0_wp ) THEN
                    line_80 = blank_80
                    WRITE( line_80, "( 'c', A, '_range : ', A, ' x', A )" )    &
                      TRIM( ADJUSTL( A_row( l ) ) ),                           &
                      TRIM( ADJUSTL( A_val( l ) ) ),                           &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  ELSE
                    line_80 = blank_80
                    WRITE( line_80, "( 'c', A, '_range : - ', A, ' x', A )" )  &
                      TRIM( ADJUSTL( A_row( l ) ) ),                           &
                      TRIM( ADJUSTL( A_val( l )( 2 : ) ) ),                    &
                      TRIM( ADJUSTL( A_col( l ) ) )
                    CALL LINE_manipulate( out_80, line_80 )
                  END IF
                END IF
                name_set = .TRUE.
              END IF
            END IF
          END DO
        END IF

        IF ( quadratic ) THEN
          DO l = CH_start( i ), CH_start( i + 1 ) - 1
            READ( CH_val( l ), * ) value
            IF ( value /= 0.0_wp ) THEN
              IF ( h_set ) THEN
                IF ( value > 0.0_wp ) THEN
                  IF ( value == 1.0_wp ) THEN
                    IF ( TRIM( ADJUSTL( CH_row( l ) ) ) ==                     &
                         TRIM( ADJUSTL( CH_col( l ) ) ) ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( '+ x', A, '^2' )" )                   &
                        TRIM( ADJUSTL( CH_row( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE
                      line_80 = blank_80
                      WRITE( line_80, "( '+ x', A, ' * x', A )" )              &
                        TRIM( ADJUSTL( CH_row( l ) ) ),                        &
                        TRIM( ADJUSTL( CH_col( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    END IF
                  ELSE
                    IF ( TRIM( ADJUSTL( CH_row( l ) ) ) ==                     &
                         TRIM( ADJUSTL( CH_col( l ) ) ) ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( '+ ', A, ' ', 'x',                    &
                      &                  A, '^2' )" )                          &
                        TRIM( CH_val( l ) ), TRIM( ADJUSTL( CH_row( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE
                      line_80 = blank_80
                      WRITE( line_80, "( '+ ', A, ' ', 'x', A,                 &
                     &                   ' * x', A ) ")                        &
                        TRIM( CH_val( l ) ), TRIM( ADJUSTL( CH_row( l ) ) ),   &
                        TRIM( ADJUSTL( CH_col( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    END IF
                  END IF
                ELSE
                  IF ( value == - 1.0_wp ) THEN
                    IF ( TRIM( ADJUSTL( CH_row( l ) ) ) ==                     &
                         TRIM( ADJUSTL( CH_col( l ) ) ) ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( '- x', A, '^2' )" )                   &
                        TRIM( ADJUSTL( CH_row( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE
                      line_80 = blank_80
                      WRITE( line_80, "( '- x', A, ' * x', A )" )              &
                        TRIM( ADJUSTL( CH_row( l ) ) ),                        &
                        TRIM( ADJUSTL( CH_col( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    END IF
                  ELSE
                    IF ( TRIM( ADJUSTL( CH_row( l ) ) ) ==                     &
                         TRIM( ADJUSTL( CH_col( l ) ) ) ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( '- ',                                 &
                     &                   A, ' ', 'x', A, '^2' )" )             &
                        TRIM( CH_val( l )( 2 : ) ),                            &
                        TRIM( ADJUSTL( CH_row( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE
                      line_80 = blank_80
                      WRITE( line_80, "( '- ', A, ' ', 'x',                    &
                     &              A, ' * x', A )") &
                        TRIM( CH_val( l )( 2 : ) ),                            &
                        TRIM( ADJUSTL( CH_row( l ) ) ),                        &
                        TRIM( ADJUSTL( CH_col( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    END IF
                  END IF
                END IF
              ELSE
                IF ( name_set ) THEN
                  IF ( TRIM( ADJUSTL( CH_row( l ) ) ) ==                       &
                       TRIM( ADJUSTL( CH_col( l ) ) ) ) THEN
                    IF ( value == 1.0_wp ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( '+ [ x', A, '^2' )" )                 &
                        TRIM( ADJUSTL( CH_row( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE IF ( value == - 1.0_wp ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( ' + [ - x', A, '^2' )" )              &
                        TRIM( ADJUSTL( CH_row( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE IF ( value > 0.0_wp ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( '+ [ ', A, ' ', 'x',                  &
                     &                   A, '^2' )" )                          &
                        TRIM( CH_val( l ) ), TRIM( ADJUSTL( CH_row( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE
                      line_80 = blank_80
                      WRITE( line_80, "( '+ [ - ', A, ' ', 'x',                &
                     &                   A, '^2' )" )   &
                        TRIM( CH_val( l )( 2 : ) ),                            &
                        TRIM( ADJUSTL( CH_row( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    END IF
                  ELSE
                    IF ( value == 1.0_wp ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( '+ [  x', A, ' * x', A )" )           &
                        TRIM( ADJUSTL( CH_row( l ) ) ),                        &
                        TRIM( ADJUSTL( CH_col( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE IF ( value == - 1.0_wp ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( '+ [ - x', A, ' * x', A )" )          &
                        TRIM( ADJUSTL( CH_row( l ) ) ),                        &
                        TRIM( ADJUSTL( CH_col( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE IF ( value > 0.0_wp ) THEN
                      line_80 = blank_80
                      WRITE( line_80, "( ' + [ ', A, ' ', 'x',                 &
                     &                     A, ' * x', A )")                    &
                        TRIM( CH_val( l ) ), TRIM( ADJUSTL( CH_row( l ) ) ),   &
                        TRIM( ADJUSTL( CH_col( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    ELSE
                      line_80 = blank_80
                      WRITE( line_80, "( '+ [ - ', A, ' ', 'x',                &
                     &             A, ' * x', A )")                            &
                        TRIM( CH_val( l )( 2 : ) ),                            &
                        TRIM( ADJUSTL( CH_row( l ) ) ),                        &
                        TRIM( ADJUSTL( CH_col( l ) ) )
                      CALL LINE_manipulate( out_80, line_80 )
                    END IF
                  END IF
                ELSE
                  IF ( j == 1 ) THEN
                    IF ( TRIM( ADJUSTL( CH_row( l ) ) ) ==                     &
                         TRIM( ADJUSTL( CH_col( l ) ) ) ) THEN
                      IF ( value == 1.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, ': [  x', A, '^2' )" )      &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_row( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE IF ( value == - 1.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, ': [ - x', A, '^2' )" )     &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_row( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE IF ( value > 0.0_wp ) THEN
                       line_80 = blank_80
                        WRITE( line_80, "( 'c', A, ': [ ', A, ' ', 'x',        &
                       &                    A, '^2' )") &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( CH_val( l ) ), TRIM( ADJUSTL( CH_row( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, ': [ - ', A, ' ', 'x',      &
                       &              A, '^2' )")                              &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( CH_val( l )( 2 : ) ),                          &
                          TRIM( ADJUSTL( CH_row( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      END IF
                    ELSE
                      IF ( value == 1.0_wp ) THEN
                         line_80 = blank_80
                        WRITE( line_80, "( 'c', A, ': [  x', A, ' * x', A )")  &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_row( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_col( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE IF ( value == - 1.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, ': [ - x', A, ' * x', A )") &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_row( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_col( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE IF ( value > 0.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, ': [ ', A, ' ', 'x',        &
                       &              A, ' * x', A ) " )                       &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( CH_val( l ) ),                                 &
                          TRIM( ADJUSTL( CH_row( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_col( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, ': [ - ', A, ' ', 'x',      &
                       &              A, ' * x', A ) " )                       &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( CH_val( l )( 2 : ) ),                          &
                          TRIM( ADJUSTL( CH_row( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_col( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      END IF
                    END IF
                  ELSE
                    IF ( TRIM( ADJUSTL( CH_row( l ) ) ) ==                     &
                         TRIM( ADJUSTL( CH_col( l ) ) ) ) THEN
                      IF ( value == 1.0_wp ) THEN
                       line_80 = blank_80
                        WRITE( line_80, "( 'c', A, '_range : [  x',            &
                       &                     A, '^2' )" )                      &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_row( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE IF ( value == - 1.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, '_range : [ - x',           &
                       &                    A, '^2' )" )   &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_row( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE IF ( value > 0.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, '_range : [ ', A, ' ', 'x', &
                       &             A, '^2' ) " )                             &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( CH_val( l ) ), TRIM( ADJUSTL( CH_row( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, '_range : [ - ', A, ' ',    &
                       &             'x', A, '^2' ) " )                        &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( CH_val( l )( 2 : ) ),                          &
                          TRIM( ADJUSTL( CH_row( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      END IF
                    ELSE
                      IF ( value == 1.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, '_range : [  x',            &
                       &                    A, ' * x', A )" ) &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_row( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_col( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE IF ( value == - 1.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, '_range : [ - x', A,        &
                       &                   ' * x', A )" )                      &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_row( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_col( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE IF ( value > 0.0_wp ) THEN
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, '_range : [ ', A, ' ',      &
                       &               'x', A, ' * x', A ) " )                 &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( CH_val( l ) ),                                 &
                          TRIM( ADJUSTL( CH_row( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_col( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      ELSE
                        line_80 = blank_80
                        WRITE( line_80, "( 'c', A, '_range : [ - ', A, ' ',    &
                       &                    'x', A, ' * x', A ) " )            &
                          TRIM( ADJUSTL( CH_ind( l ) ) ),                      &
                          TRIM( CH_val( l )( 2 : ) ),                          &
                          TRIM( ADJUSTL( CH_row( l ) ) ),                      &
                          TRIM( ADJUSTL( CH_col( l ) ) )
                        CALL LINE_manipulate( out_80, line_80 )
                      END IF
                    END IF
                  END IF
                  name_set = .TRUE.
                END IF
                h_set = .TRUE.
              END IF
            END IF
          END DO
          line_80 = blank_80
          IF ( h_set ) WRITE( line_80, "( ']' )" )
          CALL LINE_manipulate( out_80, line_80 )
        END IF
        IF ( b_type == 0 ) THEN
          line_80 = blank_80
          WRITE( line_80, "( '= ', A )" ) C_l( i )
          CALL LINE_manipulate( out_80, line_80 )
        ELSE IF ( b_type == 1 ) THEN
          line_80 = blank_80
          WRITE( line_80, "( '>= ', A )" ) C_l( i )
          CALL LINE_manipulate( out_80, line_80 )
        ELSE IF ( b_type == 2 ) THEN
          line_80 = blank_80
          WRITE( line_80, "( '<= ', A )" ) C_u( i )
          CALL LINE_manipulate( out_80, line_80 )
        ELSE
          IF ( j == 1 ) THEN
            line_80 = blank_80
            WRITE( line_80, "( '>= ', A )" ) C_l( i )
            CALL LINE_manipulate( out_80, line_80 )
          ELSE
            line_80 = blank_80
            WRITE( line_80, "( '<= ', A )" ) C_u( i )
            CALL LINE_manipulate( out_80, line_80 )
          END IF
        END IF
        CALL LINE_flush( out_80 )
      END DO
    END DO
  END IF

!  write the bounds

  IF ( no_bounds_needed ) THEN
    WRITE( 6, "( '\ default bounds [0,infinity) to be understood' )" )
  ELSE
    WRITE( 6, "( 'Bounds' )" )
    DO i = 1, n
      WRITE( index, * ) i
      index = ADJUSTL( index )
      READ( X_l( i ), * ) lower
      READ( X_u( i ), * ) upper
      IF ( lower == 0.0_wp .AND. upper >= infinity ) THEN
!  skip ... default limit
      ELSE IF ( lower <= - infinity .AND. upper >= infinity ) THEN
        WRITE( 6, "( ' x', A, ' free' )" ) TRIM( index )
      ELSE IF ( lower == upper ) THEN
        WRITE( 6, "( ' x', A, ' = ', A )" ) TRIM( index ),                     &
          TRIM( ADJUSTL( X_l( i ) ) )
      ELSE IF ( lower > - infinity .AND. upper < infinity ) THEN
        WRITE( 6, "( ' ', A, ' <= x', A, ' <= ', A )" )                        &
          TRIM( ADJUSTL( X_l( i ) ) ), TRIM( index ), TRIM( ADJUSTL( X_u( i ) ))
      ELSE IF ( lower > - infinity ) THEN
        WRITE( 6, "( ' x', A, ' >= ', A )" ) TRIM( index ),                    &
          TRIM( ADJUSTL( X_l( i ) ) )
      ELSE
        WRITE( 6, "( ' - infinity <= x', A, ' <= ', A )" ) TRIM( index ),      &
          TRIM( ADJUSTL( X_u( i ) ) )
      END IF
    END DO
  END IF
  WRITE( 6, "( 'End' )" )

  CONTAINS
    SUBROUTINE FIND_DATA_LINE( line )

!  find the next non blank, non comment line

    CHARACTER ( len = 80 ) :: line
    CHARACTER ( len = 80 ) :: blank = REPEAT( ' ', 80 )

    DO
      line = blank
      READ( 5, "( A80 )" ) line
      line = ADJUSTL( line )
      IF ( line /= blank .AND. line( 1 : 1 ) /= "!" .AND.                      &
           line( 1 : 1 ) /= "%" .AND. line( 1 : 1 ) /= "#" ) THEN
        EXIT
      END IF
    END DO
    END SUBROUTINE FIND_DATA_LINE

    SUBROUTINE LINE_manipulate( existing, new, flush )

!  merge the strings existing and new, flusing existing to unit 6 if the
!  combined string is more than 80 characters and replace it with new.
!  Optionally flush existing if requested

    CHARACTER ( len = 80 ) :: existing, new
    LOGICAL, OPTIONAL :: flush
    INTEGER :: len_existing, len_new
    CHARACTER ( len = 80 ), PARAMETER :: blank = REPEAT( ' ', 80 )

    len_new = LEN( TRIM( new ) )
    len_existing = LEN( TRIM( existing ) )
    IF ( len_new + MAX( len_existing, 1 ) <= 79 ) THEN
      existing = TRIM( existing ) // ' ' // TRIM( new )
      IF ( PRESENT( flush ) ) CALL LINE_flush( existing )
    ELSE
      WRITE( 6, "( A )" ) TRIM( existing )
      existing = blank
      existing = '  ' // TRIM( new )
      IF ( PRESENT( flush ) ) CALL LINE_flush( existing )
    END IF

    END SUBROUTINE LINE_manipulate

    SUBROUTINE LINE_flush( existing )
    CHARACTER ( len = 80 ) :: existing
    IF ( LEN( TRIM( existing ) ) > 0 ) THEN
      WRITE( 6, "( A )" ) TRIM( existing )
      existing = REPEAT( ' ', 80 )
    END IF
    END SUBROUTINE LINE_flush

!  --------------------- T R I M _ V A L U E   F U N C T I O N ----------------

      FUNCTION TRIM_VALUE( value )
      CHARACTER ( LEN = 24 ) :: TRIM_VALUE
      REAL ( KIND = wp ) :: value

!  write a real value into 24 characters trimming as much as possible
!  without losing precision

      INTEGER :: i, i_start, i_point, i_end, j, k, l, zs
      REAL ( KIND = wp ) :: minus_value
      LOGICAL :: zeros
      CHARACTER ( LEN = 22 ) :: field22
      CHARACTER ( LEN = 23 ) :: field
      CHARACTER ( LEN = 24 ) :: field24

!  cram value into 23 characters

!write(6,*) value
      IF ( value == 0.0_dp ) THEN
        field = "0.0         "
      ELSE IF ( SIGN( 1.0_dp, value ) > 0.0_dp ) THEN
        IF ( value >= ( 10.0_dp ) ** 100 ) THEN
          WRITE( field24, "( ES24.15E3 )" ) value
          field = field24( 1 : 20 ) // field24( 22 : 24 )
        ELSE IF ( value >= ( 10.0_dp ) ** 16 ) THEN
          WRITE( field24, "( ES24.15 )" ) value
          field = field24( 1 : 21 ) // field24( 23 : 24 )
        ELSE IF ( value >= ( 10.0_dp ) ** 15 ) THEN
          WRITE( field, "( F23.0 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 14 ) THEN
          WRITE( field, "( F23.1 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 13 ) THEN
          WRITE( field, "( F23.2 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 12 ) THEN
          WRITE( field, "( F23.3 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 11 ) THEN
          WRITE( field, "( F23.4 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 10 ) THEN
          WRITE( field, "( F23.5 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 9 ) THEN
          WRITE( field, "( F23.6 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 8 ) THEN
          WRITE( field, "( F23.7 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 7 ) THEN
          WRITE( field, "( F23.8 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 6 ) THEN
          WRITE( field, "( F23.9 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 5 ) THEN
          WRITE( field, "( F23.10 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 4 ) THEN
          WRITE( field, "( F23.11 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 3 ) THEN
          WRITE( field, "( F23.12 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 2 ) THEN
          WRITE( field, "( F23.13 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 1 ) THEN
          WRITE( field, "( F23.14 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** 0 ) THEN
          WRITE( field, "( F23.15 )" ) value
        ELSE IF ( value >= ( 10.0_dp ) ** ( - 1 ) ) THEN
          WRITE( field24, "( F24.16 )" ) value
          field = field24( 2 : 24 )
        ELSE IF ( value >= ( 10.0_dp ) ** ( - 2 ) ) THEN
          WRITE( field24, "( F24.17 )" ) value
          field = field24( 2 : 24 )
        ELSE IF ( value >= ( 10.0_dp ) ** ( - 3 ) ) THEN
          WRITE( field24, "( F24.18 )" ) value
          field = field24( 2 : 24 )
        ELSE IF ( value >= ( 10.0_dp ) ** ( - 4 ) ) THEN
          WRITE( field24, "( F24.16 )" ) value
          field = field24( 2 : 24 )
        ELSE IF ( value >= ( 10.0_dp ) ** ( - 9 ) ) THEN
          WRITE( field24, "( ES24.15 )" ) value
          field = field24( 1 : 22 ) // field24( 24 : 24 )
        ELSE IF ( value >= ( 10.0_dp ) ** ( - 99 ) ) THEN
          WRITE( field, "( ES23.15 )" ) value
!       ELSE IF ( value >= ( 10.0_dp ) ** ( - 999 ) ) THEN
!         WRITE( field, "( ES23.15E3 )" ) value
        ELSE
          WRITE( field, "( ES23.15E4 )" ) value
        END IF
      ELSE
        minus_value = - value
        IF ( ABS( minus_value - 1.0_dp ) <= teneps_d ) minus_value = 1.0_dp
        IF ( minus_value >= ( 10.0_dp ) ** 100 ) THEN
          WRITE( field, "( ES23.15E3 )" ) minus_value
          field22 = field( 1 : 19 ) // field( 21 : 23 )
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 16 ) THEN
          WRITE( field, "( ES23.15 )" ) minus_value
          field22 = field( 1 : 20 ) // field( 22 : 23 )
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 15 ) THEN
          WRITE( field22, "( F22.0 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 14 ) THEN
          WRITE( field22, "( F22.1 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 13 ) THEN
          WRITE( field22, "( F22.2 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 12 ) THEN
          WRITE( field22, "( F22.3 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 11 ) THEN
          WRITE( field22, "( F22.4 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 10 ) THEN
          WRITE( field22, "( F22.5 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 9 ) THEN
          WRITE( field22, "( F22.6 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 8 ) THEN
          WRITE( field22, "( F22.7 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 7 ) THEN
          WRITE( field22, "( F22.8 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 6 ) THEN
          WRITE( field22, "( F22.9 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 5 ) THEN
          WRITE( field22, "( F22.10 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 4 ) THEN
          WRITE( field22, "( F22.11 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 3 ) THEN
          WRITE( field22, "( F22.12 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 2 ) THEN
          WRITE( field22, "( F22.13 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 1 ) THEN
          WRITE( field22, "( F22.14 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** 0 ) THEN
          WRITE( field22, "( F22.15 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_dp ) ** ( - 1 ) ) THEN
          WRITE( field, "( F23.16 )" ) minus_value
          field22 = field( 2 : 23 )
        ELSE IF ( minus_value >= ( 10.0_dp ) ** ( - 2 ) ) THEN
          WRITE( field, "( F23.17 )" ) minus_value
          field22 = field( 2 : 23 )
        ELSE IF ( minus_value >= ( 10.0_dp ) ** ( - 3 ) ) THEN
          WRITE( field, "( F23.18 )" ) minus_value
          field22 = field( 2 : 23 )
        ELSE IF ( minus_value >= ( 10.0_dp ) ** ( - 4 ) ) THEN
          WRITE( field, "( F23.15 )" ) minus_value
          field22 = field( 2 : 23 )
        ELSE IF ( minus_value >= ( 10.0_dp ) ** ( - 9 ) ) THEN
          WRITE( field, "( ES23.15 )" ) minus_value
          field22 = field( 1 : 21 ) // field( 23 : 23 )
        ELSE IF ( minus_value > ( 10.0_dp ) ** ( - 99 ) ) THEN
          WRITE( field22, "( ES22.15 )" ) minus_value
!       ELSE IF ( minus_value > ( 10.0_dp ) ** ( - 999 ) ) THEN
!         WRITE( field22, "( ES22.15E3 )" ) minus_value
        ELSE
          WRITE( field22, "( ES22.15E4 )" ) minus_value
        END IF
        field = "-" //  ADJUSTL( field22 )
      END IF

!  shift the value left

      field24 = ADJUSTL( field ) // ' '

!  find the positions of the first digit in the mantissa

      IF ( field24( 1 : 1 ) == '-' ) THEN
        i_start = 2
      ELSE
        i_start = 1
      END IF

!  find the positions of the decimal point and last digit in the mantissa

      DO i = 1, 23
        IF ( field24( i : i ) == '.' ) i_point = i
        IF ( field24( i + 1 : i + 1 ) == ' ' .OR.                             &
             field24( i + 1 : i + 1 ) == 'e' .OR.                             &
             field24( i + 1 : i + 1 ) == 'E' .OR.                             &
             field24( i + 1 : i + 1 ) == 'd' .OR.                             &
             field24( i + 1 : i + 1 ) == 'D' ) THEN
          i_end = i
          EXIT
        END IF
      END DO

!     IF ( i_end - i_point >= 15 ) THEN
      IF ( i_end - i_start >= 15 ) THEN

!  round down any *01 to *00

        IF ( field24( i_end - 1 : i_end ) == '01' )                       &
          field24( i_end - 1 : i_end ) = '00'

!  round any *9r to **0r where ** = *+1

        IF ( field24( i_end - 1 : i_end ) == '99' ) THEN
          DO i = i_end, i_point + 1, - 1
            IF ( field24( i : i ) == '9' ) THEN
              field24( i : i ) = '0'
            ELSE
              READ( field24( i : i ), "( I1 )" ) l
              WRITE( field24( i : i ), "( I1 )" ) l + 1
              EXIT
            END IF
            IF ( i == i_point + 1 ) THEN
              DO j = i_point - 1, i_start, - 1
                IF ( field24( j : j ) == '9' ) THEN
                  field24( j : j ) = '0'
                ELSE
                  READ( field24( j : j ), "( I1 )" ) l
                  WRITE( field24( j : j ), "( I1 )" ) l + 1
                  EXIT
                END IF
                IF ( j == i_start ) THEN
                  DO l = i_end - 1, i_start, - 1
                    field24( l + 1 : l + 1 ) = field24( l : l )
                  END DO
                  field24( i_start : i_start ) = '1'
                END IF
              END DO
            END IF
          END DO
        END IF
      END IF

!     field24 = REPEAT( ' ', 24 )
!     IF ( value > - 10.0_wp .AND. value < 10.0_wp ) THEN
!       WRITE( field24, "( F19.16 )" ) value
!     ELSE
!       WRITE( field24, "( ES23.16 )" ) value
!     END IF

      TRIM_VALUE = field24

!  remove any leading space

!     IF ( TRIM_VALUE( 1 : 1 ) == ' ' ) THEN
!       DO i = 2, 24
!         TRIM_VALUE( i - 1 : i - 1 ) = TRIM_VALUE( i : i )
!       END DO
!     END IF

      zeros = .FALSE.
      DO i = 1, 24
        IF ( TRIM_VALUE( i : i ) == '0' ) THEN
          IF ( .NOT. zeros ) THEN
            zs = i
            zeros = .TRUE.
          END IF
        ELSE IF ( TRIM_VALUE( i : i ) == 'E' .OR.                              &
                  TRIM_VALUE( i : i ) == 'e' .OR.                              &
                  TRIM_VALUE( i : i ) == 'D' .OR.                              &
                  TRIM_VALUE( i : i ) == 'd' ) THEN
          IF ( zeros ) THEN
            DO j = zs + 1, zs + 25 - i
              k = i + ( j - zs - 1 )
              TRIM_VALUE( j : j ) = TRIM_VALUE( k : k  )
            END DO
            DO j = zs + 26 - i, 24
              TRIM_VALUE( j : j ) = ' '
            END DO
          END IF
          zeros = .FALSE.
          EXIT
        ELSE IF ( TRIM_VALUE( i : i ) == ' ' ) THEN
          IF ( zeros ) THEN
            DO j = zs + 1, i
              TRIM_VALUE( j : j ) = ' '
            END DO
          END IF
          zeros = .FALSE.
          EXIT
        ELSE
          zeros = .FALSE.
        END IF
      END DO
      IF ( zeros ) THEN
        DO j = zs + 1, i
          TRIM_VALUE( j : j ) = ' '
        END DO
      END IF

!  remove superflous 0 from the exponent

      DO i = 1, 24
        IF ( TRIM_VALUE( i : i ) == 'E' .OR.                                   &
             TRIM_VALUE( i : i ) == 'e' .OR.                                   &
             TRIM_VALUE( i : i ) == 'D' .OR.                                   &
             TRIM_VALUE( i : i ) == 'd' ) THEN
          IF ( TRIM_VALUE( i + 1 : i + 1 ) == '+' .OR.                         &
               TRIM_VALUE( i + 1 : i + 1 ) == '-' ) THEN
            IF ( TRIM_VALUE( i + 2 : i + 2 ) == '0' ) THEN
              IF ( TRIM_VALUE( i + 3 : i + 3 ) == '0' ) THEN
                IF ( TRIM_VALUE( i + 4 : i + 4 ) == ' ' ) THEN
                  TRIM_VALUE( i + 3 : i + 3 ) = '0'
                ELSE
                  TRIM_VALUE( i + 2 : i + 2 ) = TRIM_VALUE( i + 4 : i + 4 )
                  TRIM_VALUE( i + 3 : i + 4 ) = '  '
                END IF
              ELSE
                IF ( TRIM_VALUE( i + 4 : i + 4 ) == ' ' ) THEN
                ELSE
                  TRIM_VALUE( i + 2 : i + 2 ) = TRIM_VALUE( i + 3 : i + 3 )
                  TRIM_VALUE( i + 3 : i + 3 ) = TRIM_VALUE( i + 4 : i + 4 )
                  TRIM_VALUE( i + 4 : i + 4 ) = ' '
                END IF
              END IF
            END IF
          ELSE
            IF ( TRIM_VALUE( i + 1 : i + 1 ) == '0' ) THEN
              IF ( TRIM_VALUE( i + 2 : i + 2 ) == '0' ) THEN
                IF ( TRIM_VALUE( i + 3 : i + 3 ) == ' ' ) THEN
                  TRIM_VALUE( i + 2 : i + 2 ) = '0'
                ELSE
                  TRIM_VALUE( i + 1 : i + 1 ) = TRIM_VALUE( i + 3 : i + 3 )
                  TRIM_VALUE( i + 2 : i + 3 ) = '  '
                END IF
              ELSE
                IF ( TRIM_VALUE( i + 3 : i + 3 ) == ' ' ) THEN
                ELSE
                  TRIM_VALUE( i + 1 : i + 1 ) = TRIM_VALUE( i + 2 : i + 2 )
                  TRIM_VALUE( i + 2 : i + 2 ) = TRIM_VALUE( i + 3 : i + 3 )
                  TRIM_VALUE( i + 3 : i + 3 ) = ' '
                END IF
              END IF
            END IF
          END IF
          EXIT
        END IF

!  remove trailing 0 unless it is preceeded by a .

        IF ( TRIM_VALUE( i : i ) == ' ' ) THEN
          IF ( i < 3 ) EXIT
          IF ( TRIM_VALUE( i - 1 : i - 1 ) == '0' .AND.                        &
               TRIM_VALUE( i - 2 : i - 2 ) /= '.' ) THEN
               TRIM_VALUE( i - 1 : i - 1 ) = ' '
          END IF
          EXIT
        END IF

      END DO

!  if the string starts with a ., add a 0 at the front

      IF ( TRIM_VALUE( 1 : 1 ) == '.' ) THEN
        DO i = 24, 2, -1
          TRIM_VALUE( i : i ) = TRIM_VALUE( i - 1 : i - 1 )
        END DO
        TRIM_VALUE( 1 : 1 ) = '0'
      END IF

!  if the string starts with a ., add a 0 at the front

      IF ( TRIM_VALUE( 1 : 1 ) == '.' ) THEN
        DO i = 24, 2, -1
          TRIM_VALUE( i : i ) = TRIM_VALUE( i - 1 : i - 1 )
        END DO
        TRIM_VALUE( 1 : 1 ) = '0'
      END IF

!  if the string starts with a -., replace by -0. at the front

      IF ( TRIM_VALUE( 1 : 2 ) == '-.' ) THEN
        DO i = 24, 3, -1
          TRIM_VALUE( i : i ) = TRIM_VALUE( i - 1 : i - 1 )
        END DO
        TRIM_VALUE( 2 : 2 ) = '0'
      END IF
      RETURN

!  end of function TRIM_VALUE

      END FUNCTION TRIM_VALUE

  END PROGRAM QPLIB2LP
