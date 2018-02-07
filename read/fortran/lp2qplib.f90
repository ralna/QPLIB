!  lp2qplib2 - convert an ILOG lp file to qplib format
!  The ILOG lp file should be read from standard input (unit 5)
!  and the qplib format will be written to a problem-named file on unit 61
!  Nick Gould, 10th Dec 2016, revised 2nd Feb 2017

      PROGRAM LP2QPLIB

      INTEGER, PARAMETER :: wp = KIND( 1.0D+0 )
      INTEGER, PARAMETER :: dp = KIND( 1.0D+0 )
      REAL( KIND = dp ), PARAMETER :: teneps_d = 10.0_dp * EPSILON( 1.0_dp )

      CHARACTER ( len = 80 ) :: line
      CHARACTER ( len = 1 ) :: rhs
      INTEGER :: status, phase, start, first, last, len_line, range
      LOGICAL :: startc, minus, quad, objmax

      INTEGER, PARAMETER :: n_max = 100000
      INTEGER, PARAMETER :: m_max = 200000
      INTEGER, PARAMETER :: a_max = 1000000
      INTEGER, PARAMETER :: h_max = 1000000
      INTEGER, PARAMETER :: ch_max = 2000000
      REAL ( KIND = wp ), PARAMETER :: infinity = ( 10.0_wp ) ** 19
      REAL ( KIND = wp ), PARAMETER :: infinity_used = ( 10.0_wp ) ** 20

      INTEGER :: i, i_x, i_x2, ind, ip1, out, var_type, m_actual
      INTEGER :: ias, iae, ichs, iche, ipas, ipae, ipchs, ipche
      INTEGER :: a_ne = 0, h_ne = 0, ch_ne = 0, m = 0, n = 0, m_l = 0
      INTEGER :: a_new = 0, ch_new = 0, m_new = 0
      REAL ( KIND = wp ) :: val
      INTEGER, DIMENSION( n_max + 1 ) :: X_type = 0
      INTEGER, DIMENSION( m_max + 1 ) :: A_start
      INTEGER, DIMENSION( m_max + 1 ) :: CH_start
      INTEGER, DIMENSION( m_max + 1 ) :: C_actual
      REAL ( KIND = wp ), DIMENSION( n_max ) :: G = 0.0_wp
      REAL ( KIND = wp ), DIMENSION( n_max ) :: X_l = 0.0_wp
      REAL ( KIND = wp ), DIMENSION( n_max ) :: X_u = infinity_used
      REAL ( KIND = wp ), DIMENSION( m_max ) :: C_l = - infinity_used
      REAL ( KIND = wp ), DIMENSION( m_max ) :: C_u = infinity_used
      CHARACTER ( len = 34 ), DIMENSION( m_max ) :: char_c
      INTEGER, DIMENSION( a_max ) :: A_row
      INTEGER, DIMENSION( a_max ) :: A_col
      REAL ( KIND = wp ), DIMENSION( a_max ) :: A_val
      INTEGER, DIMENSION( h_max ) :: H_row
      INTEGER, DIMENSION( h_max ) :: H_col
      REAL ( KIND = wp ), DIMENSION( h_max ) :: H_val
      INTEGER, DIMENSION( ch_max ) :: CH_ind
      INTEGER, DIMENSION( ch_max ) :: CH_row
      INTEGER, DIMENSION( ch_max ) :: CH_col
      REAL ( KIND = wp ), DIMENSION( ch_max ) :: CH_val

      INTEGER, PARAMETER :: objective = 1
      INTEGER, PARAMETER :: constraints = 2
      INTEGER, PARAMETER :: bounds = 3
      INTEGER, PARAMETER :: generals = 4
      INTEGER, PARAMETER :: binaries = 5
      INTEGER, PARAMETER :: end = 6

!  extracts from qplib_main.f90 from CUTEst

      INTEGER, PARAMETER :: input = 55
      INTEGER, PARAMETER :: standard_out = 6
      INTEGER, PARAMETER :: qplib_out = 61

      INTEGER, PARAMETER :: qp = 1
      INTEGER, PARAMETER :: qcqp = 2
      INTEGER, PARAMETER :: bqp = 3
      INTEGER, PARAMETER :: lp = 4
      INTEGER, PARAMETER :: qcp = 5

      INTEGER :: int_var, bin_var, discrete_var, continuous_var
      INTEGER :: problem_type, l
      LOGICAL :: filexx
      CHARACTER ( len = 16 ) :: char_i, char_j, char_l
      CHARACTER ( len = 24 ) :: char_val
      CHARACTER ( len = 28 ) :: out_p_name = REPEAT( ' ', 28 )
      CHARACTER ( len = 34 ) :: out_p_name_qplib
      REAL ( KIND = wp ) :: mode_v

      out = standard_out

      DO
        CALL FIND_DATA_LINE( line, status )
        IF ( status /= 0 ) EXIT
        IF ( line( 2 : 9 ) == "ENCODING" ) THEN
          CYCLE

!  obtain output problem name

        ELSE IF ( line( 2 : 8 ) == "Problem" ) THEN
          len_line = len_trim( line )
          i = INDEX( line, "QPLIB" )
          out_p_name = TRIM( line( i : len_line - 3 ) )
          CYCLE
        ELSE IF ( line( 1 : 8 ) == "Minimize" ) THEN
          phase = objective
          objmax = .FALSE.
          CYCLE
        ELSE IF ( line( 1 : 8 ) == "Maximize" ) THEN
          phase = objective
          objmax = .TRUE.
          CYCLE
        ELSE IF ( line( 1 : 10 ) == "Subject To" ) THEN
          phase = constraints
          CYCLE
        ELSE IF ( line( 1 : 6 ) == "Bounds" ) THEN
          phase = bounds
          CYCLE
        ELSE IF ( line( 1 : 8 ) == "Generals" ) THEN
          phase = generals
          CYCLE
        ELSE IF ( line( 1 : 8 ) == "Binaries" ) THEN
          phase = binaries
          CYCLE
        ELSE IF ( line( 1 : 3 ) == "End" ) THEN
          phase = end
          A_start( m + 1 ) = a_ne + 1
          CH_start( m + 1 ) = ch_ne + 1
          EXIT
        ELSE
          len_line = len_trim( line )

!  disect lines within the objective function section

          IF ( phase == objective ) THEN
            start = 1
            DO
              CALL NEXT_substring( line, len_line, start, first, last )
              startc = line( last : last ) == ':'
              IF ( startc ) THEN
                minus = .FALSE.
                quad = .FALSE.
                val = 1.0_wp
              ELSE
                IF ( line(first:first) == '[' ) THEN
                  quad = .TRUE.
                  minus = .FALSE.
                  i_x = - 1 ; i_x2 = - 1
                ELSE IF ( line(first:first) == ']' ) THEN
                  EXIT
                ELSE IF ( line(first:last) == '+' ) THEN
                  minus = .FALSE.
                  val = 1.0_wp
                  IF ( quad ) THEN
                    i_x = - 1  ; i_x2 = - 1
                  END IF
                ELSE IF ( line(first:last) == '-' ) THEN
                  minus = .TRUE.
                  val = 1.0_wp
                  IF ( quad ) THEN
                    i_x = - 1  ; i_x2 = - 1
                  END IF
                ELSE IF ( line(first:first) == '*' ) THEN
                ELSE IF ( line(first:last) == '^2' ) THEN
                  h_ne = h_ne + 1
                  H_row( h_ne ) = i_x
                  H_col( h_ne ) = i_x
                  IF ( minus ) val = - val
!                 IF ( objmax ) val = - val
                  H_val( h_ne ) = val
                ELSE IF ( line(first:first) == 'x' ) THEN
                  IF ( quad ) THEN
                    IF ( i_x < 0 ) THEN
                      READ( line(first+1:last), * ) i_x
                      n = MAX( n, i_x )
                      IF ( n > n_max ) THEN
                        WRITE( 6, "( ' error: n_max = ', I0, ' too small')" )  &
                          n_max
                        STOP
                      END IF
                    ELSE
                      READ( line(first+1:last), * ) i_x2
                      n = MAX( n, i_x2 )
                      IF ( n > n_max ) THEN
                        WRITE( 6, "( ' error: n_max = ', I0, ' too small')" )  &
                          n_max
                        STOP
                      END IF
                      h_ne = h_ne + 1
                      IF ( h_ne > h_max ) THEN
                        WRITE( 6, "( ' error: h_max = ', I0, ' too small')" )  &
                          h_max
                        STOP
                      END IF
                      H_row( h_ne ) = i_x
                      H_col( h_ne ) = i_x2
                      IF ( minus ) val = - val
!                     IF ( objmax ) val = - val
                      H_val( h_ne ) = val
                    END IF
                  ELSE
                    READ( line(first+1:last), * ) i_x
                    n = MAX( n, i_x )
                    IF ( n > n_max ) THEN
                      WRITE( 6, "( ' error: n_max = ', I0, ' too small')" )    &
                        n_max
                      STOP
                    END IF
                    IF ( minus ) val = - val
!                   IF ( objmax ) val = - val
                    g( i_x ) = val
                  END IF
                ELSE
                  READ( line(first:last), * ) val
                END IF
              END IF
              IF ( last == len_line ) EXIT
              start = last + 1
            END DO

!  disect lines within the constraint section

          ELSE IF ( phase == constraints ) THEN
            start = 1
            DO
              CALL NEXT_substring( line, len_line, start, first, last )
              startc = line( first : first ) == 'c' .OR.                       &
                       line( first : first ) == 'q'
              IF ( startc ) THEN
                IF ( line( first : first ) == 'c' ) THEN
                  READ( line(first+1:last-1), * ) m
                  m_l = m
                ELSE
                  READ( line(first+1:last-1), * ) i
                  m = m_l + i
                END IF
!               write(6,*) line(first:last), m
                IF ( m > m_max ) THEN
                  WRITE( 6, "( ' error: m_max = ', I0, ' too small')" ) m_max
                  STOP
                END IF
                minus = .FALSE.
                val = 1.0_wp
                quad = .FALSE.
                A_start( m ) = a_ne + 1
                CH_start( m ) = ch_ne + 1
                rhs = ' '
              ELSE
                IF ( line(first:first) == '[' ) THEN
                  quad = .TRUE.
                  minus = .FALSE.
                  i_x = - 1 ; i_x2 = - 1
                ELSE IF ( line(first:first) == ']' ) THEN
                ELSE IF ( line(first:last) == '=' ) THEN
                  rhs = 'e'
                ELSE IF ( line(first:last) == '>=' ) THEN
                  rhs = 'l'
                ELSE IF ( line(first:last) == '<=' ) THEN
                  rhs = 'u'
                ELSE IF ( line(first:last) == '+' ) THEN
                  minus = .FALSE.
                  val = 1.0_wp
                  IF ( quad ) THEN
                    i_x = - 1  ; i_x2 = - 1
                  END IF
                ELSE IF ( line(first:last) == '-' ) THEN
                  minus = .TRUE.
                  val = 1.0_wp
                  IF ( quad ) THEN
                    i_x = - 1  ; i_x2 = - 1
                  END IF
                ELSE IF ( line(first:first) == '*' ) THEN
                ELSE IF ( line(first:last) == '^2' ) THEN
                  ch_ne = ch_ne + 1
                  IF ( ch_ne > ch_max ) THEN
                    WRITE( 6, "( ' error: ch_max = ', I0, ' too small')" ) ch_ne
                    STOP
                  END IF
                  CH_ind( ch_ne ) = m
                  CH_row( ch_ne ) = i_x
                  CH_col( ch_ne ) = i_x
                  IF ( minus ) val = - val
                  CH_val( ch_ne ) = 2.0_wp * val
                ELSE IF ( line(first:first) == 'x' ) THEN
                  IF ( quad ) THEN
                    IF ( i_x < 0 ) THEN
                      READ( line(first+1:last), * ) i_x
                      n = MAX( n, i_x )
                      IF ( n > n_max ) THEN
                        WRITE( 6, "( ' error: n_max = ', I0, ' too small')" )  &
                          n_max
                        STOP
                      END IF
                    ELSE
                      READ( line(first+1:last), * ) i_x2
                      n = MAX( n, i_x2 )
                      IF ( n > n_max ) THEN
                        WRITE( 6, "( ' error: n_max = ', I0, ' too small')" )  &
                          n_max
                        STOP
                      END IF
                      ch_ne = ch_ne + 1
                      IF ( ch_ne > ch_max ) THEN
                        WRITE( 6, "( ' error: ch_max = ', I0, ' too small')" ) &
                          ch_ne
                        STOP
                      END IF
                      CH_ind( ch_ne ) = m
                      CH_row( ch_ne ) = i_x
                      CH_col( ch_ne ) = i_x2
                      IF ( minus ) val = - val
                      CH_val( ch_ne ) = 2.0_wp * val
                    END IF
                  ELSE
                    READ( line(first+1:last), * ) i_x
                    n = MAX( n, i_x )
                    IF ( n > n_max ) THEN
                      WRITE( 6, "( ' error: n_max = ', I0, ' too small')" )    &
                        n_max
                      STOP
                    END IF
                    a_ne = a_ne + 1
                    IF ( a_ne > a_max ) THEN
                      WRITE( 6, "( ' error: a_max = ', I0, ' too small')" )    &
                        a_max
                      STOP
                    END IF
                    A_row( a_ne ) = m
                    A_col( a_ne ) = i_x
                    IF ( minus ) val = - val
                    A_val( a_ne ) = val
                  END IF
                ELSE IF ( line(last:last) == ':' ) THEN
                  WRITE( out, "( ' Error: keyword = ', A, ' in ', A )" )       &
                    line(first:last), TRIM( out_p_name )
                  STOP
                ELSE
                  READ( line(first:last), * ) val
                  IF ( rhs == 'e' ) THEN
                    C_l( m ) = val ; C_u( m ) = val
                  ELSE IF ( rhs == 'l' ) THEN
                    C_l( m ) = val
                  ELSE IF ( rhs == 'u' ) THEN
                    C_u( m ) = val
                  ELSE
                  END IF
                END IF
              END IF
              IF ( last == len_line ) EXIT
              start = last + 1
            END DO

!  disect lines within the bounds section

          ELSE IF ( phase == bounds ) THEN
            start = 1
            CALL NEXT_substring( line, len_line, start, first, last )
            IF ( line(first:first) == 'x' ) THEN
              READ( line(first+1:last), * ) i_x
              n = MAX( n, i_x )
              IF ( n > n_max ) THEN
                WRITE( 6, "( ' error: n_max = ', I0, ' too small')" ) n_max
                STOP
              END IF
              start = last + 1
              CALL NEXT_substring( line, len_line, start, first, last )
              IF ( line(first:last) == '=' ) THEN
                start = last + 1
                CALL NEXT_substring( line, len_line, start, first, last )
                READ( line(first:last), * ) val
                x_l( i_x ) = val ; x_u( i_x ) = val
              ELSE IF ( line(first:last) == '<=' ) THEN
                rhs = 'u'
                start = last + 1
                CALL NEXT_substring( line, len_line, start, first, last )
                READ( line(first:last), * ) val
                x_u( i_x ) = val
              ELSE IF ( line(first:last) == 'Free' ) THEN
                x_l( i_x ) = - infinity ; x_u( i_x ) = infinity
              END IF
            ELSE IF ( line(first:first+3) == '-Inf' ) THEN
              val = - infinity_used
              start = last + 1
              CALL NEXT_substring( line, len_line, start, first, last )
              start = last + 1
              CALL NEXT_substring( line, len_line, start, first, last )
              READ( line(first+1:last), * ) i_x
              n = MAX( n, i_x )
              IF ( n > n_max ) THEN
                WRITE( 6, "( ' error: n_max = ', I0, ' too small')" ) n_max
                STOP
              END IF
              x_l( i_x ) = val
              IF ( last < len_line ) THEN
                start = last + 1
                CALL NEXT_substring( line, len_line, start, first, last )
                start = last + 1
                CALL NEXT_substring( line, len_line, start, first, last )
                IF ( line(first:last) == 'Inf' ) THEN
                  x_u( i_x ) = infinity_used
                ELSE
                  READ( line(first:last), * ) val
                  x_u( i_x ) = val
                END IF
              END IF
            ELSE
              READ( line(first:last), * ) val
              start = last + 1
              CALL NEXT_substring( line, len_line, start, first, last )
              start = last + 1
              CALL NEXT_substring( line, len_line, start, first, last )
              READ( line(first+1:last), * ) i_x
              n = MAX( n, i_x )
              IF ( n > n_max ) THEN
                WRITE( 6, "( ' error: n_max = ', I0, ' too small')" ) n_max
                STOP
              END IF
              x_l( i_x ) = val
              IF ( last < len_line ) THEN
                start = last + 1
                CALL NEXT_substring( line, len_line, start, first, last )
                start = last + 1
                CALL NEXT_substring( line, len_line, start, first, last )
                IF ( line(first:last) == 'Inf' ) THEN
                  x_u( i_x ) = infinity_used
                ELSE
                  READ( line(first:last), * ) val
                  x_u( i_x ) = val
                END IF
              END IF
            END IF

!  disect lines within the binaries section

          ELSE IF ( phase == binaries ) THEN
            start = 1
            DO
              CALL NEXT_substring( line, len_line, start, first, last )
              READ( line(first+1:last), * ) i_x
              n = MAX( n, i_x )
              IF ( n > n_max ) THEN
                WRITE( 6, "( ' error: n_max = ', I0, ' too small')" ) n_max
                STOP
              END IF
              X_type( i_x ) = 2
              X_l( i_x ) = 0.0_wp ; X_u( i_x ) = 1.0_wp
              IF ( last == len_line ) EXIT
              start = last + 1
            END DO

!  disect lines within the generals section

          ELSE IF ( phase == generals ) THEN
            start = 1
            DO
              CALL NEXT_substring( line, len_line, start, first, last )
              READ( line(first+1:last), * ) i_x
              n = MAX( n, i_x )
              IF ( n > n_max ) THEN
                WRITE( 6, "( ' error: n_max = ', I0, ' too small')" ) n_max
                STOP
              END IF
              X_type( i_x ) = 1
              IF ( last == len_line ) EXIT
              start = last + 1
            END DO
          END IF
        END IF
      END DO

      m_q = m - m_l

!  lp file now processed

      WRITE( 6, "( A, 6( I10 ) )" )                                            &
        TRIM( out_p_name ), n, m_l, m_q, h_ne, a_ne, ch_ne

!  trivial check for superflous range constraints

      IF ( m > 0 ) THEN
        range = 0 ; i = 0 ; m_actual = 0
    90  CONTINUE
        i = i + 1
        m_actual = m_actual + 1
        c_actual( i ) = m_actual
        IF ( i <= m - 1 ) THEN
          ip1 = i + 1
          ias = A_start( i ) ; iae = A_start( ip1 ) - 1
          ipas = A_start( ip1 ) ; ipae = A_start( i + 2 ) - 1
          IF ( iae - ias /= ipae - ipas ) GO TO 90
          ichs = CH_start( i ) ; iche = CH_start( ip1 ) - 1
          ipchs = CH_start( ip1 ) ; ipche = CH_start( i + 2 ) - 1
          IF ( iche - ichs /= ipche - ipchs ) GO TO 90
          IF ( iae >= ias ) THEN
            IF ( COUNT( A_col( ias : iae ) /= A_col( ipas : ipae ) ) > 0 .OR.  &
                 COUNT( A_val( ias : iae ) /= A_val( ipas : ipae ) ) > 0 )     &
                   GO TO 90
            IF ( iche >= ichs ) THEN
              IF ( COUNT( CH_col( ichs : iche ) /= CH_col( ipchs : ipche ) )>0 &
              .OR. COUNT( CH_val( ichs : iche ) /= CH_val( ipchs : ipche ) )>0 &
                ) GO TO 90
          ELSE
            IF ( COUNT( CH_col( ichs : iche ) /= CH_col( ipchs : ipche ) )>0   &
            .OR. COUNT( CH_val( ichs : iche ) /= CH_val( ipchs : ipche ) )>0 ) &
              GO TO 90
          END IF
          range = range + 1
!         WRITE(6, "( ' whooa, constraints ', I0, ' and ',I0,' are the same')")&
!            i, ip1
            c_actual( ip1 ) = - i
            C_l( i ) = MAX( C_l( i ), C_l( ip1 ) )
            C_u( i ) = MIN( C_u( i ), C_u( ip1 ) )
            i = ip1
          END IF
          GO TO 90
        END IF

!  if there are superflous range constraints, remove them

        IF ( range > 0 ) THEN
          WRITE( 6, "( 'NB: there are ', I0, ' range constraint(s)' )" ) range
!         DO i = 1, m
!           IF ( c_actual( i ) > 0 ) THEN
!             WRITE(6,"( ' constraint ', I6, ' is now ', I6 )" ) i, c_actual( i)
!           ELSE
!             WRITE(6,"( ' constraint ', I6, ' is a copy of the new ', I6 )" ) &
!               i, c_actual( - c_actual( i ) )
!           END IF
!         END DO

!  remove superflous range constraints from CH

          IF ( problem_type == qcqp .OR. problem_type == qcp ) THEN
            ch_new = 0
            DO l = 1, ch_ne
              ind = c_actual( CH_ind( l ) )
              IF ( ind > 0 ) THEN
                ch_new = ch_new + 1
                CH_ind( ch_new ) = ind
                CH_row( ch_new ) = CH_row( l )
                CH_col( ch_new ) = CH_col( l )
                CH_val( ch_new ) = CH_val( l )
              END IF
            END DO
            ch_ne = ch_new
          END IF

!  remove superflous range constraints from A

          IF ( problem_type /= bqp ) THEN
            a_new = 0
            DO l = 1, A_ne
              i = c_actual( A_row( l ) )
              IF ( i > 0 ) THEN
                a_new = a_new + 1
                A_row( a_new ) = i
                A_col( a_new ) = A_col( l )
                A_val( a_new ) = A_val( l )
              END IF
            END DO
            A_ne = a_new

!  remove superflous lower and upper bounds from c_l and c_u

            m_new = 0
            DO i = 1, m_l
              ind = c_actual( i )
              IF ( ind > 0 ) THEN
                m_new = m_new + 1
                C_l( m_new ) = C_l( i )
                C_u( m_new ) = C_u( i )
                char_c( m_new ) = REPEAT( ' ', 34 )
                char_c( m_new ) = 'c' // TRIM( TRIM_INT( i ) )
              ELSE
                char_c( m_new ) = 'c' // TRIM( TRIM_INT( i - 1 ) ) // '&' //   &
                                  TRIM( TRIM_INT( i ) )
              END IF
            END DO
            i_x = m_new
            DO i = m_l + 1, m
              ind = c_actual( i )
              IF ( ind > 0 ) THEN
                m_new = m_new + 1
                C_l( m_new ) = C_l( i )
                C_u( m_new ) = C_u( i )
                char_c( m_new ) = REPEAT( ' ', 34 )
                char_c( m_new ) = 'q' // TRIM( TRIM_INT( i - m_l ) )
              ELSE
                char_c( m_new ) = 'q' // TRIM( TRIM_INT( i - m_l - 1 ) ) //    &
                                    '&' // TRIM( TRIM_INT( i - m_l ) )
              END IF
            END DO
            m_l = i_x
            m = m_new
            m_q = m - m_l
          END IF
        END IF
      END IF

!  ===================================
!  write the data in QPLIB file format
!  ===================================

!  open output file

      out_p_name_qplib = REPEAT( ' ', 34 )
      out_p_name_qplib = TRIM( out_p_name ) // ".qplib"
      INQUIRE( FILE = out_p_name_qplib, EXIST = filexx )
      IF ( filexx ) THEN
         OPEN( qplib_out, FILE = out_p_name_qplib, FORM = 'FORMATTED',         &
               STATUS = 'OLD', IOSTAT = status )
      ELSE
         OPEN( qplib_out, FILE = out_p_name_qplib, FORM = 'FORMATTED',         &
               STATUS = 'NEW', IOSTAT = status )
      END IF
      IF ( status /= 0 ) GO TO 900
      out = qplib_out
      REWIND ( out )

!  determine problem type: 1(QP,default),2,(QCQP),3(BQP),4(LP),5(QCP)

      IF ( ch_ne > 0 ) THEN
        IF ( h_ne > 0 ) THEN
          problem_type = 2
        ELSE
          problem_type = 5
        END IF
      ELSE
        IF ( a_ne == 0 ) THEN
          IF ( h_ne > 0 ) THEN
            problem_type = 3
          ELSE
            problem_type = 4
          END IF
        ELSE
          IF ( h_ne > 0 ) THEN
            problem_type = 1
          ELSE
            problem_type = 4
          END IF
        END IF
      END IF

!  see if the problem has integer variables

      int_var = COUNT( X_type( : n ) == 1 )
      bin_var = COUNT( X_type( : n ) == 2 )
      discrete_var = int_var + bin_var
      continuous_var = n - discrete_var

!  set header

      IF ( m > 0 ) THEN
        IF ( LEN( TRIM( out_p_name ) ) <= 24 ) THEN
          WRITE( out, "( A24, ' ILOG ', A, '.lp with n = ', I0, ', m = ',      &
         &   I0 )" ) out_p_name( 1 : 24 ), TRIM( out_p_name ), n, m
        ELSE
          WRITE( out, "( A28, /, 24X, ' ILOG ', A, '.lp with n = ', I0,        &
         & ', m = ', I0 )" ) out_p_name, TRIM( out_p_name ), n, m
        END IF
      ELSE
        WRITE( out, "( A24, ' ILOG ', A,                                       &
       &       '.lp with n = ', I0 )" ) out_p_name, TRIM( out_p_name ), n
      END IF
      IF ( discrete_var == 0 ) THEN
        SELECT CASE ( problem_type )
        CASE ( qcqp )
          WRITE( out, "( 'QCQ                      a quadratic program',       &
         &               ' with quadratic constraints' )" )
        CASE ( bqp )
          WRITE( out, "( 'QCB                      a bound-constrained',       &
         &               ' quadratic program' )" )
        CASE ( lp )
          WRITE( out, "( 'LCL                      a linear program' )" )
        CASE ( qcp )
          WRITE( out, "( 'LCQ                      a linear program',          &
         &                ' with quadratic constraints' )" )
        CASE DEFAULT
          WRITE( out, "( 'QCL                      a quadratic program' )")
        END SELECT
      ELSE IF ( bin_var == n ) THEN
        SELECT CASE ( problem_type )
        CASE ( qcqp )
          WRITE( out, "( 'QBQ                      a binary',                  &
         &    ' QP with quadratic constraints' )" )
        CASE ( bqp )
          WRITE( out, "( 'QBB                      a binary',                  &
         &    ' bound-constrained quadratic program' )" )
        CASE ( lp )
          WRITE( out, "( 'LBL                      a binary',                  &
       &     ' linear program' )" )
        CASE ( qcp )
          WRITE( out, "( 'LBQ                      a binary',                  &
         &    ' LP with quadratic constraints' )" )
        CASE DEFAULT
          WRITE( out, "( 'QBL                      a binary',                  &
         &   ' quadratic program' )")
        END SELECT
      ELSE IF ( discrete_var == n ) THEN
        SELECT CASE ( problem_type )
        CASE ( qcqp )
          WRITE( out, "( 'QIQ                      an integer',                &
         &    ' QP with quadratic constraints' )" )
        CASE ( bqp )
          WRITE( out, "( 'QIB                      an integer',                &
         &    ' bound-constrained quadratic program' )" )
        CASE ( lp )
          WRITE( out, "( 'LIL                      an integer',                &
       &     ' linear program' )" )
        CASE ( qcp )
          WRITE( out, "( 'LIQ                      an integer',                &
         &    ' LP with quadratic constraints' )" )
        CASE DEFAULT
          WRITE( out, "( 'QIL                      an integer',                &
         &   ' quadratic program' )")
        END SELECT
      ELSE
        SELECT CASE ( problem_type )
        CASE ( qcqp )
          WRITE( out, "( 'QGQ                      a mixed-integer',           &
         &    ' QP with quadratic constraints' )" )
        CASE ( bqp )
          WRITE( out, "( 'QGB                      a mixed-integer',           &
         &    ' bound-constrained quadratic program' )" )
        CASE ( lp )
          WRITE( out, "( 'QGL                      a mixed-integer',           &
       &     ' linear program' )" )
        CASE ( qcp )
          WRITE( out, "( 'LGQ                      a mixed-integer',           &
         &    ' LP with quadratic constraints' )" )
        CASE DEFAULT
          WRITE( out, "( 'QGL                      a mixed-integer',           &
         &   ' quadratic program' )")
        END SELECT
      END IF
      IF ( objmax ) THEN
        WRITE( out, "( 'Maximize                ',                             &
       &               ' # maximize the objective function' )" )
      ELSE
        WRITE( out, "( 'Minimize                ',                             &
       &               ' # minimize the objective function' )" )
      END IF

      char_l = TRIM_INT( n )
      WRITE( out, "( A16, 8X, ' # variables ' )" ) char_l
      IF ( problem_type /= bqp ) THEN
        char_l = TRIM_INT( m )
        WRITE( out, "( A16, 8X, ' # general linear constraints ' )" ) char_l
      END IF
!     IF ( objmax ) WRITE( out, "( '#', 24X, 'minus objective function',       &
!    &   ' specified for maximization' )" )

!  Hessian values

      IF ( problem_type == qp .OR. problem_type == bqp .OR.                    &
           problem_type == qcqp ) THEN
        char_l = TRIM_INT( H_ne )
        IF ( H_ne == 0 ) THEN
          WRITE( out, "( /, A16, 8X, ' # nonzeros in upper triangle of H' )" ) &
            char_l
        ELSE
          WRITE( out, "( /, A16, 8X, ' # nonzeros in upper triangle of H:',    &
         &   ' row,column,value' )" ) char_l
          DO l = 1, H_ne
            char_i = TRIM_INT( H_row( l ) ) ; char_j = TRIM_INT( H_col( l ) )
            char_val = TRIM_VALUE( H_val( l ) )
            WRITE( out, "( A16, 1X, A16, 1X, A24 )" ) char_i, char_j, char_val
          END DO
        END IF
      END IF

!  gradient values

      mode_v = MODE( n, G )
      l = COUNT( G( : n ) /= mode_v )
      char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
      WRITE( out, "( /, A24, ' default value for entries in g' )" ) char_val
      IF ( l == 0 ) THEN
        WRITE( out, "( A16, 8X, ' # non default entries in g' )" ) char_l
      ELSE
        WRITE( out, "( A16, 8X, ' # non default entries in g: index,value' )") &
          char_l
        DO i = 1, n
          IF ( G( i ) /= mode_v ) THEN
            char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( G( i ) )
            WRITE( out, "( A16, 1X, A24 )" ) char_i, char_val
          END IF
        END DO
      END IF

!  function value

      char_val = TRIM_VALUE( 0.0_wp )
      WRITE( out, "( /, A24, ' value of f' )" ) char_val

!  Hessian values for constraints

      IF ( problem_type == qcqp .OR. problem_type == qcp ) THEN

!  record the total number of constraint Hessian values

        char_l = TRIM_INT( ch_ne )
        WRITE( out, "( /, A16, 8X, ' # nonzeros in upper triangle',            &
       &  ' of the H_i')") char_l

!  append the constraint Hessian values to the qplib file

        DO l = 1, ch_ne
          char_l = TRIM_INT( CH_ind( l ) )
          char_i = TRIM_INT( CH_row( l ) ) ; char_j = TRIM_INT( CH_col( l ) )
          char_val = TRIM_VALUE( CH_val( l ) )
          WRITE( out, "( A16, 1X, A16, 1X, A16, 1X, A24 )" )                   &
             char_l, char_i, char_j, char_val
        END DO
      END IF

!  constraint Jacobian values

      IF ( problem_type /= bqp ) THEN
        char_l = TRIM_INT( A_ne )
        IF ( A_ne == 0 ) THEN
          WRITE( out, "( /, A16, 8X, ' # nonzeros in A' )" ) char_l
        ELSE
          WRITE( out, "( /, A16, 8X, ' # nonzeros in A:',                      &
         &   ' row,column,value' )" ) char_l
          DO l = 1, A_ne
            char_i = TRIM_INT( A_row( l ) ) ; char_j = TRIM_INT( A_col( l ) )
            char_val = TRIM_VALUE( A_val( l ) )
            WRITE( out, "( A16, 1X, A16, 1X, A24 )" ) char_i, char_j, char_val
          END DO
        END IF
      END IF

!  infinity

      char_val = TRIM_VALUE( infinity )
      WRITE( out, "( /, A24, ' value of infinite bounds' )" ) char_val

!  constraint lower bounds

      IF ( problem_type /= bqp ) THEN
        IF ( m > 0 ) THEN
          mode_v = MODE( m, C_l )
          l = COUNT( C_l( : m ) /= mode_v )
          char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
          WRITE( out, "( /, A24, ' default value for entries in c_l' )" )      &
            char_val
          IF ( l == 0 ) THEN
            WRITE( out, "( A16, 8X, ' # non default entries in c_l' )" ) char_l
          ELSE
            WRITE( out, "( A16, 8X, ' # non default entries in c_l:',          &
           &  ' index,value' )" ) char_l
            DO i = 1, m
              IF ( C_l( i ) /= mode_v ) THEN
                char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( C_l( i ) )
                WRITE( out, "( A16, 1X, A24 )" ) char_i, char_val
              END IF
            END DO
          END IF
        ELSE
          mode_v = - infinity_used
          l = 0
          char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
          WRITE( out, "( /, A24, ' default value for entries in c_l' )" )      &
            char_val
          WRITE( out, "( A16, 8X, ' # non default entries in c_l' )" ) char_l
        END IF

!  constraint upper bounds

        IF ( m > 0 ) THEN
          mode_v = MODE( m, C_u )
          l = COUNT( C_u( : m ) /= mode_v )
          char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
          WRITE( out, "( /, A24, ' default value for entries in c_u' )" )      &
            char_val
          IF ( l == 0 ) THEN
            WRITE( out, "( A16, 8X, ' # non default entries in c_u' )" ) char_l
          ELSE
            WRITE( out, "( A16, 8X, ' # non default entries in c_u:',          &
          &   ' index,value' )" ) char_l
            DO i = 1, m
              IF ( C_u( i ) /= mode_v ) THEN
                char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( C_u( i ) )
                WRITE( out, "( A16, 1X, A24 )" ) char_i, char_val
              END IF
            END DO
          END IF
        ELSE
          mode_v = infinity_used
          l = 0
          char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
          WRITE( out, "( /, A24, ' default value for entries in c_u' )" )      &
            char_val
          WRITE( out, "( A16, 8X, ' # non default entries in c_u' )" ) char_l
        END IF
      END IF

!  variable lower bounds

      IF ( bin_var < n ) THEN
        mode_v = MODE( n, X_l )
        l = COUNT( X_l( : n ) /= mode_v )
        char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
        WRITE( out, "( /, A24, ' default value for entries in x_l' )" ) char_val
        IF ( l == 0 ) THEN
          WRITE( out, "( A16, 8X, ' # non default entries in x_l' )" ) char_l
        ELSE
          WRITE( out, "( A16, 8X, ' # non default entries in x_l:',            &
         &  ' index,value' )" ) char_l
          DO i = 1, n
            IF ( X_l( i ) /= mode_v ) THEN
              char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( X_l( i ) )
              WRITE( out, "( A16, 1X, A24 )" ) char_i, char_val
            END IF
          END DO
        END IF

!  variable upper bounds

        mode_v = MODE( n, X_u )
        l = COUNT( X_u( : n ) /= mode_v )
        char_l = TRIM_INT( l ) ; char_val = TRIM_VALUE( mode_v )
        WRITE( out, "( /, A24, ' default value for entries in x_u' )" ) char_val
        IF ( l == 0 ) THEN
          WRITE( out, "( A16, 8X, ' # non default entries in x_u' )" ) char_l
        ELSE
          WRITE( out, "( A16, 8X, ' # non default entries in x_u:',            &
         &  ' index,value' )") char_l
          DO i = 1, n
            IF ( X_u( i ) /= mode_v ) THEN
              char_i = TRIM_INT( i ) ; char_val = TRIM_VALUE( X_u( i ) )
              WRITE( out, "( A16, 1X, A24 )" ) char_i, char_val
            END IF
          END DO
        END IF
      END IF

!  variable types

      IF ( discrete_var > 0 .AND. discrete_var < n ) THEN

!  find the dominant variable type

        IF ( continuous_var >= int_var ) THEN
          IF ( continuous_var >= bin_var ) THEN
            var_type = 0
            char_j = TRIM_INT( n - continuous_var )
          ELSE
            var_type = 2
            char_j = TRIM_INT( n - bin_var )
          END IF
        ELSE
          IF ( int_var >= bin_var ) THEN
            var_type = 1
            char_j = TRIM_INT( n - int_var )
          ELSE
            var_type = 2
            char_j = TRIM_INT( n - bin_var )
          END IF
        END IF

!  record non-dominant type variables

        char_l = TRIM_INT( var_type )
        WRITE( out, "( /, A16, 8X, ' default variable type',                   &
       &  ' (0=contin., 1=integer, 2=binary)' )" ) char_l

        WRITE( out, "( A16, 8X, ' # non default variables: index,type' )")     &
       &  char_j
        DO i = 1, n
          IF ( X_type( i ) /= var_type ) THEN
            char_i = TRIM_INT( i ) ; char_j = TRIM_INT( X_type( i ) )
            WRITE( out, "( A16, 1X, A16 )" ) char_i, char_j
          END IF
        END DO
      END IF

!  initial primal variables

      char_l = TRIM_INT( 0 )
      char_val = TRIM_VALUE( 0.0_wp )
      WRITE( out, "( /, A24, ' default value for entries in initial x' )" )    &
        char_val
      WRITE( out, "( A16, 8X, ' # non default entries in x' )" ) char_l

!  initial Lagrange multipliers

      IF ( problem_type /= bqp ) THEN
        WRITE( out, "( /, A24, ' default value for entries in initial y' )" )  &
          char_val
        WRITE( out, "( A16, 8X, ' # non default entries in y' )" ) char_l
      END IF

!  initial dual variables

      WRITE( out, "( /, A24, ' default value for entries in initial z' )" )    &
        char_val
      WRITE( out, "( A16, 8X, ' # non default entries in z' )" ) char_l

!  variable names

      WRITE( out, "( /, A16, 8X, ' # non default names for variables:',        &
     &   ' index,name' )" ) char_l

!  constraint names

      IF ( problem_type /= bqp ) THEN
        IF ( range > 0 ) THEN
          char_l = TRIM_INT( m )
          WRITE( out, "( /, A16, 8X, ' # non default names for constraints:',  &
         &   ' index,name' )" ) char_l
          DO i = 1, m
            char_i = TRIM_INT( i )
            WRITE( out, "( A16, 1X, A )" ) char_i, char_c( i )
          END DO
        ELSE
          IF ( m_q == 0 ) THEN
            WRITE( out, "( /, A16, 8X, ' # non default names for constraints:',&
           &   ' index,name' )" ) char_l
          ELSE
            char_l = TRIM_INT( m_q )
            WRITE( out, "( /, A16, 8X, ' # non default names for constraints:',&
           &   ' index,name' )" ) char_l
            DO i = 1, m_q
              char_i = TRIM_INT( m_l + i )
              WRITE( out, "( A16, 1X, A )" ) char_i, 'q' // TRIM( TRIM_INT( i ))
            END DO
          END IF
        END IF
      END IF

      CLOSE( qplib_out )

      STOP

!  error exits

 900  CONTINUE
      WRITE( out, "( ' error status = ', I0 )" ) status
      CLOSE( INPUT  )
      STOP

    CONTAINS

!  -------------- F I N D _ D A T A _ L I N E   S U B R O U T I N E -----------

      SUBROUTINE FIND_DATA_LINE( line, status )

!  find the next non blank, non comment line

      CHARACTER ( len = 80 ) :: line
      INTEGER :: status
      CHARACTER ( len = 80 ) :: blank = REPEAT( ' ', 80 )

      DO
        line = blank
        READ( 5, "( A80 )", err = 90, end = 90 ) line
        line = ADJUSTL( line )
!       IF ( line /= blank .AND. line( 1 : 1 ) /= "\" ) THEN
        IF ( line /= blank ) THEN
          EXIT
        END IF
      END DO
      status = 0
      RETURN

   90 CONTINUE
      status = 1
      RETURN

      END SUBROUTINE FIND_DATA_LINE

!  -------------- N E X T _ S U B S T R I N G   S U B R O U T I N E -----------

      SUBROUTINE NEXT_substring( line, len_line, start, first, last )

!  find the next substring starting at or beyond position start in line;
!  the substring is line(first:last)

      CHARACTER ( len = 80 ) :: line
      INTEGER, INTENT( IN ) :: len_line, start
      INTEGER, INTENT( OUT ) :: first, last
      INTEGER :: i

      DO i = start, len_line
        IF ( line( i : i ) /= ' ' ) THEN
          first = i
          EXIT
        END IF
      END DO

      IF ( first == len_line ) THEN
        last = len_line
      ELSE
        DO i = first + 1, len_line
          IF ( line( i : i ) == ' ' ) THEN
            last = i - 1
            EXIT
          END IF
          IF ( i == len_line ) last = len_line
        END DO
      END IF
      RETURN

      END SUBROUTINE NEXT_substring

!  extracts from qplib_main.f90 from CUTEst

!  ------------------------ M O D E  F U N C T I O N --------------------------

      FUNCTION MODE( n, V )
      IMPLICIT NONE
      REAL ( KIND = wp ) :: MODE
      INTEGER, INTENT( IN ) :: n
      REAL ( KIND = wp ), INTENT( IN ), DIMENSION( n ) :: V

!  find the "mode", i.e., the most commonly-occuring value, of a vector v

      INTEGER :: i, mode_start, max_len, same, len, m, inform

      REAL ( KIND = wp ), DIMENSION( n ) :: V_sorted

      IF ( n > 0 ) THEN

!  sort a copy of v into increasing order

        V_sorted = V
        CALL SORT_heapsort_build( n, V_sorted, inform ) !  build the heap
        DO i = 1, n
          m = n - i + 1
          CALL SORT_heapsort_smallest( m, V_sorted, inform ) !  reorder v
        END DO

!  run through the sorted values, finding adjacent entries that are identical

        mode_start = 1 ; max_len = 1
        same = 1 ; len = 1
        DO i = 2, n
          IF ( V_sorted( i ) /= V_sorted( same ) ) THEN
            IF ( len > max_len ) THEN
              mode_start = same
              max_len = len
            END IF
            same = i ; len = 1
          ELSE
            len = len + 1
          END IF
        END DO
        IF ( len > max_len ) THEN
          mode_start = same
          max_len = len
        END IF
!       write(6,*) max_len
!       write(6,*) V_sorted( : n )
        MODE = V_sorted( mode_start )
      ELSE
        MODE = 0.0_wp
      END IF
      RETURN

      END FUNCTION MODE

!  --------------- H E A P S O R T   S U B R O U T I N E S -------------------

!  heapsort routines extracted from GALAHAD

      SUBROUTINE SORT_heapsort_build( n, A, inform )

!  Given an array A, elements A(1), ...., A(N), subroutine SORT_heapsort_build
!  rearranges the elements to form a heap in which each parent has a smaller
!  value than either of its children.

!  Algorithm 232 of CACM (J. W. J. Williams): a combination of
!  ALGOL procedures SETHEAP and INHEAP

!  Programming: Nick Gould, January 26th 1995.

!  ------------------------- dummy arguments --------------------------
!
!  n      integer, which gives the number of values to be sorted.
!         n must be positive
!
!  A      real array of length n. On input, A must contain the
!         values which are to be sorted. On output, these values
!         will have been permuted so as to form a heap
!
!  inform integer, which informs the user of the success of SORT_heapsort_build.
!         If inform = 0 on exit, the heap has been formed.
!         If inform = 1 on exit, n was input with a value less than
!                       or equal to 0 and the heap has not been formed.
!
!  ------------------ end of dummy arguments --------------------------

      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: n
      INTEGER, INTENT( OUT ) :: inform
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( n ) :: A

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, j, k
      REAL ( KIND = wp )  :: rin

!  Add the elements to the heap one at a time

      IF ( n <= 0 ) THEN
         inform = 1
         RETURN
      ENDIF

      DO k = 2, n
        rin = A( k )

!  The cycle may be repeated log2(k) times, but on average is repeated
!  merely twice

        i = k
        DO
          IF ( i <= 1 ) EXIT
          j = i / 2
          IF ( A( j ) <= rin ) EXIT
          A( i ) = A( j )
          i = j
        END DO
        A( i ) = rin
      END DO
      inform = 0

      RETURN

!  End of subroutine SORT_heapsort_build

     END SUBROUTINE SORT_heapsort_build

     SUBROUTINE SORT_heapsort_smallest( m, A, inform )

!  Given an array A, elements A(1), ...., A(m) forming a heap,
!  SORT_heapsort_smallest assigns to rout the value of A(1), the smallest
!  member of the heap, and arranges the remaining members as elements
!  1 to m - 1 of A. rout is then placed in A(m)

!  Algorithm 232 of CACM (J. W. J. Williams): a combination of
!  ALGOL procedures OUTHEAP and SWOPHEAP

!  Programming: Nick Gould, January 26th 1995.

!  ------------------------- dummy arguments --------------------------
!
!  m      integer, which gives the number of values to be sorted.
!         m must be positive
!
!  A      real array of length m. On input, A must contain the values which
!         are to be sorted stored in a heap. On output, the smallest value
!         will have been moved into A(m) and the remaining values A(k),
!         k = 1,..., m-1 will have been restored to a heap
!
!  inform integer, which informs the user of the success of
!         SORT_heapsort_smallest.
!         If inform = 0 on exit, the smallest value has been found.
!         If inform = 1 on exit, m was input with a value less than
!                       or equal to 0 and the heap has not been formed
!
!  ------------------ end of dummy arguments --------------------------

      IMPLICIT NONE
      INTEGER, INTENT( IN ) :: m
      INTEGER, INTENT( OUT ) :: inform
      REAL ( KIND = wp ), INTENT( INOUT ), DIMENSION( m ) :: A

!-----------------------------------------------
!   L o c a l   V a r i a b l e s
!-----------------------------------------------

      INTEGER :: i, j
      REAL ( KIND = wp ) :: rin, rout

!  Add the element rin to the heap, extract and assign to rout
!  the value of the smallest member of the resulting set, and
!  leave the remaining elements in a heap of the original size.
!  In this process, elements 1 to n+1 of the array A may be disturbed

      IF ( m <= 0 ) THEN
         inform = 1
         RETURN
      ENDIF

      IF ( m > 1 ) THEN
        i = 1
        rout = A( 1 )
        rin = A( m )

!  Move from the top of the heap comparing the value of node i
!  with its two daughters. If node i is smallest, the heap has been
!  restored. If one of the children is smallest, promote this child
!  in the heap and move to the now vacated node.
!  This cycle may be repeated log2(m) times

        DO
          j = i + i
          IF ( j > m - 1 ) EXIT

!  Determine which of the two daughters is smallest

          IF ( A( j + 1 ) < A( j ) ) j = j + 1

!  Determine if the smaller daughter is less than the value from node i

          IF ( A( j ) >= rin ) EXIT
          A( i ) = A( j )
          i = j
        END DO

!  The heap has been restored

        A( i ) = rin

!  Store the smallest value in the now vacated m-th position of the list

        A( m ) = rout

      END IF
      inform = 0

      RETURN

!  End of subroutine SORT_heapsort_smallest

      END SUBROUTINE SORT_heapsort_smallest

!  ---------------------- T R I M _ I N T    F U N C T I O N ------------------

      FUNCTION TRIM_INT( i )
      CHARACTER ( LEN = 16 ) :: TRIM_INT
      INTEGER :: i

!  write integer as a left shifted length 16 character

      TRIM_INT = REPEAT( ' ', 16 )
      WRITE( TRIM_INT, "( I0 )" ) i
      RETURN

!  end of function TRIM_INT

      END FUNCTION TRIM_INT

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
      IF ( value == 0.0_wp ) THEN
        field = "0.0         "
      ELSE IF ( SIGN( 1.0_wp, value ) > 0.0_wp ) THEN
        IF ( value >= ( 10.0_wp ) ** 100 ) THEN
          WRITE( field24, "( ES24.15E3 )" ) value
          field = field24( 1 : 20 ) // field24( 22 : 24 )
        ELSE IF ( value >= ( 10.0_wp ) ** 16 ) THEN
          WRITE( field24, "( ES24.15 )" ) value
          field = field24( 1 : 21 ) // field24( 23 : 24 )
        ELSE IF ( value >= ( 10.0_wp ) ** 15 ) THEN
          WRITE( field, "( F23.0 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 14 ) THEN
          WRITE( field, "( F23.1 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 13 ) THEN
          WRITE( field, "( F23.2 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 12 ) THEN
          WRITE( field, "( F23.3 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 11 ) THEN
          WRITE( field, "( F23.4 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 10 ) THEN
          WRITE( field, "( F23.5 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 9 ) THEN
          WRITE( field, "( F23.6 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 8 ) THEN
          WRITE( field, "( F23.7 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 7 ) THEN
          WRITE( field, "( F23.8 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 6 ) THEN
          WRITE( field, "( F23.9 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 5 ) THEN
          WRITE( field, "( F23.10 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 4 ) THEN
          WRITE( field, "( F23.11 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 3 ) THEN
          WRITE( field, "( F23.12 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 2 ) THEN
          WRITE( field, "( F23.13 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 1 ) THEN
          WRITE( field, "( F23.14 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** 0 ) THEN
          WRITE( field, "( F23.15 )" ) value
        ELSE IF ( value >= ( 10.0_wp ) ** ( - 1 ) ) THEN
          WRITE( field24, "( F24.16 )" ) value
          field = field24( 2 : 24 )
        ELSE IF ( value >= ( 10.0_wp ) ** ( - 2 ) ) THEN
          WRITE( field24, "( F24.17 )" ) value
          field = field24( 2 : 24 )
        ELSE IF ( value >= ( 10.0_wp ) ** ( - 3 ) ) THEN
          WRITE( field24, "( F24.18 )" ) value
          field = field24( 2 : 24 )
        ELSE IF ( value >= ( 10.0_wp ) ** ( - 4 ) ) THEN
          WRITE( field24, "( F24.16 )" ) value
          field = field24( 2 : 24 )
        ELSE IF ( value >= ( 10.0_wp ) ** ( - 9 ) ) THEN
          WRITE( field24, "( ES24.15 )" ) value
          field = field24( 1 : 22 ) // field24( 24 : 24 )
        ELSE IF ( value >= ( 10.0_wp ) ** ( - 99 ) ) THEN
          WRITE( field, "( ES23.15 )" ) value
!       ELSE IF ( value >= ( 10.0_wp ) ** ( - 999 ) ) THEN
!         WRITE( field, "( ES23.15E3 )" ) value
        ELSE
          WRITE( field, "( ES23.15E4 )" ) value
        END IF
      ELSE
        minus_value = - value
        IF ( ABS( minus_value - 1.0_wp ) <= teneps_d ) minus_value = 1.0_wp
        IF ( minus_value >= ( 10.0_wp ) ** 100 ) THEN
          WRITE( field, "( ES23.15E3 )" ) minus_value
          field22 = field( 1 : 19 ) // field( 21 : 23 )
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 16 ) THEN
          WRITE( field, "( ES23.15 )" ) minus_value
          field22 = field( 1 : 20 ) // field( 22 : 23 )
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 15 ) THEN
          WRITE( field22, "( F22.0 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 14 ) THEN
          WRITE( field22, "( F22.1 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 13 ) THEN
          WRITE( field22, "( F22.2 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 12 ) THEN
          WRITE( field22, "( F22.3 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 11 ) THEN
          WRITE( field22, "( F22.4 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 10 ) THEN
          WRITE( field22, "( F22.5 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 9 ) THEN
          WRITE( field22, "( F22.6 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 8 ) THEN
          WRITE( field22, "( F22.7 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 7 ) THEN
          WRITE( field22, "( F22.8 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 6 ) THEN
          WRITE( field22, "( F22.9 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 5 ) THEN
          WRITE( field22, "( F22.10 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 4 ) THEN
          WRITE( field22, "( F22.11 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 3 ) THEN
          WRITE( field22, "( F22.12 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 2 ) THEN
          WRITE( field22, "( F22.13 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 1 ) THEN
          WRITE( field22, "( F22.14 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** 0 ) THEN
          WRITE( field22, "( F22.15 )" ) minus_value
        ELSE IF ( minus_value >= ( 10.0_wp ) ** ( - 1 ) ) THEN
          WRITE( field, "( F23.16 )" ) minus_value
          field22 = field( 2 : 23 )
        ELSE IF ( minus_value >= ( 10.0_wp ) ** ( - 2 ) ) THEN
          WRITE( field, "( F23.17 )" ) minus_value
          field22 = field( 2 : 23 )
        ELSE IF ( minus_value >= ( 10.0_wp ) ** ( - 3 ) ) THEN
          WRITE( field, "( F23.18 )" ) minus_value
          field22 = field( 2 : 23 )
        ELSE IF ( minus_value >= ( 10.0_wp ) ** ( - 4 ) ) THEN
          WRITE( field, "( F23.15 )" ) minus_value
          field22 = field( 2 : 23 )
        ELSE IF ( minus_value >= ( 10.0_wp ) ** ( - 9 ) ) THEN
          WRITE( field, "( ES23.15 )" ) minus_value
          field22 = field( 1 : 21 ) // field( 23 : 23 )
        ELSE IF ( minus_value > ( 10.0_wp ) ** ( - 99 ) ) THEN
          WRITE( field22, "( ES22.15 )" ) minus_value
!       ELSE IF ( minus_value > ( 10.0_wp ) ** ( - 999 ) ) THEN
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

      i_point = 24 ; i_end = 1
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

  END PROGRAM LP2QPLIB
