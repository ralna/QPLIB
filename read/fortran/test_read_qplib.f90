   PROGRAM TEST_READ_QPLIB
   USE GALAHAD_RPD_double
   INTEGER, PARAMETER :: input = 5
   TYPE ( QPT_problem_type ) :: prob
   TYPE ( RPD_inform_type ) :: inform
   CHARACTER ( LEN = 10 ) :: pname
   OPEN( input )
   CALL RPD_read_problem_data( input, prob, inform )
   CLOSE( input )
   IF ( inform%status == 0 ) THEN  !  Successful return
     pname = TRANSFER( prob%name, pname )
     WRITE( 6, "( A10, 5( I10 ) )" )                                           &
        pname, prob%n, prob%m, prob%H%ne, prob%A%ne, prob%H_c%ne
!    WRITE( 6, "( ' problem ', A, ', n = ', I0, ', m = ', I0 )" )              &
!       TRIM( pname ), prob%n, prob%m
   ELSE
     WRITE( 6, "( ' error return from qplib file read, status = ', I0 )" )     &
       inform%status
     WRITE( 6, "( ' last line read = ', I0 )" ) inform%line
   END IF
   STOP
   END PROGRAM TEST_READ_QPLIB
