*> \brief \b CTREVC3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download CTREVC3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ctrevc3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ctrevc3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ctrevc3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*     SUBROUTINE CTREVC3( SIDE, HOWMNY, SELECT, N,
*                          T, LDT, VL, LDVL, VR, LDVR,
*                          MM, M, WORK, LWORK, INFO)
*
*       .. Scalar Arguments ..
*       CHARACTER          HOWMNY, SIDE
*       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
*       ..
*       .. Array Arguments ..
*       LOGICAL            SELECT( * )
*       COMPLEX            T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> CTREVC3 computes some or all of the right and/or left eigenvectors of
*> a complex upper triangular matrix T.
*> Matrices of this type are produced by the Schur factorization of
*> a complex general matrix:  A = Q*T*Q**H, as computed by CHSEQR.
*>
*> The right eigenvector x and the left eigenvector y of T corresponding
*> to an eigenvalue w are defined by:
*>
*>              T*x = w*x,     (y**H)*T = w*(y**H)
*>
*> where y**H denotes the conjugate transpose of the vector y.
*> The eigenvalues are not input to this routine, but are read directly
*> from the diagonal of T.
*>
*> This routine returns the matrices X and/or Y of right and left
*> eigenvectors of T, or the products Q*X and/or Q*Y, where Q is an
*> input matrix. If Q is the unitary factor that reduces a matrix A to
*> Schur form T, then Q*X and Q*Y are the matrices of right and left
*> eigenvectors of A.
*>
*> This uses a blocked version of the algorithm.
*> \endverbatim
*
*  Arguments:
*  ==========
*
*> \param[in] SIDE
*> \verbatim
*>          SIDE is CHARACTER*1
*>          = 'R':  compute right eigenvectors only;
*>          = 'L':  compute left eigenvectors only;
*>          = 'B':  compute both right and left eigenvectors.
*> \endverbatim
*>
*> \param[in] HOWMNY
*> \verbatim
*>          HOWMNY is CHARACTER*1
*>          = 'A':  compute all right and/or left eigenvectors;
*>          = 'B':  compute all right and/or left eigenvectors,
*>                  backtransformed using the matrices supplied in
*>                  VR and/or VL;
*>          = 'S':  compute selected right and/or left eigenvectors,
*>                  as indicated by the logical array SELECT.
*> \endverbatim
*>
*> \param[in] SELECT
*> \verbatim
*>          SELECT is LOGICAL array, dimension (N)
*>          If HOWMNY = 'S', SELECT specifies the eigenvectors to be
*>          computed.
*>          The eigenvector corresponding to the j-th eigenvalue is
*>          computed if SELECT(j) = .TRUE..
*>          Not referenced if HOWMNY = 'A' or 'B'.
*> \endverbatim
*>
*> \param[in] N
*> \verbatim
*>          N is INTEGER
*>          The order of the matrix T. N >= 0.
*> \endverbatim
*>
*> \param[in,out] T
*> \verbatim
*>          T is COMPLEX array, dimension (LDT,N)
*>          The upper triangular matrix T.  T is modified, but restored
*>          on exit.
*> \endverbatim
*>
*> \param[in] LDT
*> \verbatim
*>          LDT is INTEGER
*>          The leading dimension of the array T. LDT >= max(1,N).
*> \endverbatim
*>
*> \param[in,out] VL
*> \verbatim
*>          VL is COMPLEX array, dimension (LDVL,MM)
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of
*>          Schur vectors returned by CHSEQR).
*>          On exit, if SIDE = 'L' or 'B', VL contains:
*>          if HOWMNY = 'A', the matrix Y of left eigenvectors of T;
*>          if HOWMNY = 'B', the matrix Q*Y;
*>          if HOWMNY = 'S', the left eigenvectors of T specified by
*>                           SELECT, stored consecutively in the columns
*>                           of VL, in the same order as their
*>                           eigenvalues.
*>          Not referenced if SIDE = 'R'.
*> \endverbatim
*>
*> \param[in] LDVL
*> \verbatim
*>          LDVL is INTEGER
*>          The leading dimension of the array VL.
*>          LDVL >= 1, and if SIDE = 'L' or 'B', LDVL >= N.
*> \endverbatim
*>
*> \param[in,out] VR
*> \verbatim
*>          VR is COMPLEX array, dimension (LDVR,MM)
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of
*>          Schur vectors returned by CHSEQR).
*>          On exit, if SIDE = 'R' or 'B', VR contains:
*>          if HOWMNY = 'A', the matrix X of right eigenvectors of T;
*>          if HOWMNY = 'B', the matrix Q*X;
*>          if HOWMNY = 'S', the right eigenvectors of T specified by
*>                           SELECT, stored consecutively in the columns
*>                           of VR, in the same order as their
*>                           eigenvalues.
*>          Not referenced if SIDE = 'L'.
*> \endverbatim
*>
*> \param[in] LDVR
*> \verbatim
*>          LDVR is INTEGER
*>          The leading dimension of the array VR.
*>          LDVR >= 1, and if SIDE = 'R' or 'B', LDVR >= N.
*> \endverbatim
*>
*> \param[in] MM
*> \verbatim
*>          MM is INTEGER
*>          The number of columns in the arrays VL and/or VR. MM >= M.
*> \endverbatim
*>
*> \param[out] M
*> \verbatim
*>          M is INTEGER
*>          The number of columns in the arrays VL and/or VR actually
*>          used to store the eigenvectors.
*>          If HOWMNY = 'A' or 'B', M is set to N.
*>          Each selected eigenvector occupies one column.
*> \endverbatim
*>
*> \param[out] WORK
*> \verbatim
*>          WORK is COMPLEX array, dimension (MAX(1,LWORK))
*> \endverbatim
*>
*> \param[in] LWORK
*> \verbatim
*>          LWORK is INTEGER
*>          The dimension of array WORK. LWORK >= max(1,2*N).
*>          For optimum performance, LWORK >= 2*N*NB, where NB is
*>          the optimal blocksize.
*>
*>          If LWORK = -1, then a workspace query is assumed; the routine
*>          only calculates the optimal size of the WORK array, returns
*>          this value as the first entry of the WORK array, and no error
*>          message related to LWORK is issued by XERBLA.
*> \endverbatim
*>
*> \param[out] INFO
*> \verbatim
*>          INFO is INTEGER
*>          = 0:  successful exit
*>          < 0:  if INFO = -i, the i-th argument had an illegal value
*> \endverbatim
*
*  Authors:
*  ========
*
*> \author Univ. of Tennessee
*> \author Univ. of California Berkeley
*> \author Univ. of Colorado Denver
*> \author NAG Ltd.
*
*> \date November 2017
*
*> \ingroup complexOTHERcomputational
*
*> \par Further Details:
*  =====================
*>
*> \verbatim
*>
*>  The algorithm used in this program is basically blocked backward
*>  (forward) substitution, with scaling to make the the code robust
*>  against possible overflow.
*>
*>  Each eigenvector is normalized so that the element of largest
*>  magnitude has magnitude 1; here the magnitude of a complex number
*>  (x,y) is taken to be |x| + |y|.
*> \endverbatim
*>
*  =====================================================================
      SUBROUTINE CTREVC3( SIDE, HOWMNY, SELECT, N,
     $                    T, LDT, VL, LDVL, VR, LDVR,
     $                    MM, M, WORK, LWORK, INFO)
     IMPLICIT NONE
*
*  -- LAPACK computational routine (version 3.8.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      CHARACTER HOWMNY, SIDE
      INTEGER   INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL   SELECT( * )
      COMPLEX   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ), WORK( N, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      REAL      ZERO, ONE, HALF
      PARAMETER ( ZERO = 0.0E+0, ONE = 1.0E+0, HALF = 5.0E-1 )
      COMPLEX   CZERO, CONE, CNONE
      PARAMETER ( CZERO = ( 0.0E+0, 0.0E+0 ),
     $            CONE  = ( 1.0E+0, 0.0E+0 ),
     $            CNONE = ( -1.0E+0, 0.0E+0 ) )
      INTEGER   NBMAX, MBMIN, MBMAX
      PARAMETER ( NBMAX = 32, MBMIN = 8, MBMAX = 128 )
*     ..
*     .. Local Scalars ..
      LOGICAL   LEFTV, RIGHTV, BACKTRANSFORM, SOMEV
      REAL      UNFL, OVFL, ULP, SMLNUM, SCALE, BOUND, TNORM, VNORM
      COMPLEX   CTEMP
      CHARACTER NORMIN
      INTEGER   I, IB, IBEND, IMAX, J, JV, JB, JOUT, JMAX, NB, MB
*     ..
*     .. Local Arrays ..
      INTEGER   JLIST( NBMAX )
      REAL      SCALES( NBMAX ), BOUNDS( NBMAX )
      COMPLEX   SHIFTS( NBMAX ), DIAG( NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL   LSAME
      INTEGER   ILAENV
      REAL      SLAMCH, SCASUM, CLANGE
      EXTERNAL  LSAME, ILAENV, ICAMAX, SLAMCH, SCASUM, CLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL  XERBLA, CSSCAL, CGEMM, CLATRS, SLABAD
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC ABS, REAL, AIMAG, CMPLX, MIN, MAX
*     ..
*     .. Statement Functions ..
      REAL      CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CTEMP ) = ABS( REAL( CTEMP ) ) + ABS( AIMAG( CTEMP ) )
*     ..
*     .. Executable Statements ..
*
*     Decode and test the input parameters
*
      RIGHTV = LSAME( SIDE, 'R' ) .OR. LSAME( SIDE, 'B' )
      LEFTV = LSAME( SIDE, 'L' ) .OR. LSAME( SIDE, 'B' )
      BACKTRANSFORM = LSAME( HOWMNY, 'B' )
      SOMEV = LSAME( HOWMNY, 'S' )
*
*     Set M to the number of columns required to store the selected
*     eigenvectors.
*
      IF( SOMEV ) THEN
         M = 0
         DO 10 J = 1, N
            IF( SELECT( J ) )
     $         M = M + 1
   10    CONTINUE
      ELSE
         M = N
      END IF
*
      INFO = 0
      MB = ILAENV( 1, 'CTREVC', SIDE // HOWMNY, N, -1, -1, -1 )
      WORK( 1, 1 ) = 2 * N * MB
      IF( .NOT.RIGHTV .AND. .NOT.LEFTV ) THEN
         INFO = -1
      ELSE IF( .NOT.LSAME( HOWMNY, 'A' )
     $         .AND. .NOT.BACKTRANSFORM
     $         .AND. .NOT.SOMEV ) THEN
         INFO = -2
      ELSE IF( N.LT.0 ) THEN
         INFO = -4
      ELSE IF( LDT.LT.MAX( 1, N ) ) THEN
         INFO = -6
      ELSE IF( LDVL.LT.1 .OR. ( LEFTV .AND. LDVL.LT.N ) ) THEN
         INFO = -8
      ELSE IF( LDVR.LT.1 .OR. ( RIGHTV .AND. LDVR.LT.N ) ) THEN
         INFO = -10
      ELSE IF( MM.LT.M ) THEN
         INFO = -11
      ELSE IF( LWORK.LT.MAX( 1, 2*N ) .AND. LWORK.NE.-1 ) THEN
         INFO = -14
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'CTREVC3', -INFO )
         RETURN
      ELSE IF( LWORK.EQ.-1 ) THEN
         RETURN
      END IF
*
*     Quick return if possible.
*
      IF( N.EQ.0 )
     $   RETURN
*
*     Use blocked version of back-transformation if sufficient workspace.
*
      IF( LWORK .GE. 2*N*MBMIN ) THEN
         JMAX = MIN( LWORK / (2*N), NBMAX )
      ELSE
         JMAX = 1
      END IF
*
*     Set the constants to control overflow.
*
      UNFL = SLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL SLABAD( UNFL, OVFL )
      ULP = SLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )

*
*     ----------------------------------------------------------
*     Compute right eigenvectors
*
      IF( RIGHTV ) THEN
         JOUT = M
         JB = 0
         DO 20 JV = N, 1, -1
*
*           ----------------------------------------------------
*           Set up workspace for current eigenvector if needed
*
            IF( SOMEV ) THEN
               SELECTV = SELECT( JV )
            ELSE
               SELECTV = .TRUE.
            END IF
            IF( SELECTV ) THEN
               IF( JB.EQ.0 ) THEN
                  JB = JMAX
                  IMAX = JV - 1
               ELSE
                  JB = JB - 1
               END IF
               DO 30 I = 1, JV - 1
                  WORK( I, JB ) = -T( I, JV )
 30            CONTINUE
               DO 40 I = JV, IMAX
                  WORK( I, JB ) = CZERO
 40            CONTINUE
               JLIST( JB ) = JV
               SHIFTS( JB ) = T( JV, JV )
               BOUNDS( JB ) = SCASUM( JV - 1, WORK( 1, JB ), 1 )
               IF( BOUNDS( JB ).GT.OVFL ) THEN
                  SCALES( JB ) = HALF * OVFL / BOUNDS( JB )
                  BOUNDS( JB ) = HALF * OVFL
                  CALL CSSCAL( JV - 1, SCALES( JB ), WORK( 1, JB ), 1 )
               ELSE
                  SCALES( JB ) = ONE
               END IF
            END IF
*
*           ----------------------------------------------------
*           Process workspace if needed
*
            IF( JB.EQ.1 .OR. ( JV.EQ.1 .AND. JB.NE.0 ) ) THEN
               MB = JMAX - JB + 1
*
*              -------------------------------------------------
*              Compute triangular eigenvectors
*     
               DO 50 IBEND = IMAX, 1, -NBMAX
                  IB = MAX( IBEND - NBMAX + 1, 1 )
                  DO 60 I = IB, IBEND
                     DIAG( I - IB + 1 ) = T( I, I )
 60               CONTINUE
                  TNORM = CLANGE( 'O', IB - 1, MB,
     $                            T( 1, IB ), LDT, CTEMP )
                  DO 70 J = JMAX, JB, -1
                     NB = MIN( JLIST( J ) - IB, IBEND - IB + 1 )
                     IF( NB.LE.0 )
     $                  GO TO 70
*                  
*                    Safeguarded solve against shifted diagonal block
*
                     BOUND = MAX( ULP * CABS1( SHIFTS( J ) ), SMLNUM )
                     DO 80 I = IB, IB + NB - 1
                        T( I, I ) = DIAG( I - IB + 1 ) - SHIFTS( J )
                        IF( CABS1( T( I, I ) ).LT.BOUND )
     $                     T( I, I ) = CMPLX( BOUND )
 80                  CONTINUE
                     NORMIN = 'Y'
                     IF( J.EQ.JMAX )
     $                  NORMIN = 'N'
                     CALL CLATRS( 'U', 'N', 'N', NORMIN, NB,
     $                            T( IB, IB ), LDT,
     $                            WORK( IB, J ), SCALE,
     $                            CNORMS, INFO )
*
*                    Rescale if needed
*
                     VNORM = SCASUM( NB, WORK( IB, J ), 1 )
                     BOUNDS( J ) = BOUNDS( J ) * SCALE
                     BOUND = OVFL - BOUNDS( J )
                     IF( VNORM.GE.ONE .AND. TNORM.GT.BOUND/VNORM ) THEN
                        SCALE2 = ( OVFL * HALF / TNORM ) / VNORM
                        SCALE = SCALE * SCALE2
                        BOUNDS( J ) = SCALE2 * BOUNDS( J ) + OVFL * HALF
                     ELSE IF( VNORM.LT.ONE
     $                        .AND. TNORM*VNORM.GT.BOUND ) THEN
                        SCALE2 = OVFL * HALF / TNORM
                        SCALE = SCALE * SCALE2
                        BOUNDS( J ) = SCALE2 * BOUNDS( J )
                        BOUNDS( J ) = BOUNDS( J ) + OVFL * HALF * VNORM
                     ELSE
                        BOUNDS( J ) = BOUNDS( J ) + TNORM * VNORM
                     END IF
                     IF( SCALE.NE.ONE ) THEN
                        CALL CSSCAL( IB - 1, SCALE, WORK( 1, J ), 1 )
                        CALL CSSCAL( JLIST( J ) - IBEND - 1, SCALE,
     $                               WORK( IBEND + 1, J ), 1 )
                        SCALES( J ) = SCALES( J ) * SCALE
                     END IF
 70               CONTINUE
                  DO 90 I = IB, IBEND
                     T( I, I ) = DIAG( I - IB + 1 )
 90               CONTINUE
*
*                 Back substitution
*
                  IF( NB.EQ.1 ) THEN
                     CALL CGEMV( 'N', IB - 1, NB,
     $                           CNONE, T( 1, IB ), LDT,
     $                           WORK( IB, JB ), 1,
     $                           CONE, WORK( 1, JB ), 1)
                  ELSE
                     CALL CGEMM( 'N', 'N', IB - 1, MB, NB,
     $                           CNONE, T( 1, IB ), LDT,
     $                           WORK( IB, JB ), N,
     $                           CONE, WORK( 1, JB ), N )
                  END IF
 50            CONTINUE
*
*              Put scale factors on diagonal
*
               DO 100 J = JB, JMAX
                  WORK( JLIST( J ), J ) = SCALES( J )
 100           CONTINUE
*
*              -------------------------------------------------
*              Copy results to output
*                  
               IF( BACKTRANSFORM ) THEN
*                  
*                 Back solve to get full eigenvectors
*
                  IF( NB.EQ.1 ) THEN
                     CALL CGEMV( 'N', N, IMAX,
     $                           CNONE, VR, LDVR, WORK( 1, 1 ), 1,
     $                           CZERO, WORK( 1, 2 ), 1)
                  ELSE
                     CALL CGEMM( 'N', 'N', N, MB, IMAX,
     $                           CONE, VR, LDVR, WORK( 1, JB ), N,
     $                           CZERO, WORK( 1, JMAX + JB ), N )
                  END IF
                  DO 110 J = JMAX, JB, -1
                     SCALE = ONE / CLANGE( 'M', N, 1,
     $                                     WORK( 1, JMAX + J ), N,
     $                                     CTEMP )
                     DO 120 I = N, 1, -1
                        VR( I, JOUT ) = SCALE * WORK( I, JMAX + J )
 120                 CONTINUE
                     JOUT = JOUT - 1
 110              CONTINUE
               ELSE
*                  
*                 Copy triangular eigenvectors to output
*
                  DO 130 J = JMAX, JB, -1
                     SCALE = ONE / CLANGE( 'M', JLIST( J ), 1,
     $                                     WORK( 1, J ), N, CTEMP )
                     DO 140 I = N, JLIST( J ), -1
                        VR( I, JOUT ) = CZERO
 140                 CONTINUE
                     DO 150 I = JLIST( J ), 1, -1
                        VR( I, JOUT ) = SCALE * WORK( I, J )
 150                 CONTINUE
                     JOUT = JOUT - 1
 130              CONTINUE
                  
               END IF
*                  
*              Workspace is now empty
*
               JB = 0
            END IF
 20      CONTINUE
      END IF


      IF( LEFTV ) THEN
*     TODO
      ENDIF

      RETURN
*
*     End of CTREVC3
*
      END
