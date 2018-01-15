*> \brief \b ZTREVC3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download ZTREVC3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/ztrevc3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/ztrevc3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/ztrevc3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE ZTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT,
*      $                    VL, LDVL, VR, LDVR, MM, M, WORK, LWORK, INFO)
*
*       .. Scalar Arguments ..
*       CHARACTER          HOWMNY, SIDE
*       INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
*       ..
*       .. Array Arguments ..
*       LOGICAL            SELECT( * )
*       COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
*      $                   WORK( * )
*       ..
*
*
*> \par Purpose:
*  =============
*>
*> \verbatim
*>
*> ZTREVC3 computes some or all of the right and/or left eigenvectors of
*> a complex upper triangular matrix T.
*> Matrices of this type are produced by the Schur factorization of
*> a complex general matrix:  A = Q*T*Q**H, as computed by ZHSEQR.
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
*>          T is COMPLEX*16 array, dimension (LDT,N)
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
*>          VL is COMPLEX*16 array, dimension (LDVL,MM)
*>          On entry, if SIDE = 'L' or 'B' and HOWMNY = 'B', VL must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of
*>          Schur vectors returned by ZHSEQR).
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
*>          VR is COMPLEX*16 array, dimension (LDVR,MM)
*>          On entry, if SIDE = 'R' or 'B' and HOWMNY = 'B', VR must
*>          contain an N-by-N matrix Q (usually the unitary matrix Q of
*>          Schur vectors returned by ZHSEQR).
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
*>          WORK is COMPLEX*16 array, dimension (MAX(1,LWORK))
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
      SUBROUTINE ZTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT,
     $                    VL, LDVL, VR, LDVR, MM, M, WORK, LWORK, INFO )
      IMPLICIT NONE
*
*  -- LAPACK computational routine (version 3.8.0) --
*  -- LAPACK is a software package provided by Univ. of Tennessee,    --
*  -- Univ. of California Berkeley, Univ. of Colorado Denver and NAG Ltd..--
*     November 2017
*
*     .. Scalar Arguments ..
      CHARACTER          HOWMNY, SIDE
      INTEGER            INFO, LDT, LDVL, LDVR, LWORK, M, MM, N
*     ..
*     .. Array Arguments ..
      LOGICAL            SELECT( * )
      COMPLEX*16         T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( N, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0, HALF = 5.0D-1 )
      COMPLEX*16         CZERO, CONE, CNONE
      PARAMETER          ( CZERO = ( 0.0D+0, 0.0D+0 ),
     $                     CONE  = ( 1.0D+0, 0.0D+0 ),
     $                     CNONE = ( -1.0D+0, 0.0D+0 ) )
      INTEGER            NBMAX, MBMIN, MBMAX
      PARAMETER          ( NBMAX = 32, MBMIN = 8, MBMAX = 128 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV, BACKTRANSFORM, SOMEV, SELECTV
      DOUBLE PRECISION   UNFL, OVFL, ULP, SMLNUM, SCALE, TNORM, VNORM,
     $                   TEMP
      COMPLEX*16         CTEMP
      CHARACTER          NORMIN
      INTEGER            I, IB, IBEND, IMIN, IMAX,
     $                   J, JV, JB, JBEND, JOUT, JMAX,
     $                   NB, MB
*     ..
*     .. Local Arrays ..
      INTEGER            JLIST( NBMAX )
      DOUBLE PRECISION   SCALES( NBMAX ), BOUNDS( NBMAX ),
     $                   CNORMS( NBMAX )
      COMPLEX*16         SHIFTS( NBMAX ), DIAG( NBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV, IZAMAX
      DOUBLE PRECISION   DLAMCH, DZASUM, ZLANGE
      EXTERNAL           LSAME, ILAENV, DLAMCH, IZAMAX, DZASUM, ZLANGE
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DLABAD, ZDSCAL, ZGEMM, ZLACPY, ZLATRS
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, AIMAG, CONJG, MIN, MAX
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   CABS1
*     ..
*     .. Statement Function definitions ..
      CABS1( CTEMP ) = ABS( DBLE( CTEMP ) ) + ABS( AIMAG( CTEMP ) )
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
      MB = ILAENV( 1, 'ZTREVC', SIDE // HOWMNY, N, -1, -1, -1 )
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
         CALL XERBLA( 'ZTREVC3', -INFO )
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
      UNFL = DLAMCH( 'Safe minimum' )
      OVFL = ONE / UNFL
      CALL DLABAD( UNFL, OVFL )
      ULP = DLAMCH( 'Precision' )
      SMLNUM = UNFL*( N / ULP )

*
*     --------------------------------------------------------------
*     Compute right eigenvectors.
*
      IF( RIGHTV ) THEN
         JOUT = M
         JB = JMAX + 1
         JBEND = JMAX
         DO 100 JV = N, 1, -1
*
*           --------------------------------------------------------
*           Add current eigenvector to workspace (if needed).
*
            IF( SOMEV ) THEN
               SELECTV = SELECT( JV )
            ELSE
               SELECTV = .TRUE.
            END IF
            IF( SELECTV ) THEN
               IF( JB.GT.JBEND ) THEN
                  IMAX = JV
               END IF
               JB = JB - 1
               WORK( 1 : JV - 1, JB ) = -T( 1 : JV - 1, JV )
               WORK( JV : N, JB ) = CZERO
               JLIST( JB ) = JV
               SHIFTS( JB ) = T( JV, JV )
               BOUNDS( JB ) = DZASUM( JV - 1, WORK( 1, JB ), 1 )
               IF( BOUNDS( JB ).GT.OVFL ) THEN
                  SCALES( JB ) = HALF * OVFL / BOUNDS( JB )
                  BOUNDS( JB ) = HALF * OVFL
                  CALL ZDSCAL( JV - 1, SCALES( JB ), WORK( 1, JB ), 1 )
               ELSE
                  SCALES( JB ) = ONE
               END IF
            END IF
*
*           --------------------------------------------------------
*           Process workspace if full or if all eigenvectors are
*           found.
*
            IF( JB.EQ.1 .OR. ( JV.EQ.1 .AND. JB.LE.JBEND ) ) THEN
               MB = JBEND - JB + 1
*
*              -----------------------------------------------------
*              Compute triangular eigenvectors with safe,
*              multi-shift, blocked back substitution.
*     
               DO 110 IBEND = IMAX, 1, -NBMAX
                  IB = MAX( IBEND - NBMAX + 1, 1 )
                  DO 120 I = IB, IBEND
                     DIAG( I - IB + 1 ) = T( I, I )
 120              CONTINUE
                  TNORM = ZLANGE( 'O', IB - 1, MB,
     $                            T( 1, IB ), LDT, CTEMP )
                  DO 130 J = JBEND, JB, -1
                     NB = MIN( JLIST( J ) - 1, IBEND) - IB + 1
                     IF( NB.LE.0 )
     $                  GO TO 130
*                  
*                    Safeguarded solve with shifted diagonal block.
*
                     TEMP = MAX( ULP * CABS1( SHIFTS( J ) ), SMLNUM )
                     DO 140 I = IB, IB + NB - 1
                        T( I, I ) = DIAG( I - IB + 1 ) - SHIFTS( J )
                        IF( CABS1( T( I, I ) ).LT.TEMP )
     $                     T( I, I ) = DCMPLX( TEMP )
 140                 CONTINUE
                     IF( J.EQ.JBEND ) THEN
                        NORMIN = 'N'
                     ELSE
                        NORMIN = 'Y'
                     END IF
                     CALL ZLATRS( 'U', 'N', 'N', NORMIN, NB,
     $                            T( IB, IB ), LDT,
     $                            WORK( IB, J ), SCALE,
     $                            CNORMS, INFO )
*
*                    Rescale solution (if needed).
*
                     VNORM = DZASUM( NB, WORK( IB, J ), 1 )
                     BOUNDS( J ) = BOUNDS( J ) * SCALE
                     TEMP = OVFL - BOUNDS( J )
                     IF( VNORM.GE.ONE .AND. TNORM.GT.TEMP/VNORM ) THEN
                        TEMP = ( OVFL * HALF / TNORM ) / VNORM
                        SCALE = TEMP * SCALE
                        BOUNDS( J ) = TEMP * BOUNDS( J ) + OVFL * HALF
                     ELSE IF( VNORM.LT.ONE
     $                        .AND. TNORM*VNORM.GT.TEMP ) THEN
                        TEMP = OVFL * HALF / TNORM
                        SCALE = TEMP * SCALE
                        BOUNDS( J ) = TEMP * BOUNDS( J )
     $                                + OVFL * HALF * VNORM
                     ELSE
                        BOUNDS( J ) = BOUNDS( J ) + TNORM * VNORM
                     END IF
                     IF( SCALE.NE.ONE ) THEN
                        CALL ZDSCAL( IB - 1, SCALE, WORK( 1, J ), 1 )
                        CALL ZDSCAL( JLIST( J ) - IBEND - 1, SCALE,
     $                               WORK( IBEND + 1, J ), 1 )
                        SCALES( J ) = SCALES( J ) * SCALE
                     END IF
 130              CONTINUE
                  DO 150 I = IB, IBEND
                     T( I, I ) = DIAG( I - IB + 1 )
 150              CONTINUE
*
*                 Back substitution with block of solution.
*
                  NB = IBEND - IB + 1
                  IF( IB.GT.1 ) THEN
                     CALL ZGEMM( 'N', 'N', IB - 1, MB, NB,
     $                           CNONE, T( 1, IB ), LDT,
     $                           WORK( IB, JB ), N,
     $                           CONE, WORK( 1, JB ), N )
                  END IF
 110           CONTINUE
*
*              Put scale factors on diagonal to get triangular
*              eigenvectors.
*
               DO 160 J = JB, JBEND
                  WORK( JLIST( J ), J ) = SCALES( J )
 160           CONTINUE
*
*              -----------------------------------------------------
*              Copy results to output.
*                  
               JOUT = JOUT - MB + 1
               IF( BACKTRANSFORM ) THEN
*                  
*                 Back transform with Schur vectors to get full
*                 eigenvectors.
*
                  CALL ZGEMM( 'N', 'N', N, MB, IMAX,
     $                        CONE, VR, LDVR, WORK( 1, JB ), N,
     $                        CZERO, WORK( 1, JMAX + JB ), N )
                  CALL ZLACPY( 'F', N, MB,
     $                         WORK( 1, JMAX + JB ), N,
     $                         VR( 1, JOUT ), LDVR )
               ELSE
                  CALL ZLACPY( 'F', N, MB,
     $                         WORK( 1, JB ), N,
     $                         VR( 1, JOUT ), LDVR )
               END IF
*              
*              Workspace is now empty.
*
               JB = JMAX + 1
            END IF
 100     CONTINUE
*
*        -----------------------------------------------------
*        Normalize eigenvectors.
*
         DO 170 J = 1, M
            I = IZAMAX( N, VR( 1, J ), 1 )
            SCALE = ONE / CABS1( VR( I, J ) )
            CALL ZDSCAL( N, SCALE, VR( 1, J ), 1 )
 170     CONTINUE
      END IF

*
*     --------------------------------------------------------------
*     Compute left eigenvectors.
*
      IF( LEFTV ) THEN
         JOUT = 1
         JB = 1
         JBEND = 0
         DO 300 JV = 1, N
*
*           --------------------------------------------------------
*           Add current eigenvector to workspace (if needed).
*
            IF( SOMEV ) THEN
               SELECTV = SELECT( JV )
            ELSE
               SELECTV = .TRUE.
            END IF
            IF( SELECTV ) THEN
               IF( JBEND.EQ.0 ) THEN
                  IMIN = JV
               END IF
               JBEND = JBEND + 1
               WORK( 1 : JV, JBEND ) = CZERO
               DO 310 I = JV + 1, N
                  WORK( I, JBEND ) = -CONJG( T( JV, I ) )
 310           CONTINUE
               JLIST( JBEND ) = JV
               SHIFTS( JBEND ) = T( JV, JV )
               BOUNDS( JBEND ) = DZASUM( N - JV,
     $                                   WORK( JV + 1, JBEND ), 1 )
               IF( BOUNDS( JBEND ).GT.OVFL ) THEN
                  SCALES( JBEND ) = HALF * OVFL / BOUNDS( JBEND )
                  BOUNDS( JBEND ) = HALF * OVFL
                  CALL ZDSCAL( N - JV, SCALES( JBEND ),
     $                         WORK( JV + 1, JBEND ), 1 )
               ELSE
                  SCALES( JBEND ) = ONE
               END IF
            END IF
*
*           --------------------------------------------------------
*           Process workspace if full or if all eigenvectors are
*           found.
*
            IF( JBEND.EQ.JMAX .OR. ( JV.EQ.N .AND. JB.LE.JBEND ) ) THEN
               MB = JBEND - JB + 1
*
*              -----------------------------------------------------
*              Compute triangular eigenvectors with safe,
*              multi-shift, blocked forward substitution.
*
               DO 320 IB = IMIN, N, NBMAX
                  IBEND = MIN( IB + NBMAX - 1, N )
                  DO 330 I = IB, IBEND
                     DIAG( I - IB + 1 ) = T( I, I )
 330              CONTINUE
                  TNORM = ZLANGE( 'O', MB, N - IBEND,
     $                            T( IB, IBEND + 1 ), LDT, CTEMP )
                  DO 340 J = JB, JBEND
                     NB = IBEND - MAX( JLIST( J ) + 1, IB ) + 1
                     IF( NB.LE.0 )
     $                    GO TO 340
*                  
*                    Safeguarded solve with shifted diagonal block.
*
                     TEMP = MAX( ULP * CABS1( SHIFTS( J ) ), SMLNUM )
                     DO 350 I = IBEND - NB + 1, IBEND
                        T( I, I ) = DIAG( I - IB + 1 ) - SHIFTS( J )
                        IF( CABS1( T( I, I ) ).LT.TEMP )
     $                     T( I, I ) = DCMPLX( TEMP )
 350                 CONTINUE
                     IF( J.EQ.JB ) THEN
                        NORMIN = 'N'
                     ELSE
                        NORMIN = 'Y'
                     END IF
                     CALL ZLATRS( 'U', 'C', 'N', NORMIN, NB,
     $                            T( IBEND - NB + 1, IBEND - NB + 1 ),
     $                            LDT,
     $                            WORK( IBEND - NB + 1, J ), SCALE,
     $                            CNORMS, INFO )
*
*                    Rescale solution (if needed).
*
                     VNORM = DZASUM( NB, WORK( IBEND - NB + 1, J ), 1 )
                     BOUNDS( J ) = BOUNDS( J ) * SCALE
                     TEMP = OVFL - BOUNDS( J )
                     IF( VNORM.GE.ONE .AND. TNORM.GT.TEMP/VNORM ) THEN
                        TEMP = ( OVFL * HALF / TNORM ) / VNORM
                        SCALE = TEMP * SCALE
                        BOUNDS( J ) = TEMP * BOUNDS( J ) + OVFL * HALF
                     ELSE IF( VNORM.LT.ONE
     $                        .AND. TNORM*VNORM.GT.TEMP ) THEN
                        TEMP = OVFL * HALF / TNORM
                        SCALE = TEMP * SCALE
                        BOUNDS( J ) = TEMP * BOUNDS( J )
     $                                + OVFL * HALF * VNORM
                     ELSE
                        BOUNDS( J ) = BOUNDS( J ) + TNORM * VNORM
                     END IF
                     IF( SCALE.NE.ONE ) THEN
                        CALL ZDSCAL( IB - JLIST( J ) - 1, SCALE,
     $                               WORK( JLIST( J ) + 1, J ), 1 )
                        CALL ZDSCAL( N - IBEND, SCALE,
     $                               WORK( IBEND + 1, J ), 1 )
                        SCALES( J ) = SCALES( J ) * SCALE
                     END IF
 340              CONTINUE
                  DO 360 I = IB, IBEND
                     T( I, I ) = DIAG( I - IB + 1 )
 360              CONTINUE
*
*                 Forward substitution with block of solution.
*
                  NB = IBEND - IB + 1
                  IF( IBEND.LT.N ) THEN
                     CALL ZGEMM( 'C', 'N', N - IBEND, MB, NB,
     $                           CNONE, T( IB, IBEND + 1 ), LDT,
     $                           WORK( IB, JB ), N,
     $                           CONE, WORK( IBEND + 1, JB ), N )
                  END IF
 320           CONTINUE
*
*              Put scale factors on diagonal to get triangular
*              eigenvectors.
*
               DO 370 J = JB, JBEND
                  WORK( JLIST( J ), J ) = SCALES( J )
 370           CONTINUE
*
*              -----------------------------------------------------
*              Copy results to output.
*                  
               IF( BACKTRANSFORM ) THEN
*                  
*                 Back transform with Schur vectors to get full
*                 eigenvectors.
*
                  CALL ZGEMM( 'N', 'N', N, MB, N - IMIN + 1,
     $                        CONE, VL( 1, IMIN ), LDVL,
     $                        WORK( IMIN, JB ), N,
     $                        CZERO, WORK( 1, JMAX + JB ), N )
                  CALL ZLACPY( 'F', N, MB,
     $                         WORK( 1, JMAX + JB ), N,
     $                         VL( 1, JOUT ), LDVL )
               ELSE
                  CALL ZLACPY( 'F', N, MB,
     $                         WORK( 1, JB ), N,
     $                         VL( 1, JOUT ), LDVL )
               END IF
               JOUT = JOUT + MB
*                  
*              Workspace is now empty.
*
               JBEND = 0
            END IF
 300     CONTINUE
*
*        -----------------------------------------------------
*        Normalize eigenvectors.
*
         DO 380 J = 1, M
            I = IZAMAX( N, VL( 1, J ), 1 )
            SCALE = ONE / CABS1( VL( I, J ) )
            CALL ZDSCAL( N, SCALE, VL( 1, J ), 1 )
 380     CONTINUE
      END IF

      RETURN
*
*     End of ZTREVC3
*
      END
