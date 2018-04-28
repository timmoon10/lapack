*> \brief \b DTREVC3
*
*  =========== DOCUMENTATION ===========
*
* Online html documentation available at
*            http://www.netlib.org/lapack/explore-html/
*
*> \htmlonly
*> Download DTREVC3 + dependencies
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.tgz?format=tgz&filename=/lapack/lapack_routine/dtrevc3.f">
*> [TGZ]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.zip?format=zip&filename=/lapack/lapack_routine/dtrevc3.f">
*> [ZIP]</a>
*> <a href="http://www.netlib.org/cgi-bin/netlibfiles.txt?format=txt&filename=/lapack/lapack_routine/dtrevc3.f">
*> [TXT]</a>
*> \endhtmlonly
*
*  Definition:
*  ===========
*
*       SUBROUTINE DTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT,
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
*> DTREVC3 computes some or all of the right and/or left eigenvectors of
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
      SUBROUTINE DTREVC3( SIDE, HOWMNY, SELECT, N, T, LDT,
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
      DOUBLE PRECISION   T( LDT, * ), VL( LDVL, * ), VR( LDVR, * ),
     $                   WORK( N, * )
*     ..
*
*  =====================================================================
*
*     .. Parameters ..
      DOUBLE PRECISION   ZERO, ONE, NEGONE, HALF
      PARAMETER          ( ZERO = 0.0D+0, ONE = 1.0D+0,
     $                     NEGONE = -1.0D+0, HALF = 5.0D-1 )
      INTEGER            NBMAX, MBMAX
      PARAMETER          ( NBMAX = 64, MBMAX = 64 )
*     ..
*     .. Local Scalars ..
      LOGICAL            LEFTV, RIGHTV, BACKTRANSFORM, SOMEV, SELECTV
      DOUBLE PRECISION   UNFL, OVFL, ULP, SMLNUM, BOUND, SCALE, SMIN,
     $                   XNORM, CNORM, TOFFNORM, TDIAGNORM, TEMP
      COMPLEX*16         SHIFT, CTEMP
      INTEGER            I, IB, IBEND, IMIN, IMAX,
     $                   J, JV, JB, JBEND, JOUT, JMAX, K,
     $                   NB, MB, DSIZE, IERR
*     ..
*     .. Local Arrays ..
      INTEGER            JLIST( MBMAX ), VDSIZES( MBMAX ),
     $                   TDSIZES( NBMAX )
      DOUBLE PRECISION   CNORMS( NBMAX ), VBOUNDS( MBMAX ),
     $                   VDIAGS( MBMAX, 2 ), VDIAG( 2, 2 ),
     $                   D( 2, 2 ), B( 2, 2 ), X( 2, 2 )
      COMPLEX*16         SHIFTS( MBMAX )
*     ..
*     .. External Functions ..
      LOGICAL            LSAME
      INTEGER            ILAENV
      DOUBLE PRECISION   DLAMCH, DASUM
      EXTERNAL           LSAME, ILAENV, DLAMCH, DASUM
*     ..
*     .. External Subroutines ..
      EXTERNAL           XERBLA, DLABAD, DAXPY, DSCAL, DLALN2, DGEMM,
     $                   DTREVC
*     ..
*     .. Intrinsic Functions ..
      INTRINSIC          ABS, DBLE, DCMPLX, DIMAG, CONJG, MIN, MAX,
     $                   MAXVAL
*     ..
*     .. Statement Functions ..
      DOUBLE PRECISION   SQ, CABS1
*     ..
*     .. Statement Function definitions ..
      SQ( TEMP ) = TEMP * TEMP
      CABS1( CTEMP ) = ABS( DBLE( CTEMP ) ) + ABS( DIMAG( CTEMP ) )
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
*     eigenvectors and standardize the array SELECT if necessary.
*
      IF( SOMEV ) THEN
         M = 0
         J = 1
         DO WHILE( J.LT.N )
            IF( T( J+1, J ).EQ.ZERO ) THEN
               IF( SELECT( J ) )
     $              M = M + 1
               J = J + 1
            ELSE
               IF( SELECT( J ).OR.SELECT( J+1 ) ) THEN
                  M = M + 2
                  SELECT( J ) = .TRUE.
                  SELECT( J+1 ) = .FALSE.
               END IF
               J = J + 2
            END IF
         END DO
         IF( J.EQ.N ) THEN
            IF( SELECT( J ) )
     $           M = M + 1
         END IF
      ELSE
         M = N
      END IF
*
      INFO = 0
      MB = ILAENV( 1, 'DTREVC', SIDE // HOWMNY, N, -1, -1, -1 )
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
      ELSE IF(LWORK.NE.-1) THEN
         IF( ( BACKTRANSFORM .AND. LWORK.LT.MAX( 1, 4*N ) )
     $       .OR. LWORK.LT.MAX( 1, 2*N ) ) THEN
            INFO = -14
         END IF
      END IF
      IF( INFO.NE.0 ) THEN
         CALL XERBLA( 'DTREVC3', -INFO )
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
*     Determine block size
*
      IF( BACKTRANSFORM ) THEN
         JMAX = MIN( LWORK / (2*N), MBMAX )
      ELSE
         JMAX = MIN( LWORK / N, MBMAX )
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
         IMIN = 1
         IMAX = 0
         JB = JMAX + 1
         JBEND = JMAX
         JOUT = M + 1
         JLIST = 0
         DO JV = N, 1, -1
*
*           --------------------------------------------------------
*           Add current eigenvector to workspace if needed.
*
            DSIZE = 1
            IF( JV.LT.N ) THEN
               IF( T( JV+1, JV ).NE.ZERO )
     $              DSIZE = 2
            END IF
            IF( JV.GT.1 ) THEN
               IF( T( JV, JV-1 ).NE.ZERO )
     $              DSIZE = 0
            END IF
            SELECTV = DSIZE.NE.0
            IF( SOMEV ) THEN
               SELECTV = SELECTV .AND. SELECT( JV )
            ENDIF
            IF( SELECTV ) THEN
               IMAX = MAX( IMAX, JV + DSIZE - 1 )
               JB = JB - DSIZE
               JLIST( JB ) = JV
               VDSIZES( JB ) = DSIZE
               SELECT CASE( DSIZE )               
               CASE( 1 )
*
*                 --------------------------------------------------
*                 Add real eigenvector to workspace.
*
                  SHIFT = DCMPLX( T( JV, JV ) )
                  X( 1, 1 ) = ONE
                  WORK( IMIN : JV - 1, JB ) = -T( IMIN : JV - 1, JV )
                  WORK( JV : IMAX, JB ) = ZERO
                  BOUND = DASUM( JV - IMIN, WORK( IMIN, JB ), 1 )
                  IF( BOUND.GE.OVFL ) THEN
                     SCALE = HALF * OVFL / BOUND
                     X( 1, 1 ) = SCALE * X( 1, 1 )
                     BOUND = HALF * OVFL
                     CALL DSCAL( JV - IMIN, SCALE, WORK( IMIN, JB ), 1 )
                  END IF
                  SHIFTS( JB ) = SHIFT
                  VDIAGS( JB , 1 ) = X( 1, 1 )
                  VBOUNDS( JB ) = BOUND
               CASE( 2 )
*     
*                 --------------------------------------------------
*                 Add complex eigenvector to workspace.
*
                  JLIST( JB+1 ) = 0
                  VDSIZES( JB+1 ) = 0
*
*                 Compute eigenvalues of 2x2 diagonal block
*                 Note: We assume 2x2 diagonal blocks of Schur matrices
*                 have complex eigenvalues/eigenvectors. DHSEQR
*                 guarantees that these blocks satisfy D(1,1)=D(2,2) and
*                 D(2,1)*D(1,2)<0, so we don't worry too much about
*                 numerical nastiness.
*
                  D = T( JV : JV+1, JV : JV+1 )
                  SHIFT = DCMPLX( HALF * ( D(1,1) + D(2,2) ) )
                  IF( D(1,1) .NE. D(2,2)
     $                .OR. D(2,1) * D(1,2) .GT. ZERO ) THEN
                     TEMP = -D(1,2) * D(2,1) - SQ( D(1,1) - D(2,2) ) / 4
                     TEMP = SQRT( MAX( TEMP, SMLNUM ) )
                  ELSE IF( D(2,1) .LT. ZERO ) THEN
                     TEMP = SQRT( -D(2,1) ) * SQRT( D(1,2) )
                  ELSE
                     TEMP = SQRT( D(2,1) ) * SQRT( -D(1,2) )
                  END IF
                  SHIFT = SHIFT + DCMPLX( ZERO, TEMP )
*
*                 Compute eigenvectors of 2x2 diagonal block
*                 Note: We apply a safe solve against a vector
*                 orthogonal to the range of D-SHIFT*I.
*     
                  SMIN = MAX( ULP * CABS1( SHIFT ), SMLNUM )
                  B( 1, 1 ) = -D( 2, 1 )
                  B( 2, 1 ) = D( 1, 1 ) - DBLE( SHIFT )
                  B( 1, 2 ) = ZERO
                  B( 2, 2 ) = DIMAG( SHIFT )
                  B = B / SUM( ABS( B ) )
                  CALL DLALN2( .FALSE., 2, 2, SMIN,
     $                         ONE, D, 2, ONE, ONE, B, 2,
     $                         DBLE( SHIFT ), DIMAG( SHIFT ), X, 2,
     $                         SCALE, XNORM, IERR )
                  X = X / XNORM
*
*                 Populate workspace and scale if needed
*                 
                  CALL DGEMM( 'N', 'N', JV - IMIN, 2, 2,
     $                        -ONE, T( IMIN, JV ), LDT, X, 2,
     $                        ZERO, WORK( IMIN, JB ), N )
                  WORK( JV : IMAX, JB : JB+1 ) = ZERO
                  BOUND = DASUM( JV - IMIN, WORK( IMIN, JB ), 1 )
     $                    + DASUM( JV - IMIN, WORK( IMIN, JB+1 ), 1 )
                  IF( BOUND.GE.OVFL ) THEN
                     SCALE = HALF * OVFL / BOUND
                     X = SCALE * X
                     BOUND = HALF * OVFL
                     CALL DSCAL( JV - IMIN, SCALE,
     $                           WORK( IMIN, JB ), 1 )
                     CALL DSCAL( JV - IMIN, SCALE,
     $                           WORK( IMIN, JB+1 ), 1 )
                  END IF
                  SHIFTS( JB ) = SHIFT
                  VDIAGS( JB:JB+1 , : ) = X
                  VBOUNDS( JB ) = BOUND
               END SELECT
            END IF
*
*           --------------------------------------------------------
*           Process workspace if full or if all eigenvectors are
*           found.
*
            IF( JB.LE.2 .OR. ( JV.EQ.1 .AND. JB.LE.JBEND ) ) THEN
               MB = JBEND - JB + 1
*
*              Blocked back substitution.
*              Note: We avoid splitting up 2x2 blocks.
*               
               IBEND = IMAX
               DO WHILE( IBEND.GT.0 )
                  IB = MAX( IBEND - NBMAX + 1, IMIN )
                  IF( IB.GT.1 ) THEN
                     IF( T( IB, IB-1 ).NE.ZERO )
     $                    IB = IB + 1
                  END IF
*
*                 Analyze diagonal block.
*                 Note: We determine whether each diagonal entry belongs
*                 to a 2x2 block, compute the 1-norms of the
*                 off-diagonal columns, and compute the entry-wise
*                 1-norm.
*     
                  TDIAGNORM = ZERO
                  DO I = IB, IBEND
                     DSIZE = 1
                     IF( I.LT.IBEND ) THEN
                        IF( T( I+1, I ).NE.ZERO )
     $                       DSIZE = 2
                     END IF
                     IF( I.GT.IB ) THEN
                        IF( T( I, I-1 ).NE.ZERO )
     $                       DSIZE = 0
                     END IF
                     SELECT CASE( DSIZE )
                     CASE( 1 )
                        BOUND = DASUM( I - IB, T( IB, I ), 1 )
                        TEMP = ABS( T( I, I ) )
                     CASE( 2 )
                        BOUND = DASUM( I - IB, T( IB, I ), 1 )
     $                          + DASUM( I - IB, T( IB, I+1 ), 1 )
                        TEMP = SUM( ABS( T( I:I+1, I:I+1 ) ) )
                     CASE( 0 )
                        BOUND = ZERO
                        TEMP = ZERO
                     END SELECT
                     TDSIZES( I - IB + 1 ) = DSIZE
                     CNORMS( I - IB + 1 ) = BOUND
                     TDIAGNORM = TDIAGNORM + BOUND + TEMP
                  END DO
*
*                 Analyze off-diagonal block.
*                 Note: Rough upper bound for infinity-norm.
*     
                  TOFFNORM = ZERO
                  DO I = IMIN, IB - 1, NBMAX
                     TEMP = ZERO
                     DO K = IB, IBEND
                        TEMP = TEMP + DASUM( MIN( IB - I, NBMAX ),
     $                                       T( I, K ), 1 )
                     END DO
                     TOFFNORM = MAX( TOFFNORM, TEMP )
                  END DO
*
*                 --------------------------------------------------
*                 Solve each RHS against diagonal block.
*                 Note: Back substitution is safeguarded to avoid
*                 overflow.
*     
                  DO J = JB, JBEND
                     NB = MIN( JLIST( J ) - 1, IBEND) - IB + 1
                     SHIFT = SHIFTS( J )
                     BOUND = VBOUNDS( J )
                     SMIN = MAX( ULP * CABS1( SHIFT ), SMLNUM )
                     IF( VDSIZES( J ).EQ.1 .AND. NB.GT.0 ) THEN
*
*                       Solve real RHS against diagonal block.
*
                        VDIAG( 1, 1 ) = VDIAGS( J, 1 )
                        DO K = IB + NB - 1, IB, -1
                           CNORM = CNORMS( K - IB + 1 )
                           SELECT CASE( TDSIZES( K - IB + 1 ) )
                           CASE( 1 )
*
*                             1x1 diagonal block with real RHS.
*
                              CALL DLALN2( .FALSE., 1, 1, SMIN,
     $                                     ONE, T( K, K ), LDT,
     $                                     ONE, ONE,
     $                                     WORK( K, J ), N,
     $                                     DBLE( SHIFT ),
     $                                     DIMAG( SHIFT ),
     $                                     X, 2, SCALE, XNORM, IERR )
                              TEMP = OVFL - BOUND
                              IF( XNORM.GE.ONE
     $                            .AND. CNORM.GT.TEMP/XNORM ) THEN
                                 TEMP = ( OVFL * HALF / CNORM ) / XNORM
                                 X( 1, 1 ) = TEMP * X( 1, 1 )
                                 SCALE = TEMP * SCALE
                              ELSE IF( XNORM.LT.ONE
     $                                 .AND. XNORM*BOUND.GT.TEMP ) THEN
                                 TEMP = OVFL * HALF / CNORM
                                 X( 1, 1 ) = TEMP * X( 1, 1 )
                                 SCALE = TEMP * SCALE
                              ELSE
                                 BOUND = XNORM * CNORM + BOUND
                              END IF
                              IF( SCALE.NE.ONE ) THEN
                                 CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                       WORK( IMIN, J ), 1 )
                                 VDIAG( 1, 1 ) = SCALE * VDIAG( 1, 1 )
                              END IF
                              WORK( K, J ) = X( 1, 1 )
                              CALL DAXPY( K - IB, -WORK( K, J ),
     $                                    T( IB, K ), 1,
     $                                    WORK( IB, J ), 1 )
                              IF( SCALE.NE.ONE )
     $                             BOUND = DASUM( K - IMIN,
     $                                            WORK( IMIN, J ), 1 )
                           CASE( 2 )
*
*                             2x2 diagonal block with real RHS.
*
                              CALL DLALN2( .FALSE., 2, 1, SMIN,
     $                                     ONE, T( K, K ), LDT,
     $                                     ONE, ONE,
     $                                     WORK( K, J ), N,
     $                                     DBLE( SHIFT ),
     $                                     DIMAG( SHIFT ),
     $                                     X, 2, SCALE, XNORM, IERR )
                              TEMP = OVFL - BOUND
                              IF( XNORM.GE.ONE
     $                            .AND. CNORM.GT.TEMP/XNORM ) THEN
                                 TEMP = ( OVFL * HALF / CNORM ) / XNORM
                                 X( : , 1 ) = TEMP * X( : , 1 )
                                 SCALE = TEMP * SCALE
                              ELSE IF( XNORM.LT.ONE
     $                                 .AND. XNORM*BOUND.GT.TEMP ) THEN
                                 TEMP = OVFL * HALF / CNORM
                                 X( : , 1 ) = TEMP * X( : , 1 )
                                 SCALE = TEMP * SCALE
                              ELSE
                                 BOUND = XNORM * CNORM + BOUND
                              END IF
                              IF( SCALE.NE.ONE ) THEN
                                 CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                       WORK( IMIN, J ), 1 )
                                 VDIAG( 1, 1 ) = SCALE * VDIAG( 1, 1 )
                              END IF
                              WORK( K:K+1, J ) = X( : , 1 )
                              CALL DGEMM( 'N', 'N', K - IB, 1, 2,
     $                                    NEGONE, T( IB, K ), LDT,
     $                                    WORK( K, J ), N,
     $                                    ONE, WORK( IB, J ), N )
                              IF( SCALE.NE.ONE )
     $                             BOUND = DASUM( K - IMIN,
     $                                            WORK( IMIN, J ), 1 )
                           END SELECT
                        END DO
*
*                       Rescale real RHS for back substitution
*     
                        SCALE = ONE
                        XNORM = DASUM( NB, WORK( IB, J ), 1 )
                        TEMP = OVFL - BOUND
                        IF( XNORM.GE.ONE
     $                      .AND. TOFFNORM.GT.TEMP/XNORM ) THEN
                           SCALE = ( OVFL * HALF / TOFFNORM ) / XNORM
                        ELSE IF( XNORM.LT.ONE
     $                           .AND. TOFFNORM*XNORM.GT.TEMP) THEN
                           SCALE = OVFL * HALF / TOFFNORM
                        END IF
                        IF( SCALE.NE.ONE ) THEN
                           CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                 WORK( IMIN, J ), 1 )
                           VDIAG( 1, 1 ) = SCALE * VDIAG( 1, 1 )
                           BOUND = DASUM( IB - IMIN,
     $                                    WORK( IMIN, J ), 1 )
                        END IF
                        VBOUNDS( J ) = BOUND
                        VDIAGS( J, 1 ) = VDIAG( 1, 1 )
                     ELSE IF( VDSIZES( J ).EQ.2 .AND. NB.GT.0 ) THEN
*
*                       Solve complex RHS against diagonal block.
*
                        VDIAG = VDIAGS( J:J+1, : )
                        DO K = IB + NB - 1, IB, -1
                           CNORM = CNORMS( K - IB + 1 )
                           SELECT CASE( TDSIZES( K - IB + 1 ) )
                           CASE( 1 )
*
*                             1x1 diagonal block with complex RHS.
*
                              CALL DLALN2( .FALSE., 1, 2, SMIN,
     $                                     ONE, T( K, K ), LDT,
     $                                     ONE, ONE,
     $                                     WORK( K, J ), N,
     $                                     DBLE( SHIFT ),
     $                                     DIMAG( SHIFT ),
     $                                     X, 2, SCALE, XNORM, IERR )
                              TEMP = OVFL - BOUND
                              IF( XNORM.GE.ONE
     $                            .AND. CNORM.GT.TEMP/XNORM ) THEN
                                 TEMP = ( OVFL * HALF / CNORM ) / XNORM
                                 X( 1, : ) = TEMP * X( 1, : )
                                 SCALE = TEMP * SCALE
                              ELSE IF( XNORM.LT.ONE
     $                                 .AND. XNORM*BOUND.GT.TEMP ) THEN
                                 TEMP = OVFL * HALF / CNORM
                                 X( 1, : ) = TEMP * X( 1, : )
                                 SCALE = TEMP * SCALE
                              ELSE
                                 BOUND = XNORM * CNORM + BOUND
                              END IF
                              IF( SCALE.NE.ONE ) THEN
                                 CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                       WORK( IMIN, J ), 1 )
                                 CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                       WORK( IMIN, J+1 ), 1 )
                                 VDIAG = SCALE * VDIAG
                              END IF
                              WORK( K, J:J+1 ) = X( 1, : )
                              CALL DGEMM( 'N', 'N', K - IB, 2, 1,
     $                                    NEGONE, T( IB, K ), LDT,
     $                                    WORK( K, J ), N,
     $                                    ONE, WORK( IB, J ), N )
                              IF( SCALE.NE.ONE )
     $                             BOUND = SUM( ABS( WORK( IMIN : K-1,
     $                                                     J : J+1 ) ) )
                           CASE( 2 )
*
*                             2x2 diagonal block with complex RHS.
*
                              CALL DLALN2( .FALSE., 2, 2, SMIN,
     $                                     ONE, T( K, K ), LDT,
     $                                     ONE, ONE,
     $                                     WORK( K, J ), N,
     $                                     DBLE( SHIFT ),
     $                                     DIMAG( SHIFT ),
     $                                     X, 2, SCALE, XNORM, IERR )
                              TEMP = OVFL - BOUND
                              IF( XNORM.GE.ONE
     $                            .AND. CNORM.GT.TEMP/XNORM ) THEN
                                 TEMP = ( OVFL * HALF / CNORM ) / XNORM
                                 X = TEMP * X
                                 SCALE = TEMP * SCALE
                              ELSE IF( XNORM.LT.ONE
     $                                 .AND. XNORM*BOUND.GT.TEMP ) THEN
                                 TEMP = OVFL * HALF / CNORM
                                 X = TEMP * X
                                 SCALE = TEMP * SCALE
                              ELSE
                                 BOUND = XNORM * CNORM + BOUND
                              END IF
                              IF( SCALE.NE.ONE ) THEN
                                 CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                       WORK( IMIN, J ), 1 )
                                 CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                       WORK( IMIN, J+1 ), 1 )
                                 VDIAG = SCALE * VDIAG
                              END IF
                              WORK( K:K+1, J:J+1 ) = X
                              CALL DGEMM( 'N', 'N', K - IB, 2, 2,
     $                                    NEGONE, T( IB, K ), LDT,
     $                                    WORK( K, J ), N,
     $                                    ONE, WORK( IB, J ), N )
                              IF( SCALE.NE.ONE )
     $                             BOUND = SUM( ABS( WORK( IMIN:K-1,
     $                                                     J:J+1 ) ) )
                           END SELECT
                        END DO
*
*                       Rescale complex RHS for back substitution.
*
                        SCALE = ONE
                        XNORM = SUM( ABS( WORK( IB:IB+NB-1,
     $                                          J:J+1 ) ))
                        TEMP = OVFL - BOUND
                        IF( XNORM.GE.ONE
     $                      .AND. TOFFNORM.GT.TEMP/XNORM ) THEN
                           SCALE = ( OVFL * HALF / TOFFNORM ) / XNORM
                        ELSE IF( XNORM.LT.ONE
     $                           .AND. TOFFNORM*XNORM.GT.TEMP) THEN
                           SCALE = OVFL * HALF / TOFFNORM
                        END IF
                        IF( SCALE.NE.ONE ) THEN
                           CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                 WORK( IMIN, J ), 1 )
                           CALL DSCAL( JLIST( J ) - IMIN, SCALE,
     $                                 WORK( IMIN, J+1 ), 1 )
                           VDIAG = SCALE * VDIAG
                           BOUND = SUM( ABS( WORK( IMIN:IB-1,
     $                                             J:J+1 ) ) )
                        END IF
                        VBOUNDS( J ) = BOUND
                        VDIAGS( J:J+1, : ) = VDIAG
                     END IF
                  END DO
*
*                 Back substitution with block of solution.
*
                  NB = IBEND - IB + 1
                  IF( IB.GT.IMIN )
     $                 CALL DGEMM( 'N', 'N', IB - IMIN, MB, NB,
     $                             NEGONE, T( IMIN, IB ), LDT,
     $                             WORK( IB, JB ), N,
     $                             ONE, WORK( IMIN, JB ), N )
                  IBEND = IB - 1
               END DO
*
*              Put diagonal blocks to get triangular eigenvectors.
*
               DO J = JB, JBEND
                  I = JLIST( J )
                  SELECT CASE( VDSIZES( J ) )
                  CASE( 1 )
                     WORK( I, J ) = VDIAGS( J, 1 )
                  CASE( 2 )
                     WORK( I:I+1, J:J+1 ) = VDIAGS( J:(J+1), : )
                  END SELECT
               END DO
*
*              -----------------------------------------------------
*              Copy results to output.
*                  
               JOUT = JOUT - MB
               IF( BACKTRANSFORM ) THEN
*                  
*                 Back transform with Schur vectors to get full
*                 eigenvectors.
*
                  CALL DGEMM( 'N', 'N', N, MB, IMAX - IMIN + 1,
     $                        ONE, VR( 1, IMIN ), LDVR,
     $                        WORK( IMIN, JB ), N,
     $                        ZERO, WORK( 1, JMAX + JB ), N )
                  VR( 1 : N, JOUT : JOUT + MB - 1 )
     $                 = WORK( 1 : N, JMAX + JB : JMAX + JBEND )
               ELSE
                  VR( 1 : IMIN - 1, JOUT : JOUT + MB - 1 ) = ZERO
                  VR( IMIN : IMAX , JOUT : JOUT + MB - 1 )
     $                 = WORK( IMIN : IMAX, JB : JBEND )
                  VR( IMAX + 1 : N, JOUT : JOUT + MB - 1 ) = ZERO
               END IF
*
*              Normalize eigenvectors.
*
               DO K = JB, JBEND
                  J = JOUT + K - JB
                  SELECT CASE( VDSIZES( K ) )
                  CASE( 1 )
                     SCALE = ONE / MAXVAL( ABS( VR( : , J ) ) )
                     CALL DSCAL( N, SCALE, VR( 1, J ), 1 )
                  CASE( 2 )
                     BOUND = SMLNUM
                     DO I = 1, N
                        BOUND = MAX( BOUND,
     $                               SUM( ABS( VR( I, J:J+1 ) ) ) )
                     END DO
                     SCALE = ONE / BOUND
                     CALL DSCAL( N, SCALE, VR( 1, J ), 1 )
                     CALL DSCAL( N, SCALE, VR( 1, J+1 ), 1 )
                  END SELECT
               END DO
*              
*              Workspace is now empty.
*
               IMIN = 1
               IMAX = 0
               JB = JMAX + 1
               JBEND = JMAX
            END IF
         END DO
      END IF

*
*     --------------------------------------------------------------
*     Compute left eigenvectors.
*     TODO
*     
      IF( LEFTV ) THEN
         CALL DTREVC( 'L', HOWMNY, SELECT, N,
     $                T, LDT, VL, LDVL, VR, LDVR, MM, M, WORK, INFO )
      END IF

      RETURN
*
*     End of DTREVC3
*
      END
