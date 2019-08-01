! -------------------------------------------------------------------------------------------------
! Utilities for MTcox R package
!
! Performs the Multi-task sparse Cox model solution path algorithm
! Author : Simon FONTAINE
! Last update : 2018-12-27
! -------------------------------------------------------------------------------------------------



! -------------------------------------------------------------------------------------------------
SUBROUTINE standardize(ntasks, ii, io, p, w, x, n, iex, &
                        wsum, xmean, xsd, flnullsd)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This subroutine standardizes the weights and the predictor matrix
! so that weights sum to 1 in each task
! and predictors have mean 0 and variance 1 within each task.
! Flsgs an error when variance is 0 in some task
! Also returns the initial means and sd to destandardize at the end
INTEGER             :: ntasks					!number of tasks
INTEGER             :: ii(ntasks)				!id of first observation per task
INTEGER             :: io(ntasks)				!id of last observation by task
INTEGER             :: p					    !number of variables
INTEGER             :: n                        !number of observations
DOUBLE PRECISION    :: w(n)		                !vector of weights
DOUBLE PRECISION    :: x(n,p)		            !matrix of variables (n rows, p columns)
INTEGER             :: iex(p)			        !indicator to exclude variable from algorithm (0=exclude)

DOUBLE PRECISION    :: wsum(ntasks)             !sum of weights
DOUBLE PRECISION    :: xmean(p, ntasks)         !mean of predictors
DOUBLE PRECISION    :: xsd(p, ntasks)           !sd of predictors
INTEGER             :: flnullsd (2)             !flag for predictor with no variation (j,k)

DOUBLE PRECISION, PARAMETER :: small = 1.0D-16  !small number of numerical checks
INTEGER             :: j					    !to cycle through variables
INTEGER             :: k					    !to cycle through tasks
DOUBLE PRECISION    :: tmp                      !general temporary variable
! - INITIALIZATION
    wsum = 0.0D0
    xmean = 0.0D0
    xsd = 0.0D0
    flnullsd = 0
    tmp = 0.0D0
! - STANDARDIZATION
    ! standardize weights to 1 across all tasks
    DO k=1,ntasks
        wsum(k) = sum(w(ii(k):io(k)))
        IF(wsum(k) == 0.0D0) THEN
            w(ii(k):io(k)) = 1.0D0/(io(k)-io(k)+1)
            wsum(k) = 1.0D0
        ENDIF
        w(ii(k):io(k)) = w(ii(k):io(k)) / wsum(k)
    ENDDO
    ! standardize variables wihtin each task, check if no variation
    DO k=1,ntasks
        DO j=1,p
            IF(iex(j) == 0) CYCLE
            xmean(j,k) = sum(x(ii(k):io(k),j) * w(ii(k):io(k)))
            x(ii(k):io(k),j) = x(ii(k):io(k),j) - xmean(j,k)
            tmp = sum(x(ii(k):io(k),j)*x(ii(k):io(k),j)*w(ii(k):io(k)))
            IF(abs(tmp)<small) THEN
                flnullsd(1) = j
                flnullsd(2) = k
                RETURN
            ENDIF
            xsd(j,k) = sqrt(tmp)
            x(ii(k):io(k),j) = x(ii(k):io(k),j) / xsd(j,k)
        ENDDO
    ENDDO
! -------------------------------------------------------------------------------------------------
END SUBROUTINE standardize
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE gradnull(ntasks, p, n, ns, iski, isko, nski, nsko, w, x, delta, d, grad)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This subroutine computes the gradient at the initial stage, when beta=0
! It is used only to compute the initial value of lambda to begin the solution path
! Varaibles are the ssame as in main routine
! - INPUT
INTEGER                     :: ntasks,p,n
INTEGER                     :: ns(ntasks)
INTEGER                     :: iski(sum(ns))
INTEGER                     :: isko(sum(ns))
INTEGER                     :: nski(ntasks)
INTEGER                     :: nsko(ntasks)
DOUBLE PRECISION            :: w(n)
DOUBLE PRECISION            :: x(n,p)
DOUBLE PRECISION            :: delta(n)
DOUBLE PRECISION            :: d(sum(ns))

! - LOCAL
INTEGER                     :: k
INTEGER                     :: j
INTEGER                     :: isk
DOUBLE PRECISION            :: SE(ntasks)
DOUBLE PRECISION            :: SR(p,ntasks)
DOUBLE PRECISION            :: IR(p,ntasks)
DOUBLE PRECISION            :: hhat

! - OUTPUT
DOUBLE PRECISION            :: grad(p,ntasks)
! - INITILIZATION
hhat = 0.0D0
grad = 0.0D0
SE = 0.0D0
SR = 0.0D0
IR = 0.0D0



    DO k=1,ntasks
        DO isk=nsko(k),nski(k),-1
        ! - UPDATE SE (eta=1 here)
            SE(k) = SE(k) + sum(w(iski(isk):isko(isk)))
        ! - COMPUTE HHAT
            hhat = 1.0D0/SE(k)
        ! - UPDATE SR
            DO j=1,p
                SR(j,k) = SR(j,k) + sum(x(iski(isk):isko(isk), j) * w(iski(isk):isko(isk)))
            ENDDO
        ! - COMPUTE IR
            DO j=1,p
                IR(j,k) = sum(x(iski(isk):isko(isk), j)* &
                                        w(iski(isk):isko(isk)) *&
                                         (1-delta(iski(isk):isko(isk))))
            ENDDO
        ! - UPDATE GRADIENT
            grad(:,k) = grad(:,k) + IR(:,k) - d(isk)*hhat*SR(:,k)
        ENDDO
    ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE gradnull
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE llk_saturated(ntasks, ns, nski, nsko, d, llk_sat)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function computes the likelihood of the saturated model
! - INPUT
INTEGER                         :: ntasks
INTEGER                         :: ns(ntasks)
INTEGER                         :: nski(ntasks)
INTEGER                         :: nsko(ntasks)
DOUBLE PRECISION                :: d(sum(ns))
! - LOCAL
INTEGER                         :: k
! - OUTPUT
DOUBLE PRECISION                :: llk_sat(ntasks)
llk_sat = 0.0D0
DO k=1,ntasks
    llk_sat(k) = -sum(d(nski(k):nsko(k)) * log(d(nski(k):nsko(k))))
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE llk_saturated
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE lambda_max(p, ntasks, reg, grad, pf, al, iex)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function computes value of lambda max such that it is the smallest with all variables excluded
! - INPUT
INTEGER                         :: p,ntasks,reg
DOUBLE PRECISION                :: grad(p,ntasks)
DOUBLE PRECISION                :: pf(p)
INTEGER                         :: iex(p)
! - LOCAL
INTEGER                         :: j
DOUBLE PRECISION                :: tmp
! - OUTPUT
DOUBLE PRECISION                :: al
al = 0.0D0
tmp = 0.0D0
DO j=1,p
    if(iex(j) ==0 ) CYCLE
    SELECT CASE (reg)
        CASE (0)
            tmp = maxval(abs(grad(j,:)))/pf(j)
        CASE (2)
            tmp = sqrt(sum(grad(j,:)**2))/pf(j)
    END SELECT
    al = max(al, tmp)
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE lambda_max
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE loglik(ntasks,n,p,ns,iski,isko,nski,nsko,w,x,delta,d,eta,llk)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function computes the log-likelihood of the model given parameters

! - INPUT
INTEGER                 :: ntasks,n,p
INTEGER                 :: ns(ntasks)
INTEGER                 :: iski(sum(ns)),isko(sum(ns))
INTEGER                 :: nski(ntasks),nsko(ntasks)
DOUBLE PRECISION        :: w(n),x(n,p)
DOUBLE PRECISION        :: delta(n)
DOUBLE PRECISION        :: d(sum(ns))
DOUBLE PRECISION        :: eta(n),mu(n)
! - LOCAL
DOUBLE PRECISION        :: SE, SI
INTEGER                 :: isk,i,k
INTEGER                 :: length                 ! for allocation
INTEGER, ALLOCATABLE    :: ind(:)                 ! will contain indices
! - OUTPUT
DOUBLE PRECISION        :: llk(ntasks)
!--------------------------------------------------------------------------------------------------
! - INITIALIZATION
!--------------------------------------------------------------------------------------------------
    llk = 0.0D0
    SE = 0.0D0
    SI = 0.0D0
    mu = exp(eta)
!--------------------------------------------------------------------------------------------------
! - COMPUTATION
!--------------------------------------------------------------------------------------------------
    DO k=1,ntasks
        SE = 0.0D0
        DO isk=nsko(k),nski(k),-1
            ! - VECTOR OF INDICES
            length = isko(isk) - iski(isk) + 1
            ALLOCATE(ind(length))
            ind = (/ (i,i=iski(isk),isko(isk),1) /)
            ! - UPDATE SE
                SE = SE + sum( w(ind)* mu(ind) )
            ! - COMPUTE SI
                SI = sum(eta(ind) * w(ind) * (1.0D0-delta(ind)) )
            ! - DEALLOCATE
            DEALLOCATE(ind)
            llk(k) = llk(k) + SI - d(isk) * log(SE)
        ENDDO
    ENDDO
! -------------------------------------------------------------------------------------------------
END SUBROUTINE loglik
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE gpg_cycle(ntasks,p,n,ii,io,iex,ns,iski,isko,nski,nsko,&
                                w,x,delta,d,beta,eta,lam,pf,reg,alpha,grad,sig,ierr,nupdates,llka)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function performs one cycle of proximal gradient descent on an active set
! either with stepsizes minorized by majorizing hessian

! - INPUT
INTEGER                     :: ntasks,p,n
INTEGER                     :: ii(ntasks)				!id of first observation per task
INTEGER                     :: io(ntasks)				!id of last observation by task
INTEGER                     :: iex(p)                   ! active variables (0=exluded), contains iex, sr and active set
INTEGER                     :: ns(ntasks)
INTEGER                     :: iski(sum(ns))
INTEGER                     :: isko(sum(ns))
INTEGER                     :: nski(ntasks)
INTEGER                     :: nsko(ntasks)
DOUBLE PRECISION            :: w(n)
DOUBLE PRECISION            :: x(n,p)
DOUBLE PRECISION            :: delta(n)
DOUBLE PRECISION            :: d(sum(ns))
DOUBLE PRECISION            :: beta(p,ntasks)
DOUBLE PRECISION            :: lam
DOUBLE PRECISION            :: pf(p)
INTEGER                     :: reg
DOUBLE PRECISION            :: alpha
DOUBLE PRECISION            :: eta(n)
DOUBLE PRECISION            :: sig(p)                      ! stepsize
! - LOCAL
DOUBLE PRECISION, PARAMETER :: small = 1.0D-16
INTEGER                     :: j,k
DOUBLE PRECISION            :: g(ntasks),h(ntasks)
DOUBLE PRECISION            :: beta_old(p,ntasks)
! - OUTPUT
INTEGER                     :: ierr(5)                  ! error flag and infos
DOUBLE PRECISION            :: grad(p,ntasks)           ! gradient
DOUBLE PRECISION            :: hess(p,ntasks)           ! hessian
INTEGER                     :: nupdates
DOUBLE PRECISION            :: llka(ntasks)
grad = 0.0D0
hess = 0.0D0
g = 0.0D0
h = 0.0D0
beta_old = beta
llka = 0.0D0
!--------------------------------------------------------------------------------------------------
! - THE CYCLE
!--------------------------------------------------------------------------------------------------
DO j=1,p
    if(iex(j) == 0 ) CYCLE
    !----------------------------------------------------------------------------------------------
    ! - SKIP INACTIVE
    !----------------------------------------------------------------------------------------------
    !IF (iex(j) == 1) CYCLE
    !----------------------------------------------------------------------------------------------
    ! - COMPUTE GRADIENT AND HESSIAN
    !----------------------------------------------------------------------------------------------
    call grad_hessj(ntasks,n,ns,iski,isko,nski,nsko,w,x(:,j),&
                    delta,d,eta,g,h)
    grad(j,:) = g
    hess(j,:) = h
    !----------------------------------------------------------------------------------------------
    ! - COMPUTE MAJORIZATION IF STEPSIZE NOT SUPPLIED
    !----------------------------------------------------------------------------------------------
    IF(sig(j) == 0.0D0) THEN
        sig(j) = maxval(h)
        IF(sig(j) < small) THEN
            ierr = (/1,4,0,j,0/)
            EXIT
        ENDIF
        sig(j) = 1.0D0 / sig(j)
    ENDIF
    !----------------------------------------------------------------------------------------------
    ! - UPDATE beta and mu
    !----------------------------------------------------------------------------------------------
    call prox(lam,pf(j),sig(j),alpha,reg,ntasks,beta(j,:),g)
    DO k=1,ntasks
        eta(ii(k):io(k)) = eta(ii(k):io(k)) + x(ii(k):io(k),j) * (beta(j,k) - beta_old(j,k))
    ENDDO
    nupdates = nupdates + 1
    !----------------------------------------------------------------------------------------------
    ! - UPDATE likelihood
    !----------------------------------------------------------------------------------------------
    call loglik(ntasks,n,p,ns,iski,isko,nski,nsko,w,x,delta,d,eta,llka)
!--------------------------------------------------------------------------------------------------
! - END CYCLE
!--------------------------------------------------------------------------------------------------
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE gpg_cycle
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE gpg_cycle_backtracking(ntasks,p,n,ii,io,iex,ns,iski,isko,nski,nsko,&
                                w,x,delta,d,beta,eta,lam,pf,reg,alpha,&
                                llka,frac,grad,hess,sig,ierr,nupdates)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function performs one cycle of proximal gradient descent on an active set
! either with stepsizes minorized by majorizing hessian

! - INPUT
INTEGER                     :: ntasks,p,n
INTEGER                     :: ii(ntasks)				!id of first observation per task
INTEGER                     :: io(ntasks)				!id of last observation by task
INTEGER                     :: iex(p)                   ! active variables (0=exluded), contains iex, sr and active set
INTEGER                     :: ns(ntasks)
INTEGER                     :: iski(sum(ns))
INTEGER                     :: isko(sum(ns))
INTEGER                     :: nski(ntasks)
INTEGER                     :: nsko(ntasks)
DOUBLE PRECISION            :: w(n)
DOUBLE PRECISION            :: x(n,p)
DOUBLE PRECISION            :: delta(n)
DOUBLE PRECISION            :: d(sum(ns))
DOUBLE PRECISION            :: beta(p,ntasks)
DOUBLE PRECISION            :: eta(n)
DOUBLE PRECISION            :: lam
DOUBLE PRECISION            :: pf(p)
INTEGER                     :: reg
DOUBLE PRECISION            :: alpha
DOUBLE PRECISION            :: llka(ntasks)             ! log-likelihood with actual beta
DOUBLE PRECISION            :: frac                     ! decreasing factor for line search
DOUBLE PRECISION            :: sig(p)                   ! stepsize
! - LOCAL
DOUBLE PRECISION, PARAMETER :: small = 1.0D-16
INTEGER                     :: j,k
DOUBLE PRECISION            :: g(ntasks),h(ntasks)
DOUBLE PRECISION            :: psg(ntasks)              ! pseudo gradient
DOUBLE PRECISION            :: llk(ntasks)              ! running log-likelihood
DOUBLE PRECISION            :: b(p,ntasks)              ! updated coefficient
DOUBLE PRECISION            :: cond
INTEGER                     :: flarmijo
DOUBLE PRECISION            :: eta_tmp(n)
INTEGER                     :: ntries
! - OUTPUT
INTEGER                     :: ierr(5)                  ! error flag and infos
DOUBLE PRECISION            :: grad(p,ntasks)           ! gradient
DOUBLE PRECISION            :: hess(p,ntasks)           ! hessian
INTEGER                     :: nupdates
! - INITIALIZATIONS
g = 0.0D0
h = 0.0D0
psg = 0.0D0
b = beta
cond = 0.0D0
eta_tmp = 0.0D0
grad = 0.0D0
hess = 0.0D0
call loglik(ntasks,n,p,ns,iski,isko,nski,nsko,w,x,delta,d,eta,llka)

DO j=1,p
    if(iex(j) == 0) CYCLE
    !--------------------------------------------------------------------------------------------------
    ! COMPUTE GRADIENT
    !--------------------------------------------------------------------------------------------------
    call grad_hessj(ntasks,n,ns,iski,isko,nski,nsko,w,x(:,j),&
                    delta,d,eta,g,h)
    grad(j,:) = g
    hess(j,:) = h
    !----------------------------------------------------------------------------------------------
    ! - COMPUTE INITIAL STEPSIZE
    ! - if not supplied (from previous cycle), then taken to be largest possible
    !----------------------------------------------------------------------------------------------
    IF(sig(j) == 0.0D0) THEN
        sig(j) = minval(h)
        IF(sig(j) < small) THEN
            ierr = (/1,4,0,j,0/)
            EXIT
        ENDIF
        sig(j) = 1.0D0 / sig(j)
    ENDIF
    !--------------------------------------------------------------------------------------------------
    ! BACKTRACKING
    !--------------------------------------------------------------------------------------------------
    flarmijo = 0
    ntries = 0
    DO WHILE (flarmijo == 0 .AND. ntries <= 20)
        eta_tmp = eta
        nupdates = nupdates + 1
        !----------------------------------------------------------------------------------------------
        ! COMPUTE UPDATE WITH PROXIMAL
        !----------------------------------------------------------------------------------------------
        b(j,:) = beta(j,:)
        call prox(lam,pf(j),sig(j),alpha,reg,ntasks,b(j,:),g)
        DO k=1,ntasks
            eta_tmp(ii(k):io(k)) = eta_tmp(ii(k):io(k)) +  x(ii(k):io(k),j) * (b(j,k) - beta(j,k))
        ENDDO
        !----------------------------------------------------------------------------------------------
        ! PSEUDO GRADIENT
        !----------------------------------------------------------------------------------------------
        psg = (b(j,:) - beta(j,:))/sig(j)
        !----------------------------------------------------------------------------------------------
        ! NEW LOG LIKELIHOOD
        !----------------------------------------------------------------------------------------------
        call loglik(ntasks,n,p,ns,iski,isko,nski,nsko,w,x,delta,d,eta_tmp,llk)
        !----------------------------------------------------------------------------------------------
        ! CHECL ARMIJO CONDITION
        !----------------------------------------------------------------------------------------------
        cond = sum(llka) - sig(j) * sum(g * psg) + sig(j) * 0.5D0 * sum(psg * psg)
        IF (cond + 1.0D-3 < sum(llk)) THEN
            ! Condition fail, reduce stepsize
            sig(j) = sig(j) * frac
        ELSE
            ! condition passes, so we exit the loop
            flarmijo = 1
        ENDIF
    ENDDO
    !--------------------------------------------------------------------------------------------------
    ! SAVE LOGLIKELIHOOD AND ESTIMATE
    !--------------------------------------------------------------------------------------------------
    beta(j,:) = b(j,:)
    eta = eta_tmp
    llka = llk
ENDDO
! -------------------------------------------------------------------------------------------------
END SUBROUTINE gpg_cycle_backtracking
! -------------------------------------------------------------------------------------------------



! -------------------------------------------------------------------------------------------------
SUBROUTINE gpg_descent(ntasks,p,n,ii,io,iex,ns,iski,isko,nski,nsko,&
                                w,x,delta,d,beta,eta,lam,pf,reg,alpha,&
                                llkb,ierr,alg,eps,frac,ncycles,nupdates,maxit)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function performs proximal gradient descent on an active set
! either with stepsizes minorized by majorizing hessian

! - INPUT
INTEGER                     :: ntasks,p,n
INTEGER                     :: ii(ntasks)				!id of first observation per task
INTEGER                     :: io(ntasks)				!id of last observation by task
INTEGER                     :: iex(p)                   !active variables (0=exluded), contains iex, sr
INTEGER                     :: ns(ntasks)
INTEGER                     :: iski(sum(ns))
INTEGER                     :: isko(sum(ns))
INTEGER                     :: nski(ntasks)
INTEGER                     :: nsko(ntasks)
DOUBLE PRECISION            :: w(n)
DOUBLE PRECISION            :: x(n,p)
DOUBLE PRECISION            :: delta(n)
DOUBLE PRECISION            :: d(sum(ns))
DOUBLE PRECISION            :: beta(p,ntasks)
DOUBLE PRECISION            :: eta(n)          !
DOUBLE PRECISION            :: lam
DOUBLE PRECISION            :: pf(p)
INTEGER                     :: reg
INTEGER                     :: maxit
DOUBLE PRECISION            :: alpha
DOUBLE PRECISION            :: llkb(ntasks)             ! log-likelihood within
DOUBLE PRECISION            :: eps                      ! convergence threshold
INTEGER                     :: alg                      ! which algorithm to use (1=majorized, 2=backtracking)
DOUBLE PRECISION            :: frac
! - LOCAL
DOUBLE PRECISION, PARAMETER :: small = 1.0D-16
INTEGER                     :: j
DOUBLE PRECISION            :: sig(p)                   ! stepsize
DOUBLE PRECISION            :: b(p,ntasks)              ! updated coefficient
INTEGER                     :: conv
DOUBLE PRECISION            :: grad(p,ntasks), hess(p,ntasks)
DOUBLE PRECISION            :: llka(ntasks)             ! llk outside
INTEGER                     :: ncycles_local
INTEGER                     :: iactive(p)
! - OUTPUT
INTEGER                     :: ierr(5)                  ! error flag and infos
INTEGER                     :: ncycles
INTEGER                     :: nupdates
! - INITIALIZATIONS
iactive = iex
conv = 0
ncycles_local = 0
b = 0.0D0
grad = 0.0D0
hess = 0.0D0
llka = 0.0D0


DO WHILE (conv == 0)
    sig = 0.0D0
    ncycles_local = ncycles_local + 1
    ncycles = ncycles + 1
    IF(ncycles > maxit) THEN
        ierr = (/ 1,6,0,0,0 /)
        EXIT
    ENDIF
    ! Save previous estimates in current variable to update
    b = beta
    llka = llkb
    SELECT CASE (alg)
        CASE (1)
            call gpg_cycle(ntasks,p,n,ii,io,iactive,ns,iski,isko,nski,nsko,&
                                w,x,delta,d,b,eta,lam,pf,reg,alpha,grad,sig,ierr,nupdates,llka)
            !call loglik(ntasks,n,p,ns,iski,isko,nski,nsko,w,x,delta,eta,d,llka)
        CASE (2)
            call gpg_cycle_backtracking(ntasks,p,n,ii,io,iactive,ns,iski,isko,nski,nsko,&
                                w,x,delta,d,b,eta,lam,pf,reg,alpha,llka,frac,grad,hess,sig,ierr,nupdates)
    END SELECT
    ! Check for fatal errors
    IF(ierr(1) == 1) EXIT
    ! Now b containes the new estimates and llka the new likelihood
    ! Compute difference in estimates and llk
    IF(maxval(abs(b-beta)) < eps)THEN
        conv = 1
    ENDIF
    beta = b
    llkb = llka
    IF(ncycles_local >= 1000) THEN
        ierr = (/ 1,5,0,0,0 /)
        EXIT
    ENDIF
    ! Update the active set
    IF(ncycles_local == 1)THEN
        DO j=1,p
            IF(maxval(abs(beta(j,:))) < small) iactive(j) = 0
        ENDDO
    ENDIF
ENDDO

! -------------------------------------------------------------------------------------------------
END SUBROUTINE gpg_descent
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE grad_hessj(ntasks,n,ns,iski,isko,nski,nsko,w,x,delta,d,eta,grad,hess)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function produces the dradient and hessian matrix (diagonal)

! - INPUT
INTEGER                     :: ntasks,n
INTEGER                     :: ns(ntasks)
INTEGER                     :: iski(sum(ns))
INTEGER                     :: isko(sum(ns))
INTEGER                     :: nski(ntasks)
INTEGER                     :: nsko(ntasks)
DOUBLE PRECISION            :: w(n)
DOUBLE PRECISION            :: x(n)
DOUBLE PRECISION            :: delta(n)
DOUBLE PRECISION            :: d(sum(ns))
DOUBLE PRECISION            :: eta(n)
! - LOCAL
INTEGER                     :: i,k
INTEGER                     :: isk
INTEGER                     :: length                 ! for allocation
INTEGER, ALLOCATABLE        :: ind(:)                 ! will contain indices
DOUBLE PRECISION            :: mu(n)
DOUBLE PRECISION            :: SE
DOUBLE PRECISION            :: SR
DOUBLE PRECISION            :: SR2
DOUBLE PRECISION            :: SI
DOUBLE PRECISION            :: hhat
! - OUTPUT
DOUBLE PRECISION            :: grad(ntasks)         ! gradient
DOUBLE PRECISION            :: hess(ntasks)         ! hessian
!--------------------------------------------------------------------------------------------------
! - INITIALIZATION
!--------------------------------------------------------------------------------------------------
grad = 0.0D0
hess = 0.0D0
mu = exp(eta)
!--------------------------------------------------------------------------------------------------
! - ALGORITHM
!--------------------------------------------------------------------------------------------------
DO k=1,ntasks
    SE = 0.0D0
    SR = 0.0D0
    SR2 = 0.0D0
    SI = 0.0D0
    hhat = 0.0D0
    DO isk=nsko(k),nski(k),-1
        ! - VECTOR OF INDICES
        length = isko(isk) - iski(isk) + 1
        ALLOCATE(ind(length))
        ind = (/ (i,i=iski(isk),isko(isk),1) /)
        ! - UPDATE SE
            SE = SE + sum( w(ind)* mu(ind) )
        ! - COMPUTE HHAT
            hhat = 1.0D0/SE
        ! - UPDATE SR and SR2
            SR = SR + sum(w(ind) * mu(ind) * x(ind) )
            SR2 = SR2 + sum(w(ind) * mu(ind) * x(ind)**2.0D0 )
        ! - COMPUTE SI
            SI = sum(x(ind) * w(ind) * (1.0D0-delta(ind)) )
        ! - UPDATE GRADIENT
            grad(k) = grad(k) - SI + d(isk)*hhat*SR
        ! - UPDATE HESSIAN
            hess(k) = hess(k) - (d(isk)*hhat*SR)**2.0D0 + d(isk)*hhat*SR2
        ! - DEALLOCATE
        DEALLOCATE(ind)
    ENDDO
ENDDO
! -------------------------------------------------------------------------------------------------
END SUBROUTINE grad_hessj
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE hhat_compute(ntasks,n,ns,iski,isko,nski,nsko,w,eta,hhat)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function produces the dradient and hessian matrix (diagonal)

! - INPUT
INTEGER                     :: ntasks,n
INTEGER                     :: ns(ntasks)
INTEGER                     :: iski(sum(ns))
INTEGER                     :: isko(sum(ns))
INTEGER                     :: nski(ntasks)
INTEGER                     :: nsko(ntasks)
DOUBLE PRECISION            :: w(n)
DOUBLE PRECISION            :: eta(n)
! - LOCAL
INTEGER                     :: i,k
INTEGER                     :: isk
INTEGER                     :: length                 ! for allocation
INTEGER, ALLOCATABLE        :: ind(:)                 ! will contain indices
DOUBLE PRECISION            :: mu(n)
DOUBLE PRECISION            :: SE
DOUBLE PRECISION            :: hhat(sum(ns))
!--------------------------------------------------------------------------------------------------
! - INITIALIZATION
!--------------------------------------------------------------------------------------------------
mu = exp(eta)
hhat = 0.0D0
!--------------------------------------------------------------------------------------------------
! - ALGORITHM
!--------------------------------------------------------------------------------------------------
DO k=1,ntasks
    SE = 0.0D0
    DO isk=nsko(k),nski(k),-1
        ! - VECTOR OF INDICES
        length = isko(isk) - iski(isk) + 1
        ALLOCATE(ind(length))
        ind = (/ (i,i=iski(isk),isko(isk),1) /)
        ! - UPDATE SE
            SE = SE + sum( w(ind)* mu(ind) )
        ! - COMPUTE HHAT
            hhat(isk) = 1.0D0/SE
        ! - DEALLOCATE
        DEALLOCATE(ind)
    ENDDO
ENDDO
! -------------------------------------------------------------------------------------------------
END SUBROUTINE hhat_compute
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE prox(lam, pf, sig, alpha, reg, ntasks, beta, grad)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! This function produces the proximal operator

! - INPUT
DOUBLE PRECISION        :: lam              ! lambda
DOUBLE PRECISION        :: pf               ! penalty factor
DOUBLE PRECISION        :: sig              ! stepsize
DOUBLE PRECISION        :: alpha            ! mixing of regularization
INTEGER                 :: reg              ! 0=infty of 2=2
INTEGER                 :: ntasks           ! dimension
DOUBLE PRECISION        :: beta(ntasks)     ! previous estimate to update
DOUBLE PRECISION        :: grad(ntasks)     ! gradient
! - LOCAL
DOUBLE PRECISION        :: s(ntasks)        ! for soft-thresholding
DOUBLE PRECISION        :: sgn(ntasks)      ! sign of s
INTEGER                 :: k                ! counter
DOUBLE PRECISION        :: proj(ntasks)     ! for l1 projection
DOUBLE PRECISION        :: tmp, norm        ! for l2 projection
! - OUTPUT
INTEGER                 :: ierr(5)          ! error flag and infos
!--------------------------------------------------------------------------------------------------
! - SOFT-THRESHOLDING (l1)
!--------------------------------------------------------------------------------------------------
s = beta - sig*grad
sgn = 0.0D0
proj = 0.0D0
tmp = 0.0D0
norm = 0.0D0
!IF (alpha > 0.0D0) THEN
    ! - RETRIEVE SIGN
    DO k=1,ntasks
        IF(s(k)>0.0D0) sgn(k) = 1.0D0
        IF(s(k)<0.0D0) sgn(k) = -1.0D0
    ENDDO
    ! - SOFT THRESHOLDING
    s = abs(s) - alpha*lam*sig*pf
    DO k=1,ntasks
        IF(s(k)<0.0D0) s(k) = 0.0D0
    ENDDO
    s=s*sgn
!ENDIF
!--------------------------------------------------------------------------------------------------
! - REGULARIZATION (lq)
!--------------------------------------------------------------------------------------------------
SELECT CASE (reg)
    !--------------------------------------------------------------------------------------------------
    ! - q=infty
    !--------------------------------------------------------------------------------------------------
    CASE (0)
        ! - SKIP ACCORDING TO RULE
        IF (sum(abs(s)) < lam * pf * sig) THEN
            beta = 0.0D0
        ! - ELSE DO THE PROJECTION
        ELSE
            call ProjB1Mich(s, ntasks, lam*pf*(1.0D0 - alpha)*sig, proj)
            beta = s - proj
        ENDIF
    !--------------------------------------------------------------------------------------------------
    ! - q=2
    !--------------------------------------------------------------------------------------------------
    CASE (2)
        norm = sqrt(dot_product(s,s))
        tmp = norm - pf * lam * (1.0D0 - alpha) * sig
        IF (tmp > 0.0D0) THEN
            beta = s * tmp / norm
        ELSE
            beta = 0.0D0
        ENDIF
END SELECT
! -------------------------------------------------------------------------------------------------
END SUBROUTINE prox
! -------------------------------------------------------------------------------------------------


! -------------------------------------------------------------------------------------------------
SUBROUTINE ProjB1(v, p, z, w)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE

! - VARIABLES DECLARATIONS -------------------------
    ! - INPUTS -
    INTEGER :: p
    DOUBLE PRECISION :: v(p)
    DOUBLE PRECISION :: z
    ! - OUTPUTS -
    DOUBLE PRECISION :: w(p)
    ! - LOCAL VARS -
    DOUBLE PRECISION :: vp(p)
    DOUBLE PRECISION :: r
    DOUBLE PRECISION :: s,ds
    DOUBLE PRECISION :: theta
    INTEGER :: U(p),G(p),L(p)
    INTEGER :: rho, drho
    INTEGER :: i,j,k,m,fl,cl
    DOUBLE PRECISION :: sgn(p)
    ! - rng tests fix seed
       integer, allocatable :: seed(:)
       integer :: size
       call random_seed(size=size)
       allocate(seed(size))
       call random_seed(put=seed)
    w = 0.0D0
    vp = 0.0D0
    r = 0.0D0
    s = 0.0D0
    ds = 0.0D0
    theta = 0.0D0

    ! - absolute components -
    vp = abs(v)
    ! - unit ball test -
    !z = 1.0D0
    ! - if already within the ball
    IF(sum(vp) <= z) THEN
        w = v
    ! - otherwise we have to project -
    ELSE
    ! - retrieve sgn -
        DO j=1,p
            IF(v(j) == 0.0D0) THEN
                sgn(j) = 0.0D0
            ELSEIF(v(j) > 0.0D0)THEN
                sgn(j) = 1.0D0
            ELSE
                sgn(j) = -1.0D0
            ENDIF
        ENDDO
    ! - initializations -
        U = 1
        s = 0.0D0
        rho = 0.0D0
    ! - while loop -
        DO
        ! - produce random index -
            !CALL RANDOM_SEED()
            CALL RANDOM_NUMBER(r)
            cl = CEILING(sum(U) * r)
            i = 0
            DO j=1,p
                IF (U(j) == 1) i = i + 1
                IF(i == cl) THEN
                    k = j
                    i = i+1
                    EXIT
                ENDIF
            ENDDO
        ! - partition -
            G = 0
            L = 0
            DO j=1,p
                IF(U(j) == 0) CYCLE
                IF(vp(j) >= vp(k)) THEN
                    G(j) = 1
                ELSE
                    L(j) = 1
                ENDIF
            ENDDO
        ! - update rule -
            drho = sum(G)
            ds = sum(vp * G)
            IF ((s+ds)-(rho+drho)*vp(k) < z ) THEN
                s = s + ds
                rho = rho + drho
                U = L
            ELSE
                U = G
                U(k) = 0
            ENDIF
        ! - while A non empty -
            IF (sum(U)==0) EXIT
        ENDDO
        ! - output -
        theta = (s-z) / rho
        DO j=1,p
            w(j) = max(vp(j) - theta, 0.0D0)*sgn(j)
        ENDDO
    ENDIF
    z=theta
! -------------------------------------------------------------------------------------------------
END SUBROUTINE ProjB1
! -------------------------------------------------------------------------------------------------




! -------------------------------------------------------------------------------------------------
SUBROUTINE ProjB1Mich(v, p, z, w)
! -------------------------------------------------------------------------------------------------
    
! - VARIABLES DECLARATIONS -------------------------
    ! - INPUTS -
    INTEGER :: p
    DOUBLE PRECISION :: v(p)
    DOUBLE PRECISION :: z
    ! - OUTPUTS -
    DOUBLE PRECISION :: w(p)
    ! - LOCAL VARS -
    DOUBLE PRECISION :: vp(p)
    DOUBLE PRECISION :: tau, rho
    INTEGER :: j,ch
    INTEGER :: sgn(p)
    INTEGER :: A(p)
! - END VARIABLE DECLARATIONS ----------------------

    ! - absolute components -
    vp = abs(v)
    ! - if already within the ball
    IF(sum(vp) <= z) THEN
        w = v
    ! - otherwise we have to project -
    ELSE
        ! - retrieve sgn -
        DO j=1,p
            IF(v(j) == 0.0D0) THEN
                sgn(j) = 0
            ELSEIF(v(j) > 0.0D0)THEN
                sgn(j) = 1
            ELSE
                sgn(j) = -1
            ENDIF
        ENDDO
        ! - cycle -
        A = 1
        rho = (sum(vp) - z) / sum(A)
        DO
            ch = 0
            DO j=1,p
                IF ( A(j) == 0 ) CYCLE
                IF (vp(j) <= rho) THEN
                    A(j) = 0
                    ch = 1
                ENDIF
            ENDDO
            rho = (sum(vp*A) - z) / sum(A)
            IF ( ch == 0 ) EXIT
        ENDDO
        tau = rho
        ! - produce projection on simplex and ajust sign -
            DO j=1,p
                w(j) = max(vp(j) - tau, 0.0D0)*sgn(j)
            ENDDO
        ENDIF
    z=tau
! -------------------------------------------------------------------------------------------------
END SUBROUTINE ProjB1Mich
! -------------------------------------------------------------------------------------------------

    
    
    