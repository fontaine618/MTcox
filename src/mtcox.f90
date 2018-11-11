! -------------------------------------------------------------------------------------------------
! MTcox subroutine for MTcox R package
!
! Performs the Multi-task sparse Cox model solution path algorithm
! Author : Simon FONTAINE
! Last update : 2018-11-11
! -------------------------------------------------------------------------------------------------

! -------------------------------------------------------------------------------------------------
SUBROUTINE mtcox_solutionpath(ntasks,ii,io,p,ns,delta,w,x,d,iski,isko,nski,nsko,n,&
                                optInt,optDbl,iex,pf,nlam,lam,hhat,llk_path,dev,pdev_path, &
                                eta_path,null_dev,ierr,beta,betanorm,kkt,kkt0,&
                                ncycles,nupdates,nbeta,entered,llk_null,llk_sat,grad,alf)
! -------------------------------------------------------------------------------------------------
IMPLICIT NONE
! -------------------------------------------------------------------------------------------------
! INPUT DESCRIPTION
! -------------------------------------------------------------------------------------------------
! - DIMENSIONS
INTEGER             :: ntasks					!number of tasks
INTEGER             :: ii(ntasks)				!id of first observation per task
INTEGER             :: io(ntasks)				!id of last observation by task
INTEGER             :: p					    !number of variables
INTEGER             :: ns(ntasks)				!number of different failure time (=S^k)
INTEGER             :: iski(sum(ns))	        !address of group first observation
INTEGER             :: isko(sum(ns))	        !address of group last observation
INTEGER             :: nski(ntasks)	            !identifies indeces for isk start
INTEGER             :: nsko(ntasks)	            !identifies indeces for isk last
INTEGER             :: n                        !number of observations
! - INPUTS
INTEGER             :: delta(n)				    !vector indicator of failure(0) or censored(1)
DOUBLE PRECISION    :: w(n)		                !vector of weights
DOUBLE PRECISION    :: x(n,p)		            !matrix of variables (n rows, p columns)
! - ALGORITHM DETAILS

INTEGER             :: nlam         		    !number of lambda iterations
DOUBLE PRECISION    :: lam(nlam)                !actual lambda values (or user supplied values)

INTEGER             :: iex(p)			        !indicator to exclude variable from algorithm (0=exclude)
DOUBLE PRECISION    :: pf(p)                    !penalty factor by variable

INTEGER             :: optInt(7)                !contains the diffrents options which are integer

INTEGER             :: str          			!indicator to perform strong rule or not
INTEGER             :: maxit        		    !maximum number of updates
INTEGER             :: earlystop        	    !whether or not to stop early if condition is met
INTEGER             :: alg                      !=0 : IRLS, =1 : majorized, =2 : backtracking
INTEGER             :: dfmax        		    !maximum number of variables in the model (all non-zero at one time)
INTEGER             :: pmax         		    !maximum number of variables to entre the model (includes those back to zero)
INTEGER             :: reg          			!regularization 0 (q=infty) or 2 (q=2)


DOUBLE PRECISION    :: optDbl(4)                !contains the diffrents options which are doubles

DOUBLE PRECISION    :: eps                      !convergence threshold
DOUBLE PRECISION    :: frac            	        !fraction for backtracing
DOUBLE PRECISION    :: alpha                    !parameter for sparse model
DOUBLE PRECISION    :: lamfrac                  !fraction of lambda_max to compute lambda_min (if =1, then user supplied)



! -------------------------------------------------------------------------------------------------
! OUTPUT DESCRIPTION
! -------------------------------------------------------------------------------------------------
INTEGER             :: nalam                    !actual number of lambda iterations performed
DOUBLE PRECISION    :: beta(p, ntasks, nlam)    !estimates along solution path
DOUBLE PRECISION    :: betanorm(p, nlam)   !norm estimates along solution path before destandardization
DOUBLE PRECISION    :: kkt(p, nlam)             !kkt condition for non zero variable along solution path
DOUBLE PRECISION    :: kkt0(p, nlam)            !kkt condition for excluded variable along solution path
DOUBLE PRECISION    :: llk_null(ntasks)         !null log-likelihood
DOUBLE PRECISION    :: llk_sat(ntasks)          !log-likelihood of saturated model
DOUBLE PRECISION    :: llk_path(ntasks,nlam)    !log-likelihood along solution path
DOUBLE PRECISION    :: pdev_path(nlam)          !percent deviance explained along solution path
DOUBLE PRECISION    :: dev(nlam)                !deviance along solution path
DOUBLE PRECISION    :: hhat(sum(ns), nlam)      !nonparametric baseline estimate along solution path
DOUBLE PRECISION    :: eta_path(n,nlam)         !linear predictor along the solution path
DOUBLE PRECISION    :: null_dev                 !null deviance
INTEGER             :: ncycles                  !number of cycles performed per lambda
INTEGER             :: nupdates                 !number of updates performed per lambda
INTEGER             :: nbeta(nlam)              !number of variables in model
INTEGER             :: entered(p)
INTEGER             :: nbetaever                !number of variables ever to enter the model
INTEGER             :: idvars(p)                !sequence in which the variables enter the model
INTEGER             :: ierr(5)                  !error informations
                                                !   1: either fatal errors(1) or warnings(0)
                                                !   2: error code
                                                !       0 0: no errors
                                                !         1: no variation in some variable
                                                !         2: no positive penalty factor
                                                !         3: negative weights
                                                !         4: max val of hessian is numerically 0
                                                !         5: no convergence in gpg_descent (>100 cycles)
                                                !         6: maximum number of cycles reached in gpg_descent
                                                !         7: lambda max too small
                                                !         8: alf too big
                                                !         9:
                                                !         10:
                                                !
                                                !       1 0: no warnings
                                                !         1: armijo convergence fails after maxcycle
                                                !         2:
                                                !   3: lambda iteration
                                                !   4: variable id
                                                !   5: task

! -------------------------------------------------------------------------------------------------
! LOCAL VARIABLES AND PARAMETERS DESCRIPTION
! -------------------------------------------------------------------------------------------------
!- PARAMETERS
DOUBLE PRECISION, PARAMETER :: small = 1.0D-16  !small number of numerical checks
DOUBLE PRECISION, PARAMETER :: big = 9.9D30     !big number of numerical checks
! - COUNTERS
INTEGER             :: i     				    !to cycle through observations
INTEGER             :: isk   				    !to remember the position of first observation for a s^k
INTEGER             :: j					    !to cycle through variables
INTEGER             :: k					    !to cycle through tasks
INTEGER             :: l					    !to cycle through lambda

! - TEMPORARY STORING
DOUBLE PRECISION    :: d(sum(ns))	            !value of d_s^k
DOUBLE PRECISION    :: eta(n)                   !current linear predictor
DOUBLE PRECISION    :: b(p, ntasks)             !current estimate
DOUBLE PRECISION    :: llk(ntasks)              !current log likelihood

DOUBLE PRECISION    :: wsum(ntasks)             !sum of weights
DOUBLE PRECISION    :: xmean(p, ntasks)         !mean of predictors
DOUBLE PRECISION    :: xsd(p, ntasks)           !sd of predictors
INTEGER             :: flnullsd (2)             !flag for predictor with no variation (j,k)
INTEGER             :: isr (p)                  !strong rule set (1= in, 0=out)

DOUBLE PRECISION    :: al                       !current lambda value
DOUBLE PRECISION    :: al0                      !previous lambda value
DOUBLE PRECISION    :: alf                      !lambda multiplication factor

INTEGER             :: flearlystop              !flag for early stop criteria (1 = will stop)
INTEGER             :: fl                       !flag general
DOUBLE PRECISION    :: grad(p,ntasks)           !gradient
DOUBLE PRECISION    :: g(ntasks)                !gradient
DOUBLE PRECISION    :: h(ntasks)                !hessian matrix
DOUBLE PRECISION    :: tmp                      !temporary storing
DOUBLE PRECISION    :: hhata(sum(ns))           !estimate of the nonparametric baselinetemporary

str=optInt(1)
maxit=optInt(2)
earlystop=optInt(3)
alg=optInt(4)
dfmax=optInt(5)
pmax=optInt(6)
reg=optInt(7)

eps=optDbl(1)
frac=optDbl(2)
alpha=optDbl(3)
lamfrac=optDbl(4)
! -------------------------------------------------------------------------------------------------
! INITIALIZATION
! -------------------------------------------------------------------------------------------------
! - INITIAL CHECKS (PUT IN R PREPARATION)
    ! need at least one positive penalty factor
    IF(maxval(pf) <= 0.0D0) THEN
        ierr = (/1,2,0,0,0/)
        RETURN
    ENDIF
    ! negative pf are set to zero
    pf=max(0.0D0,pf)
    ! check positive weights
    IF (minval(w) < 0) THEN
        ierr = (/1,3,0,0,0/)
        RETURN
    ENDIF

! - INITIALIZE VARIABLES
    grad = 0.0D0
    g = 0.0D0
    h = 0.0D0
    tmp = 0.0D0
    b = 0.0D0
    hhata = 0.0D0
    eta = 0.0D0
    llk_null = 0.0D0
    llk_sat = 0.0D0
    al0 = big
    flearlystop = 0
    flnullsd = 0
    nupdates = 0
    ncycles = 0
    al = big
    nbeta = 0
    nbetaever = 0
    kkt = 0.0D0
    kkt0 = 0.0D0
! - STANDARDIZATION
    call standardize(ntasks, ii, io, p, w, x, n, iex, wsum, xmean, xsd, flnullsd)
    IF(flnullsd(1)>0) THEN
        ierr = (/1,1,0,flnullsd(1), flnullsd(2)/)
        RETURN
    ENDIF

! ---------------------------------------------------------------------------------------------
! NULL MODEL
! ---------------------------------------------------------------------------------------------
    ! - COMPUTE LIKELIHOOD (NULL MODEL)
        call loglik(ntasks,n,p,ns,iski,isko,nski,nsko,w,x,delta,d,eta,llk_null)
    ! - COMPUTE LIKELIHOOD (SATURATED MODEL)
        call llk_saturated(ntasks, ns, nski, nsko, d, llk_sat)
    ! - COMPUTE NULL DEVIANCE
        null_dev = 2.0D0 * sum(llk_sat - llk_null)
! ---------------------------------------------------------------------------------------------
! CASE WHERE WE NEED TO COMPUTE THE FIRST LAMBDA
! ---------------------------------------------------------------------------------------------
    IF(lamfrac<1.0D0) THEN
    ! - COMPUTATION OF GRADIENT WITH BETA=0
        DO j=1,p
            IF(iex(j)==0)CYCLE
            call grad_hessj(ntasks,n,ns,iski,isko,nski,nsko,w,x(:,j),&
                            delta,d,eta,grad(j,:),h)
        ENDDO
    ! - COMPUTATION OF LAMBDA MAX
        al = 0.0D0
        call lambda_max(p, ntasks, reg, grad, pf, al, iex)
        IF(al < small)THEN
            ierr = (/1,7,1,0,0/)
            RETURN
        ENDIF
    ! - VALUE OF LOG DECREASE
        alf=lamfrac**(1.0D0/(nlam-1.0D0))
        IF(alf > big )THEN
            ierr = (/1,8,1,0,0/)
            RETURN
        ENDIF
    ! - INITIALIZE llk
        llk = llk_null
    ENDIF
! -------------------------------------------------------------------------------------------------
! BEGIN LAMBDA LOOP
! -------------------------------------------------------------------------------------------------
    DO l=1,nlam

    ! ---------------------------------------------------------------------------------------------
    ! COMPUTE (OR GET) LAMBDA VALUE
    ! ---------------------------------------------------------------------------------------------
        ! - USER DEFINED LAMBDA
        IF(lamfrac >= 1.0D0) THEN
            al = lam(l)
        ! - LOGARITHMIC SEQUENCE
        ELSE
            IF(l>1) al = al * alf
        ENDIF
    ! - FIRST LAMBDA EXCEPTION
        IF(lamfrac >= 1.0D0 .OR. l>1) THEN
    ! ---------------------------------------------------------------------------------------------
    ! INITIALIZE STRONG RULE SET
    ! ---------------------------------------------------------------------------------------------
        isr = iex
        IF(str ==1)THEN
            DO j=1,p
                IF (isr(j) == 0) CYCLE
                ! - COMPUTE GRADIENT
                call grad_hessj(ntasks,n,ns,iski,isko,nski,nsko,w,x(:,j),delta,d,eta,g,h)
                ! - COMPUTE CONDITION AND UPDATE SET
                SELECT CASE (reg)
                    CASE (0)
                        tmp = maxval(abs(g))/pf(j)
                    CASE (2)
                        tmp = sqrt(sum(g**2))/pf(j)
                END SELECT
                IF(tmp < 2.0D0*al - al0) isr(j) = 0
            ENDDO
        ENDIF
    ! ---------------------------------------------------------------------------------------------
    ! BEGIN STRONG RULE LOOP
    ! ---------------------------------------------------------------------------------------------
        DO
        ! -----------------------------------------------------------------------------------------
        ! SOLVE PENALIZED MLE WITH GPG ALGORITHM
        ! -----------------------------------------------------------------------------------------
        call gpg_descent(ntasks,p,n,ii,io,isr,ns,iski,isko,nski,nsko,&
                                w,x,delta,d,b,eta,al,pf,reg,alpha,&
                                llk,ierr,alg,eps,frac,ncycles,nupdates,maxit)
        ! - Catch error
        IF(ierr(1) == 1)THEN
            ierr(3) =l
            RETURN
        ENDIF
        ! -----------------------------------------------------------------------------------------
        ! STRONG RULE CHECK
        ! -----------------------------------------------------------------------------------------
        fl = 0
        IF(str ==1)THEN
            DO j=1,p
                IF(isr(j) == 1)CYCLE
                IF(iex(j) == 0)CYCLE
                ! - COMPUTE GRADIENT
                call grad_hessj(ntasks,n,ns,iski,isko,nski,nsko,w,x(:,j),delta,d,eta,g,h)
                ! - COMPUTE CONDITION AND UPDATE SET
                SELECT CASE (reg)
                    CASE (0)
                        tmp = maxval(abs(g(:)))/pf(j)
                    CASE (2)
                        tmp = sqrt(sum(g(:)**2))/pf(j)
                END SELECT
                IF(tmp > al) THEN
                    isr(j) = 1
                    fl = 1
                ENDIF
            ENDDO
        ENDIF
        IF(fl == 0) EXIT
    ! ---------------------------------------------------------------------------------------------
    ! END SRTONG RULE LOOP
    ! ---------------------------------------------------------------------------------------------
        ENDDO
    ! ---------------------------------------------------------------------------------------------
    ! FINAL CHECKS
    ! ---------------------------------------------------------------------------------------------

    ! - NUMBER OF VARIABLES, NORM AND ENTRY
        ! initialize counter
        DO j=1,p
            IF(iex(j) == 0) CYCLE
            ! compute norm
            SELECT CASE (reg)
                CASE (0)
                    betanorm(j,l) = maxval(abs(b(j,:)))
                CASE (2)
                    betanorm(j,l) = sqrt(sum(b(j,:)**2))
            END SELECT
            ! check if positive
            IF(betanorm(j,l) > small) THEN
                nbeta(l) = nbeta(l) + 1
                IF(entered (j) == 0 )THEN
                    entered(j) = l
                    nbetaever = nbetaever + 1
                ENDIF
            ENDIF
        ENDDO
    ! - KKT CONDITIONS
        DO j=1,p
            IF(iex(j) == 0) CYCLE
            ! compute gradient
            call grad_hessj(ntasks,n,ns,iski,isko,nski,nsko,w,x(:,j),delta,d,eta,g,h)
            ! gradient norm
            SELECT CASE (reg)
                CASE (0)
                    tmp = maxval(abs(g))
                CASE (2)
                    tmp = sqrt(sum(g**2))
            END SELECT
            ! condition value and check
            IF(betanorm(j,l) > small) THEN
                !non-zero case
                kkt(j,l) = tmp
            ELSE
                !zero case
                kkt0(j,l) = tmp
            ENDIF
        ENDDO
    ! - END FIRST LAMBDA CONDITION
        ENDIF
    ! - SAVE FINAL ESTIMATES AND OTHER OUTPUT (ETA, LLK, DEV, PDEV, LAMBDA, hhat)
        eta_path(:,l) = eta
        beta(:,:,l) = b
        llk_path(:,l) = llk
        dev(l) = 2.0D0 * sum(llk_sat - llk)
        pdev_path(l) = -(dev(l)-null_dev)/null_dev
        lam(l) = al
        call hhat_compute(ntasks,n,ns,iski,isko,nski,nsko,w,eta,hhat(:,l))
        ! store lambda for next iteration
        al0 = al
    ! - EARLY STOP
        ! CHECK KKT CONDITIONS
        ! CHECK PERCENT DEVIANCE EXPLAINED
        IF(flearlystop == 1) RETURN
! -------------------------------------------------------------------------------------------------
! END LAMBDA LOOP
! -------------------------------------------------------------------------------------------------
    ENDDO
! -------------------------------------------------------------------------------------------------
! FINAL CHECKS
! -------------------------------------------------------------------------------------------------

! - DESTANDARDIZATION
    DO l=1,nlam
        beta(:,:,l) = beta(:,:,l) / xsd
    ENDDO
! - CLEAN OUTPUT

! -------------------------------------------------------------------------------------------------
END SUBROUTINE mtcox_solutionpath
! -------------------------------------------------------------------------------------------------

