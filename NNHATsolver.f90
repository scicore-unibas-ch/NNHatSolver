program model

  implicit none

  integer i,iter,modelnum,nmodels,iselect,index

  integer,parameter::maxiter=10000000        !Maximum number of iterations/parameter set
  integer,parameter::ndata=13                !yr of actual data for phase II (yr)
  integer,parameter::maxtime1=100            !Maximum physical time for phase I  (yr)
  integer,parameter::maxtime2=maxtime1+ndata !Maximum physical time for phase II (yr)
  
  integer,parameter::maxnmodels=1            !Number of parameter sets
  integer,parameter::nvariables=26           !Number of evolved variables
  integer,parameter::nargs=6                 !Number of input arguments/model

  logical,parameter::rk2=.false.             !Runge-Kutta 2nd order
  logical,parameter::rk4=.false.             !Runge-Kutta 4th order
  logical,parameter::debug=.true.           !Set debug mode
  !If debug mode is true, the code will generate evolution.d and conservation.d
  !As it writes to hard disk every iteration, debug is considerably slower.

  double precision,parameter:: deltatimeini=1.d-7       !In yr
  double precision,parameter:: maxtimestep=0.01d0       !In yr
  double precision,parameter:: month=1.d0/12.d0         !1 month
  double precision,parameter:: dkt1=0.01d0              !Maximum % of change in variables Phase I
  double precision,parameter:: dkt2=0.002d0             !Maximum % of change in variables Phase II
  double precision,parameter:: factor=1.2d0             !Maximim change factor in timestep
  double precision,parameter:: tolerance=1.d-2          !Tolerance to keep the parameter set
  double precision,parameter:: theta=0.d0               !Integration order
                                                        !theta = 1 (1st order implicit)
                                                        !theta = 0.5 (2nd order)
                                                        !theta = 0 (1st order explicit)

  !Fixed parameters of HAT models
  double precision,parameter:: daysInYear = 365.d0
  double precision,parameter:: eta = 0.085d0 * daysInYear !31.025d0
  double precision,parameter:: epsilon = 0.7d0 
  double precision,parameter:: gamma = 0.001d0 * daysInYear !0.365d0
  
  double precision,dimension(nvariables)::y,fvar,fvar1,fvar2,fvar3,fvar4,y1,y2,y3

  double precision:: a,add
  double precision:: deltatime,time,dtmin,N2N1,VH1,VH2,r1h,r1l,r2h,r2l,ch,mugamma
  double precision:: Nl,Nh,Nal,Nah
  double precision:: Nh0,Nh10,Nh20,Nv10,Nv20
  double precision:: coef1,coef2,rPD_stg2
  double precision:: cumI1Lt,cumI1Ht,cumI2Lt,cumI2Ht
  double precision,dimension(ndata)::screened_as,pop_as,rate_as,screened_pd,pop_pd,slope,rPD_stg1

  integer len,status
  character starg*15
  double precision,dimension(nargs)::arg
  
  character number*5
  logical output


  !Read model parameters from arguments

  i=1
  do
     call get_command_argument(i,starg)
     if (len(trim(starg)).eq.0) stop 'Missing argument'
     read(starg,*) arg(i)
     if(debug)print '(a9,1x,i2,a2,es18.11)','Argument ',i,': ',arg(i)
     i=i+1
     if(i.gt.nargs) exit
  enddo
  if(i.lt.nargs) stop 'Too few arguments'

  !TO-DO
  !Now conver arg(:) into the actual parameters that we need for the model

  
  !After 100 yrs, r1l,r2l,r1h will change.
  !r1h is an annual vector/function (<20 yr)
  !r1l=r1h+raS (ras=value for 1st month and 0 for the rest 11 every year)
  !r2l=r2h+raS (r2h is fixed to the initial value)
  !Add cumulators after 100 yr
  !raS (depends on data and epsilon,N2N1=Nh/Nl fixed value)
  !implement initialization denpending pretty much on N2N1,VH1,VH2
  !Nh0=10000 fixed
  !Possible to run only second part for extrapolation  

  if(rk2.and.rk4)then
     print *,'Both methods RK2 and RK4 cannot be true at the same time!'
     stop
  endif
  
  !Initialization
  fvar(:)=0.d0
  y(:)=0.d0
  deltatime=deltatimeini
  time=0.d0
  if (debug) then
     nmodels=1
  else
     nmodels=maxnmodels
  endif
  
  !Open conservation file
  if (debug) then
     open(10,file='conservation.d')
     open(11,file='evolution.d')
  endif

  !Open accumulators file
  open(12,file='accumulators.d')
  
  !Variables number
  !10-Sl,  13-El,  11-I1l, 12-I2l,  9-Tl
  !16-Sh,  17-Eh,  14-I1h, 15-I2h, 18-Th
  !19-Sal, 21-Eal, 22-Ial, 20-Ral
  !23-Sah, 24-Eah, 26-Iah, 25-Rah
  ! 4-Svl,  3-Evl,  1-Ivl,  7-Uvl
  ! 5-Svh,  6-Evh,  2-Ivh,  8-Uvh


  global: do modelnum=1,nmodels

     !Set variable parameters
     !1-ch, 2-mugamma, 3-betavl, 4-betavh, 5-Nl, 6-Nh

     N2N1 = 0.08144843d0
     VH1  = 3.358809d0
     VH2  = 3.359019d0
     Nh0  = 10000.d0
     r1h  = 0.03285d0
     r1l  = r1h
     r2h  = 1.27895d-05
     r2l  = r2h
     ch   = 0.003330896d0
     mugamma = 0.5851145d0

     !Set initial conditions
     y=0.d0

     Nh10 = Nh0 / (1.d0 + N2N1)
     Nh20 = Nh10 * N2N1
     Nv10 = Nh10 * VH1  
     Nv20 = Nh20 * VH2


     
     y(10) = Nh10
     y(16) = Nh20
     y(19) = 4000.d0     !Depend on initial conditions
     y(23) = 4000.d0     !idem
     y(4)  = Nv10*0.95d0
     y(3)  = Nv10*0.02d0
     y(1)  = Nv10*0.02d0
     y(7)  = Nv10*0.01d0
     y(5)  = Nv20*0.95d0
     y(6)  = Nv20*0.02d0
     y(2)  = Nv20*0.02d0
     y(8)  = Nv20*0.01d0
     Nl=y(10)+y(13)+y(11)+y(12)
     Nh=y(16)+y(17)+y(14)+y(15)
     Nal=y(19)+y(21)+y(22)+y(20)
     Nah=y(23)+y(24)+y(26)+y(25)

     !Write initial conditions in evolution.d
     if(debug) then
        write(11,'(27(1x,es18.11))') time,y(10),y(13),y(11),y(12),y(9),y(16),&
             &    y(17),y(14),y(15),y(18),y(19),y(21),y(22),y(20),y(23),y(24),&
             &    y(26),y(25),y(4),y(3),y(1),y(7),y(5),y(6),y(2),y(8)
     endif

     write(number,'(I5.5)') modelnum
     open(13,file='steadystate'//number//'.d')

!********* PHASE I *********
     
     !Accumulators
     cumI1Lt=0.d0
     cumI1Ht=0.d0
     cumI2Lt=0.d0
     cumI2Ht=0.d0

     
     loop1: do iter=1,maxiter

        !Solve ODEs
        if(rk2) then
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y,fvar1,deltatime)
           y1(:)=y(:)+fvar1(:)
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y1,fvar2,deltatime)
           fvar(:)=0.5d0*fvar1(:)+0.5d0*fvar2(:)
        else if (rk4) then
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y,fvar1,deltatime*0.5d0)
           y1(:)=y(:)+fvar1(:)
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y1,fvar2,deltatime*0.5d0)
           y2(:)=y(:)+fvar2(:)
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y2,fvar3,deltatime)
           y3(:)=y(:)+fvar3(:)
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y3,fvar4,deltatime)
           fvar(:)=(fvar1(:)+2.d0*(fvar2(:)+fvar3(:))+fvar4(:))/6.d0
        else           
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,theta,y,fvar,deltatime)
        endif

        !Find new timestep
        dtmin=1.d60
        iselect=0
        var1: do i=1,nvariables
           if(y(i).ne.0.d0.and.fvar(i).ne.0.d0)dtmin=min(dtmin,abs(dkt1*y(i)/fvar(i)*deltatime))
        enddo var1

        deltatime=min(dtmin,factor*deltatime,maxtimestep)
        
        !Update variables and time
        y(:)=y(:)+fvar(:)
        time=time+deltatime

        !Accumulators
        !Cumulative number of infected people by stage and setting
        cumI1Lt=cumI1Lt+eta*y(13)*deltatime
        cumI2Lt=cumI2Lt+gamma*y(11)*deltatime
        cumI1Ht=cumI1Ht+eta*y(17)*deltatime
        cumI2Ht=cumI2Ht+gamma*y(14)*deltatime

        if(debug) then
           !Write conservation
           Nl=y(10)+y(13)+y(11)+y(12)
           Nh=y(16)+y(17)+y(14)+y(15)
           Nal=y(19)+y(21)+y(22)+y(20)
           Nah=y(23)+y(24)+y(26)+y(25)
           write(10,'(7(1x,es18.11),1x,i3,1x,es18.11)') time,deltatime,dtmin,Nl,Nh,Nal,Nah,iselect,y(12)

           !10-Sl,  13-El,  11-I1l, 12-I2l,  9-Tl
           !16-Sh,  17-Eh,  14-I1h, 15-I2h, 18-Th
           !19-Sal, 21-Eal, 22-Ial, 20-Ral
           !23-Sah, 24-Eah, 26-Iah, 25-Rah
           ! 4-Svl,  3-Evl,  1-Ivl,  7-Uvl
           ! 5-Svh,  6-Evh,  2-Ivh,  8-Uvh
     
           !Write evolution
           write(11,'(27(1x,es18.11))') time,y(10),y(13),y(11),y(12),y(9),y(16),&
                &    y(17),y(14),y(15),y(18),y(19),y(21),y(22),y(20),y(23),y(24),&
                &    y(26),y(25),y(4),y(3),y(1),y(7),y(5),y(6),y(2),y(8) 
        endif
        
        if (time.ge.maxtime1) exit loop1
        
     enddo loop1

     write(12,'(5(1x,es18.11))') time,cumI1Lt,cumI1Ht,cumI2Lt,cumI2Ht


 !********* PHASE II *********
     coef1 = 4.614399d0
     coef2 = 0.001796526d0
     output=.false.

     !Input data for human intervention
     open(1,file='data_AS.txt')
     open(2,file='data_PD.txt')
     read(1,*) !Header
     do i=1,ndata
        read(1,*) a,screened_as(i),a,a,a,pop_as(i)
     enddo
     read(2,*) !Header
     do i=1,ndata
        read(2,*) a,screened_pd(i),a,a,a,pop_pd(i)
     enddo
     close(1)
     close(2)
     rate_as(:)=-12.d0*log(1.d0-(screened_as(:)/((1.d0+N2N1)*epsilon*pop_as(:))))
     slope(:)=screened_pd(:)/(epsilon*pop_pd(:))
     rPD_stg1(:) = coef1*slope(:)
     rPD_stg2 = coef2*slope(1)

     !Initial timestep for Phase II
     deltatime=1.d-4

     loop2: do iter=1,maxiter

        !Human intervention at the beginning of every year after Phase I
        index=int(floor(time-dble(maxtime1)))+1
        
        !Only during the first month of every year of Phase II
        if(time-floor(time).gt.month) then
           add=0.d0
           if(output) then
              output=.false.
              write(12,'(5(1x,es18.11))') time,cumI1Lt,cumI1Ht,cumI2Lt,cumI2Ht
           endif
        else
           add=rate_as(index)
           output=.true.
        endif
        r1h=rPD_stg1(index)
        r2h=rPD_stg2
        r1l=rPD_stg1(index)+add
        r2l=rPD_stg2+add

        !Solve ODEs
        if(rk2) then
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y,fvar1,deltatime)
           y1(:)=y(:)+fvar1(:)
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y1,fvar2,deltatime)
           fvar(:)=0.5d0*fvar1(:)+0.5d0*fvar2(:)
        else if (rk4) then
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y,fvar1,deltatime*.5d0)
           y1(:)=y(:)+fvar1(:)
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y1,fvar2,deltatime*.5d0)
           y2(:)=y(:)+fvar2(:)
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y2,fvar3,deltatime)
           y3(:)=y(:)+fvar3(:)
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,1.d0,y3,fvar4,deltatime)
           fvar(:)=(fvar1(:)+2.d0*(fvar2(:)+fvar3(:))+fvar4(:))/6.d0
        else
           call solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,theta,y,fvar,deltatime)
        endif

        !Find new timestep
        dtmin=1.d60
        iselect=0
        var2: do i=1,nvariables
           if(y(i).ne.0.d0.and.fvar(i).ne.0.d0)dtmin=min(dtmin,abs(dkt2*y(i)/fvar(i)*deltatime))
        enddo var2

        deltatime=min(dtmin,factor*deltatime,maxtimestep)

        !Update variables and time
        y(:)=y(:)+fvar(:)
        time=time+deltatime
        Nl=y(10)+y(13)+y(11)+y(12)
        Nh=y(16)+y(17)+y(14)+y(15)
        Nal=y(19)+y(21)+y(22)+y(20)
        Nah=y(23)+y(24)+y(26)+y(25)

        !10-Sl,  13-El,  11-I1l, 12-I2l,  9-Tl
        !16-Sh,  17-Eh,  14-I1h, 15-I2h, 18-Th
        !19-Sal, 21-Eal, 22-Ial, 20-Ral
        !23-Sah, 24-Eah, 26-Iah, 25-Rah
        ! 4-Svl,  3-Evl,  1-Ivl,  7-Uvl
        ! 5-Svh,  6-Evh,  2-Ivh,  8-Uvh       

        !Accumulators
        !Cumulative number of infected people by stage and setting
        cumI1Lt=cumI1Lt+eta*y(13)*deltatime
        cumI2Lt=cumI2Lt+gamma*y(11)*deltatime
        cumI1Ht=cumI1Ht+eta*y(17)*deltatime
        cumI2Ht=cumI2Ht+gamma*y(14)*deltatime

!#
!# ## Store cummulative number of reported cases by stage and by setting
!# cumCasesPassive1t.L <- r_a * I_1L
!# cumCasesPassive2t.L <- r_b * I_2L
!# cumCasesPassive1t.H <- r_a * I_1H
!# cumCasesPassive2t.H <- r_b * I_2H
!#
!# cumCasesActive1t <- r_AS *I_1L
!# cumCasesActive2t <- r_AS *I_2L
        
        
        if(debug) then
           !Write conservation
           write(10,'(7(1x,es18.11),1x,i3,1x,es18.11)') time,deltatime,dtmin,Nl,Nh,Nal,Nah,iselect,y(12)

           !Write evolution
           write(11,'(31(1x,es18.11))') time,y(10),y(13),y(11),y(12),y(9),y(16),&
                &    y(17),y(14),y(15),y(18),y(19),y(21),y(22),y(20),y(23),y(24),&
                &    y(26),y(25),y(4),y(3),y(1),y(7),y(5),y(6),y(2),y(8),&
                &    cumI1Lt,cumI1Ht,cumI2Lt,cumI2Ht
        endif

        if(time.ge.maxtime2) exit loop2

     enddo loop2
     

     if(debug) then
        close(10)
        close(11)
     endif

     !Check if I1l and I1h latest change is below tolerance.
     !If so, write the final state as steady state.
     !Otherwise discard it
     if (abs(fvar(11)/y(11)).le.tolerance.and.abs(fvar(14)/y(14)).le.tolerance)then
        write(13,'(1x,i9,31(1x,es18.11))') modelnum,time,y(10),y(13),y(11),&
             &    y(12),y(9),y(16),y(17),y(14),y(15),y(18),y(19),y(21),y(22),&
             &    y(20),y(23),y(24),y(26),y(25),y(4),y(3),y(1),y(7),y(5),y(6),&
             &    y(2),y(8),cumI1Lt,cumI1Ht,cumI2Lt,cumI2Ht
     else
        write(13,'(1x,i9,1x,a9)') modelnum,'Discarded'
     endif
  
  enddo global
  close(12)
  close(13)


end program model


subroutine solver(r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah,theta,y,fvar,deltatime)

  implicit none

  integer i,j
  integer,parameter:: nvar = 26
  integer,parameter:: ndiagsup=3,ndiaginf=4,leftbox=6

  double precision,intent(in):: r1h,r1l,r2h,r2l,ch,mugamma,Nv10,Nv20,Nl,Nh,Nal,Nah
  double precision,intent(in)::deltatime,theta
  double precision,intent(in),dimension(nvar)::y
  double precision,intent(out),dimension(nvar)::fvar

  double precision,parameter:: daysInYear = 365.d0

  double precision,parameter:: alpha = 0.2d0*daysInYear  !73.d0
  double precision,parameter:: b = 0.8d0
  double precision,parameter:: betaal = 0.0014d0 * daysInYear * 4000.d0 !2044.d0
  double precision,parameter:: betaah = 0.0014d0 * daysInYear * 4000.d0 !2044.d0
  double precision,parameter:: cal = 0.d0
  double precision,parameter:: cah = 0.d0
  double precision,parameter:: delta = 0.006d0 * daysInYear !2.19d0
  double precision,parameter:: deltaa = 0.003d0 * daysInYear !1.095d0
  double precision,parameter:: eta = 0.085d0 * daysInYear !31.025d0
  double precision,parameter:: f = 0.333d0 * daysInYear !121.545d0
  double precision,parameter:: gamma = 0.001d0 * daysInYear !0.365d0
  double precision,parameter:: gammaah = 0.0019d0 * daysInYear !0.6935d0
  double precision,parameter:: gammaal = 0.002d0 * daysInYear !0.73d0
  double precision,parameter:: mu = 5.4795d0 * 1.d-5 * daysInYear !2.0000175d-2
  double precision,parameter:: mual = 0.0014d0 *daysInYear !0.511d0
  double precision,parameter:: muah = mual !0.511d0
  double precision,parameter:: mut = 0.d0 * daysInYear !0.d0
  double precision,parameter:: muv = 0.03d0 * daysInYear !10.95d0
  double precision,parameter:: nu = 0.037d0 * daysInYear !13.505d0
  double precision,parameter:: sigma = 0.4d0
  double precision,parameter:: sigmaal = 0.3d0
  double precision,parameter:: sigmaah = 0.3d0
  double precision,parameter:: xi = 0.62d0

  double precision betavl,betavh
  double precision factor1,theta_vlhl,theta_vlhh,theta_vlal
  double precision factor2,theta_vhah,theta_vhhh
  double precision lambdal,lambdaal,lambdaah,lambdah1,lambdah2
  double precision wll,wlh,wal,whh,wah,c1,c2
  double precision biglambdaL,biglambdaH
  
  double precision,dimension(nvar)::norm
  
  double precision,dimension(nvar,nvar+1)::phi
  
  !Assigning varying input parameters
  betavl = Nv10 * muv 
  betavh = Nv20 * muv

  !Calculating derived parameters
  !Terms for biting probabilities
  factor1 = sigma*(Nl+(1.d0-xi)*Nh)+sigmaal*Nal
  theta_vlhl = sigma*Nl/factor1
  theta_vlhh = sigma*Nh*(1.d0-xi)/factor1
  theta_vlal = sigmaal*Nal/factor1

  factor2 = sigma*xi*Nh+sigmaah*Nah
  theta_vhah = sigmaah*Nah/factor2
  theta_vhhh = sigma*xi*Nh/factor2


  lambdal = b*f*theta_vlhl/Nl
  lambdaal = b*f*theta_vlal*cal/Nal
  lambdaah = b*f*theta_vhah*cah/Nah
  lambdah1 = b*f*theta_vlhh/Nh
  lambdah2 = b*f*theta_vhhh/Nh

  wll = f*theta_vlhl/Nl
  wlh = f*theta_vlhh/Nh
  wal = f*theta_vlal*cal/Nal
  whh = f*theta_vhhh/Nh
  wah = f*theta_vhah*cah/Nah

  c1=ch
  c2=ch

  biglambdaL = wll*(c1*y(11) + c2*y(12)) + wlh*(c1*y(14) + c2*y(15)) + wal*y(22)
  biglambdaH = whh*(c1*y(14) + c2*y(15)) + wah*y(26)

  
  !10-Sl,  13-El,  11-I1l, 12-I2l,  9-Tl
  !16-Sh,  17-Eh,  14-I1h, 15-I2h, 18-Th
  !19-Sal, 21-Eal, 22-Ial, 20-Ral
  !23-Sah, 24-Eah, 26-Iah, 25-Rah
  ! 4-Svl,  3-Evl,  1-Ivl,  7-Uvl
  ! 5-Svh,  6-Evh,  2-Ivh,  8-Uvh


  !Independent vector
  fvar(1) = nu*y(3) - muv*y(1)
  fvar(2) = nu*y(6) - muv*y(2)
  fvar(3) = biglambdaL*y(4) - (muv + nu)*y(3)
  fvar(4) = betavl - (muv + biglambdaL + alpha)*y(4)
  fvar(5) = betavh - (muv + biglambdaH + alpha)*y(5)
  fvar(6) = biglambdaH*y(5) - (muv + nu)*y(6)
  fvar(7) = alpha*y(4) - muv*y(7)
  fvar(8) = alpha*y(5) - muv*y(8)
  fvar(9) = r1l*y(11) + r2l*y(12) - (mu + mut + delta)*y(9)
  fvar(10) = mu*(y(13) + y(11)) + (mu + mugamma)*y(12) + (mu + mut + delta)*y(9) - lambdal*y(1)*y(10)
  fvar(11) = eta*y(13) - (mu + gamma + r1l)*y(11)
  fvar(12) = gamma*y(11) - (mu + mugamma + r2l)*y(12)
  fvar(13) = lambdal*y(10)*y(1) - (mu + eta)*y(13)
  fvar(14) = eta*y(17) - (mu + gamma + r1h)*y(14)
  fvar(15) = gamma*y(14) - (mu + mugamma + r2h)*y(15)
  fvar(16) = mu*(y(17) + y(14)) + (mu + mugamma)*y(15) + (mu + mut + delta)*y(18) - (lambdah1*y(1) + lambdah2*y(2))*y(16)
  fvar(17) = (lambdah1*y(1)+lambdah2*y(2))*y(16) - (mu + eta)*y(17)
  fvar(18) = r1h*y(14) + r2h*y(15) - (mu + mut + delta)*y(18)
  fvar(19) = betaal + deltaa*y(20) - (mual + lambdaal*y(1))*y(19)
  fvar(20) = gammaal*y(22) - (mual + deltaa)*y(20)
  fvar(21) = lambdaal*y(1)*y(19) - (mual + eta)*y(21)
  fvar(22) = eta*y(21) - (mual + gammaal)*y(22)
  fvar(23) = betaah + deltaa*y(25) - (muah + lambdaah*y(2))*y(23)
  fvar(24) = lambdaah*y(2)*y(23) - (muah + eta)*y(24)
  fvar(25) = gammaah*y(26) - (muah + deltaa)*y(25)
  fvar(26) = eta*y(24) - (muah + gammaah)*y(26)
 
  !+1 terms in diagonal will be included later
  !so that we can multiply the whole matrix by deltatime efficiently.

  !Jacobian initialization
  phi(:,:)=0.d0

  !10-Sl,  13-El,  11-I1l, 12-I2l,  9-Tl
  !16-Sh,  17-Eh,  14-I1h, 15-I2h, 18-Th
  !19-Sal, 21-Eal, 22-Ial, 20-Ral
  !23-Sah, 24-Eah, 26-Iah, 25-Rah
  ! 4-Svl,  3-Evl,  1-Ivl,  7-Uvl
  ! 5-Svh,  6-Evh,  2-Ivh,  8-Uvh

  !Jacobian elements
  phi(1,1) = muv
  phi(1,3) = -nu
  phi(2,2) = muv
  phi(2,6) = -nu
  phi(3,3) = muv + nu
  phi(3,4) = -biglambdaL
  phi(3,11) = -wll*y(4)
  phi(3,12) = phi(3,11)
  phi(3,14) = -wlh*y(4)
  phi(3,15) = phi(3,14)
  phi(3,22) = -wal*y(4)
  phi(4,4) = muv + biglambdaL + alpha
  phi(4,11) = -phi(3,11)
  phi(4,12) = phi(4,11)
  phi(4,14) = -phi(3,14)
  phi(4,15) = phi(4,14)
  phi(4,22) = -phi(3,22)
  phi(5,5) = muv + biglambdaH + alpha
  phi(5,14) = whh*y(5)
  phi(5,15) = phi(4,14)
  phi(5,26) = wah*y(5)
  phi(6,5) = -biglambdaH
  phi(6,6) = muv + nu
  phi(6,14) = -phi(5,14)
  phi(6,15) = phi(6,14)
  phi(6,26) = -phi(5,26)
  phi(7,4) = -alpha
  phi(7,7) = muv
  phi(8,5) = -alpha
  phi(8,8) = muv
  phi(9,9) = mu + mut + delta
  phi(9,11) = -r1l
  phi(9,12) = -r2l
  phi(10,1) = lambdal*y(10)
  phi(10,9) = -(mu + mut + delta)
  phi(10,10) = lambdal*y(1)      
  phi(10,11) = -mu
  phi(10,12) = -(mu + mugamma)
  phi(10,13) = -mu
  phi(11,11) = mu + gamma + r1l
  phi(11,13) = -eta
  phi(12,11) = -gamma
  phi(12,12) = mu + mugamma + r2l
  phi(13,1) = -lambdal*y(10)
  phi(13,10) = -lambdal*y(1)
  phi(13,13) = mu + eta
  phi(14,14) = mu + gamma + r1h
  phi(14,17) = -eta
  phi(15,14) = -gamma
  phi(15,15) = mu + mugamma + r2h
  phi(16,1) = lambdah1*y(16)
  phi(16,2) = lambdah2*y(16)
  phi(16,14) = -mu
  phi(16,15) = -(mu + mugamma)
  phi(16,16) = lambdah1*y(1) + lambdah2*y(2)
  phi(16,17) = -mu
  phi(16,18) = -(mu + mut + delta)
  phi(17,1) = -phi(16,1)
  phi(17,2) = -phi(16,2)
  phi(17,17) = phi(13,13)
  phi(18,14) = -r1h
  phi(18,15) = -r2h
  phi(18,18) = phi(9,9)
  phi(19,1) = lambdaal*y(19)
  phi(19,19) = mual + lambdaal*y(1)
  phi(19,20) = -deltaa
  phi(20,20) = mual + deltaa
  phi(20,22) = -gammaal
  phi(21,1) = -phi(19,1)
  phi(21,19) = -lambdaal*y(1)
  phi(21,21) = mual + eta
  phi(22,21) = -eta
  phi(22,22) = mual + gammaal
  phi(23,2) = lambdaah*y(23)
  phi(23,23) = muah + lambdaah*y(2)
  phi(23,25) = -deltaa
  phi(24,2) = -phi(23,2)
  phi(24,23) = -lambdaah*y(2)
  phi(24,24) = muah + eta
  phi(25,25) = muah + deltaa
  phi(25,26) = -gammaah
  phi(26,24) = -eta
  phi(26,26) = muah + gammaah


  !Timestep inclussion
  phi=phi*deltatime*theta

  !Diagonal terms +1 due to self-derivative
  do i=1,nvar
     phi(i,i)=1.d0-phi(i,i)
  enddo

  
  !Normalization of the rows in the Jacobian
  do i=1,nvar
     norm(i)=abs(phi(i,1))
     do j=1,nvar
        norm(i)=max(norm(i),abs(phi(i,j)))
     enddo
     if(norm(i).eq.0.d0) then
        write(*,*) 'Zero row in the matrix...!!'
        stop
     endif
     phi(i,:)=phi(i,:)/norm(i)
     fvar(i)=fvar(i)*deltatime/norm(i)
  enddo

  !Additional column in Jacobian to pass independent variables vector to solver.
  phi(:,27) = fvar(:)

  !Sparse matrix solver. Returns solution in fvar
  call eigen(nvar,ndiagsup,ndiaginf,leftbox,phi,fvar,fvar)

end subroutine solver


! *******************************************************************           
    subroutine eigen(nisp,nds,ndi,ijd,phi,c,x)
  
      implicit none
!-----------------------------------------------------------------------        
!c    Prantzos, Arnould & Arcoragi  (ApJ 315,209 - 1987)                        
!c    inversion of the matrix (special sparse form)                             
!-----------------------------------------------------------------------        
!c    solves the equation: phi*x=c                                             
!c    warning: for computation convenience the array phi    
!c             allocates an additional column, which is equal       
!c             to vector c                                                      
!-----------------------------------------------------------------------        
!c nisp == dimension of the true matrix a (nisp x nisp), and number of          
!c         columns of array a. INPUT                             
!c c == vector containing the independent terms in the equation system. INPUT   
!c x == vector containing the solution. OUTPUT                                  
!c                                                                              
!c phi == array of dimensions (nisp x nisp+1) containing matrix a and vector c.
!c nds == width of the upper diagonal. INPUT           
!c ndi == width of the lower diagonal. INPUT             
!c ijd == number of columns in the left box + 1. INPUT    
!-----------------------------------------------------------------------        
      integer,intent(in)::nisp,nds,ndi,ijd

      DOUBLE PRECISION,dimension(nisp)::c,x
      DOUBLE PRECISION,dimension(nisp,nisp+1)::phi

      integer i,j
      integer nn,idel,jdel,jbal,jmax,imax,jmin
      DOUBLE PRECISION divj,som
      DOUBLE PRECISION,dimension(nisp)::div
      
      nn=nisp
      idel=ijd-1
      jdel=ijd-1

!---beginning of the gaussian elimination procedure                             
!-(1)-elimination of the lower diagonals                                        
                                                                               
      do 1000 jbal=ijd,nn-1
         div(jbal)=-1./phi(jbal,jbal)
         jmax=jbal+ndi
         if(jmax.gt.nn)jmax=nn
         imax=jbal+nds
         if(imax.gt.nn)imax=nn

!*vdir ignore recrdeps
         do 1001 j=jbal+1,jmax
            divj=div(jbal)*phi(j,jbal)
            do 10 i=1,ijd-1
               if(divj.ne.0.)phi(j,i)=divj*phi(jbal,i)+phi(j,i)
10          enddo
            do 20 i=jbal+1,imax
               if(divj.ne.0.)phi(j,i)=divj*phi(jbal,i)+phi(j,i)
20          enddo
1001     enddo
!*vdir ignore recrdeps
         do 1002 j=jbal+1,jmax
            if(phi(j,jbal).ne.0.)c(j)=div(jbal)*phi(j,jbal)*c(jbal)+c(j)
1002     enddo
1000  enddo
      div(nn)=-1./phi(nn,nn)

!-(2)-elimination of the upper diagonals and of the horizontal band
      do 2000 jbal=nn,ijd+1,-1
         jmin=jbal-ndi
         if(jmin.lt.ijd)jmin=ijd
         do 200 j=jmin,jbal-1
            if(phi(j,jbal).ne.0.)then
               divj=div(jbal)*phi(j,jbal)
               do 30 i=1,idel
                  phi(j,i)=divj*phi(jbal,i)+phi(j,i)
30             enddo
               c(j)=divj*c(jbal)+c(j)
            endif
200      enddo
         do 300 j=1,jdel
            if(phi(j,jbal).ne.0.)then
               divj=div(jbal)*phi(j,jbal)
               do 40 i=1,idel
                  phi(j,i)=divj*phi(jbal,i)+phi(j,i)
40             enddo
               c(j)=divj*c(jbal)+c(j)
            endif
300      enddo
2000  enddo
      do 400 j=1,jdel
         if(phi(j,ijd).ne.0.)then
            divj=div(ijd)*phi(j,ijd)
            do 50 i=1,idel
               phi(j,i)=divj*phi(ijd,i)+phi(j,i)
50          enddo
            c(j)=divj*c(ijd)+c(j)
         endif
400   enddo

!-(3)-gaussian elimination of the upper left square matrix                      
      do 3000 jbal=1,ijd-2
         div(jbal)=-1./phi(jbal,jbal)
         do 3001 j=jbal+1,ijd-1
            if(phi(j,jbal).eq.0.)exit
            divj=div(jbal)*phi(j,jbal)
            do 60 i=jbal+1,ijd-1
               phi(j,i)=divj*phi(jbal,i)+phi(j,i)
60          enddo
            c(j)=divj*c(jbal)+c(j)
3001     enddo
3000  enddo
      div(ijd-1)=-1./phi(ijd-1,ijd-1)
      x(ijd-1)=-div(ijd-1)*c(ijd-1)
      do 4000 jbal=ijd-2,1,-1
         som=0.
         do 70 i=ijd-1,jbal+1,-1
            som=som+phi(jbal,i)*x(i)
70       enddo
         x(jbal)=div(jbal)*(som-c(jbal))
4000  enddo
      som=0.
      do 80 i=1,idel
         som=som+phi(nn,i)*x(i)
80    enddo
      x(nn)=div(nn)*(som-c(nn))
      do 5000 jbal=nn-1,ijd,-1
         som=0.
         do 90 i=1,idel
            som=som+phi(jbal,i)*x(i)
90       enddo
         x(jbal)=div(jbal)*(som-c(jbal))
5000  enddo
      return
    end subroutine eigen

! **************************************************************                                   

