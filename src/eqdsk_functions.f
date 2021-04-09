c ------------------------------------------	
	subroutine equ_init
	
	USE ka2jet_mod, only: fmafl,bsig,Btfl,ieqdsk ! common definitions
	USE eqdsk_mod
	USE spline_mod
	
	implicit none
	
c       ==========================
        ! definitions
	INTEGER i,j,k,jd,ismid,idum
	REAL*8 psi_limit,BgeoRgeo
	REAL*8 dum,mindum,stp1,stp2,stp3
        REAL*8 stpx,stpz,dpp,dpx,dp2
	REAL*8,dimension(:),allocatable::fpp1,fpp2,fpp3,psi1d 
	REAL*8,dimension(:),allocatable::ps1_
	REAL*8 xp,yp,intpsi_
        REAL*8 f1,f2,f3,af,bf,cf
	
c       ==========================	

       
      ! -------------------------------------------------
      ! read info from EQDSK
      call equ_from_eqdsk
      
      ! ------------------------------------------------- 
      ! set steps in radial & vertical direction
      stpx=1e-2*(req(2)-req(1)) ! [m]
      stpz=1e-2*(zeq(2)-zeq(1))
      
      ! ------------------------------------------------- 
      ! first, initialize arrays with correct size
      call init_splines(nw,nh)  
       
      ! allocate array
      if(allocated(fpp1)) deallocate(fpp1)
      allocate(fpp1(nw))
      if(allocated(fpp2)) deallocate(fpp2)
      allocate(fpp2(nw))
      if(allocated(fpp3)) deallocate(fpp3)
      allocate(fpp3(nw))
      if(allocated(psi1d)) deallocate(psi1d)
      allocate(psi1d(nw))
      
      
      ! -------------------------------------------------
      ! PV1 ... PV9 : Spline representation os Psi(R,Z)
      
      pv1 = psirz
      call  splnrz(pv1,pv2,pv3,pv4,pv5,pv6,pv7,pv8,pv9,nw,nh)      
            
      ! make sure Psi on axis =0 and Psi-bdry=1.0
      ! If not, fix values -
      call check_psi
      ! At this point, we should have Fmafl=1.d0 (norm.)
      ! and Sibry=actual flux at the boundary.
      
      
      ! -------------------------------------------------
      ! reconstruct a 1D spline representation for
      ! Fpol(Psi) from eqdsk file (used for Bphi in
      ! regions inside the LCFS)
      ! Note: it looks like FLUSH puts 1D profiles
      ! vs. midplane Radius, instead of vs Psi=[0 ... Psiwall]
      
      ! look for Z index nearest to midplane, remap at correct Zmid
      if(ieqdsk.eq.2) then 
        if(allocated(ps1_)) deallocate(ps1_)
        allocate(ps1_(nw)) ! Psi vs Radius at midplane
	yp=zmaxis
	do k=1,nw 
	  xp=req(k)
	  call get_flux(xp,yp,intpsi_,idum) 
	  ps1_(k)= intpsi_
	enddo
      endif
      

        dpp=fmafl/(nw-1.)
        do j=1,nw
          psi1d(j)=(j-1.)*dpp
	  
	  if(ieqdsk.eq.1) then ! from standard EFIT
 	    fpp1(j)=fpol(j)
	  endif
	  
	  if(ieqdsk.eq.2) then ! from FLUSH EFIT - more involved...
	    
	    ! look for nearest point in Psi
	    idum=nw/2
	    mindum=1e12
	    do k=2,nw-1 
	      if(abs(psi1d(j)-ps1_(k)).lt.mindum) then
	       idum = k
	       mindum=abs(psi1d(j)-ps1_(k))
	      endif
	    enddo
	    
	    ! local coefficients for quadratic expansion
	    f1=fpol(idum-1)
	    f2=fpol(idum)
	    f3=fpol(idum+1)
	    stp1=ps1_(idum-1)-ps1_(idum)
	    stp2=ps1_(idum+1)-ps1_(idum)
	    
	    call get_coeffs(f1,f2,f3,stp1,stp2,af,bf,cf)
	    
	    ! infer Fpol(psi1d(j))
	    stp3=psi1d(j)-ps1_(idum)
	    fpp1(j)=af+bf*stp3+cf*stp3*stp3
	    
	  endif
        enddo
	
      ! Compute 1D spline coefficients for Fpol
      call splin1d(fpp1,fpp2,fpp3,dpp,nw)
      
c      ! debug:
c      do k=nw/2,nw
c       write(0,'(i6,f12.8)') k,fpp2(k)+2.*fpp3(k)*stp3
c      enddo
      
      
      ! -------------------------------------------------
      ! Find value of product Bgeo x Rgeo based on continuity of
      ! toroidal field across the LCFS
      BgeoRgeo=fpp1(nw)+fpp2(nw)*dpp+fpp3(nw)*dpp*dpp ! [T m]
c      write(0,'(a,f12.6)') 'BgeoRgeo=',BgeoRgeo
            
	    
      ! -------------------------------------------------
      ! compute field components Bx, Bz from Psi
      
      
      ! compute derivative
      do j=2,nw-1 ! loop over radius, height
        do k=1,nh
	  f1=psirz(j-1,k)
	  f2=psirz(j,k)
	  f3=psirz(j+1,k)
	  call get_coeffs(f1,f2,f3,-stpx,stpx,af,bf,cf)
	  bz1(j,k)=-sibry*bf/(1e-2*req(j))
	enddo
      enddo
      do j=1,nw ! loop over radius, height
	do k=2,nh-1
	  f1=psirz(j,k-1)
	  f2=psirz(j,k)
	  f3=psirz(j,k+1)
	  call get_coeffs(f1,f2,f3,-stpz,stpz,af,bf,cf)
	  bx1(j,k)=sibry*bf/(1e-2*req(j)) ! sign reversed as needed in ka2jet.f (???)
	enddo
      enddo
      
      ! fix boundaries, simple method
      bx1(1,:)  = 2.*bx1(2,:)-bx1(3,:)
      bx1(nw,:)= 2.*bx1(nw-1,:)-bx1(nw-2,:)
      bx1(:,1)  = 2.*bx1(:,2)-bx1(:,3)
      bx1(:,nh)= 2.*bx1(:,nh-1)-bx1(:,nh-2)
      bz1(1,:)  = 2.*bz1(2,:)-bz1(3,:)
      bz1(nw,:)= 2.*bz1(nw-1,:)-bz1(nw-2,:)
      bz1(:,1)  = 2.*bz1(:,2)-bz1(:,3)
      bz1(:,nh)= 2.*bz1(:,nh-1)-bz1(:,nh-2)
      
      
      ! Compute field component Bphi from Fpol(Psi) and
      ! fill in pv1 coefficient for Psi
      ! Since FPOL is a flux function, will need to 
      ! first remap (R,Z) -> Psi, then get Fpol(Psi),
      ! from which Bphi can be computed inside the LCFS.
      ! Ouside the LCFS, Bphi is set from the vacuum solution:
      ! Bphi,vac(R) = Bgeo x Rgeo / R, where Bgeo and Rgeo
      ! are the toroidal field and radius at the geometric axis.
      
      bsig = -sign(1.d0,bcentr) ! set sign of B field, used later on in Main [NEED TO CHECK VS FLUSH,EFIT,...]
      Btfl = bcentr ! field on axis, [T]
      
      do j=1,nw ! loop over radius, height
	do k=1,nh
	  if(psirz(j,k).le.fmafl) then ! inside LCFS - TBD: check wrt LCFS(R,Z) instead!
            
	    jd = psirz(j,k)*(nw-1.)/fmafl + 1
            jd = min(jd,nw-1)
            jd = max(jd,1)
            dpx = psirz(j,k)-(jd-1)*fmafl/(nw-1.)
            dp2 = dpx*dpx
	    dum=fpp1(jd)+fpp2(jd)*dpx+fpp3(jd)*dp2
	    
	    bp1(j,k)=dum/(1e-2*req(j))
	    
	  else ! outside LCFS, vacuum region
	    bp1(j,k)=BgeoRgeo/(1e-2*req(j))
	  endif
	  
	enddo
      enddo
      
      
      ! ----------------------------------------------------------
      ! Convert Bfield components from [T] to [G] and
      ! R,Z coordinates from [m] to [cm]
      bx1 = 1e4*bx1
      bz1 = 1e4*bz1
      bp1 = 1e4*bp1
      
      
      ! ----------------------------------------------------------
      ! compute spline coefficients. 
      ! Here I use a spline method adapted from ORBIT extended
      ! beyond the LCFS, work by R. B. White.
      
      call  splnrz(bx1,bx2,bx3,bx4,bx5,bx6,bx7,bx8,bx9,nw,nh)
      call  splnrz(bz1,bz2,bz3,bz4,bz5,bz6,bz7,bz8,bz9,nw,nh)
      call  splnrz(bp1,bp2,bp3,bp4,bp5,bp6,bp7,bp8,bp9,nw,nh)
      
c      ! debug:
c      call dump_bsplines(nw,nh)
c      call dump_psisplines(nw,nh)    

      
      
      ! -------------------------------------------------
      ! bookkeeping - deallocate arrays that are not used

       if(allocated(fpp1)) deallocate(fpp1)
       if(allocated(fpp2)) deallocate(fpp2)
       if(allocated(fpp3)) deallocate(fpp3)
       if(allocated(psi1d)) deallocate(psi1d)
       if(allocated(ps1_)) deallocate(ps1_)
	
       ! temp. arrays for R,Z splines
       if(allocated(bmat)) deallocate(bmat)
       if(allocated(wk2)) deallocate(wk2)
	
	
	write(0,*) '... done.'
	write(0,*) '--------------------------------------------------'
	write(0,*) ''
	
	return
	end

c -----------------------------------------------------------	
	subroutine check_psi
	
	USE ka2jet_mod, only: fmafl,ieqdsk
	USE eqdsk_mod
	USE spline_mod
	
	implicit none
        
c       ==========================
        !  definitions
	integer ier
	real*8 psi_min,psi_bdry
	real*8 xp,yp
c       ==========================	

	
	! Flux on axis
	xp = rmaxis
	yp = zmaxis
	call get_flux(xp,yp,psi_min,ier) 
	
	if(ieqdsk.eq.2) then
	if(psi_min.ne.0.d0) then ! fix it
	  
	  write(0,*)
	  if(psi_min.ne.0.d0) then
	   write(0,'(a,f8.4,a)') '*** Psi on axis is ',psi_min,
     >                        ': will be reset to zero' 
          endif
	  write(0,*)
	  
	  pv1=pv1-psi_min
	  psirz=psirz-psi_min
	  sibry=sibry-simag
          
	  ! It looks like EQDSK from FLUSH has FMAFL.ne.1,
	  ! And yet the Psi(R,Z) values are already normalized!
	  ! So, reset Fmafl and Sibry here, but don't (re)normalize Psi
	  fmafl=1.d0

	  
	  ! re-compute spline with updated values
	  call  splnrz(pv1,pv2,pv3,pv4,pv5,pv6,pv7,pv8,pv9,nw,nh)      
	  
	endif
	endif
	
	! T.B.D.
	if(ieqdsk.eq.1) then
	 write(0,*) '*** Subroutine CHECK_PSI not yet tested'
	 write(0,*) '*** for EQDSK input from standard EFIT:'
	 write(0,*) '*** Psi normalization may be inaccurate!' 
	endif
	
	return
	end

c -----------------------------------------------------------	
	subroutine equ_from_eqdsk 
	
	USE ka2jet_mod, only: ieqdsk,eq_filename,fmafl
	USE eqdsk_mod
	
	implicit none
        
c       ==========================
        !  definitions
	INTEGER nshot
	INTEGER i,j,k,jd,idum
        INTEGER nbbbs,limitr
	INTEGER neqdsk
        REAL*8 rcentr,zmid
        REAL*8 current,xdum,dum
        REAL*8,dimension(:),allocatable::press,ffprim,pprime
        REAL*8,dimension(:),allocatable::qpsi
  	
        character*12 caseid(6)
	character*24 ifrmt1,ifrmt2,ifrmt3,ifrmt4,ifrmt5
	character*24 strdum
c       ==========================	

        neqdsk=61 ! file ID
	
	
	write(0,*) ''
	write(0,'(a)') '-----------------------------------------------------'
	write(0,'(2a)') 'Init equilibrium from EQDSK file ',eq_filename
	
	
      ! Read content of EQDSK file
      open(neqdsk,file=eq_filename,status='unknown') 
      
      if(ieqdsk.eq.1) then  ! Standard EQDSK format from EFIT
         ifrmt1='(6a8,3i4)'
	 ifrmt2='(5e16.9)'
	 ifrmt3='(2i5)'
	 ifrmt4='(5e16.9)'
	 ifrmt5='(5e16.9)'
	 read (neqdsk,ifrmt1) (caseid(i),i=1,6),idum,nw,nh
      endif
      if(ieqdsk.eq.2) then  ! EQDSK format from FLUSH
         ifrmt1='(a8,a10,a4,a6,a12,3i4)'
	 ifrmt2='(5e21.12)'
	 ifrmt3='(2i4)'
	 ifrmt4='(2e21.12)'
	 ifrmt5='(5f14.6)'
         read (neqdsk,ifrmt1) (caseid(i),i=1,5),idum,nw,nh
	 read(neqdsk,*) ! skip line
      endif
      
      
      ! allocate array
      if(allocated(psirz)) deallocate(psirz)
      allocate(psirz(nw,nh))
      if(allocated(fpol)) deallocate(fpol)
      allocate(fpol(nw))
      if(allocated(press)) deallocate(press)
      allocate(press(nw))
      if(allocated(ffprim)) deallocate(ffprim)
      allocate(ffprim(nw))
      if(allocated(pprime)) deallocate(pprime)
      allocate(pprime(nw))
      if(allocated(qpsi)) deallocate(qpsi)
      allocate(qpsi(nw))
      if(allocated(req)) deallocate(req)
      allocate(req(nw))
      if(allocated(zeq)) deallocate(zeq)
      allocate(zeq(nh))
      
      
      
      ! read data
      read (neqdsk,ifrmt5) rdim,zdim,rcentr,rleft,zmid
      read (neqdsk,ifrmt5) rmaxis,zmaxis,simag,sibry,bcentr
      read (neqdsk,ifrmt5) current,simag,xdum,rmaxis,xdum
      read (neqdsk,ifrmt5) zmaxis,xdum,sibry,xdum,xdum
      
      ! debug:
c      write(0,'(5f12.6)') rdim,zdim,rcentr,rleft,zmid
c      write(0,'(5f12.6)') rmaxis,zmaxis,simag,sibry,bcentr
c      write(0,'(5f12.6)') current,simag,xdum,rmaxis,xdum
c      write(0,'(5f12.6)') zmaxis,xdum,sibry,xdum,xdum
      
      if(ieqdsk.eq.2) read(neqdsk,*)  ! skip empty lines

      
      read (neqdsk,ifrmt2) (fpol(i),i=1,nw)
      if(ieqdsk.eq.2) then
        read(neqdsk,*)  ! skip empty lines
	read(neqdsk,*)
      endif
      
      read (neqdsk,ifrmt2) (press(i),i=1,nw)
      if(ieqdsk.eq.2) then
        read(neqdsk,*)  ! skip empty lines
	read(neqdsk,*)
      endif
      
      read (neqdsk,ifrmt2) (ffprim(i),i=1,nw)
      if(ieqdsk.eq.2) then
        read(neqdsk,*)  ! skip empty lines
	read(neqdsk,*)
      endif
      
      read (neqdsk,ifrmt2) (pprime(i),i=1,nw)
      if(ieqdsk.eq.2) then
        read(neqdsk,*)  ! skip empty lines
	read(neqdsk,*)
      endif
      
      read (neqdsk,ifrmt2) ((psirz(i,j),i=1,nw),j=1,nh)
      if(ieqdsk.eq.2) then
        read(neqdsk,*)  ! skip empty lines
	read(neqdsk,*)
      endif
      
      read (neqdsk,ifrmt2) (qpsi(i),i=1,nw)
      if(ieqdsk.eq.2) then
        read(neqdsk,*)  ! skip empty lines
	read(neqdsk,*)
      endif
      
      read (neqdsk,ifrmt3) nbbbs,limitr
      if(ieqdsk.eq.2) read(neqdsk,*)  ! skip empty lines
      
      ! allocate arrays
      if(nbbbs.gt.0) then
       if(allocated(rlcfs)) deallocate(rlcfs)
       allocate(rlcfs(nbbbs))
       if(allocated(zlcfs)) deallocate(zlcfs)
       allocate(zlcfs(nbbbs))
      endif

      if(limitr.gt.0) then
       if(allocated(rlim)) deallocate(rlim)
       allocate(rlim(limitr))
       if(allocated(zlim)) deallocate(zlim)
       allocate(zlim(limitr))
      endif
      
      ! read LCFS data
      if(nbbbs.gt.0) then
        read (neqdsk,ifrmt4) (rlcfs(i),zlcfs(i),i=1,nbbbs)
      endif
      
      ! Read limiter data
      if(limitr.gt.0) then
        read (neqdsk,ifrmt4) (rlim(i),zlim(i),i=1,limitr)
      endif 
      
      ! close eqdsk file
      close(neqdsk)
      
      ! -------------------------------------------------
      ! Convert lengths from [m] to [cm]
      rmaxis=1e2*rmaxis
      zmaxis=1e2*zmaxis
      rcentr=1e2*rcentr
      rleft=1e2*rleft
      rdim=1e2*rdim
      zmid=1e2*zmid
      zdim=1e2*zdim
      rlcfs=1e2*rlcfs
      zlcfs=1e2*zlcfs
      rlim=1e2*rlim
      zlim=1e2*zlim
      
      
      ! -------------------------------------------------
      ! Rescale Psi(R,Z) from eqdsk file, make it
      ! consistent with expected Psi=0 on axis
      ! NOTE: this is further enforced after the
      ! splines are computed, along with normalizations.
c      psirz = (psirz - simag)
c      sibry = (sibry - simag)
c      simag = 0.d0 ! reset
      fmafl = sibry ! set flux at LCFS, used in main program
      
      ! -------------------------------------------------
      ! reconstruct the (R,Z) grid for data in eqdsk file
      do j=1,nw
        req(j)=rleft+rdim*(j-1.d0)/(nw-1.)
      enddo
      do j=1,nh
        zeq(j)=zmid-0.5*zdim+zdim*(j-1.d0)/(nh-1.)
      enddo
      
      
      ! -------------------------------------------------
      ! bookkeeping - deallocate arrays that are not used
      if(allocated(press)) deallocate(press)
      if(allocated(ffprim)) deallocate(ffprim)
      if(allocated(pprime)) deallocate(pprime)
      if(allocated(qpsi)) deallocate(qpsi)
      
      return
      end


c ------------------------------------------	
	subroutine init_splines(nw,nh) 
	
	USE spline_mod
	implicit none
	
	integer nw,nh,nmax
	
	! init spline coefficients based on EQDSK grid size
	if(allocated(bx1)) deallocate(bx1)
	allocate(bx1(nw,nh))
	if(allocated(bx2)) deallocate(bx2)
	allocate(bx2(nw,nh))
	if(allocated(bx3)) deallocate(bx3)
	allocate(bx3(nw,nh))
	if(allocated(bx4)) deallocate(bx4)
	allocate(bx4(nw,nh))
	if(allocated(bx5)) deallocate(bx5)
	allocate(bx5(nw,nh))
	if(allocated(bx6)) deallocate(bx6)
	allocate(bx6(nw,nh))
	if(allocated(bx7)) deallocate(bx7)
	allocate(bx7(nw,nh))
	if(allocated(bx8)) deallocate(bx8)
	allocate(bx8(nw,nh))
	if(allocated(bx9)) deallocate(bx9)
	allocate(bx9(nw,nh))
	
	if(allocated(bz1)) deallocate(bz1)
	allocate(bz1(nw,nh))
	if(allocated(bz2)) deallocate(bz2)
	allocate(bz2(nw,nh))
	if(allocated(bz3)) deallocate(bz3)
	allocate(bz3(nw,nh))
	if(allocated(bz4)) deallocate(bz4)
	allocate(bz4(nw,nh))
	if(allocated(bz5)) deallocate(bz5)
	allocate(bz5(nw,nh))
	if(allocated(bz6)) deallocate(bz6)
	allocate(bz6(nw,nh))
	if(allocated(bz7)) deallocate(bz7)
	allocate(bz7(nw,nh))
	if(allocated(bz8)) deallocate(bz8)
	allocate(bz8(nw,nh))
	if(allocated(bz9)) deallocate(bz9)
	allocate(bz9(nw,nh))
	
	if(allocated(bp1)) deallocate(bp1)
	allocate(bp1(nw,nh))
	if(allocated(bp2)) deallocate(bp2)
	allocate(bp2(nw,nh))
	if(allocated(bp3)) deallocate(bp3)
	allocate(bp3(nw,nh))
	if(allocated(bp4)) deallocate(bp4)
	allocate(bp4(nw,nh))
	if(allocated(bp5)) deallocate(bp5)
	allocate(bp5(nw,nh))
	if(allocated(bp6)) deallocate(bp6)
	allocate(bp6(nw,nh))
	if(allocated(bp7)) deallocate(bp7)
	allocate(bp7(nw,nh))
	if(allocated(bp8)) deallocate(bp8)
	allocate(bp8(nw,nh))
	if(allocated(bp9)) deallocate(bp9)
	allocate(bp9(nw,nh))
	
	if(allocated(pv1)) deallocate(pv1)
	allocate(pv1(nw,nh))
	if(allocated(pv2)) deallocate(pv2)
	allocate(pv2(nw,nh))
	if(allocated(pv3)) deallocate(pv3)
	allocate(pv3(nw,nh))
	if(allocated(pv4)) deallocate(pv4)
	allocate(pv4(nw,nh))
	if(allocated(pv5)) deallocate(pv5)
	allocate(pv5(nw,nh))
	if(allocated(pv6)) deallocate(pv6)
	allocate(pv6(nw,nh))
	if(allocated(pv7)) deallocate(pv7)
	allocate(pv7(nw,nh))
	if(allocated(pv8)) deallocate(pv8)
	allocate(pv8(nw,nh))
	if(allocated(pv9)) deallocate(pv9)
	allocate(pv9(nw,nh))
	
	if(allocated(bmat)) deallocate(bmat)
	nmax=max(nw,nh)
	allocate(bmat(nmax*nmax))
	if(allocated(wk2)) deallocate(wk2)
	allocate(wk2(nmax))
	
	return
	end
	
	
c ------------------------------------------	
	subroutine get_flux(xp,yp,psi,ier) 
	
	USE eqdsk_mod, only: nw,nh,req,zeq
	USE spline_mod
	
	implicit none
C============	
	integer ier
	integer jd,kd
	real*8 xp,yp,psi
	real*8 dx,dx2,dz,dz2
C============	
	
	ier=0
	
	jd = (xp - req(1))/(req(2) - req(1)) + 1
        jd = min(jd,nw)
        jd = max(jd,1)
        kd = (yp - zeq(1))/(zeq(2) - zeq(1)) + 1
        kd = max(1,kd)
        kd = min(nh,kd)
        dx = xp  - req(jd)
        dx2 = dx*dx
        dz = yp  - zeq(kd)
        dz2 = dz*dz
	
	psi = pv1(jd,kd) + pv2(jd,kd)*dx + pv3(jd,kd)*dx2
     a  + pv4(jd,kd)*dz + pv5(jd,kd)*dx*dz + pv6(jd,kd)*dz*dx2
     b  + pv7(jd,kd)*dz2 + pv8(jd,kd)*dz2*dx + pv9(jd,kd)*dz2*dx2
	
	return
	end
	
	
c ------------------------------------------	
c	subroutine FLUPX(IER) ! not needed
c	implicit none
c	integer ier
c	ier=0	
c	return
c	end
	
	
c ------------------------------------------	
	subroutine b_from_splines(XP,YP,BR,BZ,BT,IER) 
	
	USE ka2jet_mod ! common definitions
	USE eqdsk_mod
	USE spline_mod
	implicit none
C============	
	integer ier
	integer jd,kd
	real*8 xp,yp,psi
	real*8 dx,dx2,dz,dz2
	real*8 br,bz,bt
C============	

	ier=0
	
	jd = (xp - req(1))/(req(2) - req(1)) + 1
        jd = min(jd,nw)
        jd = max(jd,1)
        kd = (yp - zeq(1))/(zeq(2) - zeq(1)) + 1
        kd = max(1,kd)
        kd = min(nh,kd)
        dx = xp  - req(jd)
        dx2 = dx*dx
        dz = yp  - zeq(kd)
        dz2 = dz*dz
	
	! Radial
	br = bx1(jd,kd) + bx2(jd,kd)*dx + bx3(jd,kd)*dx2
     a  + bx4(jd,kd)*dz + bx5(jd,kd)*dx*dz + bx6(jd,kd)*dz*dx2
     b  + bx7(jd,kd)*dz2 + bx8(jd,kd)*dz2*dx + bx9(jd,kd)*dz2*dx2
	
	! Vertical
	bz = bz1(jd,kd) + bz2(jd,kd)*dx + bz3(jd,kd)*dx2
     a  + bz4(jd,kd)*dz + bz5(jd,kd)*dx*dz + bz6(jd,kd)*dz*dx2
     b  + bz7(jd,kd)*dz2 + bz8(jd,kd)*dz2*dx + bz9(jd,kd)*dz2*dx2
	
	! Toroidal
	bt = bp1(jd,kd) + bp2(jd,kd)*dx + bp3(jd,kd)*dx2
     a  + bp4(jd,kd)*dz + bp5(jd,kd)*dx*dz + bp6(jd,kd)*dz*dx2
     b  + bp7(jd,kd)*dz2 + bp8(jd,kd)*dz2*dx + bp9(jd,kd)*dz2*dx2
      
	return
	end
	
	
ccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine splin1d(f1,f2,f3,dpx,lst)
      
      USE ka2jet_mod, only: twopi
      
      IMPLICIT NONE
      
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER lst,lspm,lsp2,j,jm,jp,jpp
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 f1,f2,f3,dpx
C============
      dimension f1(lst),f2(lst),f3(lst)
ccc    spline calculation
ccc    f(x) = f1 + f2*dx + f3*dx**2
cccc  before calling splin1 load function f into f1(j)=f(pj)
cccc   after call,  f2(j),f3(j) will be determined
      
      lspm = lst - 1
      lsp2 = lst - 2
ccccc
cccc
cccc  set f2(1) to leave f1(2) unmoved by smoothing
      f2(1)=(10.D0*f1(2)-7.D0*f1(1)-3.D0*f1(3))/(4.D0*dpx)

cccc
      do j = 2,lspm
               jm = j - 1
               jp = j + 1
               jpp = min(j + 2,lst)
               f2(j) = -f2(jm) + 2*(f1(j)-f1(jm))/dpx
cccc  smooth f1
         if(jp.ne.lst) 
     >      f1(jp) =.4D0*dpx*f2(j)+.3D0*f1(jpp)+.7D0*f1(j) 
      enddo
       
       f2(lst) = f2(lspm)
       do j = 1,lspm
                  jp = j + 1
                  f3(j) = (f2(jp)-f2(j))/(2.d0*dpx)
       enddo
       f3(lst) = 0.d0	!(f1(lst)-f1(lspm)-f2(lspm)*dpx)/dpx**2 
ccc-f1,f2,f3-finished
ccc
      return
      end
cccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine get_coeffs(f1,f2,f3,dx1,dx3,af,bf,cf)
      
      IMPLICIT NONE

C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 f1,f2,f3,dx1,dx3
      REAL*8 v1,v2,v3
      REAL*8 k,af,bf,cf
C============
       
        bf = dx1*dx3/(dx1-dx3)*((f3-f2)/(dx3*dx3)+(f2-f1)/(dx1*dx1))
        cf = 0.5d0*((f1-f2)/(dx1*dx1)-bf/dx1+(f3-f2)/(dx3*dx3)-bf/dx3)
	af = f2
	
      return
      
      end
cccccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccccc
      subroutine splnrz(f1,f2,f3,f4,f5,f6,f7,f8,f9,nw,nh)
      
      USE spline_mod, only: bmat,wk2
      USE eqdsk_mod, only: req,zeq
      IMPLICIT NONE

cc    Spline Bx, Bz, Bphi, pol on rz grid
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER nw,nh ! grid size, width and height
      INTEGER lgw,lgmw,lgcw,lgh,lgmh,lgch
      INTEGER l,ld,lmax,lmax2,k,j,jm,jmm,jp,jpp
      INTEGER km,kmm,kp,kpp,ier,jd,lm,kmid,jsit,kn,nd
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 f1,f2,f3,f4,f5,f6,f7,f8,f9,fmax
      REAL*8 dx,dz,dum,dum1,wdum,ddf,err0
C============
      integer ispln
      dimension f1(nw,nh),f2(nw,nh),f3(nw,nh),f4(nw,nh)
     a  ,f5(nw,nh),f6(nw,nh),f7(nw,nh),f8(nw,nh),f9(nw,nh)
     b  ,wdum(nw,nh)
C============
ccc     grid given by  Req, Zeq
cccc-in domain j,k  f(x,z) = f1 + f2*dx + f3*dx**2 + f4*dz + f5*dz*dx
cccc  + f6*dz*dx**2 + f7*dz**2 + f8*dz**2*dx + f9*dz**2*dx**2
ccc-  gelg solves bmat( l,m)x(m) = wk2(m), m=1,ndim
ccc   wk2 contains x after the call
cccc  bmat stored as vector bmat(lm),lm = l + (m-1)*ndim
cccc  before calling spln load function f into f1(j,k)=f(xj,zk)
cccc   after call,  f2(j,k),f3(j,k) etc will be determined
      
      lgw = nw
      lgmw = nw-1
      lgcw = nw/2
      lgh = nh
      lgmh = nh-1
      lgch = nh/2
      
      ! grid steps
      dx = req(2) - req(1)
      dz = zeq(2) - zeq(1)
      
ccccccc  SPLINE f1, f2, f3
      
      do k = 1,lgw
         do j = 2,lgmh
            jm = j - 1
            jp = j + 1
            jpp = min(j + 2,lgw)
            f2(j,k)=-f2(jm,k) + 2.d0*(f1(j,k)-f1(jm,k))/dx
cccc  smooth f1
            if(jp.ne.lgmw) then    ! 4/23/99
              f1(jp,k)=.4D0*dx*f2(j,k)+.3D0*f1(jpp,k)+.7D0*f1(j,k)
            endif
         enddo
      enddo
      
      do k = 1,lgh
         f2(lgw,k) = f2(lgmw,k)
         do j = 1,lgmw
            jp = j + 1
            f3(j,k) = (f2(jp,k)-f2(j,k))/(2.d0*dx)
         enddo
         f3(lgmw,k)=(f1(lgw,k)-f1(lgmw,k)-f2(lgmw,k)*dx)/dx**2 ! 4/23/9
      enddo
ccc-f1,f2,f3-finished
ccc
ccc
ccc   find matrix for f4

      lmax = lgh
      do j = 1,lgmw
ccc
cccc-clear wk2
        wk2 = 0.d0
        lmax2 = lmax*lmax
	
ccc
ccccc-right hand side
        do k = 1,lgh
          kp = min(k + 1,lgh)
          wk2(k) = 2*f1(j,kp) - 2.d0*f1(j,k)
        enddo
        
	call ldbmat(lgh)
        err0 = 1.D-6
        call gelg(wk2,bmat,lmax,1,err0,ier)
	
        do k = 1,lgh
          f4(j,k) = wk2(k)
        enddo
        do k = 1,lgh
           kp = min(k + 1,lgh)
           f7(j,k) = (f4(j,kp) - f4(j,k))/(2.d0*dz)
        enddo
	
      enddo ! index j

ccc
ccc-f4-f7-finished
cccccc
ccc   find matrix for f5
      do j = 1,lgmw
ccc
cccc-clear wk2,bmat
        do l = 1,lmax
           wk2(l) = 0
        enddo
ccc
ccccc-right hand side
        do k = 1,lgh
          kp = min(k + 1,lgh)
          wk2(k) = 2.d0*f2(j,kp) - 2.d0*f2(j,k)
        enddo

        call ldbmat(lgh)
        err0 = 1.D-6
        call gelg(wk2,bmat,lmax,1,err0,ier)
	
        do k = 1,lgh
          f5(j,k) = wk2(k)
        enddo
        do k = 1,lgh
          kp = min(k + 1,lgh)
          f8(j,k) = (f5(j,kp) - f5(j,k))/(2.d0*dz)
        enddo
 
      enddo ! index j
ccc
ccc-f5-f8-finished
ccc
ccc   find matrix for f6
      do j = 1,lgmw
ccc
cccc-clear wk2,bmat
        wk2 = 0.d0
ccc
ccccc-right hand side
        do k = 1,lgh
           kp = min(k + 1,lgh)
           wk2(k) = 2.d0*f3(j,kp) - 2.d0*f3(j,k)
        enddo

        call ldbmat(lgh)
        err0 = 1.D-6
        call gelg(wk2,bmat,lmax,1,err0,ier)
	
        do k = 1,lgh
           f6(j,k) = wk2(k)
        enddo
        do k = 1,lgh
           kp = min(k + 1,lgh)
           f9(j,k) = (f6(j,kp) - f6(j,k))/(2.d0*dz)
        enddo
	
      enddo ! index j
ccc
ccc-f6-f9-finished
ccc
      return
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine ldbmat(nh)
      
      USE spline_mod, only: bmat
      USE eqdsk_mod, only: zeq
      
      IMPLICIT NONE

cccc  bmat stored as vector bmat(lm),lm = l + (m-1)*lmax
cccc
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER nh
      INTEGER lmax,lmax2,lm,l,m,mm,lmm
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 dz
C============
      
      lmax = nh
      dz = zeq(2) - zeq(1)
      
ccc
cccc-clear bmat
      lmax2 = lmax*lmax
      bmat = 0.d0
ccc
cccc  now matrix-stored-column by -column
       do l = 1,lmax
          m = l
          mm = m + 1
          lm = l + (m-1)*lmax
          lmm = l + (mm-1)*lmax
          bmat(lm) = dz
          bmat(lmm) = dz
       enddo
       
      return
      
      end
ccccccccccccccccccccccccccccccccccccccccccccccccccc
      subroutine init_reqzeq(nrp,nzp)
      
      USE ka2jet_mod
      USE eqdsk_mod, only: req,zeq
      IMPLICIT NONE
      
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER i,k
      INTEGER nrp,nzp,nv_
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 rmin_,rmax_,zmin_,zmax_
      REAL*8 stpr_,stpz_
       
C============
      
      ! read vessel boundaries
        i=1
	nv_=1
	
	open(3,file='JETWALL.DAT',status='old')
1	read(3,*,end=2) rv(i),zv(i)
        i=i+1
	go to 1
2	nv_=i-1
	close(3)
      
      ! set range
      rmin_=1e12
      rmax_=-1e12
      zmin_=1e12
      zmax_=-1e12
      do i=1,nv_
       if(rv(i).lt.rmin_) rmin_=rv(i)
       if(rv(i).gt.rmax_) rmax_=rv(i)
       if(zv(i).lt.zmin_) zmin_=zv(i)
       if(zv(i).gt.zmax_) zmax_=zv(i)
      enddo
      stpr_=(rmax_-rmin_)/(nrp-1.d0)
      stpz_=(zmax_-zmin_)/(nzp-1.d0)
      
      
      ! -------------------------------
      ! build R,Z grid
      
      if(allocated(req)) deallocate(req)
      allocate(req(nrp))
      if(allocated(zeq)) deallocate(zeq)
      allocate(zeq(nzp))
      
      do i=1,nrp
       req(i)=rmin_ + (i-1.)*stpr_
      enddo
      do i=1,nzp
       zeq(i)=zmin_ + (i-1.)*stpz_
      enddo
      
      ! set units from [m] to [cm]
      req=1e2*req
      zeq=1e2*zeq
      
      return
      end
cccccccccccccccccccccccccccccccccccccccc
      subroutine dump_bxbzbp(nrp,nzp)
      
      USE ka2jet_mod	!,only: ieqdsk,bsig
      USE spline_mod,only: pv1
      USE eqdsk_mod, only: req,zeq,psirz

      IMPLICIT NONE
      
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER k,j,ics
      INTEGER nrp,nzp
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 t_,x_,y_,phi_
      REAL*8 vdum,bdum(3),ex(3)
      REAL*8,dimension(:),allocatable::bx_,bz_,bp_,pv_
       
C============
ccccc Dump B(R,Z) values
       
       if(allocated(bx_)) deallocate(bx_)
       allocate(bx_(nzp))
       if(allocated(bz_)) deallocate(bz_)
       allocate(bz_(nzp))
       if(allocated(bp_)) deallocate(bp_)
       allocate(bp_(nzp))
       if(allocated(pv_)) deallocate(pv_)
       allocate(pv_(nzp))

       open(61,file='bxbzbp.txt',status='unknown')

       write(61,'(2i6,f12.6)') nrp,nzp,bsig

       ! Bx, Bz,Bphi fields vs (R,Z)
       write(61,'(512f12.6)') (req(j),j=1,nrp)
       write(61,'(512f12.6)') (zeq(k),k=1,nzp)
       phi_=0.d0
       t_=0.d0
       do j = 1,nrp
        do k=1,nzp       
	  x_ = 1e-2*req(j) ! [m]  
	  y_ = 1e-2*zeq(k) ! [m]
	  
	  ! B-field components
	  call bfield(t_,x_,y_,phi_,bdum,ex,ics)
	  bx_(k)=bdum(1)
	  bz_(k)=bdum(2)
	  bp_(k)=bdum(3)
	  
	  
	  
	  ! Psi
	  x_ = req(j) ! [cm]  
	  y_ = zeq(k) ! [cm]
#ifdef ISFLUSH ! JET -> use flush library  	  
	  if(ieqdsk.eq.0) then	
	   call flup(x_,y_,vdum,ics)
	  else
	   call get_flux(x_,y_,vdum,ics)
          endif
	  pv_(k) = vdum
#else
          call get_flux(x_,y_,vdum,ics)
          pv_(k) = vdum
#endif
	  
        enddo  
	  write(61,'(512f16.6)') (bx_(k),k=1,nzp) ! Bx
          write(61,'(512f16.6)') (bz_(k),k=1,nzp) ! Bz
          write(61,'(512f16.6)') (bp_(k),k=1,nzp) ! Bphi
          write(61,'(512f16.6)') (pv_(k),k=1,nzp) ! Psi
       enddo
       
       close(61)
       
       if(allocated(bx_)) deallocate(bx_)
       if(allocated(bz_)) deallocate(bz_)
       if(allocated(bp_)) deallocate(bp_)
       if(allocated(pv_)) deallocate(pv_)
       
       return
      end
ccccccccccccccccccccccccccccccc	
cccccccccccccccccccccccccccccccccccccccc
      subroutine dump_bsplines(lsp,lst)
      
      USE spline_mod
      
      IMPLICIT NONE   
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER k,j
      INTEGER lsp,lst
C============   
   
      open(61,file='bsplines.txt',status='unknown')
      write(61,'(2i6)') lsp,lst
      
      do j = 1,lsp 
       
       write(61,'(512e16.8)') (bx1(j,k),k=1,lst)
       write(61,'(512e16.8)') (bx2(j,k),k=1,lst)
       write(61,'(512e16.8)') (bx3(j,k),k=1,lst)
       write(61,'(512e16.8)') (bx4(j,k),k=1,lst)
       write(61,'(512e16.8)') (bx5(j,k),k=1,lst)
       write(61,'(512e16.8)') (bx6(j,k),k=1,lst)
       write(61,'(512e16.8)') (bx7(j,k),k=1,lst)
       write(61,'(512e16.8)') (bx8(j,k),k=1,lst)
       write(61,'(512e16.8)') (bx9(j,k),k=1,lst)
       
      enddo
      
      close(61)
      
      return
      end
cccccccccccccccccccccccccccccccccccccccc  
cccccccccccccccccccccccccccccccccccccccc
      subroutine dump_psisplines(lsp,lst)
      
      USE spline_mod
      
      IMPLICIT NONE   
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER k,j
      INTEGER lsp,lst
C============   
   
      open(61,file='psisplines.txt',status='unknown')
      write(61,'(2i6)') lsp,lst
      
      do j = 1,lsp 
       
       write(61,'(512e16.8)') (pv1(j,k),k=1,lst)
       write(61,'(512e16.8)') (pv2(j,k),k=1,lst)
       write(61,'(512e16.8)') (pv3(j,k),k=1,lst)
       write(61,'(512e16.8)') (pv4(j,k),k=1,lst)
       write(61,'(512e16.8)') (pv5(j,k),k=1,lst)
       write(61,'(512e16.8)') (pv6(j,k),k=1,lst)
       write(61,'(512e16.8)') (pv7(j,k),k=1,lst)
       write(61,'(512e16.8)') (pv8(j,k),k=1,lst)
       write(61,'(512e16.8)') (pv9(j,k),k=1,lst)
       
      enddo
      
      close(61)
      
      return
      end
cccccccccccccccccccccccccccccccccccccccc 
cccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccccc
C        SUBROUTINE GELG                                                GELG  40
C                                                                       GELG  50
C        PURPOSE                                                        GELG  60
C           TO SOLVE A GENERAL SYSTEM OF SIMULTANEOUS LINEAR EQUATIONS. GELG  70
C                                                                       GELG  80
C        USAGE                                                          GELG  90
C           CALL GELG(R,A,M,N,EPS,IER)                                  GELG 100
C                                                                       GELG 110
C        DESCRIPTION OF PARAMETERS                                      GELG 120
C           R      - THE M BY N MATRIX OF RIGHT HAND SIDES.  (DESTROYED)GELG 130
C                    ON RETURN R CONTAINS THE SOLUTION OF THE EQUATIONS.GELG 140
C           A      - THE M BY M COEFFICIENT MATRIX.  (DESTROYED)        GELG 150
C           M      - THE NUMBER OF EQUATIONS IN THE SYSTEM.             GELG 160
C           N      - THE NUMBER OF RIGHT HAND SIDE VECTORS.             GELG 170
C           EPS    - AN INPUT CONSTANT WHICH IS USED AS RELATIVE        GELG 180
C                    TOLERANCE FOR TEST ON LOSS OF SIGNIFICANCE.        GELG 190
C           IER    - RESULTING ERROR PARAMETER CODED AS FOLLOWS         GELG 200
C                    IER=0  - NO ERROR,                                 GELG 210
C                    IER=-1 - NO RESULT BECAUSE OF M LESS THAN 1 OR     GELG 220
C                             PIVOT ELEMENT AT ANY ELIMINATION STEP     GELG 230
C                             EQUAL TO 0,                               GELG 240
C                    IER=K  - WARNING DUE TO POSSIBLE LOSS OF SIGNIFI-  GELG 250
C                             CANCE INDICATED AT ELIMINATION STEP K+1,  GELG 260
C                             WHERE PIVOT ELEMENT WAS LESS THAN OR      GELG 270
C                             EQUAL TO THE INTERNAL TOLERANCE EPS TIMES GELG 280
C                             ABSOLUTELY GREATEST ELEMENT OF MATRIX A.  GELG 290
C                                                                       GELG 300
C        REMARKS                                                        GELG 310
C           INPUT MATRICES R AND A ARE ASSUMED TO BE STORED COLUMNWISE  GELG 320
C           IN M*N RESP. M*M SUCCESSIVE STORAGE LOCATIONS. ON RETURN    GELG 330
C           SOLUTION MATRIX R IS STORED COLUMNWISE TOO.                 GELG 340
C           THE PROCEDURE GIVES RESULTS IF THE NUMBER OF EQUATIONS M IS GELG 350
C           GREATER THAN 0 AND PIVOT ELEMENTS AT ALL ELIMINATION STEPS  GELG 360
C           ARE DIFFERENT FROM 0. HOWEVER WARNING IER=K - IF GIVEN -    GELG 370
C           INDICATES POSSIBLE LOSS OF SIGNIFICANCE. IN CASE OF A WELL  GELG 380
C           SCALED MATRIX A AND APPROPRIATE TOLERANCE EPS, IER=K MAY BE GELG 390
C           INTERPRETED THAT MATRIX A HAS THE RANK K. NO WARNING IS     GELG 400
C           GIVEN IN CASE M=1.                                          GELG 410
C                                                                       GELG 420
C        SUBROUTINES AND FUNCTION SUBPROGRAMS REQUIRED                  GELG 430
C           NONE                                                        GELG 440
C                                                                       GELG 450
C        METHOD                                                         GELG 460
C           SOLUTION IS DONE BY MEANS OF GAUSS-ELIMINATION WITH         GELG 470
C           COMPLETE PIVOTING.                                          GELG 480
C                                                                       GELG 490
C     ..................................................................GELG 500
C                                                                       GELG 510
      SUBROUTINE GELG(R,A,M,N,EPS,IER)
C                                                                       GELG 530
C                                                                       GELG 540
C============
C idecl:  explicitize implicit INTEGER declarations:
      IMPLICIT NONE
      INTEGER m,n,ier,mm,nm,l,i,lst,k,j,ll,lend,ii,ist
C============
C idecl:  explicitize implicit REAL declarations:
      REAL*8 a,eps,r,piv,tb,tol,pivi
C============
      DIMENSION A(1),R(1)
      IF(M)23,23,1
C                                                                       GELG 570
C     SEARCH FOR GREATEST ELEMENT IN MATRIX A                           GELG 580
    1 IER=0
      PIV=0.D0
      MM=M*M
      NM=N*M
      DO 3 L=1,MM
      TB=ABS(A(L))
      IF(TB-PIV)3,3,2
    2 PIV=TB
      I=L
    3 CONTINUE
      TOL=EPS*PIV
C     A(I) IS PIVOT ELEMENT. PIV CONTAINS THE ABSOLUTE VALUE OF A(I).   GELG 700
C                                                                       GELG 710
C                                                                       GELG 720
C     START ELIMINATION LOOP                                            GELG 730
      LST=1
      DO 17 K=1,M
C                                                                       GELG 760
C     TEST ON SINGULARITY                                               GELG 770
      IF(PIV)23,23,4
    4 IF(IER)7,5,7
    5 IF(PIV-TOL)6,6,7
    6 IER=K-1
    7 PIVI=1.D0/A(I)
      J=(I-1)/M
      I=I-J*M-K
      J=J+1-K
C     I+K IS ROW-INDEX, J+K COLUMN-INDEX OF PIVOT ELEMENT               GELG 860
C                                                                       GELG 870
C     PIVOT ROW REDUCTION AND ROW INTERCHANGE IN RIGHT HAND SIDE R      GELG 880
      DO 8 L=K,NM,M
      LL=L+I
      TB=PIVI*R(LL)
      R(LL)=R(L)
    8 R(L)=TB
C                                                                       GELG 940
C     IS ELIMINATION TERMINATED                                         GELG 950
      IF(K-M)9,18,18
C                                                                       GELG 970
C     COLUMN INTERCHANGE IN MATRIX A                                    GELG 980
    9 LEND=LST+M-K
      IF(J)12,12,10
   10 II=J*M
      DO 11 L=LST,LEND
      TB=A(L)
      LL=L+II
      A(L)=A(LL)
   11 A(LL)=TB
C                                                                       GELG1070
C     ROW INTERCHANGE AND PIVOT ROW REDUCTION IN MATRIX A               GELG1080
   12 DO 13 L=LST,MM,M
      LL=L+I
      TB=PIVI*A(LL)
      A(LL)=A(L)
   13 A(L)=TB
C                                                                       GELG1140
C     SAVE COLUMN INTERCHANGE INFORMATION                               GELG1150
      A(LST)=J
C                                                                       GELG1170
C     ELEMENT REDUCTION AND NEXT PIVOT SEARCH                           GELG1180
      PIV=0.D0
      LST=LST+1
      J=0
      DO 16 II=LST,LEND
      PIVI=-A(II)
      IST=II+M
      J=J+1
      DO 15 L=IST,MM,M
      LL=L-J
      A(L)=A(L)+PIVI*A(LL)
      TB=ABS(A(L))
      IF(TB-PIV)15,15,14
   14 PIV=TB
      I=L
   15 CONTINUE
      DO 16 L=K,NM,M
      LL=L+J
   16 R(LL)=R(LL)+PIVI*R(L)
   17 LST=LST+M
C     END OF ELIMINATION LOOP                                           GELG1380
C                                                                       GELG1390
C                                                                       GELG1400
C     BACK SUBSTITUTION AND BACK INTERCHANGE                            GELG1410
   18 IF(M-1)23,22,19
   19 IST=MM+M
      LST=M+1
      DO 21 I=2,M
      II=LST-I
      IST=IST-LST
      L=IST-M
      L=A(L)+.5D0
      DO 21 J=II,NM,M
      TB=R(J)
      LL=J
      DO 20 K=IST,MM,M
      LL=LL+1
   20 TB=TB-A(K)*R(LL)
      K=J+L
      R(J)=R(K)
   21 R(K)=TB
   22 RETURN
C                                                                       GELG1600
C                                                                       GELG1610
C     ERROR RETURN                                                      GELG1620
   23 IER=-1
      RETURN
      END
 

