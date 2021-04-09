cccccccccccccccccccccccccccccccccccccc
      subroutine read_plasma_state(froot)
      
      USE plasma_state_mod
      USE fi_constants_mod
      USE fi_ps_mod
      
      IMPLICIT NONE
      INTRINSIC sqrt,trim,abs
      !include 'fitrcom'
      
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER ierr,ids(3),ncid,mid,k
C============
C idecl:  explicitize implicit REAL declarations:
      REAL(KIND=rspec) buffer(3),temp,vdum,tdum,adum,psi_ax
      REAL(KIND=rspec) dum1,dum2,dum3
C============      
C idecl:  explicitize implicit CHARACTER declarations:      
      CHARACTER*96 froot
C============
ccc   This routine reads the information from the 
ccc   Plasma State.
      
      write(0,3) froot
3     format('	-> Open file ',A) 
      call ps_get_plasma_state(ierr,filename=trim(froot),state=ps)
      call ckerr(ierr,'ps_get_plasma_state')

      write(0,*) '	  File open, reading data...'
      if(tstart.lt.0) tstart=ps%t0  ! else take it from fi_transp.f input
      if(tstop.lt.0)  tstop =ps%t1  ! else take it from fi_transp.f input
      
      raxis=ps%R_axis ! [cm]
      zaxis=ps%Z_axis ! [cm]
      baxis=ps%B_axis ! [T]
      PsiRz=ps%id_PsiRZ ! flux
      BphiRz=ps%id_BphiRZ
      BtotRz2=ps%id_BphiRZ**2+ps%id_BRRZ**2+ps%id_BZRZ**2
      
      Rmin=ps%R_min_lcfs
      Rmax=ps%R_max_lcfs
      
      write(0,*) ''
      write(0,*) '------------------------------'
      write(0,11) trim(ps%tokamak_id),ps%shot_number,tstart,tstop
11    format(A,', shot=',i6', Macro step: ',f7.4,' -',f7.4,' s')
      write(0,*) ''            
      write(0,'(A,i5,A,i5,A)') ' Grid (R,z): ',
     >            ps%nR,' x ',ps%nZ,' points'
      write(0,'(A,i5,A,i5,A)') ' Grid (rho,th): ',
     >            ps%nrho_eq,' x ',ps%nth_eq,' points'
      write(0,'(A,i5)') 'nrho=',ps%nrho
      write(0,'(A,f6.3,A,f6.3,A)') ' Axis at (',
     >                             raxis,',',zaxis,') cm'
      write(0,'(A,f6.3,A)')' |B| = ',baxis,' [T] on axis.'
      write(0,*) ''
      
      
      ! NOTE: pitch is already in the correct form if taken from 
      !		inside NUBEAM. There's no need to check & switch sign. MP, May 2015
      ptchsign=1.D0
      
      
c      Precompute values needed for orbit classification
c	NOTE: following section mostly taken from MCGEN10.F90 by D. Liu
      ids(1) = ps%id_BRRZ
      ids(2) = ps%id_BZRZ
      ids(3) = ps%id_BphiRZ
      
      ! ion charge & mass
      write(0,'(A,f6.4)') ' Species 1 mass [1e-27 kg] = ',
     >                    ionmass
      write(0,'(A,f6.4)') ' Species 1 charge [1e-19 C] = ',
     >                    ioncharge
      mass2charge = ionmass/ioncharge
	  
      !B at magnetic axis
      call ps_intrp_2d(ps%R_axis, ps%z_axis, ids, buffer, ierr)
      call ckerr(ierr,'ps_intrp_2d')
      signB=buffer(1)/abs(buffer(1)) ! sign of B at magnetic axis
      Baxis = sqrt(buffer(1)**2 + buffer(2)**2 + buffer(3)**2)
      Brataxis = buffer(3) / Baxis
      write(0,'(A,f6.3,A,f6.3)') ' |Baxis|,ratio = ',Baxis,',',Brataxis

      !Get B at right wall
      call ps_intrp_2d(ps%R_max_lcfs,ps%z_axis,ids,buffer,ierr)
      call ckerr(ierr,'ps_intrp_2d')
      Bmin = sqrt(buffer(1)**2 + buffer(2)**2 + buffer(3)**2)
      Bratrw = buffer(3) / Bmin
      write(0,'(A,f6.3,A,f6.3)') ' |Bmin|,ratio = ',Bmin,',', Bratrw

      !Get B at left wall
      call ps_intrp_2d(ps%R_min_lcfs,ps%z_axis,ids,buffer,ierr)
      call ckerr(ierr,'ps_intrp_2d')
      Bmax = sqrt(buffer(1)**2 + buffer(2)**2 + buffer(3)**2)
      Bratlw = buffer(3) / Bmax
      write(0,'(A,f6.3,A,f6.3)') ' |Bmax|,ratio = ',Bmax,',', Bratlw

      call ps_intrp_2d(ps%R_axis,ps%z_axis,ps%id_BphiRZ,temp,ierr)
      call ckerr(ierr,'ps_intrp_2d')
      write(0,'(A,f6.3)') ' B_phi_axis = ',temp

      call ps_intrp_2d(ps%R_axis, ps%z_axis, 
     >                 ps%id_PsiRZ,psi_axis,ierr)
      call ckerr(ierr,'ps_intrp_2d')
      if(psi_axis.lt.0) then
        write(0,'(A,f12.9,A)')
     > 		'	*** Warning: Psi_axis=',psi_axis, 
     >		' will be reset to Psi_axis=0' 
	psi_axis=0.D0 ! avoid "floating invalid" exceptions
      endif
      write(0,'(A,f12.10)')  ' psi_axis = ',psi_axis
      
      ! Get max flux (at lcfs)
      call ps_intrp_2d(ps%R_min_lcfs,ps%z_axis,
     >                 ps%id_PsiRZ, temp, ierr)
      call ckerr(ierr,'ps_intrp_2d')
      call ps_intrp_2d(ps%R_max_lcfs,ps%z_axis,
     >                 ps%id_PsiRZ, psi_lcfs, ierr)
      call ckerr(ierr,'ps_intrp_2d')
      if (temp.GT.psi_lcfs) psi_lcfs = temp
      write(0,'(A,f8.6)') ' psi_lcfs = ',psi_lcfs
      write(0,'(A,f8.6)') ' psi_lcfs_2 = ',ps%psipol(ps%nrho_eq)
      
      ! normalize Psi @ magnetic axis
      ! [although Psi_axis should be =0]
      psi_axis=psi_axis/psi_lcfs
      
      ! Flux at the wall
      ! NOTE: Psi_w is normalized to B_axis, as it is
      !		usually done in ORBIT - MP, May 2015
      psiw = psi_lcfs/abs(Baxis)
      write(0,'(A,f12.10)') ' Normalized Psi_wall = ',psiw
      write(0,*) ''


      ! Define normalization factors for E, B, R, rho, ...
      Rnorm=ps%R_axis	! [m]
      Bnorm=1.D0	! [T], consistent with ORBIT [???]
      Enorm=1D-27*ionmass*9D16*1D-3*
     c           (Bnorm/(m2cf*mass2charge))**2*Rnorm**2 ! something treatable
      write(0,*) 'Normalizations:'
      write(0,'(A,f8.6)') ' R_norm = ',Rnorm
      write(0,'(A,f8.6)') ' B_norm = ',Bnorm
      write(0,'(A,f8.3)') ' E_norm = ',Enorm
      write(0,*) ''
      
      
c     fill in an array with the (Psi,theta) -> (R,z)
c     rough mapping - used later on to speed up
c     the computation of particle's position
c     based on its phase space coordinates.
      call psith2Rz_table
       
      
c     Fill in a table with the max Rho vs. poloidal angle theta.
      call lcfs_table
      
c     Init array of radial Psi_pol zones
      call init_psizones      
      
      return
      end

cccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccc
      subroutine read_amode
      
      USE plasma_state_mod
      USE fi_constants_mod
      USE fi_ps_mod
      USE fi_mode_mod
      USE fi_pdedp_mod
      USE ifport
      
      IMPLICIT NONE
      
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,k,j2
      INTEGER*4 npt,npt_max,length
      INTEGER(KIND=INT_PTR_KIND( )) handle
C============
C idecl:  explicitize implicit REAL declarations:
      REAL(KIND=rspec) vdum,tdum,adum
      REAL(KIND=rspec) ascale(npdedp_max)
      CHARACTER*12 basename,sufx
      CHARACTER*50 fname
C============
C idecl:  explicitize implicit TYPE declarations:      
      TYPE (FILE$INFO) info      
C============
ccc   This routine reads the mode amplitude data from
ccc   text file(s) in Ufile format.

      ! define common name and suffix for mode amplitude file(s)
      basename='FILEAMODE*'
      sufx='.AEP'
      npdedp=0
      
      ! allocate temp variable
      allocate(npt_mode(npdedp_max))
      npt_mode=0 ! init
      
c     Check for existing file(s) determine max number of
c     time points -> size of mode info arrays
      
      do k=1,npdedp_max ! loop to look for files
        
	handle = FILE$FIRST  
	if(k.le.9) then
	 write(fname,'(A,i1,A)') trim(basename),k,trim(sufx)
	else
	 write(fname,'(A,i2,A)') trim(basename),k,trim(sufx)
	endif
	length = GETFILEINFOQQ(fname,info,handle)
	
      if (length.ge.len(trim(fname))) then ! valid file?
	
c        check for errors	
         IF (handle.eq.FILE$ERROR) THEN
           SELECT CASE (GETLASTERRORQQ( ))
            CASE (ERR$NOMEM)
              WRITE (*,*) 'Out of memory'
            CASE (ERR$NOENT)
              EXIT
            CASE DEFAULT
              WRITE (*,*) 'Invalid file or path name'
           END SELECT
         END IF
	
	
      open(88,file=info%name,status='unknown') 

c     read header (strings)
      do j=1,3
         read(88,*)
      enddo
      
c     read scaling factor
      read(88,*) ascale(k)
      
c     read remaining header (strings)
      do j=1,3
         read(88,*)
      enddo
      
c     read number of time points
      read(88,*) npt
      
      if(npt.gt.idm) then
        write(0,'(A,i7,3A)') '*** Found ',npt,' time points for',
     >      ' mode amplitude data in ',flabel
        write(0,'(A,i7)')
     >      '	-> will clip array to maximum allowed: ',idm
	npt=idm
      endif
      
      ! update number of time points for this mode.
      ! Max number will be used for array allocation
      npt_mode(k)=npt
      
        ! update number of files
        npdedp = npdedp + 1
        
	! close file
	close(88)
	
       endif ! valid file
	     
      enddo ! loop over possible files
      
      ! update number of valid files & max number
      ! of required time points
      npt_max=maxval(npt_mode)
      
c     allocate arrays for this file
      deallocate(npt_mode)
      allocate(npt_mode(npdedp)) ! re-allocate
      allocate(tmode(npt_max,npdedp))
      allocate(amode(npt_max,npdedp))
      allocate(pbaept(npt_max,10)) ! second dimension is 'orbit type'

      
c     Re-open valid files, read data into Amode, tmode arrays   
      do k=1,npdedp ! loop over valid Amode files
        
	handle = FILE$FIRST   
	if(k.le.9) then
	 write(fname,'(A,i1,A)') trim(basename),k,trim(sufx)
	else
	 write(fname,'(A,i2,A)') trim(basename),k,trim(sufx)
	endif
	length = GETFILEINFOQQ(fname,info,handle)
	open(88,file=info%name,status='unknown') 
	!rewind(88) ! make sure it reads from the beginning
        write(0,'(2A)') '	Reading data from ',info%name
	
c     read header (strings)
      do j=1,3
         read(88,*)
      enddo
      
c     read scaling factor
      read(88,*) ascale(k)
      
c     read remaining header (strings)
      do j=1,3
         read(88,*)
      enddo
      
c     read number of time points
      read(88,*) npt
      npt_mode(k)=npt ! update
       
c     read time base
      read(88,*) (tmode(j,k), j=1,npt)
      
c     read amplitude data
      read(88,*) (amode(j,k), j=1,npt)

      write(0,11) npt
11    format(' Number of time points for mode amplitude:',i8)
      
c     scale Amode
      write(0,'(A,f8.4)') ' Scaling factor: ',ascale(k)
      do j=1,npt
        amode(j,k)=ascale(k)*amode(j,k)
      enddo      

      write(0,15) tmode(1,k),tmode(npt,k)
15    format(' Time range: ',f8.2,' - ',f8.2,' ms')      
      write(0,*)
      
c     close file      
      close(88)
      
c     ! init PBAEPT array
      do j=1,npt
        do j2=1,idp
          pbaept(j,j2)=0.D0
	enddo
      enddo
           
      enddo ! loop over possible files

      return
      end
cccccccccccccccccccccccccccccccccccccc


cccccccccccccccccccccccccccccccccccccc
      subroutine read_pdedp
      
      USE plasma_state_mod
      USE fi_constants_mod
      USE fi_ps_mod
      USE fi_pdedp_mod
      USE fi_mode_mod
      USE ifport
       
      IMPLICIT NONE
      
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER j,kf			
      INTEGER je,jp,jmu,jde,jdp
      INTEGER length
      INTEGER(KIND=INT_PTR_KIND( )) handle
C============
C idecl:  explicitize implicit REAL declarations:
      REAL(KIND=rspec) vdum,tdum,adum,prdum
      CHARACTER*12 basename,sufx
      CHARACTER*50 fname
C============
C idecl:  explicitize implicit TYPE declarations:      
      TYPE (FILE$INFO) info     
C============
ccc-this routine reads the probability distribution function data
ccc- for P(DE,DP|E,P,mu) from a text file (Ufile)
    
      ! define common name and suffix for p(DE,DP) file(s)
      basename='FILEPDEDP*'
      sufx='.AEP'

c     NOTE:
c     At this point, the number of valid files has been
c     already determined, see subroutine read_amode.

      
c     Allocate arrays.      
      allocate(nptde(npdedp))
      allocate(nptdp(npdedp))
      allocate(npte(npdedp))
      allocate(nptp(npdedp))
      allocate(nptmu(npdedp)) 
      
      allocate(pdestep(npdedp))
      allocate(pdpstep(npdedp))
      allocate(pestep(npdedp))
      allocate(ppstep(npdedp))
      
      allocate(pvar_dt(npdedp))    
      allocate(pvar_de(ibins,npdedp))
      allocate(pvar_dp(ibins,npdedp))
      allocate(pvar_e(ibins,npdedp))
      allocate(pvar_p(ibins,npdedp))
      allocate(pvar_mu(ibins,npdedp))
      
      allocate(pdedp(ibins,ibins,ibins,ibins,ibins,npdedp))
      allocate(p_sumdedp(ibins,ibins,ibins,npdedp))
      
      
      
c     Read data into p(DE,DP) arrays
      do kf=1,npdedp 
        
	handle = FILE$FIRST   
	if(kf.le.9) then
	 write(fname,'(A,i1,A)') trim(basename),kf,trim(sufx)
	else
	 write(fname,'(A,i2,A)') trim(basename),kf,trim(sufx)
	endif
	length = GETFILEINFOQQ(fname,info,handle)
	
c        check for errors	
         IF (handle.eq.FILE$ERROR) THEN
           SELECT CASE (GETLASTERRORQQ( ))
            CASE (ERR$NOMEM)
              WRITE (*,*) 'Out of memory'
            CASE (ERR$NOENT)
              EXIT
            CASE DEFAULT
              WRITE (*,*) 'Invalid file or path name'
           END SELECT
         END IF
	
	 
      write(0,'(2A)') '	-> Opening p(DE,DP) file: ',info%name
      open(77,file=info%name,status='unknown') 
      
c     read header (strings)
      do 5 j=1,9
         read(77,*)
5     continue
      
c     read time step used to produce p(DE,DP)
      read(77,*) pvar_dt(kf)
      write(0,'(A,f8.5)')
     >     ' Time step for p(DE,DP) [ms]: ',pvar_dt(kf)
            
c     continue reading header
      read(77,*)
      read(77,*)
      
      
c     read number of bins points for DE,DP,E,P,mu
      read(77,*) nptde(kf)
      read(77,*) nptdp(kf)
      read(77,*) npte(kf)
      read(77,*) nptp(kf)
      read(77,*) nptmu(kf)

            
      write(0,*) 'Number of points for p(DE,DP|E,P,mu) function:'
      write(0,*) '      E       P      mu      DE      DP'
      write(0,22) npte(kf),nptp(kf),nptmu(kf),nptde(kf),nptdp(kf)
22    format(5i8)


c     Read bins.
c     The arrays for E,P,mu indicate the center value
c     of each bin.
      read(77,23) (pvar_de(jde,kf),jde=1,nptde(kf)) ! DE
      read(77,23) (pvar_dp(jdp,kf),jdp=1,nptdp(kf)) ! DP
      read(77,23) (pvar_e(je,kf) , je =1,npte(kf)) ! E
      read(77,23) (pvar_p(jp,kf) , jp =1,nptp(kf)) ! P
      read(77,23) (pvar_mu(jmu,kf),jmu=1,nptmu(kf)) ! mu
      
      
c23    format(6e13.6) ! to be used with some old pDEDP files, before Nov. 2013
23    format(6e14.6)
      
      write(0,*) 'Binning information:'
      write(0,24) '	DE : ',pvar_de(1,kf),' - ',pvar_de(nptde(kf),kf)
      write(0,24) '	DPz: ',pvar_dp(1,kf),' - ',pvar_dp(nptdp(kf),kf)
      write(0,24) '	 E : ',pvar_e(1,kf), ' - ',pvar_e(npte(kf),kf)
      write(0,24) '	 Pz: ',pvar_p(1,kf), ' - ',pvar_p(nptp(kf),kf)
      write(0,24) '	 mu: ',pvar_mu(1,kf),' - ',pvar_mu(nptmu(kf),kf)
    
24    format(A,1e14.6,A,1e14.6)    
      			      
c     Multiple [weird...] loops to read the 5D distribution p(DE,DP|E,P,mu).
c
c     In the innermost loop(s), pre-compute the total probability function P(E,P,mu),
c     i.e. sum P over the variables (DE,DP). This will be used to skip a step if the probability
c     for DE and DP kicks is identically zero outside the central (DE,DP)=(0,0) bin.
		      
      write(0,*) 'Reading p(DE,DP) data...'
      
      do 59 je=1,npte(kf)
        do 58 jp=1,nptp(kf)
          do 57 jmu=1,nptmu(kf)
            
	    prdum = 0.D0
            do 56 jde=1,nptde(kf)
              read(77,*) (pdedp(jde,jdp,je,jp,jmu,kf),jdp=1,nptdp(kf))
	      do jdp=1,nptdp(kf)
	        prdum = prdum+pdedp(jde,jdp,je,jp,jmu,kf)
	      enddo  
56          continue ! jde  

            p_sumdedp(je,jp,jmu,kf) = 0.D0
            do jde=1,nptde(kf)
	      p_sumdedp(je,jp,jmu,kf) = p_sumdedp(je,jp,jmu,kf) +
     >		      prdum
            enddo
      
57        continue ! jmu      
58      continue ! jp      
59    continue ! je      

c     compute step size in DE, DP, E, P
      pdestep(kf)=(pvar_de(nptde(kf),kf)-pvar_de(1,kf))/(nptde(kf)-1.)
      pdpstep(kf)=(pvar_dp(nptdp(kf),kf)-pvar_dp(1,kf))/(nptdp(kf)-1.)
      pestep(kf)=(pvar_e(npte(kf),kf)-pvar_e(1,kf))/(npte(kf)-1.)
      ppstep(kf)=(pvar_p(nptp(kf),kf)-pvar_p(1,kf))/(nptp(kf)-1.)
      
c     close file      
      close(88)
      
      write(0,*) '... done.'
      write(0,*) ''
            
      enddo ! loop over valid pDEDP files


      return
      end
cccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccc
      subroutine read_fbm(files,sufx)
     
      USE plasma_state_mod
      USE fi_constants_mod
      USE fi_ps_mod
      USE fi_fbm_mod
      USE ifport
      USE netcdf
      
      IMPLICIT NONE
      INTRINSIC sqrt,trim,abs,tan,atan,modulo,dsign
      INTRINSIC sin,cos,min,max
      !include 'fitrcom'
C============
C idecl:  explicitize implicit PARAMETER declarations [used for cgs units]:      
      REAL(KIND=rspec),parameter::ZC=2.9979D10	! SPEED OF LIGHT (CM/SEC)
      REAL(KIND=rspec),parameter::ZEL=4.8032D-10 ! ELECTRON CHARGE (STATCOULOMBS)
c      REAL(KIND=rspec),parameter::ZELONZC=4.8032D0/2.9979D0 ! ZEL/ZC
c      REAL(KIND=rspec),parameter::ZMP=1.6726D-24 ! PROTON MASS (GRAMS)
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER nfiles,j,ncid,ierr,istart(2),istop(2)
      INTEGER vid,nspec,npt,ndims,dimids(10),vtype,natts,dimid
      INTEGER(KIND=INT_PTR_KIND( )) handle
      INTEGER*4 k,icnt,length,ivalid
      INTEGER ier,cnt_ax,kdum
C============
C idecl:  explicitize implicit REAL declarations:
      REAL(KIND=rspec) ndum(2),pdum(1e4,2),mudum(1e4,2),pphdum(1e4,2)
      REAL(KIND=rspec) vdum(1e4,2),rdum(1e4,2),zdum(1e4,2),wdum(1e4,2)
      REAL(KIND=rspec) xdum(1e4,2),thdum(1e4,2),tdepdum(1e4,2)
      REAL(KIND=rspec) buffer,bdum,rhodum,thdumv,gdum,bpar_out(2)
      REAL(KIND=rspec) dum1,dum2,dum3,dum4,dum5
      REAL(KIND=rspec) rnd,rrec,zrec,polrec,nptcr,psq,gout,psi_wall
      REAL(KIND=rspec) knorm,eprime
C============
C idecl:  explicitize implicit CHARACTER declarations:
      CHARACTER*80 fname
      CHARACTER*80 files,sufx
      CHARACTER*20 vname,str_dum
C============
C idecl:  explicitize implicit TYPE declarations:      
      TYPE (FILE$INFO) info
C============
ccc-this routine reads the fast ion distribution function
            
      istart(1)=1
      istart(2)=1
      istop(1) =1
      istop(2) =1
      
      ivalid=0 ! counter for 'valid' particles in Fbm file(s)
      cnt_ax=0
      Rmachine=0.D0
      Zmachine=0.D0
      
c     First, loop over Fbm files to infer
c     the number of valid particles and
c     initialize various arrays
      icnt = 0
      npt_tot=0
      handle = FILE$FIRST
      DO WHILE (.TRUE.)
        length = GETFILEINFOQQ(files, info, handle)
c       check for errors	
        IF ((handle .EQ. FILE$LAST).OR. 
     >      (handle .EQ. FILE$ERROR)) THEN
          SELECT CASE (GETLASTERRORQQ( ))
           CASE (ERR$NOMEM)
             WRITE (*,*) 'Out of memory'
           CASE (ERR$NOENT)
             EXIT
           CASE DEFAULT
             WRITE (*,*) 'Invalid file or path name'
          END SELECT
        END IF
        
	ierr=nf90_open(info%name,NF90_NOWRITE,ncid)
	IF (ierr.EQ.NF90_NOERR) THEN ! no errors, proceed
	
c	  get variable ID - number of particles
          ierr = nf90_inq_varid(ncid,'minb',vid)          
	  
c         get number of particles in this file
          ierr = nf90_get_var(ncid,vid,npt) 
	
	  ! close file
	  ierr=nf90_close(ncid)
	  
c         update total number of particles in Fbm	  
	  npt_tot = npt_tot + npt
	  
	ELSE
	 write(0,'(2A)')   '	*** Error reading data from ',info%name	 
	ENDIF
      ENDDO      

      
c     Allocate arrays for Fbm
c     and initialize arrays to 0
      write(0,'(A,i6,A)')
     >         '-> Allocating arrays for ',npt_tot,' particles'         
     
      if(allocated(rfbm)) deallocate(rfbm)
      allocate(rfbm(npt_tot))
      rfbm=0.
      if(allocated(zfbm)) deallocate(zfbm)
      allocate(zfbm(npt_tot))
      zfbm=0.
      if(allocated(efbm)) deallocate(efbm)
      allocate(efbm(npt_tot))
      efbm=0.
      if(allocated(pphfbm)) deallocate(pphfbm)
      allocate(pphfbm(npt_tot))
      pphfbm=0.
      if(allocated(mufbm)) deallocate(mufbm)
      allocate(mufbm(npt_tot))
      mufbm=0.
      if(allocated(pfbm)) deallocate(pfbm)
      allocate(pfbm(npt_tot))
      pfbm=0.
      if(allocated(wfbm)) deallocate(wfbm)
      allocate(wfbm(npt_tot))
      wfbm=0.
      if(allocated(polfbm)) deallocate(polfbm)
      allocate(polfbm(npt_tot))
      polfbm=0.
      if(allocated(thfbm)) deallocate(thfbm)
      allocate(thfbm(npt_tot))
      thfbm=0.
      if(allocated(bfbm)) deallocate(bfbm)
      allocate(bfbm(npt_tot))
      bfbm=0.
      if(allocated(tdepfbm)) deallocate(tdepfbm)
      allocate(tdepfbm(npt_tot))
      tdepfbm=0.
      if(allocated(p_skew)) deallocate(p_skew)
      allocate(p_skew(npt_tot))   
      p_skew=0. 
      if(allocated(newvalfbm)) deallocate(newvalfbm)
      allocate(newvalfbm(npt_tot))    
      newvalfbm=0
      if(allocated(otp)) deallocate(otp)
      allocate(otp(npt_tot))   
      otp=0
	 
      write(0,'(A)') '    - Arrays allocated'
         
      write(0,*) '-> Read particle distribution from Fbm file(s)'
      nfiles=0
      npt_tot=0 ! reset
c     check for existing file(s) and read data into Fbm arrays
      handle = FILE$FIRST
      DO WHILE (.TRUE.)
        length = GETFILEINFOQQ(files, info, handle)
c       check for errors	
        IF ((handle .EQ. FILE$LAST).OR. 
     >      (handle .EQ. FILE$ERROR)) THEN
          SELECT CASE (GETLASTERRORQQ( ))
           CASE (ERR$NOMEM)
             WRITE (*,*) 'Out of memory'
           CASE (ERR$NOENT)
             EXIT
           CASE DEFAULT
             WRITE (*,*) 'Invalid file or path name'
          END SELECT
        END IF
	
        nfiles = nfiles + 1 ! update number of files
	
	write(0,'(2A)') '	Reading data from ',info%name
	
c       now read data from this file
        ierr=nf90_open(info%name,NF90_NOWRITE,ncid)
	IF (ierr.EQ.NF90_NOERR) THEN ! no errors, proceed
         
c	  get variable ID - number of particles
          ierr = nf90_inq_varid(ncid,'minb',vid)          
c         get number of valid particles in this file
          ierr = nf90_get_var(ncid,vid,npt) 
	  
c         update number of particles for this file	  
	  istop(1)=npt
         
c	  ---------------------------------------------------
c	  read variables	

 	  ! pitch	  
	  ierr = nf90_inq_varid(ncid,'xksidy',vid)
	  ierr = nf90_get_var(ncid,vid,pdum,istart,istop)
	  
	  ! velocity	  
	  ierr = nf90_inq_varid(ncid,'vay',vid)
	  ierr = nf90_get_var(ncid,vid,vdum,istart,istop)
	  
	  ! mu	  
	  ierr = nf90_inq_varid(ncid,'xmuay',vid)
	  ierr = nf90_get_var(ncid,vid,mudum,istart,istop)
	  
	  ! P_phi
	  !ierr = nf90_inq_varid(ncid,'pphiay',vid)
	  ierr = nf90_inq_varid(ncid,'pmechay',vid)
	  ierr = nf90_get_var(ncid,vid,pphdum,istart,istop)
	  
	  ! R
	  ierr = nf90_inq_varid(ncid,'rmjionay',vid)
	  ierr = nf90_get_var(ncid,vid,rdum,istart,istop)
	  
	  ! x
	  ierr = nf90_inq_varid(ncid,'xiay',vid)
	  ierr = nf90_get_var(ncid,vid,xdum,istart,istop)
	  
	  ! z
	  ierr = nf90_inq_varid(ncid,'xzionay',vid)
	  ierr = nf90_get_var(ncid,vid,zdum,istart,istop)
	  
	  ! theta
	  ierr = nf90_inq_varid(ncid,'thay',vid)
	  ierr = nf90_get_var(ncid,vid,thdum,istart,istop)
	  
	  ! deposition time [wrt t0] in [s]
	  ierr = nf90_inq_varid(ncid,'tdepay',vid)
	  ierr = nf90_get_var(ncid,vid,tdepdum,istart,istop)
	  
	  ! particle's weight
	  ierr = nf90_inq_varid(ncid,'wghtay',vid)
	  ierr = nf90_get_var(ncid,vid,wdum,istart,istop)
	  
	  ierr=nf90_close(ncid)
	  
	  
	  ! store in variables
	  do 222 k=icnt+1,icnt+npt
             
	     if(rdum(k-icnt,1).gt.0.and.vdum(k-icnt,1).gt.0.and.
     >          tdepdum(k-icnt,1).lt.tstop) then ! valid particle

	       ! update counter
	       ivalid = ivalid + 1
	       
c	       Start updating Fbm arrays
	       
	       wfbm(k)=wdum(k-icnt,1)/1.D15	! some treatable number
	       rfbm(k)=rdum(k-icnt,1)/1.D2	! R [m]
	       zfbm(k)=zdum(k-icnt,1)/1.D2	! z [m]	     
	     
	       ! need to re-compute theta from (R,z) to account for
	       ! the actual position of the [magnetic] axis
	       thfbm(k)=atan2(zfbm(k)-zaxis,rfbm(k)-Raxis) ! theta=[-pi...pi]
	     	     
	       ! get machine axis (Rmachine,Zmachine)
	       if (cnt_ax.lt.1e4.and.abs(thdum(k-icnt,1)).ge.0.1) then
	         cnt_ax=cnt_ax+1
	         dum1=sqrt((rfbm(k)-Raxis)**2+(zfbm(k)-zaxis)**2)
	         dum2=sin(thfbm(k))/tan(thdum(k-icnt,1))
	         dum3=rfbm(k)-dum1*dum2
	         Rmachine=Rmachine+dum3
	         Zmachine=0.D0 ! approximation
	       endif
	      
	       ! initialize flag NEWVALFBM, which at output indicates whether particle's
	       ! variables have been modified by this script (=1) or not (=0)
	       newvalfbm(k)=0
	     
	       ! initialize deposition time for current particle
	       tdepfbm(k)=0.
	     
	       ! update deposition time for current particle
	       tdepfbm(k)=1.e3*tdepdum(k-icnt,1) ! [ms]
	       
	       ! velocity
	       vdum(k-icnt,1)=vdum(k-icnt,1)/1.D2 ! [cm/s] -> [m/s]
	       
	       pfbm(k)=pdum(k-icnt,1)	! pitch
	       efbm(k)=.5D-3*m2cf*mass2charge*vdum(k-icnt,1)*vdum(k-icnt,1) ! energy [keV]
	       
	       
	       ! get [poloidal] flux coordinate, Psi
	       call ps_intrp_2d(rfbm(k),zfbm(k),ps%id_PsiRZ,buffer,ierr)
	       polfbm(k)=buffer/psi_lcfs ! pol. flux coord. normalized to Psi @ LCFS
	       
	       ! get value of 'g' function at this flux surface
	       call gfun(polfbm(k),gout)
	       
	       ! get magnetic field at particle's location       
	       call bfieldRZ(rfbm(k),zfbm(k),bdum)		! from (R,Z)
c	       call bfield(polfbm(k),thfbm(k),bdum,bpar_out)	! from (psi,theta)
	       
	       
	       ! -----------------------------------------------------
	       ! compute magnetic moment, mu=E/B*(1-p^2) 
	       ! Note - MP, Jan. 2016: definition of mu based
	       !           on the kinetic energy=total-potential.
	       bfbm(k)=bdum	
	       mufbm(k)=efbm(k)/bdum*(1.D0-pfbm(k)*pfbm(k)) ! un-normalized definition
	       !mufbm(k)=Baxis/bfbm(k)*(1.D0-pfbm(k)*pfbm(k)) ! mu'=mu Bo/E
	       ! NOTE:  the value used for p(DE,DP) mu-axis is the
	       !	normalized mu, which is:
	       !  mu_prime=Baxis/bfbm(k)*(1.D0-pfbm(k)*pfbm(k))=mufbm(k)*Bo/efbm(k)
	       
	       ! -----------------------------------------------------
	       ! use ORBIT definition for P_zeta
	       
c	       knorm0=2.D0*1D6*(ZMP/ZQEL)/(1D2*10.D0*Baxis)**2
	       knorm=(Atnumb*1D19*ZMP/ZELONZC)/(Zcharge*Baxis)**2
	       eprime=knorm*efbm(k) ! energy in normalized units
	       rhodum=pfbm(k)*sqrt(2.D0*eprime)*Baxis/bfbm(k) 
               
	       pphfbm(k)=-polfbm(k)+
     >              (rhodum*gout/Baxis)/psiw	! toroidal angular momentum
	       ! NOTE:  In the above definition of P_phi, the value Psi_w of Psi at the LCFS
	       !	is normalized to the value of B_0 to be consistent with ORBIT.
	       
	       ! Note - MP, Jan. 2016: now move to total energy from kinetic energy.
	       call phifun(polfbm(k),dum1) ! get e.s. potential from Psi
	       efbm(k) = efbm(k) + dum1
	       
	     endif ! valid particle
	     
222       continue

	  icnt = icnt + npt
	  
c         update total number of particles in Fbm	  
	  npt_tot = icnt
	  
	endif
	IF (ierr.NE.NF90_NOERR) THEN ! errors, stop execution
	  write(0,*) trim(nf90_strerror(ierr))
          stop "Stopped"
        ENDIF
	
      END DO
c      write(0,*) ' '
      write(0,'(i5,A)') nfiles,' Fbm data file(s) found'
      write(0,'(A,i6)')
     >         '	- max no. of particles=',npt_tot
      write(0,'(A,i6)')
     >         '	- no. of valid particles=',ivalid
      write(0,*) ' '
      
      ! update value of machine axis
      if (cnt_ax.gt.0) then
        Rmachine=Rmachine/cnt_ax
	write(0,*) ''
	write(0,'(A,f6.3,A)') '	Machine axis recomputed at R=',
     >                          Rmachine,' m'
	write(0,*) ''
      endif
	 
      return
      end

cccccccccccccccccccccccccccccccccccccc

cccccccccccccccccccccccccccccccccccccc
      subroutine read_fbm_orb
     
      USE plasma_state_mod
      USE fi_constants_mod
      USE fi_ps_mod
      USE fi_fbm_mod
      USE ifport
      
      IMPLICIT NONE
      INTRINSIC sqrt,trim,abs,tan,atan,modulo,dsign
      INTRINSIC sin,cos,min,max
      !include 'fitrcom'
C============
C idecl:  explicitize implicit PARAMETER declarations [used for cgs units]:      
      REAL(KIND=rspec),parameter::ZC=2.9979D10	! SPEED OF LIGHT (CM/SEC)
      REAL(KIND=rspec),parameter::ZEL=4.8032D-10 ! ELECTRON CHARGE (STATCOULOMBS)
      REAL(KIND=rspec),parameter::ZQEL=1.6022D-19 ! ELECTRON CHARGE (COULOMBS)
c      REAL(KIND=rspec),parameter::ZELONZC=4.8032D0/2.9979D0 ! ZEL/ZC
c      REAL(KIND=rspec),parameter::ZMP=1.6726D-24 ! PROTON MASS (GRAMS)
      
C============
C idecl:  explicitize implicit INTEGER declarations:
      INTEGER nfiles,j,ncid,ierr,istart(2),istop(2)
      INTEGER vid,nspec,npt,ndims,dimids(10),vtype,natts,dimid
      INTEGER(KIND=INT_PTR_KIND( )) handle
      INTEGER*4 k,icnt,length,ivalid
      INTEGER ier,cnt_ax,ndum2
C============
C idecl:  explicitize implicit REAL declarations:
      REAL(KIND=rspec) ndum,pdum,mudum,pphdum,poldum
      REAL(KIND=rspec) vdum,rdum,zdum,edum,ptdum
      REAL(KIND=rspec) xdum,thdum,tdepdum,knorm,eprime,bdum0,knorm0
      REAL(KIND=rspec) buffer,bdum,rhodum,thdumv,gdum,bpar_out(2)
      REAL(KIND=rspec) dum1,dum2,dum3,dum4,dum5
      REAL(KIND=rspec) rnd,rrec,zrec,polrec,nptcr,psq,gout,psi_wall
C============
C idecl:  explicitize implicit CHARACTER declarations:
      CHARACTER*80 fname
      CHARACTER*80 files,sufx
      CHARACTER*20 vname,str_dum
C============
C idecl:  explicitize implicit TYPE declarations:      
      
C============
ccc-this routine reads the fast ion distribution function
      
c      files="fbm_dist_orb.dat"
      files="fbm_q30_orb.dat"
      
      istart(1)=1
      istart(2)=1
      istop(1) =1
      istop(2) =1
      
      ivalid=0 ! counter for 'valid' particles in Fbm file(s)
      cnt_ax=0
      Rmachine=0.D0
      Zmachine=0.D0
      
      
      write(0,*) '-> Read Fbm file from ORBIT'
      nfiles=1
      icnt = 0

	write(0,'(2A)') '	Reading data from ',files
	
	open(88,file=files)	!,status='unkmown')
	
	! read number of lines to skip and number of particles
	read(88,*) ndum2,npt
	
	! skip header line(s)
	do k=1,ndum2
	  read(88,*)
	enddo
	
	write(0,*) ndum2,npt
	
	
c         update number of particles for this file	  
	  istop(1)=npt
         
	  
	  ! read & store variable for all particles in the file
	  do 222 k=icnt+1,icnt+npt
	  
	     read(88,'(8e12.4)') poldum,thdum,rdum,zdum,edum, 
     >                  ptdum,pphdum,mudum
	     
             rfbm(k)=rdum/1.D2	! R [m]
	     zfbm(k)=zdum/1.D2	! z [m]	     
	     thfbm(k)=atan2(zfbm(k)-zaxis,rfbm(k)-Raxis) ! theta=[-pi...pi]
	          
	     ! get machine axis (Rmachine,Zmachine)
	     if (cnt_ax.lt.1e4.and.abs(thfbm(k)).ge.0.1) then
	       cnt_ax=cnt_ax+1
	       dum1=sqrt((rfbm(k)-Raxis)**2+(zfbm(k)-zaxis)**2)
	       dum2=sin(thfbm(k))/tan(thfbm(k))
	       dum3=rfbm(k)-dum1*dum2
	       Rmachine=Rmachine+dum3
	       Zmachine=0.D0 ! approximation
	     endif
	      
	     ! initialize flag NEWVALFBM, which at output indicates whether particle's
	     ! variables have been modified by this script (=1) or not (=0)
	     newvalfbm(k)=0
	     
	     ! initialize deposition time for current particle
	     tdepfbm(k)=0.
	     
	     !if(rfbm(k).ge.Rmin.and.rfbm(k).le.Rmax) then ! ok, particle inside plasma
	     if(rfbm(k).gt.0) then ! ok, valid particle
               
	       ! update deposition time for current particle
	       tdepfbm(k)=0.D0 ! [ms]
	       
	       pfbm(k)=ptdum	! pitch
	       efbm(k)=edum	! energy [keV]
	       vdum=ptdum*sqrt(2.D0*1.D3*edum/(m2cf*mass2charge)) ! velocity [m/s]
	       
	       
	       ! update counter for 'valid' particles
	       if(rfbm(k-icnt).gt.0.and.vdum.gt.0) then
	         ivalid=ivalid+1
	       endif
	       
	       
	       ! get [poloidal] flux coordinate, Psi
	       call ps_intrp_2d(rfbm(k),zfbm(k),ps%id_PsiRZ,buffer,ierr)
	       polfbm(k)=buffer/psi_lcfs ! pol. flux coord. normalized to Psi @ LCFS
	       
	       ! get value of 'g' function at this flux surface
	       call gfun(polfbm(k),gout)
	       
	       ! get magnetic field at particle's location       
	       call bfieldRZ(rfbm(k),zfbm(k),bdum)		! from (R,Z)
	       
	       ! compute magnetic moment, mu=E/B*(1-p^2) 
	       bfbm(k)=bdum	!*baxis
	       mufbm(k)=efbm(k)/bdum*(1-pfbm(k)*pfbm(k)) ! un-normalized definition
	       !mufbm(k)=Baxis/bfbm(k)*(1.D0-pfbm(k)*pfbm(k)) ! mu'=mu Bo/E
	       ! NOTE:  the value used for p(DE,DP) mu-axis is the
	       !	normalized mu, which is:
	       !mu_prime=Baxis/bfbm(k)*(1.D0-pfbm(k)*pfbm(k))=mufbm(k)*Bo/efbm(k)
	       
	       ! -----------------------------------------------------
	       ! use ORBIT definition for P_zeta
	       
c	       knorm0=2.D0*1D6*(ZMP/ZQEL)/(1D2*10.D0*Baxis)**2
	       knorm=Atnumb*(1D19*ZMP/ZELONZC)/(Zcharge*Baxis)**2
	       eprime=knorm*efbm(k) ! energy in normalized units
	       rhodum=pfbm(k)*sqrt(2.D0*eprime)*Baxis/bfbm(k) 
c	       rhodum=m2cf*mass2charge*pfbm(k)*abs(vdum)*
c     >	              Baxis/bfbm(k)		!*sqrt(Bnorm*ps%R_axis)
               
	       pphfbm(k)=-polfbm(k)+
     >              (rhodum*gout/Baxis)/psiw	! toroidal angular momentum
	       ! NOTE:  In the above definition of P_phi, the value Psi_w of Psi at the LCFS
	       !	is normalized to the value of B_0 to be consistent with ORBIT.
	       
	       
	       ! -----------------------------------------------------
	       
	      
	       ! DEBUG
	       if (k.le.0) then
	         write(0,'(4f12.8)') rdum/1.e2,zdum/1.e2,ptdum,edum/1.e3
		 write(0,'(4f12.8)') rfbm(k),zfbm(k),pfbm(k),efbm(k)
		 write(0,*) ''
	       endif
	       
	       ! DEBUG
	       if (k.le.0) then
	         write(0,'(12f12.8)') pphdum+polfbm(k),mudum,
     >      ptdum	!,gout/Baxis,knorm0,Baxis,bfbm(k)
		 write(0,'(12f12.8)') pphfbm(k)+polfbm(k),
     >      mufbm(k),pfbm(k),eprime,gout/Baxis,rhodum,psi_lcfs
		 write(0,*) ''
	       endif
	       
	       
	       ! DEBUG:
	       if (k.eq.-1) then
	         write(0,'(5f12.6)') efbm(k),pphfbm(k),
     >                   mufbm(k),thfbm(k),bfbm(k)
                 call orb_update(k,polfbm(k),ierr)
		 write(0,*) ierr
		 write(0,*)
	       endif
	       
	       
	       ! DEBUG - check (Psi,thet) -> (R,z) conversion
	       if (k.le.0) then
		 
	         call psithet2RZ(polfbm(k),thfbm(k),dum1,dum2)
                 write(0,'(i6,4f12.8)') k,rfbm(k),dum1,zfbm(k),dum2
                 write(0,*) ''
 
	       endif
	     
	     !else 
	     !    write(0,'(i6,2f12.2)') k,vdum(k-icnt,1),rfbm(k)
		   
	     endif
	     
222       continue
	  
	  icnt = icnt + npt
	  
c         update total number of particles in Fbm	  
	  npt_tot = icnt

      write(0,'(i5,A)') nfiles,' Fbm data file(s) found'
      write(0,'(A,i6)')
     >         '	- max no. of particles=',npt_tot
      write(0,'(A,i6)')
     >         '	- no. of valid particles=',ivalid
      
      write(0,*) ' '
      
      ! update value of machine axis
      if (cnt_ax.gt.0) then
        Rmachine=Rmachine/cnt_ax
	write(0,*) ''
	write(0,'(A,f6.3,A)') '	Machine axis recomputed at R=',
     >                          Rmachine,' m'
	write(0,*) ''
      endif
	 
      return
      end

cccccccccccccccccccccccccccccccccccccc
cccccccccccccccccccccccccccccccccccccc
      SUBROUTINE ckerr(ierr,sbrtn)
      integer, intent(in) :: ierr
      character*(*), intent(in) :: sbrtn

      if(ierr.NE.0) then
        write(6,*) ' *** Plasma State, error in call: '//trim(sbrtn)
        stop
      endif
      end
cccccccccccccccccccccccccccccccccccccc


