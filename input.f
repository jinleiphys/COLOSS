        module INPUT
            use mesh
            use system
            use precision
            use constants
            use pot_class
            use rot_pot

            implicit none
           
            contains
                subroutine  read_input()

                namelist /general/ ctheta,alpha,nr,rmax,numgauss,rmaxgauss,backrot,
     &                                 matgauss,thetah,thetamax,
     &                                 bgauss,cwftype,method

                namelist /system/ zp,massp,namep,zt,masst,namet,elab,jmax,jmin,sp

                namelist /pot/
     &                                  vv,rv,av,wv,rw,aw,
     &                                  vs,rvs,avs,ws,rws,aws,
     &                                  vsov,rsov,asov,vsow,rsow,asow,
     &                                  rc
                namelist /nonlocal_pot/ nonlocal,nonlocal_beta


                write(*,*)'************************************************************'
                write(*,*)'*                  COLOSS: Complex-scaled                  *'
                write(*,*)'*          Optical and couLOmb Scattering Solver           *'
                write(*,*)'************************************************************'

c------------------general namelist-------------------------------------
                ctheta = 5.d0; nr = 20; alpha = 0
                rmaxgauss = 200.d0; numgauss = 200
                matgauss = .false.
                bgauss = .false.
                thetah=1.d0; thetamax=90.d0
                cwftype=1
                backrot=.false.
                method=1
            
                read(5, nml=general)

                theta_rad = ctheta/180.d0*pi!convert degree to rad
                eitheta = exp(iu*theta_rad)
                theta_n_max = int(thetamax/thetah)

                
                write(*,*)''
                write(*,100) ctheta
100             format(' Rotation Angle for Complex Scaling:',f6.2,' degrees')
                if(backrot) then
                        write(*,*)'Rotation Operation: Applied to Lagrange-Laguerre basis'
                else
                        write(*,*)'Rotation Operation: Applied to Potential'
                end if
                write(*,*)''

c----------------end of general namelist---------------------------------   


c------------------system namelist--------------------------------------
                massp = 0.d0; masst = 0.d0
                zp = 0.d0; zt = 0.d0
                elab = 0.d0
                sp = 0d0
                jmin = 0
                jmax = 0

                read(5, nml=system)

                ecm = masst*elab/(masst + massp)
                mu = amu * masst*massp/(massp + masst)
                k=sqrt( 2.*mu*ecm/(hbarc**2) )
                !k=sqrt( 0.5d0*ecm/20.736d0 )
                eta=zt*zp*e2*mu/hbarc/hbarc/k

                write(*,*) '-------------------Reaction system-------------------'
                write(*,200) ' Projectile: ',namep,massp,zp
                write(*,200) ' Target: ',namet,masst,zt
200             FORMAT(a15,a6,' (A: ', f8.2,' Z:', f8.2,')')
                write(*,202) elab
202             FORMAT(' Lab  Frame Energy:', f8.2, ' MeV')
                write(*,203) ecm
203             FORMAT(' C.M. Frame Energy:', f8.2, ' MeV')
                write(*,204) jmin,jmax
204             FORMAT(' Total J Range:', I3,' <= J <= ',I3)


                write(*,*) ''


c-----------------end of system namelist-------------------------------- 


c------------------potential namelist----------------------------------
                vv=0d0; rv=0d0; av=0d0
                wv=0d0; rw=0d0; aw=0d0
                vs=0d0; rvs=0d0; avs=0d0
                ws=0d0; rws=0d0; aws=0d0
                vsov=0d0; rsov=0d0; asov=0d0
                vsow=0d0; rsow=0d0; asow=0d0
                rc=0d0

c
                read(5,nml=pot)
            
c------------------end of potential namelist---------------------------

c------------------nonlocal_pot namelist------------------------------

                nonlocal = .false.
                nonlocal_beta = 0.d0

                read(5,nml=nonlocal_pot)

c------------------end of nonlocal_pot namelist-----------------------

                write(1,nml=general)
                write(1,nml=system)
                write(1,nml=pot)
                write(1,nml=nonlocal_pot) 


            end subroutine

            subroutine outinfo()
            character*45 :: outfile(100)
            logical :: writeouot(100)
            integer :: writf(100)
            integer :: i
            integer :: nwrit

            writeouot(:) = .false.
            outfile(:) = ' '

            outfile(1) = 'local copy of input'
            writeouot(1) = .TRUE.
            outfile(10) = 'scaled Laguerre mesh points and weights'
            writeouot(10) = .TRUE.
            outfile(60) = 'cross section LSJ distribution'
            writeouot(60) = .TRUE.
            outfile(61) = 'scat amplitude LSJ distribution'
            writeouot(61) = .TRUE.
            outfile(67) = 'cross section angular distribution'
            writeouot(67) = .TRUE.
            write(*,*) ''
            write(*,989)
            nwrit = 0
            do i=1,100
              if(writeouot(i)) then
                outfile(i) = trim(outfile(i))//'.'
                nwrit = nwrit+1
                writf(nwrit) = i
              endif
            enddo
            write(*,990) (writf(i),outfile(writf(i)),i=1,nwrit)
989         FORMAT('----------------- Files Created -----------------')
990         FORMAT((i3,': ',a45))

            end subroutine


            subroutine write_verison()
c               This subroutine is used to print out the compling information
        

                write(*,*) "----------------------------------------------------------"
                write(*,*) "The following is some info for compilation and git version"
                write(*,*) "----------------------------------------------------------"
#ifdef BASE
        print *, 'Base directory: ', BASE
#endif
            
#ifdef VERDATE
        print *, 'Version date: ', VERDATE
#endif
            
#ifdef VERREV
        print *, 'Version revision: ', VERREV
#endif
            
#ifdef COMPDATE
        print *, 'Compilation date: ', COMPDATE
#endif 
                end subroutine

        end module input 