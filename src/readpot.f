        module readpot

        use precision
        use constants

        implicit none

        complex*16,allocatable,dimension(:) :: Uread
        real*8,allocatable,dimension(:) :: xv
        integer :: readlenth
        real*8 :: rstep

        contains
c-----------------------------------------------------------------------
c           Read external potential
c           Format:
c           1st line: header
c           2nd line: npoints, rstep, rfirst
c           Next lines: real, imag
        subroutine extpot(filename)
            
            implicit none
            logical uu  ! logical variable, used to check whether a file exists
            integer npt,n
            integer :: ios
            integer,parameter:: kpot=50
            character*40 header
            character*7 filename
            real*8:: r,rfirst,x,y1,y2
            real*8,allocatable::faux(:),haux(:)
            real*8,parameter:: alpha=0

            uu = .false.
            inquire(file=filename,exist=uu)  ! check if file exists
            if (.not.uu) then
                write(0,*) 'Potential file:', filename,' not found!'
                stop
            endif

            write(*,'(2x, "Reading from Potential file:",a20)') filename
            open(kpot,file=filename,status='old')
            rewind(kpot)

            read(kpot,*) header  ! read headers
            npt=0
            do
                read(kpot, *, iostat=ios) x, y1, y2
                if (ios /= 0) exit  ! 当读取到文件的末尾或遇到错误时，退出循环
                npt = npt + 1
            end do
            if (npt < 1) then
                write(*,*)'error reading ext. potential !'
                stop
            end if

            rewind(kpot)
            read(kpot,*) header
            read(kpot,*) npt,rstep,rfirst
            
            readlenth = npt
            npt = npt - 1
            allocate(faux(0:npt),haux(0:npt),xv(0:npt))
            allocate(Uread(0:npt))
            faux=0d0; haux(:)=0d0; xv(:)=0d0; Uread(:)=0d0
            do n=0,npt
                read(kpot,*)  y1,y2
                xv(n)=rfirst+rstep*n
                faux(n)=y1
                haux(n)=y2
                Uread(n) = complex(y1,y2)
            enddo

            write(*,'(3x,"=> read:",i4," points")') npt+1
            write(*,250) xv(0),xv(npt),xv(2)-xv(1)
250         format(/,2x,"[Radial grid: Rmin=",1f7.3," fm,", " Rmax=",1f7.1," fm,", " Step=",1f7.3," fm]",/)

            deallocate(faux,haux)


        end subroutine

        end module readpot

