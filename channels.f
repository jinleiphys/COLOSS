        module channels

        use system
        use mesh 
        use precision
        use constants

        IMPLICIT NONE


        type :: channel2B

            integer :: ch_numbers

            integer, allocatable, dimension(:) :: L 
            real*8, allocatable, dimension(:) :: S,J

        end type channel2B

        type(channel2B) :: channel_index

        contains

        subroutine generate_channels()

            integer :: il,iJ,ich

            write(*,*) "Generating the channel index"
            
            lmin = 0


            lmax = nint(jmax - sp) 

            channel_index%ch_numbers = 0
            do il = lmin, lmax
                do iJ = nint(2 * abs(il - sp)), nint(2 * (il + sp)), 2
                    if(iJ .lt. 2*jmin) cycle
                    if(iJ .gt. 2*jmax) cycle
                    channel_index%ch_numbers = channel_index%ch_numbers + 1
                end do
            end do

            allocate(channel_index%L(1:channel_index%ch_numbers))
            allocate(channel_index%S(1:channel_index%ch_numbers))
            allocate(channel_index%J(1:channel_index%ch_numbers))

            write(11, '(A5, 2X, A5, 2X, A5, 2X, A5)') 'ich', 'L', 'S', 'J'
            ich = 1
            do il = lmin, lmax
                do iJ = nint(2 * abs(il - sp)), nint(2 * (il + sp)), 2
                    channel_index%L(ich) = il
                    channel_index%S(ich) = sp
                    channel_index%J(ich) = iJ/2d0
                    write(*,1000) ich, channel_index%L(ich), channel_index%S(ich), channel_index%J(ich)
                    write(11,1000) ich, channel_index%L(ich), channel_index%S(ich), channel_index%J(ich)
                    ich = ich + 1
                end do
            end do

            write(*,204) lmin,lmax
            write(*,205) sp, (lmax+sp)
204         FORMAT('       L Range:', I5,' <= L <= ',I5)
205         FORMAT(' Total J Range:', f5.1,' <= J <= ',f5.1)
            write(*,*)


1000  format (I5, 2X, I5, 2X, F6.1, 2X, F6.1)

        end subroutine

            
        end module