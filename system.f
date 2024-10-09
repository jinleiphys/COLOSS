        module system
c       this module define all the variable of the 3B system
            character(len=5) :: namep, namet
            real*8 :: zp, massp !define the charge, mass and spin of a
            real*8 :: zt, masst !define the charge, mass and spin of A
            real*8 :: sp

            real*8 :: elab
            real*8 :: ecm !energy of the reaction system
            real*8 :: mu
            complex*16 :: eta ! the sommerfeld parameter of the 2B system
            real*8 :: k !the wave vector

            integer :: jmin
            integer :: jmax
        end module