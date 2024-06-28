        module pot_class



            type :: pot_para

                real*8 :: a1,a2!mass for 2B system
                real*8 :: z1,z2!charge for 2B system

                real*8 :: vv,rvv,avv!real volume term
                real*8 :: wv,rw,aw!imag volume term
                real*8 :: vs,rvs,avs!real surface term
                real*8 :: ws,rws,aws!imag surface term
                real*8 :: vsov,rsov,asov!real SO term
                real*8 :: vsow,rsow,asow!imag SO term
                real*8 :: rc!coulomb interaction parameter

            end type pot_para


            type(pot_para) :: input_pot

            !variables of OMP in the input namelists
            real*8 :: vv,rv,av
            real*8 :: wv,rw,aw
            real*8 :: vs,rvs,avs
            real*8 :: ws,rws,aws
            real*8 :: rc

        end module