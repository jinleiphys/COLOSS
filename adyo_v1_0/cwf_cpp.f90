module cwf_cpp
    implicit none
    interface
        subroutine compute_coulomb_wave_functions(is_normalized, eta, z, l, Nl, F, DF, G, DG) bind(C, name="compute_coulomb_wave_functions")
        use iso_c_binding !this should be used inside the interface!
            character(kind=c_char), intent(in) :: is_normalized(*)
            complex(kind=c_double_complex), intent(in) :: eta
            complex(kind=c_double_complex), intent(in) :: z
            complex(kind=c_double_complex), intent(in) :: l
            integer(kind=c_int), intent(in) :: Nl
            complex(kind=c_double_complex), intent(out) :: F(Nl), DF(Nl), G(Nl), DG(Nl)
        end subroutine compute_coulomb_wave_functions
    end interface
end module cwf_cpp