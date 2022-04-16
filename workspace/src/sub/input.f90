module nrtype
    implicit none
    !Symbolic names for kind types of 4-, 2-, and 1-byte integers:
    INTEGER, PARAMETER :: I4B = SELECTED_INT_KIND(9)
    INTEGER, PARAMETER :: I2B = SELECTED_INT_KIND(4)
    INTEGER, PARAMETER :: I1B = SELECTED_INT_KIND(2)
    !Symbolic names for kind types of single- and double-precision reals:
    INTEGER, PARAMETER :: SP = 4!KIND(1.0)
    INTEGER, PARAMETER :: DP = 8!KIND(1.0D0)
    !Symbolic names for kind types of single- and double-precision complex:
    INTEGER, PARAMETER :: SPC = KIND((1.0, 1.0))
    INTEGER, PARAMETER :: DPC = KIND((1.0D0, 1.0D0))
    !Symbolic name for kind type of default logical:
    INTEGER, PARAMETER :: LGT = KIND(.true.)

    integer, PARAMETER :: IP = I4B
    integer, PARAMETER :: RP = DP
end module nrtype
module input
    use nrtype
    implicit none
    character(len=99) ::m_input_file, m_data_directory, m_working_directory, m_mixture_name, m_iNewtonhistory, m_gsi_mechism
    integer(I4B) :: m_maxiteration, m_subiteration
    real(DP) :: m_pressure_init, m_temperature_init, m_distance , m_pert_m, m_pert_T, m_eps

end module input

subroutine setDefaultInput
    use input
    implicit none 
    m_input_file="input.xml"
    m_data_directory=""
    m_working_directory="/D/gitee/Mutation/Mutationpp/workspace/GSI_mutation++/data"
    m_working_directory=""
    m_mixture_name = ""
    m_pressure_init = 101325;
    m_temperature_init = 3000;
    m_distance = 1e-3
    m_maxiteration = 50
    m_subiteration = 5
    m_pert_m = 1e-2
    m_pert_T = 1.0
    m_eps = 1.0e-12
    m_iNewtonhistory = "false"
end subroutine

subroutine setInput 
    use input
    implicit none
    call setDefaultInput
    !data directory can be set automaticly
    !call mpp_loadstringfromfile(m_input_file,"data_directory",m_data_directory)
    call mpp_loadstringfromfile(m_input_file,"working_directory",m_working_directory)
    call mpp_loadstringfromfile(m_input_file,"mixture_name",m_mixture_name)
    call mpp_loadstringfromfile(m_input_file,"iNewtonhistory",m_iNewtonhistory)
    call mpp_loadintfromfile(m_input_file,"maxiterations",m_maxiteration)
    call mpp_loadintfromfile(m_input_file,"subiterations",m_subiteration)
    call mpp_loaddoublefromfile(m_input_file,"pressure_init",m_pressure_init)
    call mpp_loaddoublefromfile(m_input_file,"temperature_init",m_temperature_init)
    call mpp_loaddoublefromfile(m_input_file,"distance",m_distance)
    call mpp_loaddoublefromfile(m_input_file,"pert_m",m_pert_m)
    call mpp_loaddoublefromfile(m_input_file,"pert_T",m_pert_T)
    call mpp_loaddoublefromfile(m_input_file,"resnorm",m_eps)
    

endsubroutine
