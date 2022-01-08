program main
    use mutationpp
    use input
    implicit none
    integer(kind=4) :: i, j, ne, ns, nt, neq
    integer (IP) :: set_state_with_rhoi_T, pos_T_trans

    real(kind=8) :: T, P, rho, E, mblow, solid_heat
    character(len=12), allocatable, DIMENSION(:) :: species_name
    real(RP), dimension(:), allocatable :: rhoi_s, wdot, xi_s, xi_e, mwi, dxidx, vdi, v_hi, v_h
    real(RP), dimension(:), allocatable :: T_e, T_s, dtdx, lambda, q_srad, q_gcon, v_hi_rhoi_vi, Residual

    set_state_with_rhoi_T = 1
    pos_T_trans = 1

    call setInput

    call MPP_globaloptionsworkingdir(m_working_directory)
    call mpp_initializegsi(m_mixture_name)

    !set solver options
    call mpp_setiterationssurfacebalance(m_maxstep)
    call mpp_setiterationspert_m(m_pert_m)
    call mpp_setiterationspert_t(m_pert_T)
    call mpp_setiterationseps(m_tol)
    call mpp_setiterationshistory(m_iNewtonhistory)
    call mpp_get_gsi_mechism(m_gsi_mechism)

    ns = mpp_nspecies()
    nt = mpp_n_energy_eqns()
    ne = mpp_nelements()
    neq = ns + nt

    allocate (xi_e(ns), rhoi_s(ns), xi_s(ns), dxidx(ns), vdi(ns), wdot(ns), Residual(ns+1))
    allocate (T_e(nt), T_s(nt))
    if (m_gsi_mechism == "phenomenological_mass_energy") then
        allocate (dtdx(nt), lambda(nt), q_srad(nt), q_gcon(nt))
    end if
    allocate (v_hi(ns*nt), v_h(nt), v_hi_rhoi_vi(nt))
    allocate (species_name(ns))

    !obtain species name list
    do j = 1, ns
        call mpp_species_name(j, species_name(j))
    end do

    !set diffusion model use equilibrium composition of the mixture with given T,P
    call mpp_equilibrate(m_temperature_init, m_pressure_init)
    call mpp_x(xi_e)
    !boundary edge composition can be given by input
    !call mpp_get_species_composition("Gas",xi_e)
    call mpp_set_diffusion_model(xi_e, m_distance)

    !set heat gradient if needed
    T_e = m_temperature_init
    if (m_gsi_mechism == "phenomenological_mass_energy") then
        call mpp_set_cond_heat_flux(T_e, m_distance)
    end if

    !Initial conditions of the surface are the ones in the first physical cell
    call mpp_species_densities(rhoi_s)
    T_s = m_temperature_init/2.0
    call mpp_set_surface_state(rhoi_s, T_s, set_state_with_rhoi_T)

    ! Solve balance and request solution
    call mpp_solve_surface_balance()
    call mpp_get_surface_state(rhoi_s, T_s, set_state_with_rhoi_T)
    rho = sum(rhoi_s)

    ! Verifying the solution gives low residual in the balance equations
    call mpp_set_state(rhoi_s, T_s, set_state_with_rhoi_T)
    call mpp_x(xi_s)

    ! compute diffusion velocities
    dxidx = (xi_s - xi_e)/m_distance; 
    E = 0.0
    call mpp_stefan_maxwell(dxidx, vdi, E)

    ! compute heat flux
    if (m_gsi_mechism == "phenomenological_mass_energy") then
        !dtdx = (T_s - T_e)/m_distance
        !call mpp_frozen_thermal_conductivity(lambda)
        !q_gcon = lambda*dtdx
        !gas heat flux can be calculated by (-lambda*dtdx=q_gcon>0) or function below, just for validation
        call mpp_compute_gas_heat_flux(T_s,q_gcon)
    end if

    ! Get surface production rates
    call mpp_set_surface_state(rhoi_s, T_s, set_state_with_rhoi_T)
    call mpp_surface_production_rates(wdot)

    ! Blowing flux (should be zero for catalysis)
    call mpp_mass_blowing_rate(mblow)

    ! Species and mixture enthalpies
    call mpp_species_h_mass(v_hi)
    v_h(pos_T_trans) = DOT_PRODUCT(rhoi_s(1:ns)/rho, v_hi(1:ns))

    ! Surface radiation
    if (m_gsi_mechism == "phenomenological_mass_energy") then
        q_srad = 0.0
        q_srad(pos_T_trans) = mpp_getsurfaceradiativeheatflux()
    end if

    ! Chemical Energy Contribution
    v_hi_rhoi_vi(pos_T_trans) = dot_product((v_hi(1:ns)), (rhoi_s*vdi))

    ! solid conduction
    solid_heat=mblow*mpp_get_solid_heat()
    ! MixtureOptions graphiteopt("graphite.xml");
    ! Mixture graphite(graphiteopt);
    ! const int set_state_PT = 1;
    ! graphite.setState(&P_init,&T_s[0],1);
    ! double hcp=graphite.mixtureHMass(T_s[0]);

    ! Building balance functions
    Residual(1:ns) = (rhoi_s/rho)*mblow + (rhoi_s*vdi) - wdot; 
    if (m_gsi_mechism == "phenomenological_mass_energy") then
        Residual(ns + 1) = sum(q_gcon) - q_srad(pos_T_trans) + mblow*v_h(pos_T_trans) + v_hi_rhoi_vi(pos_T_trans) + solid_heat
    end if

    ! write results
    write (*, *) "The residual of gas surface interaction:"
    write (*, "(4A18)") "Res.AVG(mass)", "Res.MAX(mass)", "Res.AVG(energy)", "Res.MAX(energy)"
    write (*, "(4E18.7)") sum(Residual(1:ns))/ns, maxval(Residual(1:ns)), sum(Residual(ns + 1:ns+1)), maxval(Residual(ns + 1:ns+1))

    write (*, *) "Surface balance temperature[K]"
    do i = 1, nt
        write (*, '(f18.5)', advance='no') T_s(i)
    end do
    write (*, *)

    write (*, *) "Surface blowing flux[kg/m^2-s]"
    write (*, "(f18.7)") mblow

    write (*, "(A30)", advance='no') "Surface mass properties"
    do i = 1, ns
        write (*, '(5X,A13)', advance='no') species_name(i)
    end do
    write (*, *)

    ! species density
    write (*, "(A30)", advance='no') "Species density[kg/m^3]"
    do i = 1, ns
        write (*, "(E18.6)", advance='no') rhoi_s(i)
    end do
    write (*, *)

    ! chemical source
    write (*, "(A30)", advance='no') "Chemical production[kg/m^2-s]"
    do i = 1, ns
        write (*, "(E18.6)", advance='no') wdot(i)
    end do
    write (*, *)

    write (*, "(A30)", advance='no') "Surface species mole fraction"
    do i = 1, ns
        write (*, "(E18.6)", advance='no') xi_s(i)
    end do
    write (*, *)

    write (*, "(A30)", advance='no') "Edge species mole fraction"
    do i = 1, ns
        write (*, "(E18.6)", advance='no') xi_e(i)
    end do
    write (*, *)

    write (*, *) "Surface energy properties"
    write (*, "(5A25)") "Conductive heat[J]", "Radiative heat[J]", "Blowing heat[J]", "Diffusion heat[J]", "Solid heat(J)"
    write (*, "(5E25.7)", advance='no') sum(q_gcon), -q_srad(pos_T_trans) , &
    mblow*v_h(pos_T_trans), v_hi_rhoi_vi(pos_T_trans), solid_heat
    write (*, *)

    deallocate (xi_e, rhoi_s, xi_s, dxidx, vdi, wdot, Residual)
    deallocate (T_e, T_s)
    if (m_gsi_mechism == "phenomenological_mass_energy") then
        deallocate (dtdx, lambda, q_srad, q_gcon)
    end if
    deallocate (v_hi, v_h, v_hi_rhoi_vi)
    deallocate (species_name)
    call mpp_destroy()
101 format(1X, F18.7)

end program
