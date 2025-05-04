program main
    use mol2
    implicit none

    integer, parameter :: pop_size = 100, max_gen = 500
    real, parameter :: crossover_rate = 0.95, mutation_rate = 0.05
    type(Individual), allocatable :: population(:), new_population(:)
    type(Atom), allocatable :: ligand_atoms(:), site_atoms(:)
    integer :: i, j, generation, n_ligand, n_site, parent1, parent2
    real :: r
    character(len=256) :: ligand_file, site_file, output_name
    integer :: unit_csv

    integer :: best_idx, p1, p2
    real :: best_fit, mean_fit, min_fit

    ! Variables for timing
    real :: t0, t1
    integer :: start_clock, end_clock, clock_rate
    real :: elapsed_time

    call random_seed()

    ligand_file = "ligand_NO_BOND.mol2"
    site_file = "site.mol2"

    call read_mol2(ligand_file, ligand_atoms, n_ligand)
    call read_mol2(site_file, site_atoms, n_site)

    allocate(population(pop_size))
    allocate(new_population(pop_size))

    do i = 1, pop_size
        call random_number(r)
        population(i)%tx = (r-0.5) * 10.0
        call random_number(r)
        population(i)%ty = (r-0.5) * 10.0
        call random_number(r)
        population(i)%tz = (r-0.5) * 10.0
        call random_number(r)
        population(i)%rx = (r-0.5) * 360.0
        call random_number(r)
        population(i)%ry = (r-0.5) * 360.0
        call random_number(r)
        population(i)%rz = (r-0.5) * 360.0
        call evaluate_individual(population(i), ligand_atoms, n_ligand, site_atoms, n_site)
    end do

    open(newunit=unit_csv, file="evolution.csv", status='replace', action='write')
    write(unit_csv,*) "Generation,Min,Max,Mean"

    !----------------------------------------------------------------
    ! Boucle principale de l'algorithme génétique (Version random sampling)
    !----------------------------------------------------------------
    ! do generation = 1, max_gen
    !     !$OMP PARALLEL DO PRIVATE(i, r, parent1, parent2)
    !     do i = 1, pop_size
    !         call random_number(r)
    !         if (r < crossover_rate) then
    !             call random_number(r)
    !             parent1 = int(r * pop_size) + 1
    !             call random_number(r)
    !             parent2 = int(r * pop_size) + 1
    !             call crossover(population(parent1), population(parent2), new_population(i))
    !         else
    !             new_population(i) = population(i)
    !             call random_number(r)
    !             if (r < mutation_rate) call mutate(new_population(i))
    !         end if
    !         call evaluate_individual(new_population(i), ligand_atoms, n_ligand, site_atoms, n_site)
    !     end do
    !     !$OMP END PARALLEL DO

    !     population = new_population

    !     if (mod(generation,10) == 0) then
    !         call save_statistics(population, pop_size, generation, unit_csv)
    !     end if
    ! end do

    !----------------------------------------------------------------
    ! Boucle principale de l'algorithme génétique (Version tournoi)
    !----------------------------------------------------------------
    call system_clock(start_clock, clock_rate)
    call cpu_time(t0)
    do generation = 1, max_gen
        ! On garde le champion (élitisme)
        call tournament_selection(population, best_idx)
        new_population(1) = population(best_idx)

        ! Génération des enfants
        !$OMP PARALLEL DO PRIVATE(i,r,p1,p2) SHARED(population,new_population)
        do i = 1, (pop_size-1)/2
            ! sélection
            call tournament_selection(population, p1)
            call tournament_selection(population, p2)

            ! crossover
            call crossover(population(p1), population(p2), new_population(2*i))
            call crossover(population(p2), population(p1), new_population(2*i+1))

            ! mutation
            call random_number(r)
            if (r < mutation_rate) call mutate(new_population(2*i))
            call random_number(r)
            if (r < mutation_rate) call mutate(new_population(2*i+1))

            ! évaluation
            call evaluate_individual(new_population(2*i), ligand_atoms, n_ligand, site_atoms, n_site)
            call evaluate_individual(new_population(2*i+1), ligand_atoms, n_ligand, site_atoms, n_site)
        end do
        !$OMP END PARALLEL DO

        ! Remplacement
        population = new_population

        ! Statistiques
        if (mod(generation,10) == 0) then
            call save_statistics(population, pop_size, generation, unit_csv)
            call population_statistics(population, min_fit, best_fit, mean_fit)
            ! print *, "Generation", generation, "Min", min_fit, "Max", best_fit, "Mean", mean_fit
            print '(A,I4,A,F8.2,A,F8.2)', "Gen=", generation, "  Best=", best_fit, "  Avg=", mean_fit
        end if
    end do
    call cpu_time(t1)
    call system_clock(end_clock, clock_rate)
    elapsed_time = real(end_clock - start_clock) / clock_rate
    print *, "Elapsed time (s) -", elapsed_time
    print *, "CPU time (s) -", t1 - t0
    print *, "Speedup -", (t1 - t0) / elapsed_time

    close(unit_csv)

    do i = 1, pop_size
        write(output_name, '("results/individual_",i3.3,".mol2")') i
        call save_individual(population(i), ligand_file, output_name)
    end do

    deallocate(population, new_population, ligand_atoms, site_atoms)

contains

    subroutine save_statistics(pop, size, gen, unit)
        type(Individual), intent(in) :: pop(:)
        integer, intent(in) :: size, gen, unit
        real :: minf, maxf, meanf
        integer :: i

        minf = pop(1)%fitness
        maxf = pop(1)%fitness
        meanf = 0.0
        do i = 1, size
            minf = min(minf, pop(i)%fitness)
            maxf = max(maxf, pop(i)%fitness)
            meanf = meanf + pop(i)%fitness
        end do
        meanf = meanf / size
        write(unit,'(I0,",",F8.3,",",F8.3,",",F8.3)') gen, minf, maxf, meanf
    end subroutine save_statistics

    subroutine tournament_selection(pop, idx_selected)
        implicit none
        type(Individual), intent(in ) :: pop(:)
        integer, intent(out) :: idx_selected
        integer :: i1, i2
        real :: r1, r2
        call random_number(r1)
        call random_number(r2)

        i1 = int(r1 * size(pop)) + 1
        i2 = int(r2 * size(pop)) + 1
        if (pop(i1)%fitness > pop(i2)%fitness) then
            idx_selected = i1
        else
            idx_selected = i2
        end if
    end subroutine tournament_selection

    subroutine population_statistics(pop, min_fit, max_fit, mean_fit)
        implicit none
        type(Individual), intent(in)  :: pop(:)
        real, intent(out) :: min_fit, max_fit, mean_fit
        integer :: i, n

        n = size(pop)
        min_fit = pop(1)%fitness
        max_fit = pop(1)%fitness
        mean_fit = 0.0
        do i = 1, n
            min_fit = min(min_fit, pop(i)%fitness)
            max_fit = max(max_fit, pop(i)%fitness)
            mean_fit = mean_fit + pop(i)%fitness
        end do
        mean_fit = mean_fit / real(n)
    end subroutine population_statistics

end program main
