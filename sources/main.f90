program main
    use mol2
    implicit none

    integer, parameter :: pop_size = 100, max_gen = 500
    real, parameter :: crossover_rate = 0.95, mutation_rate = 0.05
    type(Individual), allocatable :: population(:), new_population(:)
    type(Atom), allocatable :: ligand_atoms(:), site_atoms(:)
    integer :: idx_individual, idx_generation, n_ligand, n_site
    integer :: idx_parent1, idx_parent2
    real :: random_value
    character(len=256) :: ligand_file, site_file, output_name
    integer :: unit_csv

    integer :: best_idx, idx_tournament1, idx_tournament2
    real :: best_fit, mean_fit, min_fit

    ! Variables for timing
    real :: start_time, end_time
    integer :: start_clock, end_clock, clock_rate
    real :: elapsed_time

    call random_seed()

    ligand_file = "data/ligand_NO_BOND.mol2"
    site_file = "data/site.mol2"

    call read_mol2(ligand_file, ligand_atoms, n_ligand)
    call read_mol2(site_file, site_atoms, n_site)

    allocate(population(pop_size))
    allocate(new_population(pop_size))

    do idx_individual = 1, pop_size
        call random_number(random_value)
        population(idx_individual)%tx = (random_value-0.5) * 10.0
        call random_number(random_value)
        population(idx_individual)%ty = (random_value-0.5) * 10.0
        call random_number(random_value)
        population(idx_individual)%tz = (random_value-0.5) * 10.0
        call random_number(random_value)
        population(idx_individual)%rx = (random_value-0.5) * 360.0
        call random_number(random_value)
        population(idx_individual)%ry = (random_value-0.5) * 360.0
        call random_number(random_value)
        population(idx_individual)%rz = (random_value-0.5) * 360.0
        call evaluate_individual(population(idx_individual), ligand_atoms, n_ligand, site_atoms, n_site)
    end do

    open(newunit=unit_csv, file="evolution.csv", status='replace', action='write')
    write(unit_csv,*) "Generation,Min,Max,Mean"

    call system_clock(start_clock, clock_rate)
    call cpu_time(start_time)
    do idx_generation = 1, max_gen
        call tournament_selection(population, best_idx)
        new_population(1) = population(best_idx)

        !$OMP PARALLEL DO PRIVATE(idx_individual,random_value,idx_tournament1,idx_tournament2) SHARED(population,new_population)
        do idx_individual = 1, (pop_size-1)/2
            call tournament_selection(population, idx_tournament1)
            call tournament_selection(population, idx_tournament2)

            call crossover(population(idx_tournament1), population(idx_tournament2), new_population(2*idx_individual))
            call crossover(population(idx_tournament2), population(idx_tournament1), new_population(2*idx_individual+1))

            call random_number(random_value)
            if (random_value < mutation_rate) call mutate(new_population(2*idx_individual))
            call random_number(random_value)
            if (random_value < mutation_rate) call mutate(new_population(2*idx_individual+1))

            call evaluate_individual(new_population(2*idx_individual), ligand_atoms, n_ligand, site_atoms, n_site)
            call evaluate_individual(new_population(2*idx_individual+1), ligand_atoms, n_ligand, site_atoms, n_site)
        end do
        !$OMP END PARALLEL DO

        population = new_population

        if (mod(idx_generation,10) == 0) then
            call save_statistics(population, pop_size, idx_generation, unit_csv)
            call population_statistics(population, min_fit, best_fit, mean_fit)
            print '(A,I4,A,F8.2,A,F8.2)', "Gen=", idx_generation, "  Best=", best_fit, "  Avg=", mean_fit
        end if
    end do
    call cpu_time(end_time)
    call system_clock(end_clock, clock_rate)
    elapsed_time = real(end_clock - start_clock) / clock_rate
    print *, "Real time (s) -", elapsed_time
    print *, "CPU time (s) -", end_time - start_time
    print *, "Speedup (cpu time / real time) -", (end_time - start_time) / elapsed_time

    close(unit_csv)

    do idx_individual = 1, pop_size
        write(output_name, '("results/individual_",i3.3,".mol2")') idx_individual
        call save_individual(population(idx_individual), ligand_file, output_name)
    end do

    deallocate(population, new_population, ligand_atoms, site_atoms)

contains

    subroutine save_statistics(pop, size, gen, unit)
        type(Individual), intent(in) :: pop(:)
        integer, intent(in) :: size, gen, unit
        real :: minf, maxf, meanf
        integer :: idx

        minf = pop(1)%fitness
        maxf = pop(1)%fitness
        meanf = 0.0
        do idx = 1, size
            minf = min(minf, pop(idx)%fitness)
            maxf = max(maxf, pop(idx)%fitness)
            meanf = meanf + pop(idx)%fitness
        end do
        meanf = meanf / size
        write(unit,'(I0,",",F8.3,",",F8.3,",",F8.3)') gen, minf, maxf, meanf
    end subroutine save_statistics

    subroutine tournament_selection(pop, idx_selected)
        implicit none
        type(Individual), intent(in ) :: pop(:)
        integer, intent(out) :: idx_selected
        integer :: idx1, idx2
        real :: rand1, rand2
        call random_number(rand1)
        call random_number(rand2)

        idx1 = int(rand1 * size(pop)) + 1
        idx2 = int(rand2 * size(pop)) + 1
        if (pop(idx1)%fitness > pop(idx2)%fitness) then
            idx_selected = idx1
        else
            idx_selected = idx2
        end if
    end subroutine tournament_selection

    subroutine population_statistics(pop, min_fit, max_fit, mean_fit)
        implicit none
        type(Individual), intent(in)  :: pop(:)
        real, intent(out) :: min_fit, max_fit, mean_fit
        integer :: idx, n

        n = size(pop)
        min_fit = pop(1)%fitness
        max_fit = pop(1)%fitness
        mean_fit = 0.0
        do idx = 1, n
            min_fit = min(min_fit, pop(idx)%fitness)
            max_fit = max(max_fit, pop(idx)%fitness)
            mean_fit = mean_fit + pop(idx)%fitness
        end do
        mean_fit = mean_fit / real(n)
    end subroutine population_statistics

end program main
