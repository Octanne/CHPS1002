program genetic_docking
  use, intrinsic :: iso_fortran_env, only: wp => real64
  implicit none
  integer, parameter :: n_atoms = 20, pop_size = 100, n_generations = 500
  integer :: gen, i, j, parent1, parent2, child
  real(wp), dimension(n_atoms, 3, pop_size) :: population
  real(wp), dimension(n_atoms, 3, pop_size) :: new_population
  real(wp), dimension(pop_size) :: fitness
  real(wp), dimension(n_atoms, 3) :: tmp
  character(len=32) :: filename
  real(wp) :: best, worst, avg
  integer :: unit_csv

  call random_seed()

  ! === Initialisation population ===
  do i = 1, pop_size
    do j = 1, n_atoms
      call random_number(population(j, :, i))
      population(j, :, i) = population(j, :, i) * 10.0_wp ! Coord aléatoires [0,10]
    end do
  end do

  ! === Ouvrir CSV pour suivi des scores ===
  open(newunit=unit_csv, file="results.csv", status="replace")
  write(unit_csv, *) "Generation,Min,Max,Mean"

  do gen = 1, n_generations
    ! === Évaluer la population ===
    do i = 1, pop_size
      fitness(i) = evaluate_dummy(population(:, :, i))
    end do

    ! === Suivi des statistiques ===
    best = minval(fitness)
    worst = maxval(fitness)
    avg = sum(fitness) / pop_size
    if (mod(gen, 10) == 0) write(unit_csv, *) gen, best, worst, avg

    ! === Nouvelle génération ===
    do child = 1, pop_size
      parent1 = tournament_selection(fitness)
      parent2 = tournament_selection(fitness)
      tmp = crossover(population(:, :, parent1), population(:, :, parent2))
      call mutate(tmp)
      new_population(:, :, child) = tmp
    end do

    ! === Remplacer population actuelle ===
    population = new_population
  end do

  close(unit_csv)

  ! === Écrire les individus finaux ===
  do i = 1, pop_size
    write(filename, '("individual_", I3.3, ".mol2")') i
    call write_dummy_mol2(population(:, :, i), filename)
  end do

  print *, "Docking terminé. Résultats dans results.csv et fichiers MOL2."

contains

  function evaluate_dummy(conf) result(score)
    real(wp), dimension(n_atoms, 3), intent(in) :: conf
    real(wp) :: score
    ! Score fictif : plus le centre est proche de (5,5,5), meilleur c'est
    real(wp), dimension(3) :: center
    center = sum(conf, dim=1) / n_atoms
    score = -norm2(center - [5.0_wp, 5.0_wp, 5.0_wp]) ! Plus c'est proche, mieux c'est
  end function

  function tournament_selection(fitness) result(idx)
    real(wp), dimension(pop_size), intent(in) :: fitness
    integer :: idx, i1, i2
    call random_number(i1); i1 = int(i1 * pop_size) + 1
    call random_number(i2); i2 = int(i2 * pop_size) + 1
    if (fitness(i1) > fitness(i2)) then
      idx = i1
    else
      idx = i2
    end if
  end function

  function crossover(p1, p2) result(child)
    real(wp), dimension(n_atoms, 3), intent(in) :: p1, p2
    real(wp), dimension(n_atoms, 3) :: child
    child = 0.5_wp * (p1 + p2)
  end function

  subroutine mutate(conf)
    real(wp), dimension(n_atoms, 3), intent(inout) :: conf
    real(wp), dimension(3) :: shift
    call random_number(shift)
    shift = (shift - 0.5_wp) * 2.0_wp ! [-1,1] aléatoire
    conf = conf + shift
  end subroutine

  subroutine write_dummy_mol2(conf, filename)
    real(wp), dimension(n_atoms, 3), intent(in) :: conf
    character(len=*), intent(in) :: filename
    integer :: i, unit_id

    open(newunit=unit_id, file=filename, status="replace")
    write(unit_id, '(A)') "@<TRIPOS>MOLECULE"
    write(unit_id, '(A)') "Ligand"
    write(unit_id, '(I5,1x,I5)') n_atoms, 0
    write(unit_id, '(A)') "SMALL"
    write(unit_id, '(A)') "NO_CHARGES"
    write(unit_id, '(A)') "@<TRIPOS>ATOM"
    do i = 1, n_atoms
      write(unit_id, '(I5,1X,A4,3F10.4,1X,A3,1X,I2)') i, "C", conf(i,1), conf(i,2), conf(i,3), "LIG", 1
    end do
    close(unit_id)
  end subroutine

end program genetic_docking
