program genetic_docking
    use docking_utils
    implicit none

    type(atominfo), allocatable :: ligand(:), receptor(:)
    type(atomdata), dimension(118) :: covalent_radii
    type(bond), allocatable :: ligand_bonds(:), receptor_bonds(:)
    integer :: ligand_bond_count, receptor_bond_count
    character(len=255) :: ligand_file, receptor_file, covalent_file

    ! Lecture des arguments
    if (command_argument_count() /= 3) then
        print *, 'Usage: genetic_docking <ligand.mol2> <receptor.mol2> <covalent_radii.txt>'
        stop
    end if
    call get_command_argument(1, ligand_file)
    call get_command_argument(2, receptor_file)
    call get_command_argument(3, covalent_file)

    ! Chargement des rayons covalents
    covalent_radii = init_covalent_to_zero()
    covalent_radii = load_atoms_covalents(covalent_file)

    ! Chargement des fichiers mol2
    ligand = load_mol2_file(ligand_file)
    receptor = load_mol2_file(receptor_file)

    ! Association des rayons covalents aux atomes
    call add_covalent_into_atoms(ligand, covalent_radii)
    call add_covalent_into_atoms(receptor, covalent_radii)

    ! Détermination des topologies
    call determine_topology(ligand, covalent_radii, ligand_bonds, ligand_bond_count)
    call determine_topology(receptor, covalent_radii, receptor_bonds, receptor_bond_count)

    ! Détection des liaisons hydrogène
    call detect_hydrogen_bonds(ligand, receptor)

    ! === Implémentation de l'algorithme génétique ===
    integer :: recombination_rate, mutation_rate
    real :: hbond_score

    recombination_rate = 95  ! Taux de recombinaison en pourcentage
    mutation_rate = 5        ! Taux de mutation en pourcentage

    do gen = 1, n_generations
        ! === Évaluer la population (parallélisé avec OpenMP) ===
        !$omp parallel do private(i, hbond_score) shared(population, fitness, receptor)
        do i = 1, pop_size
            hbond_score = evaluate_hbond(population(:, :, i), receptor)
            fitness(i) = hbond_score
        end do
        !$omp end parallel do

        ! === Suivi des statistiques ===
        best = minval(fitness)
        worst = maxval(fitness)
        avg = sum(fitness) / pop_size
        if (mod(gen, 10) == 0) write(unit_csv, *) gen, best, worst, avg

        ! === Nouvelle génération ===
        do child = 1, pop_size
            parent1 = tournament_selection(fitness)
            parent2 = tournament_selection(fitness)
            if (random_number() * 100 <= recombination_rate) then
                tmp = crossover(population(:, :, parent1), population(:, :, parent2))
            else
                tmp = population(:, :, parent1)
            end if
            if (random_number() * 100 <= mutation_rate) then
                call mutate(tmp)
            end if
            new_population(:, :, child) = tmp
        end do

        ! === Remplacer population actuelle ===
        population = new_population
    end do

end program genetic_docking

module docking_utils
    implicit none
    private
    public :: atomdata, atominfo, bond
    public :: load_mol2_file, load_atoms_covalents, init_covalent_to_zero
    public :: add_covalent_into_atoms, determine_topology
    public :: detect_hydrogen_bonds

    type atomdata
        real, dimension(3) :: radii
        character(len=2) :: atom_symbol
        integer :: atom_number
    end type atomdata

    type atominfo
        character(len=10) :: atom_symbol
        real :: x, y, z
        real, dimension(3) :: radii
    end type atominfo

    type bond
        integer :: atom1, atom2
        integer :: bond_type
        real :: distance
    end type bond

contains

    ! (Réutilise ici tes fonctions existantes: load_mol2_file, load_atoms_covalents, etc.)

    subroutine detect_hydrogen_bonds(ligand, receptor)
        type(atominfo), intent(in) :: ligand(:), receptor(:)
        real :: d, max_hbond_dist
        integer :: i, j

        max_hbond_dist = 3.5  ! Distance typique pour une liaison H

        print *, "Hydrogen bond detection (distance < 3.5 Å):"
        do i = 1, size(ligand)
            if (trim(ligand(i)%atom_symbol) /= 'H') cycle
            do j = 1, size(receptor)
                if (trim(receptor(j)%atom_symbol) /= 'O' .and. trim(receptor(j)%atom_symbol) /= 'N') cycle
                d = distance(ligand(i), receptor(j))
                if (d < max_hbond_dist) then
                    print *, 'Possible H-bond: LIGAND H(', i, ') - RECEPTOR ', trim(receptor(j)%atom_symbol), '(', j, ') ->', d, 'Å'
                end if
            end do
        end do
    end subroutine detect_hydrogen_bonds

    real function distance(atom1, atom2)
        type(atominfo), intent(in) :: atom1, atom2
        distance = sqrt((atom1%x - atom2%x)**2 + (atom1%y - atom2%y)**2 + (atom1%z - atom2%z)**2)
    end function distance

end module docking_utils

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
  function evaluate_hbond(ligand, receptor) result(score)
    type(atominfo), dimension(:), intent(in) :: ligand, receptor
    real :: score
    integer :: i, j, k
    real :: d, angle, max_hbond_dist, min_hbond_dist
    logical :: is_hbond

    max_hbond_dist = 4.0
    min_hbond_dist = 2.2
    score = 0.0

    do i = 1, size(ligand)
        if (trim(ligand(i)%atom_symbol) /= 'H') cycle
        do j = 1, size(receptor)
            if (trim(receptor(j)%atom_symbol) /= 'O' .and. trim(receptor(j)%atom_symbol) /= 'N') cycle
            d = distance(ligand(i), receptor(j))
            if (d >= min_hbond_dist .and. d <= max_hbond_dist) then
                do k = 1, size(receptor)
                    if (k == j) cycle
                    angle = calculate_angle(ligand(i), receptor(j), receptor(k))
                    if (angle >= 90.0 .and. angle <= 150.0) then
                        score = score + 1.0
                        exit
                    end if
                end do
            end if
        end do
    end do
  end function

  real function calculate_angle(h, electroneg, neighbor)
    type(atominfo), intent(in) :: h, electroneg, neighbor
    real :: v1(3), v2(3)
    v1 = [h%x - electroneg%x, h%y - electroneg%y, h%z - electroneg%z]
    v2 = [neighbor%x - electroneg%x, neighbor%y - electroneg%y, neighbor%z - electroneg%z]
    calculate_angle = acos(dot_product(v1, v2) / (norm2(v1) * norm2(v2))) * 180.0 / acos(-1.0)
  end function
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

program affiche_topologie
    implicit none

    ! There is the AtomData type
    type atomdata
        real, dimension(3) :: radii
        character(len=2) :: atom_symbol
        integer :: atom_number
    end type atomdata
    type :: atominfo
        character(len=10) :: atom_symbol
        real :: x, y, z
        real, dimension(3) :: radii
    end type atominfo

    ! There is the AtomInfo array
    type(atominfo), allocatable :: atoms(:)
    ! There are AtomData arrays and simple AtomData
    type(atomdata), dimension(118) :: atoms_arrays
    type(atomdata) :: atom_data
    ! We put Covalent radii in a file
    character(len=255) :: filenameCovalent
    ! We put the MOL2 file in a file
    character(len=255) :: filenameMol2
    ! We put the MOL2 output file
    character(len=255) :: new_filename
    ! There are the bonds
    type bond
        integer :: atom1, atom2
        integer :: bond_type
        real :: distance
    end type bond
    type(bond), allocatable :: bonds(:)
    integer :: bond_count
    
    ! We ask the user for the filename of the MOL2 file in first argument and the covalent radii in second argument
    if (command_argument_count() /= 2) then
        print *, 'Usage: complete_mol2 <filenameMol2> <filenameCovalent>'
        stop
    end if
    call get_command_argument(1, filenameMol2)
    call get_command_argument(2, filenameCovalent)

    ! We load the CoValent radii
    atoms_arrays = init_covalent_to_zero()
    atoms_arrays = load_atoms_covalents(filenameCovalent)
    print *, 'Covalent Data loaded'
    ! We load the MOL2 file
    atoms = load_mol2_file(filenameMol2)
    print *, 'MOL2 Data loaded'

    ! We load the covalent radii into the atoms
    call add_covalent_into_atoms(atoms, atoms_arrays)

    ! We print the MOL2 data
    !call print_mol2_data(atoms)

    ! We determine the topology
    call determine_topology(atoms, atoms_arrays, bonds, bond_count)

    ! Génération du nouveau fichier MOL2
    new_filename = generate_output_filename(filenameMol2)
    call write_complete_mol2(filenameMol2, new_filename, bonds, bond_count)
    print *, 'Fichier MOL2 complété créé : ', trim(new_filename)

contains
    ! Fonction pour générer le nom du fichier de sortie
    function generate_output_filename(input_filename) result(output_filename)
        character(len=*), intent(in) :: input_filename
        character(len=:), allocatable :: output_filename
        integer :: dot_pos

        dot_pos = index(input_filename, '.', back=.true.)
        if (dot_pos == 0) then
            output_filename = trim(input_filename) // '_WITH_BONDS'
        else
            output_filename = input_filename(1:dot_pos-1) // '_WITH_BONDS' // input_filename(dot_pos:)
        end if
    end function generate_output_filename

    ! Subroutine pour écrire le fichier MOL2 complet
    subroutine write_complete_mol2(orig_file, new_file, bonds, bond_count)
        character(len=*), intent(in) :: orig_file, new_file
        type(bond), intent(in) :: bonds(:)
        integer, intent(in) :: bond_count
        integer :: u_in, u_out, ios, i, bond_id
        character(len=1000) :: line
        logical :: in_bond_section, after_atom

        open(newunit=u_in, file=orig_file, status='old', action='read')
        open(newunit=u_out, file=new_file, status='replace', action='write')

        after_atom = .false.
        in_bond_section = .false.

        do
            read(u_in, '(A)', iostat=ios) line
            if (ios /= 0) exit

            ! Gestion des sections spéciales
            if (index(line, '@<TRIPOS>ATOM') /= 0) then
                after_atom = .true.
                in_bond_section = .false.
                write(u_out, '(A)') trim(line)
            else if (index(line, '@<TRIPOS>BOND') /= 0) then
                in_bond_section = .true.
                cycle  ! Ignore la section BOND originale
            else if (line(1:1) == '@' .and. index(line, '<TRIPOS>') /= 0) then
                if (after_atom) then
                    call write_bond_section(u_out, bonds, bond_count)
                    after_atom = .false.
                end if
                in_bond_section = .false.
                write(u_out, '(A)') trim(line)
            else
                if (in_bond_section) then
                    cycle  ! Ignore les lignes de liaison existantes
                else
                    write(u_out, '(A)') trim(line)
                end if
            end if
        end do

        ! Écriture finale si nécessaire
        if (after_atom) call write_bond_section(u_out, bonds, bond_count)

        close(u_in)
        close(u_out)
    end subroutine write_complete_mol2

    ! Subroutine pour écrire la section BOND
    subroutine write_bond_section(u_out, bonds, bond_count)
        integer, intent(in) :: u_out, bond_count
        type(bond), intent(in) :: bonds(:)
        integer :: i
        
        write(u_out, '(A)') '@<TRIPOS>BOND'
        do i = 1, bond_count
            write(u_out, '(I6,1X,I5,1X,I5,1X,A1)') i, bonds(i)%atom1, bonds(i)%atom2, &
                char(ichar('0') + bonds(i)%bond_type)
        end do
    end subroutine write_bond_section

    ! Fonction pour calculer la distance entre deux atomes
    real function distance(atom1, atom2)
        type(atominfo), intent(in) :: atom1, atom2
        distance = sqrt( (atom1%x - atom2%x)**2 + &
                        (atom1%y - atom2%y)**2 + &
                        (atom1%z - atom2%z)**2 )
    end function distance

    ! Subroutine pour déterminer la topologie moléculaire
    subroutine determine_topology(atoms, atoms_arrays, bonds, bond_count)
        implicit none
        type(atominfo), intent(in) :: atoms(:)
        type(atomdata), dimension(118), intent(in) :: atoms_arrays
        real :: delta
        integer :: i, j, k, bond_type, natoms
        logical, allocatable :: bonded(:)
        real :: dist, sum_radii, lower, upper
        character(len=3) :: bond_type_str
        type(bond), allocatable :: bonds(:)
        integer :: bond_count

        natoms = size(atoms)
        allocate(bonded(natoms))
        delta = 0.1  ! Delta initial de 10%

        do while (delta <= 0.35)  ! Boucle jusqu'à delta=35%
            bonded = .false.
            bond_count = 0
            if (allocated(bonds)) deallocate(bonds)
            allocate(bonds(natoms*(natoms-1)/2))  ! Allocation maximale

            ! Vérification de toutes les paires d'atomes
            do i = 1, natoms
                do j = i+1, natoms
                    dist = distance(atoms(i), atoms(j))

                    ! Vérification des types de liaison (triple -> simple)
                    do bond_type = 3, 1, -1
                        sum_radii = atoms(i)%radii(bond_type) + atoms(j)%radii(bond_type)
                        lower = sum_radii * (1.0 - delta)
                        upper = sum_radii * (1.0 + delta)

                        if (dist >= lower .and. dist <= upper) then
                            bond_count = bond_count + 1
                            bonds(bond_count)%atom1 = i
                            bonds(bond_count)%atom2 = j
                            bonds(bond_count)%bond_type = bond_type
                            bonds(bond_count)%distance = dist
                            bonded(i) = .true.
                            bonded(j) = .true.
                            exit  ! On prend le type de liaison le plus fort possible
                        end if
                    end do
                end do
            end do

            ! Vérification si tous les atomes sont liés
            if (all(bonded)) exit
            delta = delta + 0.05  ! Augmentation du delta par paliers de 5%
        end do

        ! Affichage des résultats
        if (delta > 0.35) then
            print *, 'Impossible de déterminer la topologie avec un delta maximal de 35%'
        else
            print '(A, F5.1, A)', 'Topologie déterminée avec un delta de', delta*100, '%'
            print *, 'Liaisons détectées:'
            do k = 1, bond_count
                select case (bonds(k)%bond_type)
                    case (1)
                        bond_type_str = 'SIMPLE'
                    case (2)
                        bond_type_str = 'DOUBLE'
                    case (3)
                        bond_type_str = 'TRIPLE'
                end select
                print '(A, I0, A, I0, A, A)', 'Atome ', bonds(k)%atom1, ' - Atome ', bonds(k)%atom2, &
                                                ' : ', trim(bond_type_str)
            end do
        end if
        deallocate(bonded)
    end subroutine determine_topology

    subroutine print_mol2_data(atoms)
        implicit none
        type(atominfo), intent(in) :: atoms(:)
        integer :: i

        do i = 1, size(atoms)
            print *, atoms(i)%atom_symbol, atoms(i)%x, atoms(i)%y, atoms(i)%z, atoms(i)%radii
        end do
    end subroutine print_mol2_data

    function load_mol2_file(filename) result(atoms)
        implicit none
        character(len=255), intent(in) :: filename
        type(atominfo), allocatable :: atoms(:)
        character(len=100) :: line, atom_name
        integer :: ios, unit, i, nb_atoms, num_atom
        character(len=10) :: atom_symbol

        unit = 10
        open(unit=unit, file=trim(filename), status='old', action='read')

        ! We go to the molecule section
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '@<TRIPOS>MOLECULE') /= 0) then
            exit
            end if
        end do
        ! Read molecule name
        read(unit, '(A)', iostat=ios) atom_name
        ! Read number of atoms to allocate
        read(unit, '(A)', iostat=ios) line
        read(line, *) nb_atoms
        ! We allocate the array
        allocate(atoms(nb_atoms))

        ! We go the the atoms section
        i = 1
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '@<TRIPOS>ATOM') /= 0) then
                exit
            end if
        end do
        ! We read the atoms
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '@<TRIPOS>BOND') /= 0) then
                exit
            end if
            read(line, *) num_atom, atom_symbol, atoms(i)%x, atoms(i)%y, atoms(i)%z, atoms(i)%atom_symbol
            !atoms(i)%atom_symbol = adjustl(trim(atoms(i)%atom_symbol))
            i = i + 1
        end do
        close(unit)
    end function load_mol2_file

    subroutine add_covalent_into_atoms(atoms, atoms_arrays)
        implicit none
        type(atominfo), intent(inout) :: atoms(:)
        type(atomdata), dimension(118), intent(in) :: atoms_arrays
        integer :: i, j
        character(len=2) :: atom_symbol

        do i = 1, size(atoms)
            do j = 1, size(atoms_arrays)
                atoms(i)%radii = (/0, 0, 0/)
                ! We calc the atoms type from the symbol
                atom_symbol = atoms(i)%atom_symbol(1:1)
                if (len_trim(atoms(i)%atom_symbol) > 1) then
                    if (atoms(i)%atom_symbol(2:2) /= '.' .and. &
                        .not. (atoms(i)%atom_symbol(2:2) >= '0' .and. atoms(i)%atom_symbol(2:2) <= '9')) then
                        atom_symbol = atoms(i)%atom_symbol(1:2)
                    end if
                end if
                if (atom_symbol == atoms_arrays(j)%atom_symbol) then
                    atoms(i)%radii = atoms_arrays(j)%radii
                    exit
                end if
            end do
        end do
    end subroutine add_covalent_into_atoms

    function load_atoms_covalents(filename) result(atoms_arrays)
        implicit none
        character(len=255), intent(in) :: filename
        character(len=100) :: line
        type(atomdata), dimension(118) :: atoms_arrays
        type(atomdata) :: atom_data
        integer :: i, ios, unit
        integer :: a, b, c

        ! We read the file
        unit = 10
        i = 0
        open(unit=unit, file=trim(filename), status='old', action='read')
        do
            read(unit, '(A)', iostat=ios) line
            if (ios /= 0) exit
            if (index(line, '#') > 0) then
                cycle
            end if
            i = i + 1
            ! We extract data from the line
            read(line, '(I3, A2, 3I3)') atom_data%atom_number, atom_data%atom_symbol, a, b, c
            ! We convert from picometers to angstroms ( 1 pm <=> 0.01 angrstroms)
            atom_data%radii = (/a, b, c/)
            atom_data%radii = atom_data%radii * 0.01
            atoms_arrays(i) = atom_data
        end do
        close(unit)
    end function load_atoms_covalents

    function get_covalent_values(atoms_arrays, atom_symbol) result(atom_data)
        implicit none
        type(atomdata), dimension(118), intent(in) :: atoms_arrays

        character(len=2), intent(in) :: atom_symbol
        type(atomdata) :: atom_data
        integer :: i, i1

        ! We search for the atom number corresponding to the symbol
        i1 = -1
        atom_data = atomdata(radii=(/0, 0, 0/), atom_symbol='', atom_number=-1)
        do i = 1, size(atoms_arrays)
            if (atoms_arrays(i)%atom_symbol == atom_symbol) then
                i1 = i
                atom_data = atoms_arrays(i)
                exit
            end if
        end do

        ! We print the data
        if (i1 == -1) then
            print *, 'Symbol ', atom_symbol, ' not found'
        else
            print *, 'Symbol ', atom_data%atom_symbol, ' found'
            print *, 'Atomic Number is ', atom_data%atom_number
            print *, 'Covalent radii for ', atom_data%atom_symbol, ' are ', atom_data%radii
        end if
    end function get_covalent_values

    function init_covalent_to_zero() result(atoms_arrays)
        implicit none
        type(atomdata), dimension(118) :: atoms_arrays
        integer :: i

        ! We initialize the arrays to zero
        do i = 1, size(atoms_arrays)
            atoms_arrays(i) = atomdata(radii=(/0, 0, 0/), atom_symbol='', atom_number=-1)
        end do
    end function init_covalent_to_zero

end program affiche_topologie