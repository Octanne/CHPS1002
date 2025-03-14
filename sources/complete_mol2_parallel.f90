program affiche_topologie
    implicit none
    ! Include OpenMP module
    include 'omp_lib.h'

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
            !$omp parallel do private(j, dist, bond_type, sum_radii, lower, upper) shared(bonded, bond_count, bonds)
            do i = 1, natoms
                do j = i+1, natoms
                    dist = distance(atoms(i), atoms(j))

                    ! Vérification des types de liaison (triple -> simple)
                    do bond_type = 3, 1, -1
                        sum_radii = atoms(i)%radii(bond_type) + atoms(j)%radii(bond_type)
                        lower = sum_radii * (1.0 - delta)
                        upper = sum_radii * (1.0 + delta)

                        if (dist >= lower .and. dist <= upper) then
                            !$omp critical
                            bond_count = bond_count + 1
                            bonds(bond_count)%atom1 = i
                            bonds(bond_count)%atom2 = j
                            bonds(bond_count)%bond_type = bond_type
                            bonds(bond_count)%distance = dist
                            bonded(i) = .true.
                            bonded(j) = .true.
                            !$omp end critical
                            exit  ! On prend le type de liaison le plus fort possible
                        end if
                    end do
                end do
            end do
            !$omp end parallel do

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

        !$omp parallel do private(j, atom_symbol) shared(atoms, atoms_arrays)
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
        !$omp end parallel do
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