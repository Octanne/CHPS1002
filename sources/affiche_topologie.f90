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
    ! There are the bonds
    type bond
        character(len=2) :: atom1, atom2
        real :: distance
    end type bond
    type(bond), allocatable :: bonds(:)
    integer :: bond_count
    
    ! We ask the user for the filename of the MOL2 file in first argument and the covalent radii in second argument
    if (command_argument_count() /= 2) then
        print *, 'Usage: affiche_topologie <filenameMol2> <filenameCovalent>'
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
    call print_mol2_data(atoms)


    !determine_topology(atoms, atoms_arrays)

contains
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
            read(line, *) num_atom, atoms(i)%atom_symbol, atoms(i)%x, atoms(i)%y, atoms(i)%z
            i = i + 1
        end do
        close(unit)
    end function load_mol2_file    

    subroutine add_covalent_into_atoms(atoms, atoms_arrays)
        implicit none
        type(atominfo), intent(inout) :: atoms(:)
        type(atomdata), dimension(118), intent(in) :: atoms_arrays
        integer :: i, j

        do i = 1, size(atoms)
            do j = 1, size(atoms_arrays)
                atoms(i)%radii = (/-1, -1, -1/)
                if (atoms(i)%atom_symbol == atoms_arrays(j)%atom_symbol) then
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