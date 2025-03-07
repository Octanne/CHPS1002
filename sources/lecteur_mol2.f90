program read_mol2
    implicit none
    type :: atominfo
        character(len=10) :: atom_symbol
        real :: x, y, z
    end type atominfo

    type(atominfo), allocatable :: atoms(:)
    character(len=255) :: filename

    ! We ask the user for the filename
    if (iargc() < 1) then
        print *, 'Usage: read_atoms <filename>'
        stop
    end if
    call getarg(1, filename)

    atoms = load_mol2_file(filename)
    call print_mol2_data(atoms)

contains
    subroutine print_mol2_data(atoms)
        implicit none
        type(atominfo), intent(in) :: atoms(:)
        integer :: i

        do i = 1, size(atoms)
            print *, atoms(i)%atom_symbol, atoms(i)%x, atoms(i)%y, atoms(i)%z
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

end program read_mol2