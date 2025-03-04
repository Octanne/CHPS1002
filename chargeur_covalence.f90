program read_covalance
    implicit none

    ! There is the AtomData type
    type atomdata
        real, dimension(3) :: radii
        character(len=2) :: atom_symbol
        integer :: atom_number
    end type atomdata    
    
    ! There are AtomData arrays and simple AtomData
    type(atomdata), dimension(118) :: atoms_arrays
    type(atomdata) :: atom_data

    ! We initialize the arrays to store the data => int (atomic_number) which give 3 real (one_bond, double_bond, triple_bond)
    character(len=100) :: filename
    character(len=2) :: atom_symbol
    integer :: i

    ! We ask the user for the filename
    if (iargc() < 1) then
        print *, 'Usage: read_covalance <filename>'
        stop
    end if
    call getarg(1, filename)

    atoms_arrays = init_covalent_to_zero()
    atoms_arrays = load_atoms_covalents(filename)
    print *, 'Data loaded'

    ! We ask the user for the symbol three times with a loop
    do
        atom_symbol = ask_symbol()
        atom_data = get_covalent_values(atoms_arrays, atom_symbol)
        i = i + 1
        if (i == 3) exit
    end do

contains
    function init_covalent_to_zero() result(atoms_arrays)
        implicit none
        type(atomdata), dimension(118) :: atoms_arrays
        integer :: i

        ! We initialize the arrays to zero
        do i = 1, size(atoms_arrays)
            atoms_arrays(i) = atomdata(radii=(/0, 0, 0/), atom_symbol='', atom_number=-1)
        end do
    end function init_covalent_to_zero

    function ask_symbol() result(symbol)
        implicit none
        character(len=2) :: symbol

        ! We ask the user for the symbol
        print *, '> Enter the symbol of the atom'
        read(*, '(A2)') symbol
    end function ask_symbol

    function load_atoms_covalents(file) result(atoms_arrays)
        implicit none
        character(len=100), intent(in) :: file
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
            print *, 'Atom number ', atom_data%atom_number, ' symbol ', atom_data%atom_symbol, ' radii ', atom_data%radii
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

end program read_covalance