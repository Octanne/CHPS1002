program read_mol2
    implicit none
    character(len=100) :: line
    character(len=10) :: atom_symbol
    real :: x, y, z
    integer :: ios, unit, atom_number

    character(len=255) :: filename
    print *, 'Enter the name of the file to read:'
    read *, filename

    unit = 10
    open(unit=unit, file=trim(filename), status='old', action='read')
    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, '@<TRIPOS>ATOM') /= 0) then
            exit
        end if
    end do

    do
        read(unit, '(A)', iostat=ios) line
        if (ios /= 0) exit
        if (index(line, '@<TRIPOS>BOND') /= 0) then
            exit
        end if
        read(line, *) atom_number, atom_symbol, x, y, z
        print *, 'Atom number: ', atom_number  
        print *, 'Atom symbol: ', atom_symbol      
        print *, 'Coordinates: ', x, y, z
    end do

    close(unit)
end program read_mol2