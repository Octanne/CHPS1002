module mol2
    implicit none

    type :: Atom
        integer :: atom_number
        character(len=4) :: atom_name
        real :: x, y, z
        character(len=5) :: atom_type
    end type Atom

    type :: Individual
        real :: tx, ty, tz    ! translation
        real :: rx, ry, rz    ! rotation
        real :: fitness       ! quality (number of H-bonds)
    end type Individual

contains

    subroutine read_mol2(filename, atoms, n_atoms)
        character(len=*), intent(in) :: filename
        type(Atom), allocatable, intent(out) :: atoms(:)
        integer, intent(out) :: n_atoms
        character(len=256) :: line
        integer :: unit, idx_line, io_status

        open(newunit=unit, file=filename, status='old', action='read')

        do
            read(unit, '(a)', iostat=io_status) line
            if (io_status /= 0) stop "Error reading file"
            if (trim(line) == '@<TRIPOS>ATOM') exit
        end do

        n_atoms = 0
        do
            read(unit, '(a)', iostat=io_status) line
            if (io_status /= 0) exit
            if (trim(line) == '@<TRIPOS>BOND') exit
            n_atoms = n_atoms + 1
        end do

        rewind(unit)

        do
            read(unit, '(a)', iostat=io_status) line
            if (io_status /= 0) stop "Error reading file"
            if (trim(line) == '@<TRIPOS>ATOM') exit
        end do

        allocate(atoms(n_atoms))
        do idx_line = 1, n_atoms
            read(unit, *) atoms(idx_line)%atom_number, atoms(idx_line)%atom_name, &
                atoms(idx_line)%x, atoms(idx_line)%y, atoms(idx_line)%z, atoms(idx_line)%atom_type
        end do

        close(unit)
    end subroutine read_mol2

    subroutine transform_atoms(atoms, n_atoms, ind)
        type(Atom), intent(inout) :: atoms(:)
        integer, intent(in) :: n_atoms
        type(Individual), intent(in) :: ind
        integer :: idx_atom

        do idx_atom = 1, n_atoms
            atoms(idx_atom)%x = atoms(idx_atom)%x + ind%tx
            atoms(idx_atom)%y = atoms(idx_atom)%y + ind%ty
            atoms(idx_atom)%z = atoms(idx_atom)%z + ind%tz
        end do
    end subroutine transform_atoms

    subroutine compute_bonds(atoms, n_atoms, radii_file, bonds, bond_types, nbonds)
        implicit none
        character(len=*), intent(in)    :: radii_file
        type(Atom),    intent(in)       :: atoms(:)
        integer,       intent(in)       :: n_atoms
        integer,       allocatable, intent(out) :: bonds(:,:)
        character(len=2), allocatable, intent(out) :: bond_types(:)
        integer,       intent(out)      :: nbonds

        integer :: idx_atom1, idx_atom2, bond_count, io_status, unit
        real, allocatable :: cov_radii(:)
        real :: x1, y1, z1, x2, y2, z2
        real :: thr2, dist2
        real, parameter :: delta = 0.1
        integer :: atomn
        character(len=2) :: aname
        integer :: tmp1, tmp2, tmp3

        open(newunit=unit, file=trim(radii_file), status='old', action='read')
        do idx_atom1 = 1, 10
            read(unit, '(A)', iostat=io_status) aname
        end do
        allocate(cov_radii(n_atoms))
        do idx_atom1 = 1, n_atoms
            read(unit, '(I3, A2, I3, I3, I3)', iostat=io_status) atomn, aname, tmp1, tmp2, tmp3
            if (io_status /= 0) then
                cov_radii(idx_atom1) = 0.0
            else
                cov_radii(idx_atom1) = real(tmp1) / 100.0
            end if
        end do
        close(unit)

        bond_count = 0
        do idx_atom1 = 1, n_atoms - 1
            x1 = atoms(idx_atom1)%x; y1 = atoms(idx_atom1)%y; z1 = atoms(idx_atom1)%z
            do idx_atom2 = idx_atom1 + 1, n_atoms
                x2 = atoms(idx_atom2)%x; y2 = atoms(idx_atom2)%y; z2 = atoms(idx_atom2)%z
                dist2 = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
                thr2 = ((cov_radii(idx_atom1) + cov_radii(idx_atom2)) * (1.0 + delta))**2
                if (dist2 <= thr2) bond_count = bond_count + 1
            end do
        end do

        nbonds = bond_count
        if (nbonds == 0) then
            allocate(bonds(0, 2))
            allocate(bond_types(0))
            deallocate(cov_radii)
            return
        end if

        allocate(bonds(nbonds, 2))
        allocate(bond_types(nbonds))

        bond_count = 0
        do idx_atom1 = 1, n_atoms - 1
            x1 = atoms(idx_atom1)%x; y1 = atoms(idx_atom1)%y; z1 = atoms(idx_atom1)%z
            do idx_atom2 = idx_atom1 + 1, n_atoms
                x2 = atoms(idx_atom2)%x; y2 = atoms(idx_atom2)%y; z2 = atoms(idx_atom2)%z
                dist2 = (x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2
                thr2 = ((cov_radii(idx_atom1) + cov_radii(idx_atom2)) * (1.0 + delta))**2
                if (dist2 <= thr2) then
                    bond_count = bond_count + 1
                    bonds(bond_count, 1) = atoms(idx_atom1)%atom_number
                    bonds(bond_count, 2) = atoms(idx_atom2)%atom_number
                    bond_types(bond_count) = '1'
                end if
            end do
        end do

        deallocate(cov_radii)
    end subroutine compute_bonds

    subroutine save_individual(ind, ligand_file, output_file)
        implicit none
        type(Individual), intent(in) :: ind
        character(len=*), intent(in) :: ligand_file, output_file

        type(Atom), allocatable :: atoms(:)
        integer :: n_atoms, idx_atom, unit
        integer, allocatable :: bonds(:,:)
        character(len=2), allocatable :: bond_types(:)
        integer :: nbonds, idx_bond

        call read_mol2(ligand_file, atoms, n_atoms)
        call transform_atoms(atoms, n_atoms, ind)

        call compute_bonds(atoms, n_atoms, 'data/CoV_radii', bonds, bond_types, nbonds)

        open(newunit=unit, file=trim(output_file), status='replace', action='write')

        write(unit, '(A)') '@<TRIPOS>ATOM'
        do idx_atom = 1, n_atoms
            write(unit, '(I7,1X,A4,3F10.4,1X,A5)') &
                atoms(idx_atom)%atom_number, trim(atoms(idx_atom)%atom_name), &
                atoms(idx_atom)%x, atoms(idx_atom)%y, atoms(idx_atom)%z, &
                trim(atoms(idx_atom)%atom_type)
        end do

        write(unit, '(A)') '@<TRIPOS>BOND'
        do idx_bond = 1, nbonds
            write(unit, '(I7,2I6,1X,A)') idx_bond, bonds(idx_bond, 1), &
                bonds(idx_bond, 2), trim(bond_types(idx_bond))
        end do

        close(unit)

        deallocate(atoms)
        if (nbonds > 0) then
            deallocate(bonds, bond_types)
        else
            deallocate(bonds, bond_types)
        end if
    end subroutine save_individual

    subroutine mutate(ind)
        type(Individual), intent(inout) :: ind
        real :: mutation_strength

        call random_number(mutation_strength)
        ind%tx = ind%tx + (mutation_strength - 0.5) * 2.0
        call random_number(mutation_strength)
        ind%ty = ind%ty + (mutation_strength - 0.5) * 2.0
        call random_number(mutation_strength)
        ind%tz = ind%tz + (mutation_strength - 0.5) * 2.0
        call random_number(mutation_strength)
        ind%rx = ind%rx + (mutation_strength - 0.5) * 10.0
        call random_number(mutation_strength)
        ind%ry = ind%ry + (mutation_strength - 0.5) * 10.0
        call random_number(mutation_strength)
        ind%rz = ind%rz + (mutation_strength - 0.5) * 10.0
    end subroutine mutate

    subroutine crossover(parent1, parent2, child)
        type(Individual), intent(in) :: parent1, parent2
        type(Individual), intent(out) :: child
        real :: alpha

        call random_number(alpha)

        child%tx = alpha * parent1%tx + (1 - alpha) * parent2%tx
        child%ty = alpha * parent1%ty + (1 - alpha) * parent2%ty
        child%tz = alpha * parent1%tz + (1 - alpha) * parent2%tz
        child%rx = alpha * parent1%rx + (1 - alpha) * parent2%rx
        child%ry = alpha * parent1%ry + (1 - alpha) * parent2%ry
        child%rz = alpha * parent1%rz + (1 - alpha) * parent2%rz
    end subroutine crossover

    subroutine evaluate_individual(ind, ligand_atoms, n_ligand, site_atoms, n_site)
        implicit none
        type(Individual), intent(inout) :: ind
        type(Atom), intent(in) :: ligand_atoms(:)
        type(Atom), intent(in) :: site_atoms(:)
        integer, intent(in) :: n_ligand, n_site
        integer :: idx_ligand, idx_site, idx_neighbor
        real :: dist, angle
        integer :: count_hbonds
        type(Atom), allocatable :: new_atoms(:)

        allocate(new_atoms(n_ligand))
        new_atoms = ligand_atoms

        call transform_atoms(new_atoms, n_ligand, ind)

        count_hbonds = 0
        do idx_ligand = 1, n_ligand
            if (.not. is_hydrogen(ligand_atoms(idx_ligand))) cycle
            do idx_site = 1, n_site
                if (.not. is_acceptor(site_atoms(idx_site))) cycle
                dist = distance(new_atoms(idx_ligand), site_atoms(idx_site))
                if (dist >= 2.2 .and. dist <= 4.0) then
                    do idx_neighbor = 1, n_site
                        if (idx_neighbor == idx_site) cycle
                        angle = compute_angle(new_atoms(idx_ligand), site_atoms(idx_site), &
                                              site_atoms(idx_neighbor))
                        if (angle >= 90.0 .and. angle <= 150.0) then
                            count_hbonds = count_hbonds + 1
                            exit
                        end if
                    end do
                end if
            end do
        end do

        ind%fitness = real(count_hbonds)

        deallocate(new_atoms)
    end subroutine evaluate_individual

    logical function is_hydrogen(a) result(is_h)
        type(Atom), intent(in) :: a
        character(len=len_trim(a%atom_name)) :: nm
        character(len=1) :: sym

        nm = adjustl(trim(a%atom_name))
        sym = nm(1:1)
        is_h = (sym == 'H')
    end function is_hydrogen

    logical function is_acceptor(a) result(is_acc)
        type(Atom), intent(in) :: a
        character(len=len_trim(a%atom_name)) :: nm
        character(len=1) :: sym

        nm = adjustl(trim(a%atom_name))
        sym = nm(1:1)
        is_acc = (sym == 'O' .or. sym == 'N')
    end function is_acceptor

    real function distance(atom1, atom2)
        type(Atom), intent(in) :: atom1, atom2
        distance = sqrt((atom1%x - atom2%x)**2 + (atom1%y - atom2%y)**2 + (atom1%z - atom2%z)**2)
    end function distance

    real function compute_angle(h, acc, neighbor)
        type(Atom), intent(in) :: h, acc, neighbor
        real :: vec1(3), vec2(3), dot_product, norm1, norm2

        vec1 = [h%x - acc%x, h%y - acc%y, h%z - acc%z]
        vec2 = [neighbor%x - acc%x, neighbor%y - acc%y, neighbor%z - acc%z]

        dot_product = sum(vec1 * vec2)
        norm1 = sqrt(sum(vec1**2))
        norm2 = sqrt(sum(vec2**2))

        if (norm1 > 0 .and. norm2 > 0) then
            compute_angle = acos(dot_product / (norm1 * norm2)) * 180.0 / acos(-1.0)
        else
            compute_angle = 0.0
        end if
    end function compute_angle

end module mol2