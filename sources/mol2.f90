module mol2
    !use individuals
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
        real :: fitness       ! qualité (nombre de ponts H)
    end type Individual

contains

    subroutine read_mol2(filename, atoms, n_atoms)
        character(len=*), intent(in) :: filename
        type(Atom), allocatable, intent(out) :: atoms(:)
        integer, intent(out) :: n_atoms
        character(len=256) :: line
        integer :: unit, i, ok

        open(newunit=unit, file=filename, status='old', action='read')

        do
            read(unit, '(a)', iostat=ok) line
            if (ok /= 0) stop "Error reading file"
            if (trim(line) == '@<TRIPOS>ATOM') exit
        end do

        n_atoms = 0
        do
            read(unit, '(a)', iostat=ok) line
            if (ok /= 0) exit
            if (trim(line) == '@<TRIPOS>BOND') exit
            n_atoms = n_atoms + 1
        end do

        rewind(unit)

        do
            read(unit, '(a)', iostat=ok) line
            if (ok /= 0) stop "Error reading file"
            if (trim(line) == '@<TRIPOS>ATOM') exit
        end do

        allocate(atoms(n_atoms))
        do i = 1, n_atoms
            read(unit, *) atoms(i)%atom_number, atoms(i)%atom_name, atoms(i)%x, atoms(i)%y, atoms(i)%z, atoms(i)%atom_type
        end do

        close(unit)
    end subroutine read_mol2

    subroutine transform_atoms(atoms, n_atoms, ind)
        type(Atom), intent(inout) :: atoms(:)
        integer, intent(in) :: n_atoms
        type(Individual), intent(in) :: ind
        integer :: i

        do i = 1, n_atoms
            atoms(i)%x = atoms(i)%x + ind%tx
            atoms(i)%y = atoms(i)%y + ind%ty
            atoms(i)%z = atoms(i)%z + ind%tz
            ! Ici, on pourrait ajouter les rotations rx, ry, rz (facultatif si tu veux)
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

        ! Locals
        integer :: i, j, count, ios, unit
        real, allocatable :: cov_radii(:)
        real :: x1,y1,z1, x2,y2,z2
        real :: thr2, dist2
        real, parameter :: delta = 0.1
        integer :: atomn
        character(len=2) :: aname
        integer :: tmp1,tmp2,tmp3

        !----------------------------------------------------------------
        ! 1) Lire les rayons covalents (on ne garde que le 1er rayon)
        !----------------------------------------------------------------
        open(newunit=unit, file=trim(radii_file), status='old', action='read')
        ! sauter l'en-tête (10 lignes)
        do i = 1,10
            read(unit,'(A)',iostat=ios) aname
        end do
        allocate(cov_radii(n_atoms))
        do i = 1, n_atoms
            read(unit,'(I3, A2, I3, I3, I3)', iostat=ios) atomn, aname, tmp1, tmp2, tmp3
            if (ios/=0) then
                cov_radii(i) = 0.0
            else
                cov_radii(i) = real(tmp1) / 100.0
            end if
        end do
        close(unit)

        !----------------------------------------------------------------
        ! 2) première passe : compter le nombre de liaisons à créer
        !----------------------------------------------------------------
        count = 0
        do i = 1, n_atoms-1
            x1 = atoms(i)%x; y1 = atoms(i)%y; z1 = atoms(i)%z
            do j = i+1, n_atoms
                x2 = atoms(j)%x; y2 = atoms(j)%y; z2 = atoms(j)%z
                dist2 = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
                thr2  = ((cov_radii(i) + cov_radii(j)) * (1.0 + delta))**2
                if (dist2 <= thr2) count = count + 1
            end do
        end do

        nbonds = count
        if (nbonds == 0) then
            allocate(bonds(0,2))
            allocate(bond_types(0))
            deallocate(cov_radii)
            return
        end if

        !----------------------------------------------------------------
        ! 3) allouer les tableaux à la taille exacte
        !----------------------------------------------------------------
        allocate(bonds(nbonds,2))
        allocate(bond_types(nbonds))

        !----------------------------------------------------------------
        ! 4) seconde passe : remplir bonds et bond_types
        !----------------------------------------------------------------
        count = 0
        do i = 1, n_atoms-1
            x1 = atoms(i)%x; y1 = atoms(i)%y; z1 = atoms(i)%z
            do j = i+1, n_atoms
                x2 = atoms(j)%x; y2 = atoms(j)%y; z2 = atoms(j)%z
                dist2 = (x1-x2)**2 + (y1-y2)**2 + (z1-z2)**2
                thr2  = ((cov_radii(i) + cov_radii(j)) * (1.0 + delta))**2
                if (dist2 <= thr2) then
                    count = count + 1
                    bonds(count,1) = atoms(i)%atom_number
                    bonds(count,2) = atoms(j)%atom_number
                    bond_types(count) = '1'    ! ici liaison simple par défaut
                end if
            end do
        end do

        ! Cleanup
        deallocate(cov_radii)
    end subroutine compute_bonds


    !----------------------------------------------------------------
    ! Save a transformed ligand with ATOM and BOND sections
    !----------------------------------------------------------------
    subroutine save_individual(ind, ligand_file, output_file)
        implicit none
        type(Individual), intent(in) :: ind
        character(len=*), intent(in) :: ligand_file, output_file

        type(Atom), allocatable :: atoms(:)
        integer :: n_atoms, i, unit
        integer, allocatable :: bonds(:,:)
        character(len=2), allocatable :: bond_types(:)
        integer :: nbonds

        ! 1) lire et transformer
        call read_mol2(ligand_file, atoms, n_atoms)
        call transform_atoms(atoms, n_atoms, ind)

        ! 2) calcul des liaisons covalentes
        call compute_bonds(atoms, n_atoms, 'CoV_radii', bonds, bond_types, nbonds)

        ! 3) écriture du fichier .mol2
        open(newunit=unit, file=trim(output_file), status='replace', action='write')

        write(unit,'(A)') '@<TRIPOS>ATOM'
        do i = 1, n_atoms
            write(unit,'(I7,1X,A4,3F10.4,1X,A5)') &
                atoms(i)%atom_number, trim(atoms(i)%atom_name), &
                atoms(i)%x, atoms(i)%y, atoms(i)%z, trim(atoms(i)%atom_type)
        end do

        write(unit,'(A)') '@<TRIPOS>BOND'
        do i = 1, nbonds
            write(unit,'(I7,2I6,1X,A)') i, bonds(i,1), bonds(i,2), trim(bond_types(i))
        end do

        close(unit)

        ! clean
        deallocate(atoms)
        if (nbonds > 0) then
            deallocate(bonds, bond_types)
        else
            deallocate(bonds, bond_types)  ! même si length 0
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

        child%tx = alpha * parent1%tx + (1-alpha) * parent2%tx
        child%ty = alpha * parent1%ty + (1-alpha) * parent2%ty
        child%tz = alpha * parent1%tz + (1-alpha) * parent2%tz
        child%rx = alpha * parent1%rx + (1-alpha) * parent2%rx
        child%ry = alpha * parent1%ry + (1-alpha) * parent2%ry
        child%rz = alpha * parent1%rz + (1-alpha) * parent2%rz
    end subroutine crossover

    subroutine evaluate_individual(ind, ligand_atoms, n_ligand, site_atoms, n_site)
        implicit none
        type(Individual), intent(inout) :: ind
        type(Atom), intent(in) :: ligand_atoms(:)
        type(Atom), intent(in) :: site_atoms(:)
        integer, intent(in) :: n_ligand, n_site
        integer :: i, j, k
        real :: dist, angle
        integer :: count_hbonds
        type(Atom), allocatable :: new_atoms(:)

        allocate(new_atoms(n_ligand))
        new_atoms = ligand_atoms

        call transform_atoms(new_atoms, n_ligand, ind)

        count_hbonds = 0
        do i = 1, n_ligand
            if (.not. is_hydrogen(ligand_atoms(i))) cycle
            do j = 1, n_site
                if (.not. is_acceptor(site_atoms(j))) cycle
                dist = distance(new_atoms(i), site_atoms(j))
                if (dist >= 2.2 .and. dist <= 4.0) then
                    do k = 1, n_site
                        if (k == j) cycle
                        angle = compute_angle(new_atoms(i), site_atoms(j), site_atoms(k))
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

        nm  = adjustl(trim(a%atom_name))
        sym = nm(1:1)
        is_h = (sym == 'H')
    end function is_hydrogen

    logical function is_acceptor(a) result(is_acc)
        type(Atom), intent(in) :: a
        character(len=len_trim(a%atom_name)) :: nm
        character(len=1) :: sym

        nm  = adjustl(trim(a%atom_name))
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
            compute_angle = acos( dot_product / (norm1 * norm2) ) * 180.0 / acos(-1.0)
        else
            compute_angle = 0.0
        end if
    end function compute_angle

end module mol2
