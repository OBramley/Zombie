program test 

    use stuff 
    use vals 


    implicit none 

    type(zombiest), dimension(3):: zstore
    type(elecintrgl)::elecs
    type(oprts)::an_cr,an2_cr2,an2_cr2_diff
    integer::j,k,ierr,l
    real(kind=8)::r
    integer,allocatable,dimension(:,:,:)::occupancy_an
    integer,allocatable,dimension(:,:,:,:)::occupancy_2an
    complex(kind=8),allocatable, dimension(:,:,:,:)::z1jk
    complex(kind=8),allocatable, dimension(:,:,:)::z2l
    real(kind=8),dimension(10,10)::haml_diff 
    complex(kind=8),dimension(10)::zs1sin,zs1cos,zs2sin,zs2cos


    call init()

    allocate(occupancy_an(norb,2,norb),stat=ierr)
    allocate(occupancy_2an(norb,norb,2,norb),stat=ierr)
    allocate(z1jk(norb,norb,2,norb),stat=ierr)
    allocate(z2l(norb,2,norb),stat=ierr)

    

    occupancy_an=1
    occupancy_2an=1

        
       
    do j=1,norb
        occupancy_an(j,1,j)=0
        do l=j-1, 1, -1
            occupancy_an(j,1,l)=-1
        end do
    end do
    
    do j=1,norb
        do k=1, norb
            occupancy_2an(j,k,:,:)=occupancy_an(j,:,:)
            occupancy_2an(j,k,2,k)=occupancy_2an(j,k,1,k)
            occupancy_2an(j,k,1,k)=0
            do l=k-1,1,-1
                occupancy_2an(j,k,1,l)=occupancy_2an(j,k,1,l)*(-1)
            end do
        end do
    end do
      


    allocate(elecs%h1ei_old(10,10))
    allocate(elecs%h2ei_old(10,10,10,10))
        

    call electronintegrals(elecs,an_cr,an2_cr2,an2_cr2_diff)
    call electronintegrals_old(elecs)

    zstore(1)%phi(1:nel)=0.5*pirl
    zstore(1)%phi(nel+1:)=0
    zstore(1)%sin=0
    zstore(1)%cos=1
    zstore(1)%sin(1:nel)=1
    zstore(1)%cos(1:nel)=0
    zstore(1)%val(nel:norb)=zstore(1)%sin
    zstore(1)%val(norb+1:)=zstore(1)%cos
    do j=2,3
        do k=1,norb
            zstore(j)%phi(k)=2*pirl*r 
        end do 
        zstore(j)%sin=sin(zstore(j)%phi)
        zstore(j)%cos=cos(zstore(j)%phi)
        zstore(j)%val(1:)=zstore(j)%sin
        zstore(j)%val(norb+1:)=zstore(j)%cos
    end do


    ! do l=1, norb
        z2l(l,1,:)=zstore(2)%sin(:)
        z2l(l,2,:)=zstore(2)%cos(:)
        z2l(l,2,l)=z2l(l,1,l)
        z2l(l,1,l)=cmplx(0.0,0.0)
    ! end do

    do j=1, norb
        do k=1, norb
            z1jk(j,k,:,:)=z2l(j,:,:)
            z1jk(j,k,2,k)=z1jk(j,k,1,k)
            z1jk(j,k,1,k)=cmplx(0.0,0.0)
        end do
    end do

    z2l(l,1,:)=zstore(1)%sin(:)
    z2l(l,2,:)=zstore(1)%cos(:)
    z2l(l,2,l)=z2l(l,1,l)
    z2l(l,1,l)=cmplx(0.0,0.0)

    !!$omp end do simd
    !$omp end parallel
    z1jk=z1jk*occupancy_2an

    call haml_grad(haml_diff,zstore,elecs,an_cr,an2_cr2,an2_cr2_diff,2)

    zs1sin=cmplx(zstore(2)%sin,0.0d0,kind=8)
    zs1cos=cmplx(zstore(2)%cos,0.0d0,kind=8)
    zs2sin=cmplx(zstore(1)%sin,0.0d0,kind=8)
    zs2cos=cmplx(zstore(1)%cos,0.0d0,kind=8)
    
    call two_elec_part_grad(zs1sin,zs1cos,zs2sin,zs2cos,z2l,z1jk,elecs%h2ei_old,occupancy_2an,occupancy_an) 
    


end program test