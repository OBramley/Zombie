MODULE gradient_descent

    use globvars
    use alarrays
    use ham
    use imgtp

    contains

subroutine zombie_alter(zstore,grad_fin,haml,elect,en,dvecs)

    implicit none
    type(zombiest),dimension(:),intent(inout)::zstore
    type(grad),intent(in)::grad_fin
    type(elecintrgl),intent(in)::elect
    type(dvector),dimension(:),intent(inout)::dvecs
    type(energy),intent(inout)::en
    type(hamiltonian),intent(inout)::haml
    type(zombiest),dimension(:),allocatable::temp_zom
    real(kind=8)::gamma,alpha,b,t,fxtdk,gradtd 
    integer::j,break

    if (errorflag .ne. 0) return

    break=0
    alpha=0.0001
    b=0.8
    t=1
    gradtd=0

    GDflg='n'
    call alloczs(temp_zom,ndet)
    do j=1, ndet
        temp_zom(j)%phi=zstore(j)%phi(:)-(t*(grad_fin%vars(j,:)))
        temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
        temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
        gradtd=gradtd+dot_product(grad_fin%vars(j,:),((-1)*grad_fin%vars(j,:)))
    end do
    call hamgen(haml,temp_zom,elect,ndet,0)
    dvecs(1)%d=cmplx(0.0,0.0)
    dvecs(1)%d(1)=cmplx(1.0,0.0)
    en%erg=0
    en%t=0
    call imgtime_prop(dvecs,en,haml)
    fxtdk=real(en%erg(1,timesteps+1)) !ergcalc(haml%hjk,dvecs(1)%d)
    if(fxtdk.lt.(grad_fin%prev_erg+(alpha*t*gradtd)))then 
        zstore=temp_zom
    else 
        do while(break.eq.0)
            t=b*t
            ! dvecs(1)%d=0
            ! en%erg=(0.0,0.0)
            ! en%t=0
            do j=1, ndet
                temp_zom(j)%phi=zstore(j)%phi(:)-(t*(grad_fin%vars(j,:)))
                temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
                temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
            end do
            call hamgen(haml,temp_zom,elect,ndet,0)
            call imgtime_prop(dvecs,en,haml)
            fxtdk= real(en%erg(1,timesteps+1)) !ergcalc(haml%hjk,dvecs(1)%d)
            ! print*,fxtdk
            if(fxtdk.lt.(grad_fin%prev_erg+(alpha*t*gradtd)))then
                break=5
            end if
        end do
    end if
    zstore=temp_zom
    ! print*,t
    call dealloczs(temp_zom) 
    GDflg='y'
    ! stop

end subroutine zombie_alter


END MODUle gradient_descent