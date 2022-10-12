MODULE gradient_descent

    use globvars
    use alarrays
    use ham
    use imgtp

    contains

subroutine zombie_alter(zstore,grad_fin,haml,elect,en,dvecs,lralt,lralt2,step)

    implicit none
    type(zombiest),dimension(:),intent(inout)::zstore
    type(grad),intent(inout)::grad_fin
    type(elecintrgl),intent(in)::elect
    type(dvector),dimension(:),intent(inout)::dvecs
    type(energy),intent(inout)::en
    type(hamiltonian),intent(inout)::haml
    type(zombiest),dimension(:),allocatable::temp_zom
    integer,intent(in)::step
    integer,intent(inout)::lralt,lralt2
    real(kind=8)::gamma,alpha,b,t,fxtdk,gradtd,picker,mmntm,rpropa,c 
    real(kind=8),dimension(ndet,norb)::mmntmmx
    integer::j,break,pick,k,sign
    integer(kind=8)::temp, temp2
    
    if (errorflag .ne. 0) return

    ! break=0
    alpha=0.0001
    b=1.0D-2
    c=1.1
    mmntm=0.6
    !gradtd=0
    if(step.eq.1)then 
        ! do j=1,ndet
        !     grad_fin%momentum(j,:)=zstore(j)%phi(:) !for momentum
        ! end do
        mmntmmx=0
        grad_fin%rpropaevious=0
        ! grad_fin%rprop=b
    else
        if(grad_fin%current_erg.gt.grad_fin%prev_erg)then 
            ! print*,'reducing learning rate'
            ! lralt=lralt+1
        else
            ! print*,'increasing learning rate'
            ! lralt2=lralt2+1
        end if
        do j=1,ndet
            mmntmmx(j,:)=zstore(j)%phi(:)-grad_fin%momentum(j,:)
        end do
        do j=1,ndet
            grad_fin%momentum(j,:)=zstore(j)%phi(:)
        end do
    end if
    
    t=b*(alpha**(lralt-1))*(c**(lralt2-1))
    call alloczs(temp_zom,ndet)

    !irprop+ with momentum
    ! do j=1,ndet
    !     do k=1, norb
    !         rpropa=grad_fin%rpropaevious(j,k)*grad_fin%vars(j,k)
    !         if(rpropa.gt.0)then
    !             ! if(grad_fin%vars(j,k).gt.0) sign=1
    !             ! if(grad_fin%vars(j,k).lt.0) sign=-1
    !             grad_fin%rprop(j,k)=1.2*grad_fin%rprop(j,k)
    !             temp_zom(j)%phi(k)=zstore(j)%phi(k)-(grad_fin%rprop(j,k)*(grad_fin%vars(j,k)))+b*mmntmmx(j,k)
    !         else if (rpropa.lt.0)then
    !             ! if(grad_fin%vars(j,k).gt.0) sign=1
    !             ! if(grad_fin%vars(j,k).lt.0) sign=-1 
    !             if(grad_fin%current_erg.gt.grad_fin%prev_erg)then 
    !                 grad_fin%vars(j,k)= grad_fin%rpropaevious(j,k)*(-1)
    !                 temp_zom(j)%phi(k)=grad_fin%momentum(j,k)
    !                 grad_fin%vars(j,k)=0
    !             else 
    !                 grad_fin%rprop(j,k)=0.5*grad_fin%rprop(j,k)
    !                 grad_fin%vars(j,k)=0
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)
    !             end if
    !         else if(rpropa.eq.0)then 
    !             temp_zom(j)%phi(k)=zstore(j)%phi(k)-(grad_fin%rprop(j,k)*(grad_fin%vars(j,k)))+b*mmntmmx(j,k)
    !         end if
    !         grad_fin%rpropaevious(j,k)=grad_fin%vars(j,k)
    !     end do
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     grad_fin%momentum(j,:)=zstore(j)%phi(:)
    ! end do


    !irprop+
    ! do j=1,ndet
    !     do k=1, norb
    !         rpropa=grad_fin%rpropaevious(j,k)*grad_fin%vars(j,k)
    !         if(rpropa.gt.0)then
    !             ! if(grad_fin%vars(j,k).gt.0) sign=1
    !             ! if(grad_fin%vars(j,k).lt.0) sign=-1
    !             grad_fin%rprop(j,k)=1.2*grad_fin%rprop(j,k)
    !             temp_zom(j)%phi(k)=zstore(j)%phi(k)-(grad_fin%rprop(j,k)*(grad_fin%vars(j,k)))
    !         else if (rpropa.lt.0)then
    !             ! if(grad_fin%vars(j,k).gt.0) sign=1
    !             ! if(grad_fin%vars(j,k).lt.0) sign=-1 
    !             if(grad_fin%current_erg.gt.grad_fin%prev_erg)then 
    !                 grad_fin%vars(j,k)= grad_fin%rpropaevious(j,k)*(-1)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)-(grad_fin%rprop(j,k)*(grad_fin%vars(j,k)))
    !                 grad_fin%vars(j,k)=0
    !             else 
    !                 grad_fin%rprop(j,k)=0.5*grad_fin%rprop(j,k)
    !                 grad_fin%vars(j,k)=0
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)-(grad_fin%rprop(j,k)*(grad_fin%vars(j,k)))
    !             end if
    !         else if(rpropa.eq.0)then 
    !             temp_zom(j)%phi(k)=zstore(j)%phi(k)-(grad_fin%rprop(j,k)*(grad_fin%vars(j,k)))
    !         end if
    !         grad_fin%rpropaevious(j,k)=grad_fin%vars(j,k)
    !     end do
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    ! end do

    !irprop-
    ! do j=1,ndet
    !     do k=1, norb
    !         rpropa=grad_fin%rpropaevious(j,k)*grad_fin%vars(j,k)
    !         if(rpropa.gt.0)then
    !             grad_fin%rprop(j,k)=1.2*grad_fin%rprop(j,k)
    !             temp_zom(j)%phi(k)=zstore(j)%phi(k)-(grad_fin%rprop(j,k)*(grad_fin%vars(j,k)))
    !         else if (rpropa.lt.0)then
    !             grad_fin%rprop(j,k)=0.5*grad_fin%rprop(j,k)
    !             grad_fin%vars(j,k)=0
    !             temp_zom(j)%phi(k)=zstore(j)%phi(k)-(grad_fin%rprop(j,k)*(grad_fin%vars(j,k)))
    !         end if
    !         grad_fin%rpropaevious(j,k)=grad_fin%vars(j,k)
    !     end do
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    ! end do

    !momentum type 1
    do j=1, ndet
        temp_zom(j)%phi=zstore(j)%phi(:)-(t*(grad_fin%vars(j,:)))+mmntm*mmntmmx(j,:)
        temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
        temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
       ! gradtd=gradtd+dot_product(grad_fin%vars(j,:),((-1)*grad_fin%vars(j,:)))
    end do

    !momentum type 2
    ! do j=1, ndet
    !     grad_fin%momentum(j,:)=(mmntm*grad_fin%rpropaevious(j,:))+grad_fin%vars(j,:)
    !     temp_zom(j)%phi=zstore(j)%phi(:)-t*(grad_fin%momentum(j,:))
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     grad_fin%rpropaevious(j,:)=grad_fin%momentum(j,:)
    !    ! gradtd=gradtd+dot_product(grad_fin%vars(j,:),((-1)*grad_fin%vars(j,:)))
    ! end do



   ! Normal GD
    ! do j=1, ndet
    !     temp_zom(j)%phi(:)=zstore(j)%phi(:)-(grad_fin%rprop(j,:)*(grad_fin%vars(j,:)))
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    ! end do


    GDflg='y'
    
    ! call random_number(picker)
    ! pick = anint(picker*ndet)
    ! temp_zom=zstore
    ! temp_zom(pick)%phi=zstore(pick)%phi(:)-(t*(grad_fin%vars(pick,:)))
    ! temp_zom(pick)%sin=sin(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
    ! temp_zom(pick)%cos=cos(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
    ! t=b*(alpha**(lralt-1))
    ! print*,t

    

    

   

    ! print*,gradtd
    ! call hamgen(haml,temp_zom,elect,ndet,0)
    ! dvecs(1)%d=cmplx(0.0,0.0)
    ! dvecs(1)%d(1)=cmplx(1.0,0.0)
    ! en%erg=0
    ! en%t=0
    ! call imgtime_prop(dvecs,en,haml)
    ! temp=anint(real(en%erg(1,timesteps+1))*10d12)
    ! fxtdk=real(en%erg(1,timesteps+1))      ! temp/10d12
    ! temp=fxtdk*10d12
    ! temp2=grad_fin%current_erg*10d12!+(alpha*t*gradtd)
    ! grad_fin%prev_erg=real(en%erg(1,timesteps+1))
    ! fxtdk=real(en%erg(1,timesteps+1)) !ergcalc(haml%hjk,dvecs(1)%d)
    ! print*,temp
    ! print*,temp2
    ! break=0
    ! if(temp.gt.(temp2))then!+(alpha*t*gradtd)))then !+(alpha*t*gradtd)
    !     do while(break.eq.0)
    !         lralt=lralt+1
    !         t=b*(alpha**(lralt-1))
    !         dvecs(1)%d=cmplx(0.0,0.0)
    !         dvecs(1)%d(1)=cmplx(1.0,0.0)
    !         en%erg=0
    !         en%t=0
            
            ! temp_zom=zstore
            ! temp_zom(pick)%phi=zstore(pick)%phi(:)-(t*(grad_fin%vars(pick,:)))
            ! temp_zom(pick)%sin=sin(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
            ! temp_zom(pick)%cos=cos(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
            ! do j=1, ndet
            !     temp_zom(j)%phi=zstore(j)%phi(:)-(t*(grad_fin%vars(j,:)))+mmntm*mmntmmx(j,:)
            !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
            !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
            !    ! gradtd=gradtd+dot_product(grad_fin%vars(j,:),((-1)*grad_fin%vars(j,:)))
            ! end do
        
            ! do j=1, ndet
            !     temp_zom(j)%phi=zstore(j)%phi(:)-(t*(grad_fin%vars(j,:)))
            !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
            !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
            ! end do
    !         call hamgen(haml,temp_zom,elect,ndet,0)
    !         call imgtime_prop(dvecs,en,haml)
    !         ! temp=anint(real(en%erg(1,timesteps+1))*10d12)
    !         ! fxtdk=temp/10d12
    !         fxtdk=real(en%erg(1,timesteps+1))
    !         temp=int(fxtdk*10d12)
    !         ! fxtdk= real(en%erg(1,timesteps+1)) !ergcalc(haml%hjk,dvecs(1)%d)
    !         print*,temp
    !         if((temp.lt.(temp2)))then
    !             print*,fxtdk
    !             break=5
    !         end if
    !     end do
    ! end if
    grad_fin%prev_erg=grad_fin%current_erg
    ! print*,temp_zom(2)%phi
    zstore=temp_zom
    ! grad_fin%momentum=grad_fin%vars
    ! print*,t
    call dealloczs(temp_zom) 
    GDflg='y'
    ! stop

end subroutine zombie_alter


END MODUle gradient_descent