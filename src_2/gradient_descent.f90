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
    real,intent(inout)::lralt,lralt2
    real(kind=8)::gamma,alpha,b,t,fxtdk,gradtd,picker,mmntm,rpropa,c
    real(kind=16)::delmax,delmin,del0,nablaplus,nablaminu
    real(kind=8),dimension(ndet,norb)::mmntmmx
    integer::j,break,pick,k,sign,power
    integer(kind=8)::temp, temp2
    
    if (errorflag .ne. 0) return

    ! break=0
    alpha=0.1
    b=1D-1
    c=1.001
    ! mmntm=0.1
    ! delmax=
    ! delmin=
    ! del0=0.00000999999999
    ! nablaplus=1.009999999999 !old1.09999999999
    ! nablaminu=0.249999999999 !old0.49999999999

    !gradtd=0
    power=1+int(lralt)
    if(step.eq.1)then 
        ! do j=2,ndet
        !     grad_fin%prev_phi(j,:)=zstore(j)%phi(:) !for momentum
        ! end do
        ! mmntmmx=0
        ! grad_fin%rpropaevious=0
        ! grad_fin%rprop=del0
    else
        ! if(modulo(step,100).eq.0)then
        !     lralt=lralt+1
        ! end if
        if(grad_fin%current_erg.gt.grad_fin%prev_erg)then 
            print*,'Bad move'
            lralt=lralt+0.1
        else
            ! print*,'increasing learning rate'
            ! lralt2=lralt2+1
        end if
        ! heavy ball
        ! do j=2,ndet
        !     mmntmmx(j,:)=zstore(j)%phi(:)-grad_fin%prev_phi(j,:)
        !     grad_fin%prev_phi(j,:)=zstore(j)%phi(:)
        ! end do
        !Momentum
        ! do j=2,ndet
        !     mmntmmx(j,:)=mmntm*grad_fin%prev_phi(j,:)
        ! end do
    end if
    
    !decay by 0.5 when energy increases
    t=b*(alpha**(power-1))*(c**(lralt2-1))
    ! t=b*exp(-0.01*(step-1))
    call alloczs(temp_zom,ndet)
    temp_zom(1)=zstore(1)

    
    !irprop+
    ! temp_zom=zstore
    ! do j=2,ndet
    !     do k=1, norb
    !         !To compare if gradients were the same sign
    !         rpropa=grad_fin%rpropaevious(j,k)*grad_fin%vars(j,k)
    !         if(grad_fin%vars(j,k).gt.0)then !positive current gradient
    !             if(rpropa.gt.0)then !Same adjacent gradients
    !                 grad_fin%rprop(j,k)=(-1)*nablaplus*grad_fin%rprop(j,k)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)+grad_fin%rprop(j,k)
    !                 grad_fin%rpropaevious(j,k)=1
    !             else if (rpropa.lt.0)then !different adjacent gradients
    !                 if(grad_fin%prev_erg.lt.grad_fin%current_erg)then
    !                     temp_zom(j)%phi(k)=zstore(j)%phi(k)-grad_fin%rprop(j,k)
    !                 end if
    !                 grad_fin%rprop(j,k)=nablaminu*grad_fin%rprop(j,k)
    !                 grad_fin%rpropaevious(j,k)=0
    !             else
    !                 grad_fin%rprop(j,k)=(-1)*nablaplus*grad_fin%rprop(j,k)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)+grad_fin%rprop(j,k)
    !             end if
    !         else if(grad_fin%vars(j,k).lt.0)then !negative current gradient
    !             if(rpropa.gt.0)then !Same adjacent gradients
    !                 grad_fin%rprop(j,k)=nablaplus*grad_fin%rprop(j,k)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)+grad_fin%rprop(j,k)
    !                 grad_fin%rpropaevious(j,k)=-1
    !             else if (rpropa.lt.0)then !different adjacent gradients
    !                 if(grad_fin%prev_erg.lt.grad_fin%current_erg)then
    !                     temp_zom(j)%phi(k)=zstore(j)%phi(k)-grad_fin%rprop(j,k)
    !                 end if
    !                 grad_fin%rprop(j,k)= nablaminu*grad_fin%rprop(j,k)
    !                 grad_fin%rpropaevious(j,k)=0
    !             else
    !                 grad_fin%rprop(j,k)=nablaplus*grad_fin%rprop(j,k)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)+grad_fin%rprop(j,k)
    !             end if
    !         end if
    !     end do
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    ! end do

    !irprop-
 
    ! temp_zom=zstore
    ! do j=2,ndet
    !     do k=1, norb
    !         !To compare if gradients were the same sign
    !         rpropa=grad_fin%rpropaevious(j,k)*grad_fin%vars(j,k)
    !         if(grad_fin%vars(j,k).gt.0)then !positive current gradient
    !             if(rpropa.gt.0)then !Same adjacent gradients
    !                 grad_fin%rprop(j,k)=nablaplus*grad_fin%rprop(j,k)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)-grad_fin%rprop(j,k)
    !                 grad_fin%rpropaevious(j,k)=1
    !             else if (rpropa.lt.0)then !different adjacent gradients
    !                 grad_fin%rprop(j,k)=nablaminu*grad_fin%rprop(j,k)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)-grad_fin%rprop(j,k)
    !                 grad_fin%rpropaevious(j,k)=0
    !             end if
    !         else if(grad_fin%vars(j,k).lt.0)then !negative current gradient
    !             if(rpropa.gt.0)then !Same adjacent gradients
    !                 grad_fin%rprop(j,k)=nablaplus*grad_fin%rprop(j,k)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)+grad_fin%rprop(j,k)
    !                 grad_fin%rpropaevious(j,k)=-1
    !             else if (rpropa.lt.0)then !different adjacent gradients
    !                 grad_fin%rprop(j,k)= nablaminu*grad_fin%rprop(j,k)
    !                 temp_zom(j)%phi(k)=zstore(j)%phi(k)+grad_fin%rprop(j,k)
    !                 grad_fin%rpropaevious(j,k)=0
    !             end if
    !         end if
    !     end do
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    ! end do
    ! if(step.eq.1)then 
    !     do j=2,ndet
    !         do k=1, norb
    !             if(grad_fin%vars(j,k).gt.0)then
    !                 grad_fin%rpropaevious=1
    !             else if(grad_fin%vars(j,k).lt.0)then
    !                 grad_fin%rpropaevious=-1
    !             end if 
    !         end do
    !     end do
    ! end if
   
    !momentum Heavy ball
    ! do j=2, ndet
    !     temp_zom(j)%phi=zstore(j)%phi(:)-(t*(grad_fin%vars(j,:)))+mmntm*mmntmmx(j,:)
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !    ! gradtd=gradtd+dot_product(grad_fin%vars(j,:),((-1)*grad_fin%vars(j,:)))
    ! end do

    !momentum type 2

    ! do j=2, ndet
    !     ! print*, maxval(abs(grad_fin%vars(j,:))),minval(abs(grad_fin%vars(j,:))),sum(abs(grad_fin%vars(j,:)))/norb
    !     grad_fin%prev_phi(j,:)=mmntmmx(j,:)+((1-mmntm)*grad_fin%vars(j,:))
    !     temp_zom(j)%phi=zstore(j)%phi(:)-t*(grad_fin%prev_phi(j,:))
    !     temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !     temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    !    ! gradtd=gradtd+dot_product(grad_fin%vars(j,:),((-1)*grad_fin%vars(j,:)))
    ! end do



   ! Normal GD
    do j=2, ndet
        temp_zom(j)%phi(:)=zstore(j)%phi(:)-(t*(grad_fin%vars(j,:)))
        temp_zom(j)%sin=sin(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
        temp_zom(j)%cos=cos(cmplx(temp_zom(j)%phi,0.0d0,kind=8))
    end do


    ! GDflg='y'
    
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
    ! GDflg='y'
    ! stop

end subroutine zombie_alter


END MODUle gradient_descent