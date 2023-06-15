MODULE grad_calc

    use globvars
    use alarrays
    use ham 
    
    contains

    !##############################################################################################################################

    !Level 1 routine to calcualate hamiltonian and overlap gradient and hessian matrix
    !subroutine to calcualte gradient of w.r.t to a specific zombie state
    subroutine gradient_zs(haml,zstore,elecs,an_cr,an2_cr2,state,orb,cmplt,strt)

        implicit none

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::state,orb,strt
        integer,dimension(:),intent(in)::cmplt
        integer,dimension(ndet)::cmplt_2
        integer::j,loop_num
        
        if (errorflag .ne. 0) return

        loop_num=0
        do j=1,ndet
            if(cmplt(j).eq.0)then
                loop_num=loop_num+1
                cmplt_2(loop_num)=j 
            end if 
        end do 

        if (loop_num.eq.0) return

        if(orb.eq.0)then
            call haml_ovrlp_grad_comb(haml,zstore,elecs,an_cr,an2_cr2,state,loop_num,cmplt_2,strt)
        else
            call haml_ovrlp_grad_comb_one_elec(haml,zstore,elecs,an_cr,an2_cr2,state,orb,loop_num,cmplt_2)
        end if 
       
        return

    end subroutine gradient_zs

    !##############################################################################################################################


    subroutine haml_ovrlp_grad_comb(haml,zstore,elecs,an_cr,an2_cr2,state,loops,cmplt,strt)

        implicit none

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::state,loops,strt
        integer,dimension(:),intent(in)::cmplt
        real(kind=8)::h1etot,h2etot
        real(kind=8),dimension(0:2*norb)::z1d
        integer::j,k,bgn,end

        if (errorflag .ne. 0) return

        ! if(strt.eq.0)then
            bgn=1
            end=norb
        ! else if(strt.eq.1)then 
        !     bgn=1
        !     end=norb/2
        ! elseif(strt.eq.2)then  
        !     bgn=(norb/2)+1
        !     end=norb
        ! end if
        
        !$omp parallel do collapse(2) private(j,k,h1etot,h2etot,z1d) shared(cmplt,zstore,haml,state,an_cr,an2_cr2,elecs)
        do j=bgn,end 
            do k=1,loops
                z1d(0:2*norb)=zstore(state)%val(0:2*norb)
                z1d(j)=zstore(state)%cos(j)
                z1d(j+norb)=zstore(state)%sin(j)*(-1)
                if(cmplt(k).ne.state)then
                    haml%diff_ovrlp(state,j,cmplt(k))=overlap_1(z1d,zstore(cmplt(k))%val)
                    h1etot = haml_vals(z1d,zstore(cmplt(k))%val,an_cr%ham,elecs%h1ei,elecs%h1_num)
                    h2etot = haml_vals(z1d,zstore(cmplt(k))%val,an2_cr2%ham,elecs%h2ei,elecs%h2_num)
                    haml%diff_hjk(state,j,cmplt(k))=h1etot+(0.5*h2etot)+(haml%diff_ovrlp(state,j,cmplt(k))*elecs%hnuc)
                else 
                    haml%diff_ovrlp(state,j,state)=0.0 
                    h1etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an_cr%diff(j),elecs%h1ei,an_cr%dcnt(0:,j))
                    h2etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an2_cr2%diff(j),elecs%h2ei,an2_cr2%dcnt(0:,j))
                    haml%diff_hjk(state,j,state)=h1etot+(0.5*h2etot)
                end if 
            end do 
        end do 
        !$omp end parallel do  
        
        return 


    end subroutine haml_ovrlp_grad_comb

    subroutine haml_ovrlp_grad_comb_one_elec(haml,zstore,elecs,an_cr,an2_cr2,state,orb,loops,cmplt)

        implicit none

        type(hamiltonian), intent(inout)::haml 
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::state,loops,orb
        integer,dimension(:),intent(in)::cmplt
        real(kind=8)::h1etot,h2etot
        real(kind=8),dimension(0:2*norb)::z1d
        integer::k

        if (errorflag .ne. 0) return

        z1d(0:2*norb)=zstore(state)%val(0:2*norb)
        z1d(orb)=zstore(state)%cos(orb)
        z1d(orb+norb)=zstore(state)%sin(orb)*(-1)
        !$omp parallel do private(k,h1etot,h2etot) shared(cmplt,zstore,haml,state,an_cr,an2_cr2,elecs,orb,z1d)
        do k=1,loops
            if(cmplt(k).ne.state)then 
                haml%diff_ovrlp(state,orb,cmplt(k))=overlap_1(z1d,zstore(cmplt(k))%val)
            
                h1etot = haml_vals(z1d,zstore(cmplt(k))%val,an_cr%ham,elecs%h1ei,elecs%h1_num)
                h2etot = haml_vals(z1d,zstore(cmplt(k))%val,an2_cr2%ham,elecs%h2ei,elecs%h2_num)
                haml%diff_hjk(state,orb,cmplt(k))=h1etot+(0.5*h2etot)+(haml%diff_ovrlp(state,orb,cmplt(k))*elecs%hnuc)
                
            else 
                haml%diff_ovrlp(state,orb,state)=0.0 
                h1etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an_cr%diff(orb),elecs%h1ei,an_cr%dcnt(0:,orb))
                h2etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an2_cr2%diff(orb),elecs%h2ei,an2_cr2%dcnt(0:,orb))
                haml%diff_hjk(state,orb,state)=h1etot+(0.5*h2etot)
            end if 
        end do 
        !$omp end parallel do 

        return 


    end subroutine haml_ovrlp_grad_comb_one_elec

    subroutine sub_matrices(haml,state,strt)

        implicit none 
        type(hamiltonian), intent(inout)::haml 
        integer,intent(in)::state,strt!orb
        real(kind=8),allocatable,dimension(:,:,:)::temp2 
        integer::ierr,k,l,j,p,bgn,end!,orbsrt,orblim

        if (errorflag .ne. 0) return

        ierr=0
        if(strt.eq.0)then
            bgn=1
            end=norb
        else if(strt.eq.1)then 
            bgn=1
            end=norb/2
        elseif(strt.eq.2)then  
            bgn=(norb/2)+1
            end=norb
        end if
        

        ! if(orb.eq.0)then
        !     orbsrt=1
        !     orblim=norb
        !     haml%diff_ov_dov(state,:,:,:)=0.0
        !     haml%diff_in_dhjk(state,:,:,:)=0.0
        !     haml%diff_invh(state,:,:,:)=0.0
        ! else 
        !     orbsrt=orb
        !     orblim=orb
        !     haml%diff_ov_dov(state,:,orb,:)=0.0
        !     haml%diff_in_dhjk(state,:,orb,:)=0.0
        !     haml%diff_invh(state,:,orb,:)=0.0
        ! end if

        haml%diff_ov_dov(state,:,:,:)=0.0
        haml%diff_in_dhjk(state,:,:,:)=0.0
        haml%diff_invh(state,:,:,:)=0.0

        allocate(temp2(ndet,norb,ndet),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in occupancy diff_inverse temp array allocation . ierr had value ", ierr
            errorflag=1
            return
        end if
        temp2=0.0

        do j=1,norb
        ! do j=orbsrt,orblim
            do k=1,ndet
                do l=1, ndet
                    if(l.ne.state)then
                        haml%diff_ov_dov(state,k,j,l)=(haml%ovrlp(k,state)*&
                            haml%diff_ovrlp(state,j,l))/abs(haml%ovrlp(k,state))

                        haml%diff_in_dhjk(state,k,j,l)=(haml%inv(k,state)*haml%diff_hjk(state,j,l))

                        temp2(k,j,l)=(haml%inv(k,l))*haml%diff_ovrlp(state,j,l)
                    else
                        do p=1,ndet
                            haml%diff_ov_dov(state,k,j,l)=haml%diff_ov_dov(state,k,j,l)+&
                            (haml%ovrlp(k,p)*haml%diff_ovrlp(state,j,p))/abs(haml%ovrlp(k,p))

                            haml%diff_in_dhjk(state,k,j,l)=haml%diff_in_dhjk(state,k,j,l)+&
                            (haml%inv(k,p)*haml%diff_hjk(state,j,p))

                            temp2(k,j,l)= temp2(k,p,l)+(haml%inv(k,p)*haml%diff_ovrlp(state,j,p))
                        end do 
                    end if 
                end do
            end do 
    
            do k=1, ndet
                do l=1, ndet
                    do p=1,ndet
                        haml%diff_invh(state,k,j,l)=haml%diff_invh(state,k,j,l)+(temp2(k,j,p)*haml%kinvh(p,l))
                    end do 
                    haml%diff_invh(state,k,j,l)=haml%diff_invh(state,k,j,l)*(-1)
                end do  
            end do
        end do


        deallocate(temp2)
        
        return

    end subroutine sub_matrices 


    !##############################################################################################################################

    !##############################################################################################################################

    !Level 4 routine used to calcualte individual hamiltonian elements for gradient and hessian matrices

    real(kind=8) function haml_vals_mod(z1d,z2d,ops,el,el_num)

        implicit none 
        real(kind=8),dimension(0:),intent(in)::z1d,z2d
        real(kind=8),dimension(:),intent(in)::el
        integer,dimension(0:),intent(in)::el_num
        type(oprts_2),intent(in)::ops
        real(kind=8)::ov
        integer::j,k,len

        if (errorflag .ne. 0) return
   
        len=int(el_num(0))
     
        haml_vals_mod=0.0
        !!$omp parallel do reduction(+:haml_vals_mod) private(j,k,ov) shared(ops,z1d,z2d,el) 
        !!$omp do  
        do j=1,len
            ov=1.0
            !!$omp do simd reduction(*:ov)
            do k=1, norb
                ov=ov*((z1d(k)*z2d(ops%alive(k,j))*ops%neg_alive(k,j))+(z1d(k+norb)*z2d(ops%dead(k,j))*ops%neg_dead(k,j))) 
            end do
            !!$omp end do simd
            haml_vals_mod=haml_vals_mod+(ov*el(el_num(j)))
        end do
        !!$omp end do 
        !!$omp end parallel do
      
        return 
      
    end function haml_vals_mod

    !##############################################################################################################################

  
    
    subroutine grad_zom_setup(zstore,grad_fin,elect,an_cr,an2_cr2,state)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(oprts_2),intent(in)::an_cr,an2_cr2
        integer,intent(in)::state
        integer::j,k,l
        real(kind=8)::ov

        if (errorflag .ne. 0) return
        
        !$omp parallel do private(j,k,l,ov) shared(zstore,grad_fin,elect,an_cr,an2_cr2,state)
        do l=1,ndet
            do j=1,elect%h1_num
                ov=1.0
                do k=1, norb
                    ov=ov*((zstore(state)%val(k)*zstore(l)%val(an_cr%alive(k,j))*an_cr%neg_alive(k,j))+&
                    (zstore(state)%val(k+norb)*zstore(l)%val(an_cr%dead(k,j))*an_cr%neg_dead(k,j))) 
                end do
                grad_fin%one_elec(l,1,j)=(ov*elect%h1ei(j))
            end do
            
            do j=1,elect%h2_num
                ov=1.0
                do k=1, norb
                    ov=ov*((zstore(state)%val(k)*zstore(l)%val(an2_cr2%alive(k,j))*an2_cr2%neg_alive(k,j))+&
                    (zstore(state)%val(k+norb)*zstore(l)%val(an2_cr2%dead(k,j))*an2_cr2%neg_dead(k,j))) 
                end do
                grad_fin%two_elec(l,1,j)=(ov*elect%h2ei(j))
            end do
          
        end do
        !$omp end parallel do
        return 

    end subroutine grad_zom_setup

    subroutine grad_zom_div(zstore,grad_fin,haml,elect,an_cr,an2_cr2,state,orb)

        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(hamiltonian),intent(inout)::haml
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        integer,intent(in)::state,orb
        type(oprts_2),intent(in)::an_cr,an2_cr2
        integer::j,l,k
        real(kind=8)::ov

        if (errorflag .ne. 0) return
        !$omp parallel do private(j,k,l,ov) shared(zstore,grad_fin,haml,elect,an_cr,an2_cr2,state,orb)
        do l=1,ndet
            grad_fin%ovrlp_div(l)=haml%ovrlp(state,l)/(zstore(state)%val(orb)*zstore(l)%val(orb)+&
            zstore(state)%val(orb+norb)*zstore(l)%val(orb+norb))
            do j=1,elect%h1_num
                if(grad_fin%one_elec(l,1,j).ne.0)then
                    grad_fin%one_elec(l,2,j)=grad_fin%one_elec(l,1,j)/((zstore(state)%val(orb)*zstore(l)%val(an_cr%alive(orb,j))&
                    *an_cr%neg_alive(orb,j))+(zstore(state)%val(orb+norb)*zstore(l)%val(an_cr%dead(orb,j))*an_cr%neg_dead(orb,j)))
                else
                    ov=1.0
                    do k=1, norb
                        if(k.eq.orb) cycle
                        ov=ov*((zstore(state)%val(k)*zstore(l)%val(an_cr%alive(k,j))*an_cr%neg_alive(k,j))+&
                        (zstore(state)%val(k+norb)*zstore(l)%val(an_cr%dead(k,j))*an_cr%neg_dead(k,j))) 
                    end do
                    grad_fin%one_elec(l,2,j)=ov*elect%h1ei(j)
                end if
            end do
            
            do j=1,elect%h2_num
                if(grad_fin%two_elec(l,1,j).ne.0)then
                    grad_fin%two_elec(l,2,j)=grad_fin%two_elec(l,1,j)/((zstore(state)%val(orb)*zstore(l)%val(an2_cr2%alive(orb,j))&
                    *an2_cr2%neg_alive(orb,j))+(zstore(state)%val(orb+norb)*zstore(l)%val(an2_cr2%dead(orb,j))&
                    *an2_cr2%neg_dead(orb,j)))
                else
                    ov=1.0
                    do k=1, norb
                        if(k.eq.orb) cycle
                        ov=ov*((zstore(state)%val(k)*zstore(l)%val(an2_cr2%alive(k,j))*an2_cr2%neg_alive(k,j))+&
                        (zstore(state)%val(k+norb)*zstore(l)%val(an2_cr2%dead(k,j))*an2_cr2%neg_dead(k,j)))
                    end do
                    grad_fin%two_elec(l,2,j)=ov*elect%h2ei(j)
                end if 
            end do
        end do 
        !$omp end parallel do

        return 

    end subroutine grad_zom_div

    subroutine grad_zom_grad(zstore,grad_fin,haml,elect,state,orb,an_cr,an2_cr2)
            
        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(grad),intent(inout)::grad_fin
        type(hamiltonian),intent(inout)::haml
        type(elecintrgl),intent(in)::elect
        integer,intent(in)::state,orb
        type(oprts),intent(in)::an_cr,an2_cr2
        real(kind=8)::dead,alive,h1etot,h2etot
        integer::j,l
        if (errorflag .ne. 0) return
    
        alive=zstore(state)%cos(orb)
        dead=zstore(state)%sin(orb)*(-1)

        !$omp do private(h1etot,h2etot,j,l) shared(grad_fin,haml,elect,an_cr,an2_cr2,zstore,alive,dead)
        do l=1,ndet
            if(l.ne.state)then 
                haml%diff_ovrlp(state,orb,l)=grad_fin%ovrlp_div(l)*(alive*zstore(l)%val(orb)+dead*zstore(l)%val(orb+norb))
                h1etot=0.0
                do j=1,elect%h1_num
                    h1etot=h1etot+grad_fin%one_elec(l,2,j)*((alive*zstore(l)%val(an_cr%ham%alive(orb,j))*&
                    an_cr%ham%neg_alive(orb,j))+(dead*zstore(l)%val(an_cr%ham%dead(orb,j))*an_cr%ham%neg_dead(orb,j)))
                end do
                h2etot=0.0
                do j=1,elect%h2_num
                    h2etot=h2etot+grad_fin%two_elec(l,2,j)*((alive*zstore(l)%val(an2_cr2%ham%alive(orb,j))*&
                    an2_cr2%ham%neg_alive(orb,j))+(dead*zstore(l)%val(an2_cr2%ham%dead(orb,j))*an2_cr2%ham%neg_dead(orb,j)))  
                end do
                haml%diff_hjk(state,orb,l)=h1etot+0.5*h2etot+haml%diff_ovrlp(state,orb,l)*elect%hnuc
              
            else 
                haml%diff_ovrlp(state,orb,l)=0.0
                h1etot=0.0
                h2etot=0.0
                h1etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an_cr%diff(orb),elect%h1ei,an_cr%dcnt(0:,orb))
                h2etot = haml_vals_mod(zstore(state)%val,zstore(state)%val,an2_cr2%diff(orb),elect%h2ei,an2_cr2%dcnt(0:,orb))
                haml%diff_hjk(state,orb,state)=h1etot+(0.5*h2etot)
            end if
            
            
        end do 
        !$omp end do 
        return 
    
    end subroutine grad_zom_grad

    subroutine grad_zom_new_vals(zstore,zs_diff,elect,grad_fin,haml,state,orb,an_cr,an2_cr2,store_one,store_two)
            
        implicit none
        type(zombiest),dimension(:),intent(in)::zstore
        type(zombiest),intent(in)::zs_diff
        type(grad),intent(inout)::grad_fin
        type(hamiltonian),intent(inout)::haml
        type(elecintrgl),intent(in)::elect
        integer,intent(in)::state,orb
        type(oprts_2),intent(in)::an_cr,an2_cr2
        real(kind=8),dimension(:,:),intent(inout)::store_one,store_two
        real(kind=8)::dead,alive,h1etot,h2etot
        integer::j,l
        if (errorflag .ne. 0) return
     
        alive=zs_diff%sin(orb)
        dead=zs_diff%cos(orb)

        do l=1,ndet
            if(state.ne.l)then
                haml%ovrlp(state,l)=grad_fin%ovrlp_div(l)*(alive*zstore(l)%val(orb)+dead*zstore(l)%val(orb+norb))
                haml%ovrlp(l,state)=haml%ovrlp(state,l)
                h1etot=0.0
                do j=1,elect%h1_num
                    store_one(l,j)=grad_fin%one_elec(l,2,j)*((alive*zstore(l)%val(an_cr%alive(orb,j))*&
                    an_cr%neg_alive(orb,j))+(dead*zstore(l)%val(an_cr%dead(orb,j))*an_cr%neg_dead(orb,j)))
                    h1etot=h1etot+store_one(l,j)
                end do
                h2etot=0.0
                do j=1,elect%h2_num
                    store_two(l,j)=grad_fin%two_elec(l,2,j)*((alive*zstore(l)%val(an2_cr2%alive(orb,j))*&
                    an2_cr2%neg_alive(orb,j))+(dead*zstore(l)%val(an2_cr2%dead(orb,j))*an2_cr2%neg_dead(orb,j)))  
                    h2etot=h2etot+store_two(l,j)
                end do
            else 
                haml%ovrlp(state,l)=1.0
                h1etot=0.0
                do j=1,elect%h1_num
                    store_one(l,j)=grad_fin%one_elec(l,2,j)*((alive*zs_diff%val(an_cr%alive(orb,j))*&
                    an_cr%neg_alive(orb,j))+(dead*zs_diff%val(an_cr%dead(orb,j))*an_cr%neg_dead(orb,j)))
                    h1etot=h1etot+store_one(l,j)
                end do
                h2etot=0.0
                do j=1,elect%h2_num
                    store_two(l,j)=grad_fin%two_elec(l,2,j)*((alive*zs_diff%val(an2_cr2%alive(orb,j))*&
                    an2_cr2%neg_alive(orb,j))+(dead*zs_diff%val(an2_cr2%dead(orb,j))*an2_cr2%neg_dead(orb,j)))  
                    h2etot=h2etot+store_two(l,j)
                end do
            end if
            
            haml%hjk(state,l)=h1etot+0.5*h2etot+(haml%ovrlp(state,l)*elect%hnuc)
            haml%hjk(l,state)=haml%hjk(state,l)
        end do 

        return 
    
    end subroutine grad_zom_new_vals

END MODULE grad_calc