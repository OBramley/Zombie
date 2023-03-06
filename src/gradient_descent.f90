MODULE gradient_descent

    use globvars
    use alarrays
    use ham
    use grad_d
    use imgtp
    use outputs
    use infnan_mod
    ! use ieee_arithmetic
    contains

    ! Subroutine to calcualte Hamiltonian elements combines 1st and 2nd electron integral calcualations so remove double 
    ! calcualiton of certain results. This minimises the (slow) applicaiton of the creation and annihilaiton operators

    ! Subroutine that controls and calcualtes all of the hamiltonian variables 
    subroutine he_full_row(haml,zstore,elecs,size,an_cr,an2_cr2,diff_state)

        implicit none 

        type(hamiltonian), intent(inout)::haml
        type(zombiest),dimension(:),intent(in)::zstore
        type(elecintrgl),intent(in)::elecs
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::size,diff_state
        integer, allocatable,dimension(:)::IPIV1
        real(kind=8),allocatable,dimension(:)::WORK1
        
        integer::ierr


        if (errorflag .ne. 0) return
        ierr=0
        
        haml%ovrlp(:,diff_state)=ovrlp_column(zstore(diff_state)%val,zstore,diff_state)
        haml%ovrlp(diff_state,:)=haml%ovrlp(:,diff_state)
       
        haml%hjk(:,diff_state)=haml%ovrlp(:,diff_state)*elecs%hnuc 
        call haml_column(haml%hjk(:,diff_state),zstore(diff_state)%val,zstore,an_cr%ham,an2_cr2%ham,elecs,1)
        haml%hjk(diff_state,:)=haml%hjk(:,diff_state)
        haml%inv=haml%ovrlp
        allocate(WORK1(size),IPIV1(size),stat=ierr)
        if (ierr/=0) then
            write(0,"(a,i0)") "Error in IPIV or WORK1 vector allocation . ierr had value ", ierr
            errorflag=1
        end if 

        Call dgetrf(size, size, haml%inv, size, IPIV1, ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if
        if (ierr==0) call dgetri(size,haml%inv,size,IPIV1,WORK1,size,ierr)
        if (ierr/=0) then
            write(0,"(a,i0)")"Error in DGETRF",ierr
        end if

        deallocate(WORK1,IPIV1)

        call DGEMM("N","N",size,size,size,1.d0,haml%inv,size,haml%hjk,size,0.d0,haml%kinvh,size)


    end subroutine he_full_row


    subroutine grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvec,grad_fin,en,orb)

        implicit none 

        type(zombiest),dimension(:),intent(in)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvec
        type(hamiltonian),intent(inout)::haml
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(in)::pick,orb
        DOUBLE PRECISION, external::ZBQLU01
        type(energy),intent(inout)::en

        if (errorflag .ne. 0) return
        
        if(grad_fin%grad_avlb(pick).eq.0)then
            dvec(1)%d_diff(:,pick,:)=0
            call gradient_zs(haml,zstore,elect,an_cr,an2_cr2,pick,orb)
        else if(grad_fin%grad_avlb(pick).eq.1)then
            dvec(1)%d_diff(:,pick,:)=0
            call imgtime_prop(dvec,en,haml,pick,0)
        end if
      
        call final_grad(dvec(1),haml,grad_fin,pick,orb,grad_fin%grad_avlb(pick))
           
        en%erg=0
        en%t=0
        if( grad_fin%grad_avlb(pick).eq.1)then 
            grad_fin%grad_avlb(pick)=2
        else if( grad_fin%grad_avlb(pick).eq.0)then 
            grad_fin%grad_avlb(pick)=1
        end if 
       

        return

    end subroutine grad_calc

    subroutine orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,temp_zom,&
        epoc_cnt,alphain,b,picker,maxloop,an_cr,an2_cr2,rjct_cnt_in) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs,temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(inout)::epoc_cnt,rjct_cnt_in
        integer,intent(in)::maxloop
        real(kind=8),intent(in)::b,alphain
        integer,dimension(ndet-1),intent(inout)::picker
        integer::rjct_cnt,next,acpt_cnt,pick,pickorb,rjct_cnt2,loops,d_diff_flg,lralt_zs
        real(kind=8)::t,fxtdk,l2_rglrstn
        integer,dimension(ndet-1)::rsrtpass
        integer,dimension(norb)::pickerorb
        integer,dimension(norb)::chng_trk2
        logical::nanchk
        character(len=4)::ergerr
        integer::j,k,l,n
        
    
        lralt_zs=0    ! power alpha is raised to 
        chng_trk2=0 !stores which orbitals in the ZS have changed 
        rjct_cnt=0 !tracks how many rejections 
        acpt_cnt=0  !counts how many ZS have been changed
        ! l2_rglrstn=1! !L2 regularisation lambda paramter
        rjct_cnt2=0
        loops=0
        ! d_diff_flg=1
        pickerorb=scramble_norb(norb)
        do while(rjct_cnt2.lt.(norb*100))
            loops=loops+1
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent steps &
                        &   | Orbitals altered '
            do j=1,(ndet-1)
               
                grad_fin%current_erg=grad_fin%prev_erg
                pick=picker(j)
                ! pickerorb=scramble_norb(norb)
                !pickerorb=(/1,2,3,4,5,6,7,8,9,10/)
                chng_trk2=0
                acpt_cnt=0
                do n=1,norb
                   
                    rjct_cnt=0
                    pickorb=pickerorb(n)
                    lralt_zs=0
                    t=b*(alphain**lralt_zs)
                    
                    do while(t.gt.(1.0d-13))
                    
                        nanchk=.false.
                        if(is_nan(grad_fin%vars(pick,pickorb)).eqv..true.)then
                            call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,pickorb)
                        end if 
                    
                        ! Setup temporary zombie state
                        temp_zom=zstore
                        temp_zom(pick)%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*(grad_fin%vars(pick,pickorb)))
                        ! temp_zom(pick)%phi(pickorb)=temp_zom(pick)%phi(pickorb)+&
                                            ! l2_rglrstn*((grad_fin%vars(pick,pickorb))*grad_fin%vars(pick,pickorb))
                        if(is_nan(temp_zom(pick)%phi(pickorb)).eqv..true.)then
                            t=(1.0d-14)
                        end if
                        ! temp_zom(pick)%sin(pickorb)=sin(cmplx(temp_zom(pick)%phi(pickorb),0.0d0,kind=8))
                        ! temp_zom(pick)%cos(pickorb)=cos(cmplx(temp_zom(pick)%phi(pickorb),0.0d0,kind=8))
                        temp_zom(pick)%sin(pickorb)=sin(temp_zom(pick)%phi(pickorb))
                        temp_zom(pick)%cos(pickorb)=cos(temp_zom(pick)%phi(pickorb))
                        temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                        temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos
                        temp_ham=haml
                        call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                        ! call he_full_row(temp_ham,temp_zom,elect,pick,ndet,occupancy_2an,occupancy_an_cr,occupancy_an)
                        ! Imaginary time propagation for back tracing
                        
                        en%erg=0
                        en%t=0
                        call imgtime_prop(temp_dvecs,en,temp_ham,0,0)
                        ! fxtdk=real(en%erg(1,timesteps+1))
                        fxtdk=en%erg(1,timesteps+1)
                
                        if(is_nan(fxtdk).eqv..true.)then 
                            ergerr='NaN ' 
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb 
                            nanchk=.true.
                            call emergency(zstore,haml,grad_fin,dvecs,temp_zom,temp_ham,temp_dvecs,en,an_cr,an2_cr2)
                            call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)
                            temp_zom(pick)%sin(:)=sin(temp_zom(pick)%phi(:))
                            temp_zom(pick)%cos(:)=cos(temp_zom(pick)%phi(:))
                            temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                            temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos
                            temp_ham=haml
                            call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                            fxtdk=en%erg(1,timesteps+1)
                            if(is_nan(fxtdk).eqv..true.)then 
                                nanchk=.true.
                                t=(1.0d-14)
                                rjct_cnt=1
                            else 
                                nanchk=.False.
                                write(0,"(a)") "Error corrected"
                            end if
                        else if(is_posinf(fxtdk).eqv..true.)then
                            ergerr='+NaN'  
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb  
                            nanchk=.true.
                            call emergency(zstore,haml,grad_fin,dvecs,temp_zom,temp_ham,temp_dvecs,en,an_cr,an2_cr2)
                            call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)
                            temp_zom(pick)%sin(:)=sin(temp_zom(pick)%phi(:))
                            temp_zom(pick)%cos(:)=cos(temp_zom(pick)%phi(:))
                            temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                            temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos
                            temp_ham=haml
                            call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                            fxtdk=en%erg(1,timesteps+1)
                            if(is_nan(fxtdk).eqv..true.)then 
                                nanchk=.true.
                                t=(1.0d-14)
                                rjct_cnt=1
                            else 
                                nanchk=.False.
                                write(0,"(a)") "Error corrected"
                            end if
                        else if(is_neginf(fxtdk).eqv..true.)then
                            ergerr='-NaN'  
                            write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                            " for zombie state ", pick, ",orbital ", pickorb 
                            nanchk=.true.
                            call emergency(zstore,haml,grad_fin,dvecs,temp_zom,temp_ham,temp_dvecs,en,an_cr,an2_cr2)
                            call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)
                            temp_zom(pick)%sin(:)=sin(temp_zom(pick)%phi(:))
                            temp_zom(pick)%cos(:)=cos(temp_zom(pick)%phi(:))
                            temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                            temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos
                            temp_ham=haml
                            call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                            fxtdk=en%erg(1,timesteps+1)
                            if(is_nan(fxtdk).eqv..true.)then 
                                nanchk=.true.
                                t=(1.0d-14)
                                rjct_cnt=1
                            else 
                                nanchk=.False.
                                write(0,"(a)") "Error corrected"
                            end if
                        end if
                        
                        if(nanchk.eqv..false.)then
                            ! Check if energy is lower and accept or reject
                            if(fxtdk.lt.grad_fin%prev_erg)then
                                t=(1.0d-14)
                                acpt_cnt=acpt_cnt+1
                                zstore=temp_zom
                                dvecs=temp_dvecs
                                chng_trk2(acpt_cnt)=pickorb
                                rjct_cnt=0
                                rjct_cnt2=0
                                haml=temp_ham
                                grad_fin%grad_avlb=0
                                grad_fin%vars=0.0
                                haml%diff_hjk=0
                                haml%diff_invh=0
                                haml%diff_ovrlp=0
                                dvecs(1)%d_diff=0
                                grad_fin%prev_erg=fxtdk
                                rjct_cnt_in=0
                                ! grad_fin%prev_mmntm(pick,:)=zstore(pick)%phi(:)
                            else 
                                lralt_zs=lralt_zs+1
                                t=b*(alphain**lralt_zs)
                                rjct_cnt=rjct_cnt+1
                            end if
                        end if
                    end do
                    
                    if((rjct_cnt.eq.0).and.(n.ne.norb))then
                        grad_fin%grad_avlb=0
                        call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,pickerorb(n))
                    end if
                end do

                if(j.eq.(ndet-1))then 
                    epoc_cnt=epoc_cnt+1
                    picker=scramble(ndet-1)
                    next=picker(1)
                    pickerorb=scramble_norb(norb)
                    lralt_zs=0
                    !Every 100 epoc brings phi values back within the normal 0-2pi range
                    if(modulo(epoc_cnt,100).eq.0)then 
                        do k=2,ndet 
                            do l=1,norb 
                                if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
                                    zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
                                    if((zstore(k)%phi(l).gt.2*pirl))then
                                        zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
                                    else if((zstore(k)%phi(l).lt.0))then
                                        zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
                                    end if
                                end if
                            end do
                        end do
                    end if
                else 
                    next=picker(j+1)
                    pickerorb=scramble_norb(norb)
                end if

                if(grad_fin%grad_avlb(next).eq.1)then  
                    do k=1,norb
                        if(is_nan(grad_fin%vars(next,k)).eqv..true.)then
                            grad_fin%grad_avlb(next)=0
                            exit 
                        end if 
                    end do
                end if

                if(loops.lt.maxloop)then
                    grad_fin%grad_avlb=0
                    ! if(grad_fin%grad_avlb(next).eq.0)then 
                        call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,pickerorb(1))
                    ! end if
                else
                    if(rjct_cnt_in.eq.0)then
                        grad_fin%grad_avlb=0
                    end if
                    !Set up gradients for next pass
                    ! if(grad_fin%grad_avlb(next).eq.0)then 
                    call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,0)
                    ! end if
                end if

                if(acpt_cnt.gt.0)then
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,*(i0:','))") '       ', pick,'              ', &
                    grad_fin%current_erg,'               ',grad_fin%prev_erg,'            ',chng_trk2(1:acpt_cnt) 
                    rsrtpass(pick)=0
                    zstore(pick)%update_num=zstore(pick)%update_num+1
                    call zombiewriter(zstore(pick),pick,rsrtpass(pick))
                else 
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,i0)") '       ', pick,'              ', &
                    grad_fin%current_erg,'               ',grad_fin%prev_erg,'            ',0
                    rjct_cnt2=rjct_cnt2+1

                end if
                
            end do

            if(rjct_cnt2.gt.(ndet-1)*2)then
                grad_fin%grad_avlb=0
                grad_fin%vars=0.0
                haml%diff_hjk=0
                haml%diff_ovrlp=0
                haml%diff_invh=0
                haml%diff_ov_dov=0
                haml%diff_in_dhjk=0
                dvecs(1)%d_diff=0
            end if

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,pick,t,0)
            write(6,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg
            acpt_cnt=0
            if(loops.ge.maxloop)then 
                exit 
            end if
        end do

    end subroutine orbital_gd

    ! subroutine orbital_gd(zstore,grad_fin,elect,dvecs,dvecs_temp,an_cr,an2_cr2,en,haml,temp_ham,temp_zom,&
    !     epoc_cnt,alphain,b,picker,maxloop,rjct_cnt,epoc_max) 

    !     implicit none 

    !     type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
    !     type(grad),intent(inout)::grad_fin
    !     type(elecintrgl),intent(in)::elect
    !     type(dvector),dimension(:),intent(inout)::dvecs,dvecs_temp
    !     type(energy),intent(inout)::en
    !     type(hamiltonian),intent(inout)::haml,temp_ham
    !     type(oprts),intent(in)::an_cr,an2_cr2
    !     integer,intent(inout)::epoc_cnt,rjct_cnt
    !     integer,intent(in)::maxloop,epoc_max
    !     real(kind=8),intent(in)::b,alphain
    !     integer,dimension(ndet-1),intent(inout)::picker
    !     integer::next,acpt_cnt,pick,pickorb,acpt_cnt2,bloops,lralt,lralt_temp,loop_max
    !     real(kind=8)::t,fxtdk,btl,gradient_norm
    !     integer,dimension(norb)::pickerorb
    !     integer,dimension(norb)::chng_trk2
    !     integer::j,k,l,n
        
    
    !     lralt=0    ! power alpha is raised to 
    !     chng_trk2=0 !stores which orbitals in the ZS have changed 
    !     acpt_cnt2=0
    !     loop_max=20
    !     bloops=0
    !     btl=1.0d-3
    !     pickerorb=scramble_norb(norb)
        
    !     do while(rjct_cnt.lt.(ndet-1)*10+4)
    !         bloops=bloops+1
    !         write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent steps &
    !                     &   | Orbitals altered '
    !         do j=1,(ndet-1)
               
    !             pick=picker(j)
    !             chng_trk2=0
    !             acpt_cnt=0
    !             grad_fin%current_erg=grad_fin%prev_erg
    !             do n=1,norb
    !                 pickorb=pickerorb(n)
    !                 lralt_temp=lralt
    !                 gradient_norm=sqrt((grad_fin%vars(pick,pickorb))*grad_fin%vars(pick,pickorb))
    !                 do while(lralt_temp.lt.loop_max)

    !                     t=b*(alphain**lralt_temp)
                    
    !                     ! Setup temporary zombie state
                        
    !                     temp_zom=zstore
    !                     ! temp_zom(pick)%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*(grad_fin%vars(pick,pickorb)))
    !                     temp_zom(pick)%phi(pickorb)=zstore(pick)%phi(pickorb)-(t*grad_fin%vars_hess(pick,pickorb))
                       
    !                     temp_zom(pick)%sin(pickorb)=sin(temp_zom(pick)%phi(pickorb))
    !                     temp_zom(pick)%cos(pickorb)=cos(temp_zom(pick)%phi(pickorb))
    !                     temp_zom(j)%val(1:norb)=temp_zom(j)%sin
    !                     temp_zom(j)%val(norb+1:)=temp_zom(j)%cos
                       
    !                     temp_ham=haml
    !                     call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
        
    !                     ! Imaginary time propagation for back tracing
    !                     en%erg=0
    !                     en%t=0
    !                     call imgtime_prop(dvecs_temp,en,temp_ham,0,0)
    !                     fxtdk=en%erg(1,timesteps+1)
                       
    !                     ! print*,fxtdk
    !                     ! Check if energy is lower and accept or reject
    !                     ! if((fxtdk.lt.(grad_fin%prev_erg+(btl*t*(gradient_norm)))))then
    !                     if((fxtdk.lt.grad_fin%prev_erg))then
    !                         ! if((sqrt((temp_zom(pick)%phi(pickorb)*temp_zom(pick)%phi(pickorb))).le.1.0d-4).or.&
    !                         ! (sqrt((grad_fin%vars(pick,pickorb)*grad_fin%vars(pick,pickorb))).le.1.0d-4).or.&
    !                         ! (abs(fxtdk-grad_fin%prev_erg).le.1.0d-4))then
                            
    !                             dvecs=dvecs_temp
    !                             acpt_cnt=acpt_cnt+1
    !                             acpt_cnt2=acpt_cnt2+1
    !                             zstore(pick)=temp_zom(pick)
    !                             haml=temp_ham
    !                             zstore(pick)=temp_zom(pick)
    !                             chng_trk2(acpt_cnt)=pickorb
    !                             haml=temp_ham
    !                             grad_fin%grad_avlb=0
    !                             grad_fin%vars=0.0
    !                             haml%diff_hjk=0
    !                             haml%diff_ovrlp=0
    !                             haml%diff_invh=0
    !                             haml%diff_ov_dov=0
    !                             haml%diff_in_dhjk=0
    !                             dvecs(1)%d_diff=0
    !                             grad_fin%prev_erg=fxtdk
    !                             Exit      
    !                         ! end if
    !                     else 
    !                         lralt_temp=lralt_temp+1
    !                         t=b*(alphain**lralt_temp)
    !                         ! rjct_cnt=rjct_cnt+1
    !                     end if
    !                     ! end if
                            
    !                     ! lralt_temp=lralt_temp+1
                        
                        
    !                 end do
                   
    !                 if((n.ne.norb))then
    !                     grad_fin%grad_avlb(pick)=0
    !                     call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)
    !                     ! call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,pickerorb(n+1))
    !                 end if
            
    !             end do

    !             if(acpt_cnt.gt.0)then
    !                 call epoc_writer(grad_fin%prev_erg,epoc_cnt,pick,t,0)
    !                 call zombiewriter(zstore(pick),pick,0)
    !                 zstore(pick)%update_num=zstore(pick)%update_num+1
    !                 write(6,"(a,i0,a,f21.16,a,f21.16,a,*(i0:','))") '       ', pick,'              ', &
    !                 grad_fin%current_erg,'               ',grad_fin%prev_erg,'            ',chng_trk2(1:acpt_cnt) 
                    
    !             else
    !                 write(6,"(a,i0,a,f21.16,a,f21.16,a,i0)") '       ', pick,'              ', &
    !                 grad_fin%current_erg,'               ',grad_fin%prev_erg,'            ',0
    !             end if 
          
                
    !             if((j.eq.(ndet-1)))then!.or.(acpt_cnt.gt.0))then 
    !                 epoc_cnt=epoc_cnt+1
    !                 picker=scramble(ndet-1)
    !                 next=picker(1)
    !                 !Every 100 epoc brings phi values back within the normal 0-2pi range
    !                 if(modulo(epoc_cnt,100).eq.0)then 
    !                     do k=2,ndet 
    !                         do l=1,norb 
    !                             if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
    !                                 zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
    !                                 if((zstore(k)%phi(l).gt.2*pirl))then
    !                                     zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
    !                                 else if((zstore(k)%phi(l).lt.0))then
    !                                     zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
    !                                 end if
    !                             end if
    !                         end do
    !                     end do
    !                 end if
    !             else 
    !                 pickerorb=scramble_norb(norb)
    !                 next=picker(j+1)
    !                 if(grad_fin%grad_avlb(next).eq.0)then
    !                     ! call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,pickerorb(1))
    !                     call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,0)
    !                 end if 
    !             end if

                
    !         end do
    
    !         write(6,"(a,i0,a,f21.16)") "Energy after epoch no. ",epoc_cnt,": ",grad_fin%prev_erg
    !         if(acpt_cnt2.gt.0)then
    !             rjct_cnt=0
    !         else 
    !             rjct_cnt=rjct_cnt+1
    !         end if 

    !         if(bloops.ge.maxloop)then 
    !             call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,0)
    !             exit 
    !         else 
    !             pickerorb=scramble_norb(norb)
    !             call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,pickerorb(1))
    !         end if
    !         if(epoc_cnt.gt.epoc_max)then 
    !             exit 
    !         end if
    !     end do

    !     return 

    ! end subroutine orbital_gd

    subroutine full_zs_gd(zstore,grad_fin,elect,dvecs,en,haml,temp_ham,temp_zom,&
      epoc_cnt,alphain,b,picker,an_cr,an2_cr2,epoc_max) 

        implicit none 

        type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(dvector),dimension(:),intent(inout)::dvecs
        type(dvector), dimension(:), allocatable:: temp_dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        integer,intent(inout)::epoc_cnt
        integer,intent(in)::epoc_max
        real(kind=8),intent(in)::b,alphain
        integer,dimension(ndet-1),intent(inout)::picker
        integer::lralt,rjct_cnt,next,acpt_cnt,pick,orbitcnt,d_diff_flg,lralt_temp,mmntmflg,loop_max,loop_dwn
        real(kind=8)::newb,t,fxtdk,l2_rglrstn,alpha,mmnmtb
        real(kind=8),dimension(norb)::gradient_norm
        real(kind=8),dimension(ndet)::mmntm,mmntma
        integer,dimension(ndet-1)::rsrtpass
        DOUBLE PRECISION, external::ZBQLU01
        ! real::r
        logical::nanchk
        character(len=4)::ergerr
        integer::j,k,l

    
        if (errorflag .ne. 0) return


        alpha=alphain  ! learning rate reduction
        lralt=0    ! power alpha is raised to  
        newb=b
        rjct_cnt=0 !tracks how many rejections 
        t=newb*(alpha**lralt) !learning rate
        acpt_cnt=0  !counts how many ZS have been changed
        ! l2_rglrstn=1! !L2 regularisation lambda paramter
        ! orbitcnt=0
        ! mmntm=0
        ! mmntma=1
        ! mmntmflg=0
        loop_max=120
        ! loop_dwn=0 

        call allocdv(temp_dvecs,1,ndet,norb)
        temp_dvecs=dvecs

        do while(rjct_cnt.lt.200)
            t=newb*(alpha**lralt)
            write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
            &   | Learning rate | Acceptance count | Rejection count'
            
            do j=1,(ndet-1)
                
                pick=picker(j)
                lralt_temp=lralt
                do k=1,norb
                    if(is_nan(grad_fin%vars(pick,k)).eqv..true.)then
                        call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)
                        nanchk=.false. 
                        exit 
                    end if 
                end do
                ! gradient_norm=((grad_fin%vars(pick,:))*grad_fin%vars(pick,:))
                do while(lralt_temp.lt.(loop_max))
                    t=newb*(alpha**lralt_temp)
                    nanchk=.false.
                    ! mmnmtb=(t*mmntm(pick))/mmntma(pick)
                    ! Setup temporary zombie state
                    temp_zom=zstore
                    temp_zom(pick)%phi(:)=zstore(pick)%phi(:)-(t*(grad_fin%vars(pick,:)))!+&
                    ! mmnmtb*(zstore(pick)%phi(:)-grad_fin%prev_mmntm(pick,:))
                    ! temp_zom(pick)%phi(:)=temp_zom(pick)%phi(:)+l2_rglrstn*gradient_norm
                
                    do k=1,norb
                        if(is_nan(temp_zom(pick)%phi(k)).eqv..true.)then
                            temp_zom(pick)%phi=zstore(pick)%phi
                            exit
                        end if
                    end do
                    ! temp_zom(pick)%sin=sin(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
                    ! temp_zom(pick)%cos=cos(cmplx(temp_zom(pick)%phi,0.0d0,kind=8))
                    temp_zom(pick)%sin(:)=sin(temp_zom(pick)%phi(:))
                    temp_zom(pick)%cos(:)=cos(temp_zom(pick)%phi(:))
                    temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                    temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos
                    temp_ham=haml
                    call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                    
                    ! Imaginary time propagation for back tracing
                    en%erg=0
                    en%t=0
                    call imgtime_prop(temp_dvecs,en,temp_ham,0,0)
                
                    ! fxtdk=real(en%erg(1,timesteps+1))
                    fxtdk=en%erg(1,timesteps+1)
                    ! print*,fxtdk 
                    ! print*,grad_fin%vars(pick,:)
                    if(is_nan(fxtdk).eqv..true.)then 
                        nanchk=.true.
                        ergerr='NaN '
                        grad_fin%vars(pick,:)=0.0
                        grad_fin%grad_avlb(pick)=0
                        call emergency(zstore,haml,grad_fin,dvecs,temp_zom,temp_ham,temp_dvecs,en,an_cr,an2_cr2)
                        call grad_calc(haml,zstore,elect,an_cr,an2_cr2,pick,dvecs,grad_fin,en,0)
                        temp_zom(pick)%sin(:)=sin(temp_zom(pick)%phi(:))
                        temp_zom(pick)%cos(:)=cos(temp_zom(pick)%phi(:))
                        temp_zom(pick)%val(1:)=temp_zom(pick)%sin
                        temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos
                        temp_ham=haml
                        call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
                        fxtdk=en%erg(1,timesteps+1)
                        if(is_nan(fxtdk).eqv..true.)then 
                            nanchk=.true.
                        else 
                            nanchk=.False.
                        end if

                    else if(is_posinf(fxtdk).eqv..true.)then
                        nanchk=.true.
                        ergerr='+NaN'
                        grad_fin%vars(pick,:)=0.0
                        grad_fin%grad_avlb(pick)=0
                        call emergency(zstore,haml,grad_fin,dvecs,temp_zom,temp_ham,temp_dvecs,en,an_cr,an2_cr2)  
                    else if(is_neginf(fxtdk).eqv..true.)then
                        nanchk=.true.
                        ergerr='-NaN'
                        grad_fin%vars(pick,:)=0.0
                        grad_fin%grad_avlb(pick)=0 
                        call emergency(zstore,haml,grad_fin,dvecs,temp_zom,temp_ham,temp_dvecs,en,an_cr,an2_cr2)
                    end if
            
                    if(nanchk.eqv..false.)then
                        ! Check if energy is lower and accept or reject
                        if(fxtdk.lt.grad_fin%prev_erg)then
                            acpt_cnt=acpt_cnt+1
                            zstore=temp_zom
                            zstore(pick)%update_num=zstore(pick)%update_num+1
                            call zombiewriter(zstore(pick),pick,rsrtpass(pick))
                            rsrtpass(pick)=0
                            haml=temp_ham
                            dvecs=temp_dvecs
                            rjct_cnt=0
                            grad_fin%grad_avlb=0
                            grad_fin%vars=0.0
                            haml%diff_hjk=0
                            haml%diff_invh=0
                            haml%diff_ovrlp=0
                            grad_fin%prev_erg=fxtdk
                            ! d_diff_flg=0
                            ! grad_fin%prev_mmntm(pick,:)=zstore(pick)%phi(:)
                            dvecs(1)%d_diff=0
                            ! mmntma(pick)=t
                            ! mmntm(pick)=mmnmtb
                            ! loop_dwn=loop_dwn+1
                            ! if(orbitcnt.lt.0)then
                                ! orbitcnt=orbitcnt+1 
                            ! else
                                ! orbitcnt=0
                            ! end if
                            ! if((loop_max.gt.20).and.(loop_dwn.eq.5))then
                                ! loop_max=loop_max-1
                                ! loop_dwn=0
                            ! end if
                            write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',fxtdk,'            ',t,'          ',acpt_cnt,'                 ',rjct_cnt
                            Exit
                        
                        else 
                            lralt_temp=lralt_temp+1
                        end if
                    else 
                        write(0,"(a,a,a,i0,a,i0)") "Error in energy calculation which took value ",ergerr, &
                                                " for zombie state ", pick, ", on epoc ", epoc_cnt 
                        Exit
                    end if
                end do

                if(lralt_temp.ge.loop_max)then
                    rjct_cnt=rjct_cnt+1
                    write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
                grad_fin%prev_erg,'               ',fxtdk,'            ',0.0,'          ',acpt_cnt,'                 ',rjct_cnt
                    ! if(orbitcnt.ne.0)then
                        ! orbitcnt=orbitcnt+1
                    ! end if
                    ! if(loop_max.lt.35)then 
                        ! loop_max=loop_max+1
                    ! end if
                    ! loop_dwn=0 
                end if
                
            
                t=newb*(alpha**lralt)
                if(j.eq.(ndet-1))then 
                    epoc_cnt=epoc_cnt+1
                    ! if(acpt_cnt.gt.0)then 
                        picker=scramble(ndet-1)
                    ! end if
                    next=picker(1)
                    !Every 100 epoc brings phi values back within the normal 0-2pi range
                    if(modulo(epoc_cnt,100).eq.0)then 
                        do k=2,ndet 
                            do l=1,norb 
                                if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
                                    zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
                                    if((zstore(k)%phi(l).gt.2*pirl))then
                                        zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
                                    else if((zstore(k)%phi(l).lt.0))then
                                        zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
                                    end if
                                end if
                            end do
                        end do
                    end if
                else 
                    next=picker(j+1)
                end if
        
                ! if(grad_fin%grad_avlb(next).eq.1)then 
                !     d_diff_flg=1 
                ! else 
                !     d_diff_flg=0
                !     if(ZBQLU01(1).lt.0.3)then 
                !         d_diff_flg=1
                !     end if
                ! end if

                !Set up gradients for next pass
                if(grad_fin%grad_avlb(next).lt.2)then
                    call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,0)
                end if

            end do

            call epoc_writer(grad_fin%prev_erg,epoc_cnt,pick,t,0)
            write(6,"(a,i0,a,f21.16,a,f12.10,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
                grad_fin%prev_erg, ". The current learning rate is: ",t, ". ", acpt_cnt, " Zombie state(s) altered."

            
            ! If only a single ZS is being altered the learning rate is lowered to allow others to be changed
            ! if(acpt_cnt.eq.1)then
            !     if((orbitcnt.ge.0).and.(epoc_cnt.gt.100))then  
            !         call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,temp_zom,&
            !         epoc_cnt,alphain,b,picker,1,an_cr,an2_cr2)

            !         orbitcnt=-(10*ndet)
            !         lralt=0
            !     end if   

            if(rjct_cnt.gt.((ndet-1)*2)+1)then
                ! if((epoc_cnt.gt.50))then  
                    call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,temp_zom,&
                    epoc_cnt,alphain,b,picker,1,an_cr,an2_cr2,rjct_cnt)

                    ! grad_fin%grad_avlb=0
                    ! orbitcnt=-(10*ndet)
                ! end if
            end if

            ! if((epoc_cnt.gt.500).and.(orbitcnt.ge.0).and.(modulo(epoc_cnt,100).eq.0))then
            !     call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,temp_zom,&
            !         epoc_cnt,alphain,b,picker,1,an_cr,an2_cr2)
               
            !     orbitcnt=-(10*ndet)
            !     lralt=0
            ! end if

            if(rjct_cnt.gt.(ndet-1)*2)then
                grad_fin%grad_avlb=0
                grad_fin%vars=0.0
                haml%diff_hjk=0
                haml%diff_ovrlp=0
                haml%diff_invh=0
                haml%diff_ov_dov=0
                haml%diff_in_dhjk=0
                dvecs(1)%d_diff=0
            end if

          
            t=newb*(alpha**lralt)
            ! if(mmntmflg.eq.0)then 
            !     mmntm=0.4
            !     mmntmflg=1
            ! end if 
            acpt_cnt=0
            if(epoc_cnt.gt.epoc_max)then 
                exit 
            end if
        end do

        rjct_cnt=0
        call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,en,haml,temp_ham,temp_zom,&
        epoc_cnt,alphain,b,picker,epoc_max-epoc_cnt,an_cr,an2_cr2,rjct_cnt)

        return

    end subroutine full_zs_gd


    ! subroutine full_zs_gd(zstore,grad_fin,elect,dvecs,an_cr,an2_cr2,en,&
    !     haml,temp_ham,temp_zom,epoc_cnt,alphain,b,picker,epoc_max) 

    !     implicit none 

    !     type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
    !     type(grad),intent(inout)::grad_fin
    !     type(elecintrgl),intent(in)::elect
    !     type(dvector),dimension(:),intent(inout)::dvecs
    !     type(dvector),dimension(:),allocatable::dvecs_temp
    !     type(energy),intent(inout)::en
    !     type(hamiltonian),intent(inout)::haml,temp_ham
    !     type(oprts),intent(in)::an_cr,an2_cr2
    !     integer,intent(inout)::epoc_cnt,epoc_max
    !     real(kind=8),intent(in)::b,alphain
    !     integer,dimension(ndet-1),intent(inout)::picker
    !     integer::lralt,rjct_cnt,next,acpt_cnt,pick,lralt_temp,loop_max,d_cnt,btl_cnt,dvals
    !     real(kind=8)::newb,t,fxtdk,alpha,btl,d_avrg,d_chck_num,epsilon_1,epsilon_2,epsilon_3
    !     real(kind=8)::gradient_norm
    !     integer,dimension(ndet-1)::rsrtpass
    !     DOUBLE PRECISION, external::ZBQLU01
    !     integer::j,k,l,p,q
       
    
    
    !     if (errorflag .ne. 0) return

    !     if(rstrtflg.eq.'y')then
    !         rsrtpass=1
    !         rstrtflg='n'
    !     else
    !         rsrtpass=0
    !     end if 

    !     alpha=alphain  ! learning rate reduction
    !     lralt=0    ! power alpha is raised to  
    !     newb=b
    !     rjct_cnt=0 !tracks how many rejections 
    !     t=newb*(alpha**lralt) !learning rate
    !     acpt_cnt=0  !counts how many ZS have been changed
    !     loop_max=170
    !     d_cnt=0 
    !     btl_cnt=0
    !     btl=1.0d-4 !li2
    !     ! btl=1.0d-10
      
    !     epsilon_1=1.0d-2 !li2 1
    !     epsilon_2=1.0d-2 !li2 1
    !     epsilon_3=1.0d-2 !li2 0.01
    !     dvals=0
    !     call allocdv(dvecs_temp,1,ndet,norb)
    !     dvecs_temp=dvecs
       
      
    !     do while(rjct_cnt.lt.(ndet-1)*20)
           
    !         acpt_cnt=0
    !         d_avrg=0
    !         write(6,"(a)") '    Zombie state    |     Previous Energy     |    Energy after Gradient Descent step &
    !         &   | Learning rate | Acceptance count | Rejection count   |     Difference'
            
    !         do j=1,(ndet-1)
                
    !             pick=picker(j)
    !             lralt_temp=lralt
                
    !             gradient_norm=sqrt(sum((grad_fin%vars(pick,:))*grad_fin%vars(pick,:)))
    !             ! if(grad_fin%hess_sum(pick).lt.0)then 
    !             !     newb=0.1
    !             ! else 
    !             !     newb=1
    !             ! end if 
    !             do while(lralt_temp.lt.(loop_max))
                  
    !                 t=newb*(alpha**lralt_temp)
    !                 temp_zom=zstore
    !                 ! if(grad_fin%grad_avlb(pick).eq.2)then
    !                 !     temp_zom(pick)%phi(:)=zstore(pick)%phi(:)-(t*(grad_fin%vars(pick,:)))
    !                 ! else 
    !                     temp_zom(pick)%phi(:)=zstore(pick)%phi(:)-(t*(grad_fin%vars_hess(pick,:)))
    !                 ! end if

    !                 do l=1,norb 
    !                     if(is_nan( temp_zom(pick)%phi(l)).eqv..true.)then
    !                         grad_fin%grad_avlb=0
    !                         call emergency(zstore,haml,grad_fin,dvecs,temp_zom,temp_ham,dvecs_temp,en,an_cr,an2_cr2)
    !                         call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,0)
    !                         temp_zom(pick)%phi(:)=zstore(pick)%phi(:)-(t*(grad_fin%vars_hess(pick,:)))
    !                     end if
    !                 end do 
    !                 temp_zom(pick)%sin=sin(temp_zom(pick)%phi)
    !                 temp_zom(pick)%cos=cos(temp_zom(pick)%phi)
    !                 temp_zom(pick)%val(1:)=temp_zom(pick)%sin
    !                 temp_zom(pick)%val(norb+1:)=temp_zom(pick)%cos

    !                 temp_ham=haml
                   
    !                 call he_full_row(temp_ham,temp_zom,elect,ndet,an_cr,an2_cr2,pick)
        
    !                 ! Imaginary time propagation for back tracing
    !                 en%erg=0
    !                 en%t=0
    !                 call imgtime_prop(dvecs_temp,en,temp_ham,0,0)
                
    !                 fxtdk=en%erg(1,timesteps+1)
    !                 if(is_nan(fxtdk).eqv..true.)then
    !                     exit 
    !                 end if
                 
    !                 ! if((fxtdk.lt.(grad_fin%prev_erg+(btl*t*(gradient_norm)))))then
    !                     if((fxtdk.lt.(grad_fin%prev_erg)))then
                   
    !                     ! if((sqrt(sum(temp_zom(pick)%phi(:)*temp_zom(pick)%phi(:))).le.epsilon_1).or.&
    !                     ! (sqrt(sum(grad_fin%vars(pick,:)*grad_fin%vars(pick,:))).le.epsilon_2).or.&
    !                     ! (abs(fxtdk-grad_fin%prev_erg).le.epsilon_3))then
    !                         ! if((abs(fxtdk-grad_fin%prev_erg).le.0.1))then
    !                         if(fxtdk.lt.-14.871913783299298)then  !-14.871913783299298 -25.20649144705521
    !                             print*,fxtdk
    !                             print*,grad_fin%vars_hess(pick,:)
    !                             print*,grad_fin%vars(pick,:)
    !                             stop
    !                         end if 
    !                         d_avrg=d_avrg+(grad_fin%prev_erg-fxtdk)
    !                         dvecs=dvecs_temp
    !                         acpt_cnt=acpt_cnt+1
    !                         zstore=temp_zom
    !                         haml=temp_ham
    !                         zstore(pick)=temp_zom(pick)
    !                         zstore(pick)%update_num=zstore(pick)%update_num+1
    !                         call zombiewriter(zstore(pick),pick,rsrtpass(pick))
    !                         rsrtpass(pick)=0
    !                         haml=temp_ham
    !                         rjct_cnt=0
    !                         grad_fin%grad_avlb=0
    !                         grad_fin%vars=0.0
    !                         haml%diff_hjk=0
    !                         haml%diff_ovrlp=0
    !                         haml%diff_invh=0
    !                         haml%diff_ov_dov=0
    !                         haml%diff_in_dhjk=0
    !                         haml%hess_hjk=0
    !                         haml%hess_ovrlp=0
    !                         grad_fin%vars_hess=0
    !                         dvecs(1)%d_diff=0
                        
    !                         write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0,a,f21.16)") '       ', pick,'              ', &
    !                     grad_fin%prev_erg,'               ',fxtdk,'            ',t,'          ',acpt_cnt,&
    !                     '                 ',rjct_cnt, &
    !                     '          ',(fxtdk-grad_fin%prev_erg)
    !                         grad_fin%prev_erg=fxtdk
    !                         call epoc_writer(grad_fin%prev_erg,epoc_cnt,pick,t,0)
                            
    !                         Exit
    !                         ! end if 
    !                     ! end if
    !                 ! end if
    !                 else
    !                 lralt_temp=lralt_temp+1
    !                 t=newb*(alpha**lralt_temp)
    !                 end if
    !             end do

    !             if(lralt_temp.ge.loop_max)then
    !                 rjct_cnt=rjct_cnt+1
    !                 write(6,"(a,i0,a,f21.16,a,f21.16,a,f12.10,a,i0,a,i0)") '       ', pick,'              ', &
    !             grad_fin%prev_erg,'               ',fxtdk,'            ',0.0,'          ',acpt_cnt,'                 ',rjct_cnt
    !             end if
            
                
    !             if(j.eq.(ndet-1))then 
    !                 ! if(acpt_cnt.gt.0)then 
    !                     picker=scramble(ndet-1)
    !                 ! end if
    !                 next=picker(1)
    
    !                 !Every 100 epoc brings phi values back within the normal 0-2pi range
    !                 if(modulo(epoc_cnt,100).eq.0)then 
    !                     do k=2,ndet 
    !                         do l=1,norb 
    !                             if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
    !                                 zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
    !                                 if((zstore(k)%phi(l).gt.2*pirl))then
    !                                     zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
    !                                 else if((zstore(k)%phi(l).lt.0))then
    !                                     zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
    !                                 end if
    !                             end if
    !                         end do
    !                     end do
    !                 end if
    !             else 
    !                 next=picker(j+1)
    !             end if

    !             if(rjct_cnt.gt.(ndet-1)*2)then
    !                 grad_fin%grad_avlb=0
    !                 grad_fin%vars=0.0
    !                 haml%diff_hjk=0
    !                 haml%diff_ovrlp=0
    !                 haml%diff_invh=0
    !                 haml%diff_ov_dov=0
    !                 haml%diff_in_dhjk=0
    !                 dvecs(1)%d_diff=0
    !                 grad_fin%vars_hess=0
    !             end if

    !             !Set up gradients for next pass
    !             if(grad_fin%grad_avlb(next).lt.2)then
    !                 call grad_calc(haml,zstore,elect,an_cr,an2_cr2,next,dvecs,grad_fin,en,0)
    !             end if
               

    !         end do
            
           
    !         write(6,"(a,i0,a,f21.16,a,i0,a)") "Energy after epoc no. ",epoc_cnt,": ", &
    !             grad_fin%prev_erg, ".  ", acpt_cnt, " Zombie state(s) altered."
            
    !             print*,grad_fin%hess_sum
    !         epoc_cnt=epoc_cnt+1

    !         ! if(rjct_cnt.ge.(ndet-1)*10)then
    !         !     btl=0.0d0
    !         ! end if
            
    !         ! if((epoc_cnt.gt.1500).and.(modulo(epoc_cnt,250).eq.0))then!.and.(all()then
    !         ! !     ! print*,'here'
    !         !     newb=newb*0.8
    !         !     ! btl=btl*0.8
    !         !     ! epsilon1=epsilon1*0.8
    !         !     ! epsilon2=epsilon2*0.8
    !         !     ! epsilon3=epsilon3*0.8
    !         ! end if
           
    !         ! if((acpt_cnt.le.1).or.(modulo(epoc_cnt,250).eq.0))then
    !             if(acpt_cnt.eq.0)then
    !             call orbital_gd(zstore,grad_fin,elect,dvecs,dvecs_temp,an_cr,an2_cr2,en,haml,temp_ham,temp_zom,&
    !             epoc_cnt,alphain,1.0d0,picker,1,rjct_cnt,epoc_max) 
        
    !         end if 

            

    !         if(epoc_cnt.gt.epoc_max)then
    !             exit
    !         end if
    
    !     end do

    !     return 

    ! end subroutine full_zs_gd


    subroutine zombie_alter(zstore,grad_fin,haml,elect,en,dvecs,an_cr,an2_cr2)

        implicit none

        type(zombiest),dimension(:),intent(inout)::zstore
        type(grad),intent(inout)::grad_fin
        type(elecintrgl),intent(in)::elect
        type(oprts),intent(in)::an_cr,an2_cr2
        type(dvector),dimension(:),intent(inout)::dvecs
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml
        type(zombiest),dimension(:),allocatable::temp_zom
        type(hamiltonian)::temp_ham
        integer::epoc_cnt,epoc_max
        real(kind=8)::alpha,b
        integer,dimension(ndet-1)::picker

    
        integer::ierr,j,k,l
        
        if (errorflag .ne. 0) return

        epoc_cnt=0 !epoc counter
        if(rstrtflg.eq.'y')then 
            open(unit=450,file='epoc.csv',status="old",iostat=ierr)
            if(ierr/=0)then
                write(0,"(a,i0)") "Error in opening epoc file to read in. ierr had value ", ierr
                errorflag=1
                return
            end if
            do 
                read(450,*,iostat=ierr)
                if(ierr<0)then
                    close(450)
                    exit
                else if (ierr/=0) then
                    write(0,"(a,i0)") "Error in counting epocs. ierr had value ", ierr
                    errorflag=1
                    return
                end if
                epoc_cnt=epoc_cnt+1
            end do
            close(450) 

            ! call epoc_writer(grad_fin%prev_erg,epoc_cnt,chng_trk,0.0d0,1)
            call epoc_writer(grad_fin%prev_erg,epoc_cnt,0,0.0d0,1)
            ierr=0
        else
            epoc_cnt=1
        end if

        alpha=0.8  ! learning rate reduction
        b=8.0D0 !starting learning rate
        
        epoc_max=1000

        do j=1, ndet-1
            picker(j)=j+1
        end do

        call alloczs(temp_zom,ndet)
        call allocham(temp_ham,ndet,norb)
    
        if(epoc_cnt.lt.epoc_max)then
            call full_zs_gd(zstore,grad_fin,elect,dvecs,en,haml,temp_ham,temp_zom,&
            epoc_cnt,alpha,b,picker,an_cr,an2_cr2,epoc_max)
            ! call full_zs_gd(zstore,grad_fin,elect,dvecs,an_cr,an2_cr2,en,haml,temp_ham,&
            ! temp_zom,epoc_cnt,alpha,b,picker,epoc_max) 
        
            ! call orbital_gd(zstore,grad_fin,elect,dvecs,temp_dvecs,an_cr,an2_cr2,en,&
            ! haml,temp_ham,chng_trk,temp_zom,epoc_cnt,alpha,b,picker,(epoc_max-epoc_cnt))
        end if 

        !Brings phi values back within the normal 0-2pi range
        do k=2,ndet 
            do l=1,norb 
                if((zstore(k)%phi(l).gt.2*pirl).or.(zstore(k)%phi(l).lt.0))then 
                    zstore(k)%phi(l)=asin(real(zstore(k)%sin(l)))
                    if((zstore(k)%phi(l).gt.2*pirl))then
                        zstore(k)%phi(l)=zstore(k)%phi(l)-2*pirl
                    else if((zstore(k)%phi(l).lt.0))then
                        zstore(k)%phi(l)=zstore(k)%phi(l)+2*pirl
                    end if
                end if
            end do
        end do

        call deallocham(temp_ham) 
        call dealloczs(temp_zom)
        
        return

    end subroutine zombie_alter


    subroutine emergency(zstore,haml,grad_fin,dvecs,temp_zom,temp_ham,dvecs_temp,en,an_cr,an2_cr2)

        implicit none
        type(zombiest),dimension(:),intent(inout)::zstore,temp_zom
        type(grad),intent(inout)::grad_fin
        type(dvector),dimension(:),intent(inout)::dvecs,dvecs_temp
        type(energy),intent(inout)::en
        type(hamiltonian),intent(inout)::haml,temp_ham
        type(oprts),intent(in)::an_cr,an2_cr2
        real(kind=8)::r
        integer::j,k,l
        logical::checker

        ! do j=2,ndet
        !     do k=1,norb 
        !         if(is_nan(zstore(j)%phi(k)).eqv..true.)then
        !             call random_number(r)
        !             zstore(j)%phi(k)=2*pirl*r 
        !             zstore(j)%sin=sin(zstore(j)%phi(k))
        !             zstore(j)%cos=cos(zstore(j)%phi(k))
        !             zstore(j)%val(k)=zstore(j)%sin(k)
        !             zstore(j)%val(norb+k)=zstore(j)%cos(k)
        !             write(0,"(a,i0,i0)") "Error in ZS value number and orbital ", j,k
        !         end if
        !     end do 
        ! end do
        
        checker=.true.
        do j=1,ndet 
            if(is_nan(dvecs(1)%d(j)).eqv..true.)then
                dvecs(1)%d=0.0d0
                dvecs(1)%d_diff=0.0d0
                dvecs(1)%norm=0.0d0
                dvecs_temp(1)%d=0.0d0
                dvecs_temp(1)%d_diff=0.0d0
                dvecs_temp(1)%norm=0.0d0
                checker=.false.
                exit
            end if 
            do k=1,norb
                do l=1,ndet 
                    if(is_nan(dvecs(1)%d_diff(j,k,l)).eqv..true.)then
                        dvecs(1)%d=0.0d0
                        dvecs(1)%d_diff=0.0d0
                        dvecs(1)%norm=0.0d0
                        dvecs_temp(1)%d=0.0d0
                        dvecs_temp(1)%d_diff=0.0d0
                        dvecs_temp(1)%norm=0.0d0
                        checker=.false.
                        exit
                    end if
                end do 
                if(checker.eqv..false.)then
                    exit 
                end if 
            end do 
            if(checker.eqv..false.)then
                exit 
            end if 
        end do 

        call imgtime_prop(dvecs,en,haml,0,0)

        ! call  dealloc(temp_zom)

        ! call  deallocham(haml)
        ! call  deallocham(haml_temp)
        ! call  deallocgrad(grad_fin)
        ! call  deallocdv(dvecs)
        ! call  deallocdv(temp_dvecs)
        ! call  deallocerg(en)

        ! call allocham(haml,ndet,norb)
        ! call alloc(dvecs,1)
        ! call allocdv(dvecs,1,ndet,norb)
        ! call allocerg(en,1)
        ! call  allocgrad(grad_fin,ndet,norb)

        ! call hamgen(haml,zstore,elect,ndet,an_cr,an2_cr2,1)
        ! call imgtime_prop(dvecs,en,haml,diff_state,0)


        return 

    end subroutine emergency

    ! Produces a random order for the ZS to be posisbly changed
    function scramble( number_of_values ) result(out)
        
        implicit none

        integer,intent(in)::number_of_values
        integer::out(number_of_values),array(number_of_values+1)
        integer::n,m,k,j,l,jtemp
        !real::r
        DOUBLE PRECISION, external::ZBQLU01

        out=[(j,j=1,number_of_values)]
        array=[(j,j=1,number_of_values+1)]
        n=1; m=number_of_values
        do k=1,2
            do j=1,m+1
                l = n + FLOOR((m+1-n)*ZBQLU01(1))
                jtemp=array(l)
                array(l)=array(j)
                array(j)=jtemp
            end do
        end do
        n=1
        do j=1,m+1
            if(array(j).eq.1)then 
                cycle
            end if
            out(n)=array(j)
            n=n+1
        end do
        return
     end function scramble

    ! Produces a random order for the ZS to be posisbly changed
    function scramble_norb( number_of_values ) result(out)
    
        implicit none

        integer,intent(in)::number_of_values
        integer::out(number_of_values)
        integer::n,m,k,j,l,jtemp
        DOUBLE PRECISION, external::ZBQLU01
        !real::r
        out=[(j,j=1,number_of_values)]
        
        n=1; m=number_of_values
        do k=1,2
            do j=1,m
                l = n + FLOOR((m-n)*ZBQLU01(1))
                jtemp=out(l)
                out(l)=out(j)
                out(j)=jtemp

            end do
        end do

        return
    end function scramble_norb



END MODULE gradient_descent