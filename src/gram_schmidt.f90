MODULE gram_schmidt

    use mod_types
    use globvars
    use alarrays
    use readpars
    use ham 
    use outputs
    use zom 
    use imgtp
    use gradient_descent
    contains


    subroutine gram_schmidt_control(elect)

        implicit none

        type(elecintrgl),intent(in)::elect
        type(gram),dimension(:),allocatable::gramstore
        real(wp), dimension(:,:),allocatable::erg
        integer::ierr,j,k,l,epoc_max_store,p
        character(LEN=20)::filenm
        logical :: file_exists
        character(len=2)::gst_num
        real(wp),dimension(1:norb)::parts
        real(wp)::norm,tot,temp
        if(errorflag.ne.0) return

        ierr=0
        if(gramwave.eq.1)then
            write(stdout,"(a)") "The Gram Schmidt in this routine requires multiple wavefunctions to be generated."
            errorflag=1
            return
        end if

        allocate(gramstore(gramnum+1),stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error allocating zstore array ierr had value ",ierr
            errorflag=1
            return
        end if
        call allocgram(gramstore(1),1,ndet,norb)
        call allocgram(gramstore(2),2,ndet,norb)
        call allocgram(gramstore(3),3,ndet*2,norb)
        ! do j=2,gramnum+1
        !     call allocgram(gramstore(j),j,ndet_max,norb)
        ! end do
        write(stdout,"(a)") "Gram Schmidt orthogonalisation components allocated"
        if(zomgflg=='y')then
            call genzf(gramstore(1)%zstore,ndet) 
            do k=1,ndet
                call zombiewriter(gramstore(1)%zstore(k),k,gramstore(1)%zstore(k)%gram_num)
                gramstore(1)%zstore(k)%gram_num=1
            end do
            rhf_1='n'
            call genzf(gramstore(2)%zstore,ndet) 
            ! call gen_ran_zs(gramstore(2)%zstore,ndet_max) 
            gramstore(3)%zstore(1:ndet)=gramstore(2)%zstore(1:ndet)
            gramstore(3)%zstore(1+ndet:2*ndet)=gramstore(1)%zstore(1:ndet)
            do j=2,gramnum+1
                ! call gen_ran_zs(gramstore(j)%zstore,ndet_max) 
                do k=1,ndet_max
                    gramstore(j)%zstore(k)%gram_num=j
                    call zombiewriter(gramstore(j)%zstore(k),k,j)
                end do
            end do
            write(stdout,"(a)") "Zombie states generated"
        else if (zomgflg=='n') then
            call read_zombie(gramstore(1)%zstore,1)
            rhf_1='n'
            do j=2,gramnum+1
                call gen_ran_zs(gramstore(j)%zstore,ndet_max) 
                do k=1,ndet_max
                    gramstore(j)%zstore(k)%gram_num=j
                    call zombiewriter(gramstore(j)%zstore(k),k,j)
        
                end do
            end do
            ! do j=1,gramnum+1
            !     call read_zombie(gramstore(j)%zstore,gramstore(j)%zstore(k)%gram_num)
            ! end do
            write(stdout,"(a)") "Zombie read in"
        end if
 
        if(hamgflg=='y')then
            write(stdout,"(a)") "To hamiltonian gen"
           
            call hamgen(gramstore(1)%haml,gramstore(1)%zstore,elect,ndet,1)
            write(gst_num,"(i1.1)")1
            call matrixwriter(gramstore(1)%haml%hjk,ndet,"data/ham_"//trim(gst_num)//".csv")
            call matrixwriter(gramstore(1)%haml%ovrlp,ndet,"data/ovlp_"//trim(gst_num)//".csv")
            write(stdout,"(a)") "Hamiltonian 1 successfully generated"
            do j=2,2 !gramnum+1
                call hamgen(gramstore(j)%haml,gramstore(j)%zstore,elect,ndet_max,1)
                write(gst_num,"(i1.1)")j
                call matrixwriter(gramstore(j)%haml%hjk,ndet_max,"data/ham_"//trim(gst_num)//".csv")
                call matrixwriter(gramstore(j)%haml%ovrlp,ndet_max,"data/ovlp_"//trim(gst_num)//".csv")
                write(stdout,"(a,i0,a)") "Hamiltonian ",j," successfully generated"
            end do
            do j=3,3!gramnum+1
                call hamgen(gramstore(j)%haml,gramstore(j)%zstore,elect,ndet*2,1)
                write(gst_num,"(i1.1)")j
                call matrixwriter(gramstore(j)%haml%hjk,ndet*2,"data/ham_"//trim(gst_num)//".csv")
                call matrixwriter(gramstore(j)%haml%ovrlp,ndet*2,"data/ovlp_"//trim(gst_num)//".csv")
                write(stdout,"(a,i0,a)") "Hamiltonian ",j," successfully generated"
            end do
            write(stdout,"(a)") "Hamiltonian successfully generated"
        else if (hamgflg=='n') then
            do j=1,gramnum+1
                write(gst_num,"(i1.1)")j
                call read_ham(gramstore(j)%haml,ndet,"ham_"//trim(gst_num)//".csv","ovlp_"//trim(gst_num)//".csv")
            end do
            write(stdout,"(a)") "Hamiltonian read in"
        end if
        write(stdout,"(a)") "Setting up  wavefunction overlap matrices "
       
        allocate(erg(gramnum+1,timesteps+1),stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error allocating erg array ierr had value ",ierr
            errorflag=1
            return
        end if
        

        call imaginary_time(gramstore(1)%dvecs,erg(1,:),gramstore(1)%haml, ndet)
        write(stdout,"(a,i1,a,f21.16)") "State, ",1," energy: ", erg(1,timesteps+1)
        call imaginary_time(gramstore(2)%dvecs,erg(2,:),gramstore(2)%haml, ndet)
        write(stdout,"(a,i1,a,f21.16)") "State, ",2," energy: ", erg(2,timesteps+1)
        call gram_ovrlp_fill(gramstore,2)
        gramstore(3)%dvecs%d(1:ndet)=gramstore(2)%dvecs%d(1:ndet)
        gramstore(3)%dvecs%d(1+ndet:2*ndet)=(-1)*gramstore(1)%dvecs%d(1:ndet)*&
        (dot_product(gramstore(1)%dvecs%d,matmul(gramstore(2)%wf_ovrlp(1,:,:),gramstore(2)%dvecs%d))/&
        (dot_product(gramstore(1)%dvecs%d,matmul(gramstore(1)%haml%ovrlp,gramstore(1)%dvecs%d))))
        erg(3,timesteps+1)=dot_product(gramstore(3)%dvecs%d,matmul(gramstore(3)%haml%hjk,gramstore(3)%dvecs%d))
        
        ! do j=2,gramnum+1
        !     call imaginary_time(gramstore(j)%dvecs,erg(j,:),gramstore(j)%haml, ndet_max)
        !     write(stdout,"(a,i1,a,f21.16)") "State, ",j," energy: ", erg(j,timesteps+1)
        ! end do
    
        call energywriter(erg,"energy.csv",0)
        write(stdout,"(a)") "Imaginary time propagation complete"
        
        if(GDflg.eq."y")then
            call zombie_alter(gramstore(1)%zstore,gramstore(1)%haml,elect,gramstore(1)%dvecs)
        end if 
       
        GDflg='n'
        timesteps=100000
        beta=30000
        do j=1,gramnum+1
            do k=1,ndet
                call zombiewriter(gramstore(j)%zstore(k),k,j)
            end do
        end do
        do j=2,gramnum+1
            call gram_ovrlp_fill(gramstore,j)
        end do
        do j=1,gramnum+1
            call hamgen(gramstore(j)%haml,gramstore(j)%zstore,elect,ndet,1)
            write(gst_num,"(i1.1)")j
            call matrixwriter(gramstore(j)%haml%hjk,ndet,"data/ham_final_"//trim(gst_num)//".csv")
            call matrixwriter(gramstore(j)%haml%ovrlp,ndet,"data/ovlp_final_"//trim(gst_num)//".csv")
        end do
        ! do j=2,gramnum+1
        !     call gram_ovrlp_fill(gramstore,j)
        ! end do

        call imaginary_time(gramstore(1)%dvecs,erg(1,:),gramstore(1)%haml, ndet)
        write(stdout,"(a,i1,a,f21.16)") "State, ",1," energy: ", erg(1,timesteps+1)
        call imaginary_time(gramstore(2)%dvecs,erg(2,:),gramstore(2)%haml, ndet)
        write(stdout,"(a,i1,a,f21.16)") "State, ",2," energy: ", erg(2,timesteps+1)
        call gram_ovrlp_fill(gramstore,2)
        gramstore(3)%dvecs%d(1:ndet)=gramstore(2)%dvecs%d(1:ndet)
        gramstore(3)%dvecs%d(1+ndet:2*ndet)=(-1)*gramstore(1)%dvecs%d(1:ndet)*&
        (dot_product(gramstore(1)%dvecs%d,matmul(gramstore(2)%wf_ovrlp(1,:,:),gramstore(2)%dvecs%d))/&
        (dot_product(gramstore(1)%dvecs%d,matmul(gramstore(1)%haml%ovrlp,gramstore(1)%dvecs%d))))
        erg(3,timesteps+1)=dot_product(gramstore(3)%dvecs%d,matmul(gramstore(3)%haml%hjk,gramstore(3)%dvecs%d))
        write(stdout,"(a,i1,a,f21.16)") "State, ",3," energy: ", erg(3,timesteps+1)
        ! call imaginary_time(gramstore(1)%dvecs,erg(1,:),gramstore(1)%haml, ndet)
        ! write(stdout,"(a,i1,a,f21.16)") "State, ",1," energy: ", erg(1,timesteps+1)
        ! gramstore(1)%d_ovrlp_d=dot_product(gramstore(1)%dvecs%d,matmul(gramstore(1)%haml%ovrlp,gramstore(1)%dvecs%d))
        ! do j=2,gramnum+1
        !     call imaginary_time(gramstore,erg(j,:),ndet,j)
        !     ! call imaginary_time(gramstore(j)%dvecs,erg(j,:),gramstore(j)%haml, ndet)
        !     write(stdout,"(a,i1,a,f21.16)") "State, ",j," energy: ", erg(j,timesteps+1)
        !     gramstore(j)%d_ovrlp_d=dot_product(gramstore(j)%dvecs%d,matmul(gramstore(j)%haml%ovrlp,gramstore(j)%dvecs%d))
        ! end do
            
        call energywriter(erg,"energy_final.csv",0)
        
        deallocate(erg,stat=ierr)
        if(ierr/=0)then
            write(stderr,"(a,i0)") "Error deallocating erg array ierr had value ",ierr
            errorflag=1
            return
        end if

        do j=1,gramnum+1
            call deallocgram(gramstore(j))
        end do
        ! deallocate(gramstore,stat=ierr)
        ! if(ierr/=0)then
        !     write(stderr,"(a,i0)") "Error deallocating gramstore array ierr had value ",ierr
        !     errorflag=1
        !     return
        ! end if
    end subroutine gram_schmidt_control

    function overlap_sum(gramstore,state,vector)
        implicit none
        type(gram),dimension(:),allocatable::gramstore
        integer::state,vector
        real(wp)::overlap_sum
        integer::j,k
        overlap_sum=0

        do j=1,state-1
            do k=1, ndet     
            overlap_sum=overlap_sum+abs(product(gramstore(state)%zstore(vector)%val(1:norb)*gramstore(j)%zstore(k)%val(1:norb)+&
            gramstore(state)%zstore(vector)%val(1+norb:2*norb)*gramstore(j)%zstore(k)%val(1+norb:2*norb)))
            end do 
        end do

    end function overlap_sum

    function overlap_sum2(gramstore,state,vector,orb)
        implicit none
        type(gram),dimension(:),allocatable::gramstore
        integer::state,vector,orb
        real(wp)::overlap_sum2
        integer::j,k
        overlap_sum2=0

        do j=1,state-1
            do k=1, ndet     
            overlap_sum2=overlap_sum2+abs((gramstore(state)%zstore(vector)%val(orb)*gramstore(j)%zstore(k)%val(orb)+&
            gramstore(state)%zstore(vector)%val(orb+norb)*gramstore(j)%zstore(k)%val(orb+norb)))
            end do 
        end do

    end function overlap_sum2
    
End Module gram_schmidt