MODULE gram_schmidt

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
      
        real(kind=8), dimension(:,:),allocatable::erg
        integer::ierr,j,k
        character(LEN=20)::filenm
        logical :: file_exists
        character(len=2)::gst_num


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

        do j=1,gramnum+1
            call allocgram(gramstore(j),j,ndet,norb)
        end do
            
     

    end subroutine gram_schmidt_control

    ! subroutine gram_ovrlp_fill(zstore)
        
    !     implicit none

    !     type(zombiest),intent(in),dimension(:,:)::zstore
    !     type(hamiltonian),intent(inout),dimension(:)::haml
    !     integer::j,k,l,plc,state

    !     plc=1

    !     do state=2,gramnum+1
    !         do j=1, (state-1)
    !             haml(state)%gs_ovrlp(j,:,:)=haml(j)%ovrlp
    !             haml(state)%gs_kinvh(j,:,:)=haml(j)%kinvh
    !             do k=1,ndet
    !                 do l=k,ndet
    !                     haml(state)%gs_ovrlp_self(j,k,l)=overlap_1(zstore(state,k)%val,zstore(j,l)%val)
    !                     haml(state)%gs_ovrlp_self(j,l,k)=haml(state)%gs_ovrlp_self(j,k,l)
    !                 end do 
    !             end do
    !         end do
    !     end do 

    !     return

    ! end subroutine gram_ovrlp_fill


    ! !Gram-Schmidt orthogonalisation 
    ! subroutine gs(dvecs,haml)

    !     implicit none
    !     type(dvector), intent(inout)::dvecs
    !     type(hamiltonian),intent(in)::haml
    !     real(kind=8)::numer,den
    !     real(kind=8),dimension(:),allocatable::dvecs_copy
    !     integer::states,j,k,ierr

    !     if (errorflag .ne. 0) return
    !     ierr=0
    !     allocate(dvecs_copy(ndet),stat=ierr)
    !     if(ierr/=0)then
    !         write(stderr,"(a,i0)") "Error allocating dvecs_copy array ierr had value ",ierr
    !         errorflag=1
    !         return
    !     end if
       
       
    !     do j=2,gramnum+1
    !         dvecs_copy=dvecs(j)%d
    !         do k=1, j-1
    !             numer=dot_product(dvecs(j)%d,matmul(gram_ovrlp(k)%ovrlps_numer((j-k),:,:),dvecs(k)%d))
    !             den=dot_product(dvecs(k)%d,matmul(haml(k)%ovrlp,dvecs(k)%d))
    !             dvecs_copy = dvecs_copy - (dvecs(k)%d*(numer/den))
    !         end do
    !         dvecs(j)%d=dvecs_copy
    !     end do

    !     return

    ! end subroutine gs

    







End Module gram_schmidt