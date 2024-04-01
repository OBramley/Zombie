Module gram_funcs

    use mod_types
    use globvars

contains 

subroutine gram_ovrlp_fill(gramstore,state)
        
    implicit none
    type(gram),dimension(:)::gramstore
    integer::state
    integer::j,k,l

    if(errorflag.ne.0) return

    do j=1,state-1
        do k=1,ndet
            do l=k,ndet
                gramstore(state)%wf_ovrlp(j,k,l)=product(gramstore(state)%zstore(k)%val(1:norb)*&
                                                        gramstore(j)%zstore(l)%val(1:norb)+&
                                                        gramstore(state)%zstore(k)%val(1+norb:2*norb)*&
                                                        gramstore(j)%zstore(l)%val(1+norb:2*norb))
                gramstore(state)%wf_ovrlp(j,l,k)=gramstore(state)%wf_ovrlp(j,k,l)
            end do 
            
        end do
        
    end do
   
    return

end subroutine gram_ovrlp_fill

subroutine gram_ovrlp_var(gramstore,state1,state2,var)
    
    implicit none
    type(gram),dimension(:)::gramstore
    integer::state1,state2,var
    integer::k

    if(errorflag.ne.0) return

    do k=1,ndet
        gramstore(state1)%wf_ovrlp(state2,k,var)=product(gramstore(state1)%zstore(k)%val(1:norb)*&
                                                gramstore(state2)%zstore(var)%val(1:norb)+&
                                                gramstore(state1)%zstore(k)%val(1+norb:2*norb)*&
                                                gramstore(state2)%zstore(var)%val(1+norb:2*norb))
        gramstore(state1)%wf_ovrlp(state2,var,k)=gramstore(state1)%wf_ovrlp(state2,k,var)
    end do
        
    return

end subroutine gram_ovrlp_var

subroutine gs_wavefunction(gramstore,state)
    implicit none
    type(gram),dimension(:)::gramstore
    integer::state
    real(wp)::numer,denom
    real(wp),dimension(state-1)::gs
    integer::j,k
    if(errorflag.ne.0) return

    do j=1,state-1
        numer=dot_product(gramstore(state)%dvecs%d,matmul(gramstore(state)%wf_ovrlp(j,:,:),gramstore(j)%dvecs%d))
        denom=dot_product(gramstore(j)%dvecs%d,matmul(gramstore(j)%haml%ovrlp,gramstore(j)%dvecs%d))
        gs(j)=numer/denom
    end do

    do j=1,state-1
        gramstore(state)%dvecs%d=gramstore(state)%dvecs%d-(gramstore(j)%dvecs%d*gs(j))
        ! do k=1,ndet
        !     gramstore(state)%zstore(k)%phi(:)=gramstore(state)%zstore(k)%phi(:)-(gramstore(j)%zstore(k)%phi(:)*gs(j))
        ! end do 
    end do

    ! do k=1,ndet
    !     call val_set(gramstore(state)%zstore(k))
    ! end do 

    return 
end subroutine gs_wavefunction



End Module gram_funcs