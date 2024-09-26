module neural_network
    
    use mod_types
    use globvars
    use alarrays
   
    
    implicit none
    
    ! Parameters
    integer :: input_size = 5 !(norb*2)+5  ! Number of input features
    integer :: hidden_size = 1  ! Number of neurons in the hidden layer
    
contains

subroutine neural_network_control(size,elecs,haml,zstore,neural_net)
    implicit none
    integer,intent(in)::size
    type(elecintrgl),intent(in)::elecs
    type(trial_data),dimension(:),allocatable::trial_data_array
    type(neural_network_layer),allocatable,dimension(:),intent(inout)::neural_net
    type(hamiltonian),intent(in)::haml
    type(zombiest),dimension(:),intent(in)::zstore
    integer::num_epochs,j,layers,data_size
    integer::ierr=0
    real(wp)::learning_rate,trend

    trend=10
    data_size=0
    layers=2
    num_epochs=100
    do j=1,ndet
        data_size=data_size+j
    end do 
    data_size = data_size*norb*size
    input_size=5!+(norb*2)
    hidden_size=2 !norb*2
    learning_rate = 0.000001
    allocate(trial_data_array(data_size),stat=ierr)
    do j=1, data_size
        allocate(trial_data_array(j)%input_features(input_size))
    end do 

    allocate(neural_net(layers+1),stat=ierr)
    allocate(neural_net(1)%weights(input_size,hidden_size),stat=ierr)
    allocate(neural_net(1)%biases(hidden_size),stat=ierr)
    allocate(neural_net(1)%hidden_layer(hidden_size),stat=ierr)

    if(layers.gt.1)then
        do j=2,layers
            allocate(neural_net(j)%weights(hidden_size,hidden_size),stat=ierr)
            allocate(neural_net(j)%biases(hidden_size),stat=ierr)
            allocate(neural_net(j)%hidden_layer(hidden_size),stat=ierr)
        end do
    end if 

    allocate(neural_net(layers+1)%weights(1,hidden_size),stat=ierr)
    allocate(neural_net(layers+1)%biases(1),stat=ierr)
  
    if(ierr/=0)then
        print*,'Error allocating trial_data_array'
        errorflag=1
        return
    end if
    call trial_data_creation(trial_data_array,elecs,haml,zstore,data_size)
    print*,'setup data complete'
  
    call randomise_data_set(data_size,trial_data_array)
    do j=1, layers+1
        call layer_intialise(neural_net(j))
    end do 
   
    do j=1,100
        call train_neural_network(neural_net,trial_data_array,learning_rate,num_epochs,data_size,trend)
        call predict(neural_net,elecs)
        call randomise_data_set(data_size,trial_data_array)
        ! learning_rate =  learning_rate*0.1
    end do 
  
    deallocate(trial_data_array,stat=ierr)
    if(ierr/=0)then
        print*,'Error deallocating trial_data_array'
        errorflag=1
        return
    end if
    
end subroutine neural_network_control


! Train the neural network
subroutine train_neural_network(neural_net, trial_data_array, learning_rate, num_epochs, data_size,trend)
    implicit none
    type(neural_network_layer),dimension(:),intent(inout)::neural_net
    type(trial_data),intent(in)::trial_data_array(:)
    real(wp),dimension(data_size)::errors
    real(wp),intent(inout)::trend,learning_rate
    integer,intent(in)::num_epochs,data_size
    integer::j,k

    do k=1,num_epochs
        do j=1,data_size
            call forward_pass(neural_net,trial_data_array(j))
            call backpropagation(trial_data_array(j),neural_net,learning_rate,errors(j))
        end do
        print*,sum(errors)/data_size,learning_rate
    end do 
    ! if(ABS((sum(errors)/data_size)-trend).lt.learning_rate)then
    !    learning_rate=learning_rate*0.1
    ! end if
    ! trend=sum(errors)/data_size
    return 
end subroutine train_neural_network

! Forward pass through the neural network
subroutine forward_pass(neural_net,input)
    implicit none 
    type(neural_network_layer),dimension(:),intent(inout)::neural_net
    type(trial_data),intent(in)::input
    real(wp):: weighted_sum_hidden(hidden_size),weighted_sum_output
    integer::j,k

    do concurrent(j=1:hidden_size)
        weighted_sum_hidden(j)=sum(input%input_features*neural_net(1)%weights(:,j))+neural_net(1)%biases(j)
        neural_net(1)%hidden_layer(j)=sigmoid(weighted_sum_hidden(j))
    end do

    if(size(neural_net).gt.2)then
        do k=2,size(neural_net)-1
            do concurrent(j=1:hidden_size)
                weighted_sum_hidden(j)=sum(neural_net(k-1)%hidden_layer*neural_net(k)%weights(:,j))+neural_net(k)%biases(j)
                neural_net(k)%hidden_layer(j)=sigmoid(weighted_sum_hidden(j))
            end do
        end do
    end if
    weighted_sum_output = sum(neural_net(size(neural_net)-1)%hidden_layer*&
                                neural_net(size(neural_net))%weights(1,:))+&
                                neural_net(size(neural_net))%biases(1)
    neural_net(size(neural_net))%output = sigmoid(weighted_sum_output)*input%old_ham

    return 
end subroutine forward_pass

function forward_result(neural_net,input) result(output)
    implicit none 
    type(neural_network_layer),dimension(:),intent(inout)::neural_net
    real(wp),dimension(:),intent(in)::input
    real(wp):: output
    real(wp):: weighted_sum_hidden(hidden_size),weighted_sum_output
    integer::j,k

    do concurrent(j=1:hidden_size)
        weighted_sum_hidden(j)=sum(input*neural_net(1)%weights(:,j))+neural_net(1)%biases(j)
        neural_net(1)%hidden_layer(j)=sigmoid(weighted_sum_hidden(j))
    end do
    if(size(neural_net).gt.2)then
        do k=2,size(neural_net)-1
            do concurrent(j=1:hidden_size)
                weighted_sum_hidden(j)=sum(neural_net(k-1)%hidden_layer*neural_net(k)%weights(:,j))+neural_net(k)%biases(j)
                neural_net(k)%hidden_layer(j)=sigmoid(weighted_sum_hidden(j))
            end do
        end do
    end if
    weighted_sum_output = sum(neural_net(size(neural_net)-1)%hidden_layer*&
                neural_net(size(neural_net))%weights(1,:))+&
                neural_net(size(neural_net))%biases(1)
    output = sigmoid(weighted_sum_output)*input(3)
   
    return 
end function forward_result

! Backpropagation to update weights and biases
subroutine backpropagation(input,neural_net,learning_rate,mse)
    implicit none
    type(neural_network_layer),dimension(:),intent(inout)::neural_net
    type(trial_data),intent(in)::input
    real(wp),intent(in)::learning_rate
    real(wp),intent(inout)::mse
    real(wp)::output
    real(wp) :: error_output,delta_output
    real(wp),dimension(:,:),allocatable::errors,deltas
    integer :: j,length,k,val
    
    length=size(neural_net)
    output=neural_net(length)%output
  
    ! Compute error for output layer
    error_output = output * (1.0 - output) * (input%new_ham - output)
    mse=abs((input%new_ham - output))!*(input%new_ham - output)) !/abs(input%new_ham)
  
    ! Update weights and biases for output layer
    delta_output = learning_rate * error_output
   
    neural_net(length)%weights(1,:) = neural_net(length)%weights(1,:) + (delta_output * neural_net(length-1)%hidden_layer)
    neural_net(length)%biases(1) =  neural_net(length)%biases(1) + delta_output
    
    allocate(errors(length-1,hidden_size),deltas(length-1,hidden_size))
     ! Compute error for hidden layer
    
    if(length.gt.2)then
        do concurrent (j = 1:hidden_size)
            errors(1,j) =  neural_net(length-1)%hidden_layer(j) * (1.0 - neural_net(length-1)%hidden_layer(j)) * &
            (neural_net(length)%weights(1,j) * error_output)
        end do

        do concurrent(j=1:hidden_size)
            deltas(1,j) = learning_rate * errors(1,j)
            neural_net(length-1)%weights(:, j) = neural_net(length-1)%weights(:, j) + & 
                deltas(1,j) * neural_net(length-2)%hidden_layer

            neural_net(length-1)%biases(j) = neural_net(length-1)%biases(j) + deltas(1,j)
        end do
    end if 
    if(length.gt.4)then 
        do k=2,(length-2)
            val=length-k
            do concurrent(j = 1:hidden_size)
                errors(k,j) = neural_net(val)%hidden_layer(j)  * (1.0 - neural_net(val)%hidden_layer(j)) * &
                sum(neural_net(val+1)%weights(j,:) * errors(k-1,:))
            end do
            do concurrent(j=1:hidden_size)
                deltas(k,j) = learning_rate * errors(k,j)
                neural_net(val)%weights(:, j) = neural_net(val)%weights(:, j) + & 
                    deltas(k,j) * neural_net(val-1)%hidden_layer
        
                neural_net(val)%biases(j) = neural_net(val)%biases(j) + deltas(k,j)
            end do
        end do
    end if 
  
    do concurrent(j = 1:hidden_size)
        errors(length-1,j) = neural_net(1)%hidden_layer(j)  * (1.0 - neural_net(1)%hidden_layer(j)) * &
        sum(neural_net(2)%weights(j,:) * errors(length-2,:))
    end do 
   
    do concurrent(j=1:hidden_size)
        deltas(length-1,j) = learning_rate * errors(length-1,j)
        neural_net(1)%weights(:, j) = neural_net(1)%weights(:, j) + & 
            deltas(length-1,j) * input%input_features

        neural_net(1)%biases(j) = neural_net(1)%biases(j) + deltas(length-1,j)
    end do

    deallocate(errors,deltas)
   
    return 
end subroutine backpropagation

subroutine predict(neural_net,elecs)
    type(neural_network_layer),dimension(:),intent(inout)::neural_net
    type(elecintrgl),intent(in)::elecs
    type(zombiest)::z1,z2
    real(wp)::output,true_output,ovrlp  
    real(wp),dimension(input_size)::input_features
    integer::k
    call alloczf(z1)
    call alloczf(z2)
  
    do k=1,norb
        call trial_data_zom(z1)
        call trial_data_zom(z2)
        ovrlp=product((z1%val(1:norb)*z2%val(1:norb))+((z1%val(1+norb:2*norb)*z2%val(1+norb:2*norb))))
        input_features(3)=haml_val_trial_data(z1%val,z2%val,elecs,ovrlp,0)
        input_features(4)=k
        input_features(5)=0
        call trial_data_zom_alt(z2,k)
        input_features(1)=z2%val(k+norb)
        input_features(2)=z2%val(k)
        ! input_features(6:(6+norb))=z1%val(1:norb)
        ! input_features((7+norb):)=z2%val(1:norb)
        output=forward_result(neural_net,input_features)
        ovrlp=product((z1%val(1:norb)*z2%val(1:norb))+((z1%val(1+norb:2*norb)*z2%val(1+norb:2*norb))))
        true_output=haml_val_trial_data(z1%val,z2%val,elecs,ovrlp,0)
        
        print*,output,true_output,(100*abs((output-true_output)/true_output))
    end do

    do k=1,norb
        call trial_data_zom(z1)
        z2=z1
        ovrlp=1.0d0
        input_features(3)=haml_val_trial_data(z1%val,z2%val,elecs,ovrlp,1)
        input_features(4)=k
        input_features(5)=1
        call trial_data_zom_alt(z1,k)
        input_features(1)=z1%val(k+norb)
        input_features(2)=z1%val(k)
        z2=z1
        ! input_features(6:(6+norb))=z1%val(1:norb)
        ! input_features((7+norb):)=z2%val(1:norb)
        output=forward_result(neural_net,input_features)
        true_output=haml_val_trial_data(z1%val,z2%val,elecs,ovrlp,1)
        
        print*,output,true_output,(100*abs((output-true_output)/true_output))
    end do

end subroutine predict


subroutine trial_data_creation(trial_data_array,elecs,haml,zstore,data_size)
    implicit none
    type(trial_data),intent(inout) :: trial_data_array(:)
    type(elecintrgl),intent(in)::elecs
    type(hamiltonian),intent(in)::haml
    type(zombiest),dimension(:)::zstore
    integer,intent(in) :: data_size
    integer :: j,l,k,ind
    real(wp)::ovrlp
    type(zombiest)::z1,z2

    call alloczf(z1)
    call alloczf(z2)

    l=1;j=2;k=1 
    do ind=1, data_size
        trial_data_array(ind)%orbital=k
        trial_data_array(ind)%z1=l
        trial_data_array(ind)%z2=j
        trial_data_array(ind)%old_ham=(haml%hjk(l,j)-(haml%ovrlp(l,j)*elecs%hnuc))
        trial_data_array(ind)%new_ham=haml%ovrlp(l,j)
        if(l==j)then 
            trial_data_array(ind)%diag=1
        else
            trial_data_array(ind)%diag=0
        end if 
        k=k+1
        if(k.gt.norb)then
            k=1
            j=j+1
            if(j.gt.ndet)then
                l=l+1
                if(l.gt.ndet)then
                    l=1
                    j=2
                end if
                j=l
            end if
        end if
    end do
    !$omp parallel do private(z1,z2,ovrlp) shared(trial_data_array,zstore,elecs)
    do ind=1, data_size
        z1=zstore(trial_data_array(ind)%z1)
        z2=zstore(trial_data_array(ind)%z2)
      
        call trial_data_zom_alt(z2, trial_data_array(ind)%orbital)
        if(trial_data_array(ind)%diag==1)then 
            z1=z2
        end if
        ovrlp=trial_data_array(ind)%new_ham
        trial_data_array(ind)%new_ham=(haml_val_trial_data(z1%val,z2%val,elecs,ovrlp,trial_data_array(ind)%diag))

        trial_data_array(ind)%input_features(1)=trial_data_array(ind)%dead
        trial_data_array(ind)%input_features(2)=trial_data_array(ind)%alive
        trial_data_array(ind)%input_features(3)=trial_data_array(ind)%old_ham
        trial_data_array(ind)%input_features(4)=trial_data_array(ind)%orbital
        trial_data_array(ind)%input_features(5)=trial_data_array(ind)%diag
        ! trial_data_array(ind)%input_features(6:(6+norb))=z1%val(1:norb)
        ! trial_data_array(ind)%input_features((7+norb):)=z2%val(1:norb)
    end do 
    !$omp end parallel do 
    call dealloczf(z1)
    call dealloczf(z2)
  
    return 

end subroutine trial_data_creation

function haml_val_trial_data(z1d,z2d,elecs,ovrlp,sm) result(ham_tot)
    implicit none 
    real(wp),dimension(0:),intent(in)::z1d,z2d
    real(wp),intent(in)::ovrlp
    real(wp)::ham_tot
    type(elecintrgl),intent(in)::elecs
    integer,intent(in)::sm
    real(wp),dimension(4,norb)::perts
    real(wp),dimension(norb)::div,bth
    real(wp)::ov,aa,dd,ad,da
    integer::j,k,l
    real(wp),dimension(elecs%num)::ovrlp_vec
    
    ov=1.0d0
    if(sm==1)then
        do j=1,norb
            dd=z1d(j+norb)*z2d(norb+j)
            bth(j)=z1d(j)*z2d(j)
            perts(1,j)=z1d(j+norb)*z2d(j)
            perts(2,j)=z1d(j+norb)*z2d(j)
            perts(3,j)=z1d(j)*z2d(j)
            perts(4,j)=(-bth(j)+dd)
        end do
    else 
        do j=1,norb
            aa=z1d(j)*z2d(j)
            dd=z1d(j+norb)*z2d(norb+j)
            ad=z1d(j)*z2d(norb+j)
            da=z1d(j+norb)*z2d(j)
            div(j)=aa+dd
            perts(1,j)=da/div(j)
            perts(2,j)=ad/div(j) 
            perts(3,j)=aa/div(j)
            perts(4,j)=(-aa+dd)/div(j)
        end do
    end if 
    
    ovrlp_vec=ovrlp
    
    do k=1,norb
        do l=1,elecs%orbital_choice2(0,k)
            ov=ovrlp_vec(elecs%orbital_choice2(k,(l*2)-1))*&
            perts(elecs%orbital_choice(k,elecs%orbital_choice2(k,(l*2)-1)),elecs%orbital_choice3(k))
            do j=elecs%orbital_choice2(k,(l*2)-1),elecs%orbital_choice2(k,l*2)
                ovrlp_vec(j)=ov
            end do
        end do
    end do
    
    do j=1,elecs%num
        ham_tot=ham_tot+(ovrlp_vec(j)*elecs%integrals(j))
    end do
    
    return 
    
end function haml_val_trial_data


subroutine trial_data_zom(z1)
    implicit none
    type(zombiest)::z1
    integer::k
    DOUBLE PRECISION, external::ZBQLU01

    z1%phi=0.0001
    if(nel.gt.4)then
        z1%phi(1:4)=0.5*pirl
        if(modulo(nel,2)==0) then
            do k=5,nel+4
                !$omp critical
                z1%phi(k)=0.5*pirl*ZBQLU01(1)
                !$omp end critical
            end do
        else
            do k=5,nel+5
                !$omp critical
                z1%phi(k)=0.5*pirl*ZBQLU01(1)
                !$omp end critical
            end do
        end if
    else if(nel.eq.4)then
        z1%phi(1:2)=0.5*pirl
        if(modulo(nel,2)==0) then
            do k=3,nel+4
                !$omp critical
                z1%phi(k)=0.5*pirl*ZBQLU01(1)
                !$omp end critical
            end do
        else
            do k=3,nel+5
                    !$omp critical
                z1%phi(k)=0.5*pirl*ZBQLU01(1)
                !$omp end critical
            end do
        end if
    else 
        if(modulo(nel,2)==0) then
            do k=1,nel+4
                !$omp critical
                z1%phi(k)=0.5*pirl*ZBQLU01(1)
                !$omp end critical
            end do
        else
            do k=1,nel+5
                    !$omp critical
                z1%phi(k)=0.5*pirl*ZBQLU01(1)
                !$omp end critical
            end do
        end if
    end if
    z1%val(1:norb)=sin(z1%phi)
    z1%val(1+norb:2*norb)=cos(z1%phi)
 
    return 

end subroutine trial_data_zom

subroutine trial_data_zom_alt(z1,orb)
    implicit none
    type(zombiest)::z1
    integer,intent(in)::orb
    real(wp)::pert
    integer::posneg 
    DOUBLE PRECISION, external::ZBQLU01

     !$omp critical
    if(ZBQLU01(1).gt.0.3)then
        posneg=-1
    else
        posneg=1
    end if
    pert=0.1*ZBQLU01(1)
    !$omp end critical
    z1%phi(orb)=z1%phi(orb)+pert*posneg

    z1%val(orb)=sin(z1%phi(orb))
    z1%val(orb+norb)=cos(z1%phi(orb))
    return 
end subroutine trial_data_zom_alt


subroutine layer_intialise(layer)
    implicit none
    type(neural_network_layer),intent(inout)::layer
    DOUBLE PRECISION, external::ZBQLU01
    integer::j,k,l1,l2
    l1=size(layer%weights(:,1))
    l2=size(layer%weights(1,:))

    do k=1,l1
        do j=1,l2
            layer%weights(k,j)=sqrt(-2.0 * log(ZBQLU01(1))) * cos(2.0 * acos(-1.0) * ZBQLU01(1))
            layer%weights(k,j)=layer%weights(k,j)*0.001
        end do
    end do
   
    layer%biases = 0.0d0

    return 
end subroutine layer_intialise


subroutine randomise_data_set(N,input)
    implicit none
    type(trial_data),dimension(:),intent(inout)::input
    type(trial_data),dimension(N)::temp
    integer::N
    integer::j,k,l,jtemp
    integer,dimension(N)::ind
    DOUBLE PRECISION, external::ZBQLU01

    ind = (/ (j, j=1, N) /)
   
    do k=1,2
        do j=1,N
            l = 1 + FLOOR((N-1)*ZBQLU01(1))
            jtemp=ind(l)
            ind(l)=ind(j)
            ind(j)=jtemp
        end do
    end do
    
    do j=1,N
        temp(j)=input(ind(j))
    end do
    input = temp
   
    return 

end subroutine randomise_data_set

! Activation function: Sigmoid
real(wp) elemental function sigmoid(x)
    implicit none
    real(wp), intent(in) :: x
    sigmoid = 1.0d0 / (1.0d0 + exp(-x))
end function sigmoid

end module neural_network


