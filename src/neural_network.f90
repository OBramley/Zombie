module neural_network
    use mod_types
    use globvars
    use alarrays
    use ham
    
    implicit none
    
    ! Parameters
    integer, parameter :: input_size = 6   ! Number of input features
    integer, parameter :: hidden_size = 1  ! Number of neurons in the hidden layer
    
   
    
    type trial_data
        real(wp)::dead
        real(wp)::alive
        real(wp)::new_ham
        real(wp)::old_ham
        integer::orbital
        integer::diag
        integer::rhf
        real(wp),dimension(6)::input_features
    end type trial_data

    type neural_network_layer
        real(wp)::weights_input_hidden(input_size, hidden_size)
        real(wp)::biases_hidden(hidden_size)
        real(wp)::weights_hidden_output(hidden_size)
        real(wp)::bias_output
    end type neural_network_layer


contains

subroutine neural_network_control(size,elecs)
    implicit none
    integer,intent(in)::size
    type(elecintrgl),intent(in)::elecs
    type(trial_data),dimension(:),allocatable::trial_data_array
    type(neural_network_layer)::layer
    integer::num_epochs
    integer::ierr=0
    real(wp)::learning_rate
    
    num_epochs = (norb*size*size)+(norb*size)
   
    learning_rate = 0.0001
    allocate(trial_data_array(num_epochs),stat=ierr)
    if(ierr/=0)then
        print*,'Error allocating trial_data_array'
        errorflag=1
        return
    end if
    call trial_data_creation(trial_data_array,elecs,size,num_epochs)
    print*,'setup data complete'
    call randomise_data_set(num_epochs,trial_data_array)
    call layer_intialise(layer)
    call train_neural_network(layer,trial_data_array,learning_rate,num_epochs)
    call predict(layer,elecs)
    deallocate(trial_data_array,stat=ierr)
    if(ierr/=0)then
        print*,'Error deallocating trial_data_array'
        errorflag=1
        return
    end if
    stop
end subroutine neural_network_control


! Train the neural network
subroutine train_neural_network(layer, trial_data_array, learning_rate, num_epochs)
    implicit none
    type(neural_network_layer),intent(inout)::layer
    type(trial_data),intent(in)::trial_data_array(:)
    real(wp),intent(in)::learning_rate
    integer,intent(in)::num_epochs
    real(wp)::output
    real(wp),dimension(hidden_size)::hidden_layer
    integer::j

    do j=1,num_epochs
        call forward_pass(layer,trial_data_array(j),output,hidden_layer)
        call backpropagation(trial_data_array(j),layer,hidden_layer,output,learning_rate)
    end do

    return 
end subroutine train_neural_network

! Forward pass through the neural network
subroutine forward_pass(layer,input,output,hidden_layer)
    implicit none 
    type(neural_network_layer),intent(in)::layer
    type(trial_data),intent(in)::input
    real(wp),intent(out)::output
    real(wp),intent(out):: hidden_layer(hidden_size)
    real(wp):: weighted_sum_hidden,weighted_sum_output
   
    integer::j

    do j=1,hidden_size
        weighted_sum_hidden=sum(input%input_features*layer%weights_input_hidden(:,j))+layer%biases_hidden(j)
        hidden_layer(j)=sigmoid(weighted_sum_hidden)
    end do

    weighted_sum_output = sum(hidden_layer*layer%weights_hidden_output)+layer%bias_output
    output = sigmoid(weighted_sum_output)*abs(input%old_ham)

    return 
end subroutine forward_pass

! Backpropagation to update weights and biases
subroutine backpropagation(input,layer, hidden_layer, output, learning_rate)
    implicit none
    type(neural_network_layer),intent(inout)::layer
    type(trial_data),intent(in)::input
    real(wp),intent(in)::output,learning_rate
    real(wp),intent(in):: hidden_layer(hidden_size)
    real(wp) :: error_output, error_hidden(hidden_size)
    real(wp) :: delta_output, delta_hidden(hidden_size)
    integer :: j
    
    ! Compute error for output layer
    error_output = output * (1.0 - output) * (input%new_ham - output)
    
    ! Update weights and biases for output layer
    delta_output = learning_rate * error_output
    layer%weights_hidden_output = layer%weights_hidden_output + (delta_output * hidden_layer)
    layer%bias_output = layer%bias_output + delta_output

    ! Compute error for hidden layer
    do concurrent (j = 1:hidden_size)
        error_hidden(j) = hidden_layer(j) * (1.0 - hidden_layer(j)) * (layer%weights_hidden_output(j) * error_output)
    end do
    
    ! Update weights and biases for hidden layer
    do concurrent (j = 1:hidden_size)
        delta_hidden(j) = learning_rate * error_hidden(j)
        layer%weights_input_hidden(:, j) = layer%weights_input_hidden(:, j) + delta_hidden(j) * input%input_features
        layer%biases_hidden(j) = layer%biases_hidden(j) + delta_hidden(j)
    end do
    return 
end subroutine backpropagation

subroutine predict(layer,elecs)
    type(neural_network_layer),intent(in)::layer
    type(elecintrgl),intent(in)::elecs
    type(zombiest)::z1,z2
    real(wp)::output,true_output
    real(wp):: hidden_layer(hidden_size)
    real(wp):: weighted_sum_hidden,weighted_sum_output
    real(wp),dimension(6)::input_features
    integer::j,k
    call alloczf(z1)
    call alloczf(z2)
    print*,layer%weights_input_hidden
    print*,layer%biases_hidden
    print*,layer%weights_hidden_output
    print*,layer%bias_output
    do k=1,norb
        call trial_data_zom(z1)
        call trial_data_zom(z2)
        input_features(3)=haml_val_trial_data(z1%val,z2%val,elecs,1)
        input_features(4)=k
        input_features(5)=0
        input_features(6)=0
        call trial_data_zom_alt(z1,k)
        input_features(1)=z1%val(k+norb)
        input_features(2)=z1%val(k)
        true_output=haml_val_trial_data(z1%val,z2%val,elecs,1)
        do j=1,hidden_size
            weighted_sum_hidden=sum(input_features*layer%weights_input_hidden(:,j))+layer%biases_hidden(j)
            hidden_layer(j)=sigmoid(weighted_sum_hidden)
        end do
        weighted_sum_output = sum(hidden_layer*layer%weights_hidden_output)+layer%bias_output
        output = sigmoid(weighted_sum_output)*input_features(3)
        print*,output,true_output
    end do

    do k=1,norb
        call trial_data_zom(z1)
        z2=z1
        input_features(3)=haml_val_trial_data(z1%val,z2%val,elecs,0)
        input_features(4)=k
        input_features(5)=1
        input_features(6)=0
        call trial_data_zom_alt(z1,k)
        input_features(1)=z1%val(k+norb)
        input_features(2)=z1%val(k)
        z2=z1
        true_output=haml_val_trial_data(z1%val,z2%val,elecs,0)
        do j=1,hidden_size
            weighted_sum_hidden=sum(input_features*layer%weights_input_hidden(:,j))+layer%biases_hidden(j)
            hidden_layer(j)=sigmoid(weighted_sum_hidden)
        end do
        weighted_sum_output = sum(hidden_layer*layer%weights_hidden_output)+layer%bias_output
        output = sigmoid(weighted_sum_output)*input_features(3)
        print*,output,true_output
    end do

end subroutine predict


subroutine trial_data_creation(trial_data_array,elecs,size,num_epochs)
    implicit none
    type(trial_data),intent(inout) :: trial_data_array(:)
    type(elecintrgl),intent(in)::elecs
    integer,intent(in) :: size,num_epochs
    integer :: j,l,k,ind
    type(zombiest)::z1,z2

    call alloczf(z1)
    call alloczf(z2)

    !$omp parallel do private(z1,z2,ind,j,l,k) shared(trial_data_array,elecs,size,norb)
    do l=1,size 
        do j=1,norb 
            do k=1, size
                ind = (l - 1) * size * norb + (j - 1) * size + k
                print*,ind
                if(l==k)then 
                    trial_data_array(ind)%orbital=j
                    trial_data_array(ind)%diag=1
                    trial_data_array(ind)%rhf=0

                    call trial_data_zom(z1)
                    z2=z1
                    trial_data_array(ind)%old_ham=haml_val_trial_data(z1%val,z2%val,elecs,0)
                    call trial_data_zom_alt(z1,j)
        
                    trial_data_array(ind)%dead=z1%val(j+norb)
                    trial_data_array(ind)%alive=z1%val(j)

                    z2=z1
                    
                    trial_data_array(ind)%new_ham=haml_val_trial_data(z1%val,z2%val,elecs,0)
                else 
                    trial_data_array(ind)%orbital=j
                    trial_data_array(ind)%diag=0
                    trial_data_array(ind)%rhf=0
                    call trial_data_zom(z1)
                    call trial_data_zom(z2)

                    trial_data_array(ind)%old_ham=haml_val_trial_data(z1%val,z2%val,elecs,1)
                    call trial_data_zom_alt(z1,j)
        
                    trial_data_array(ind)%dead=z1%val(j+norb)
                    trial_data_array(ind)%alive=z1%val(j)

                    trial_data_array(ind)%new_ham=haml_val_trial_data(z1%val,z2%val,elecs,1)
                end if 
                trial_data_array(ind)%input_features(1)=trial_data_array(ind)%dead
                trial_data_array(ind)%input_features(2)=trial_data_array(ind)%alive
                trial_data_array(ind)%input_features(3)=trial_data_array(ind)%old_ham
                trial_data_array(ind)%input_features(4)=trial_data_array(ind)%orbital
                trial_data_array(ind)%input_features(5)=trial_data_array(ind)%diag
                trial_data_array(ind)%input_features(6)=trial_data_array(ind)%rhf
            end do 
        end do 
    end do 
    !$omp end parallel do

    z1%phi(1:nel)=0.5*pirl
    z1%phi(nel+1:)=0
    z1%val(1:norb)=0 
    z1%val(norb+1:)=1 
    z1%val(1:nel)=1
    z1%val(norb+1:norb+nel)=0
    ind=1 
    do j=(norb*size*size+1), (num_epochs)
        trial_data_array(j)%orbital=ind
        trial_data_array(j)%diag=0
        trial_data_array(j)%rhf=1
        call trial_data_zom(z2)

        trial_data_array(j)%old_ham=haml_val_trial_data(z1%val,z2%val,elecs,1)
       
        call trial_data_zom_alt(z1,ind)

        trial_data_array(j)%dead=z2%val(ind+norb)
        trial_data_array(j)%alive=z2%val(ind)

        trial_data_array(j)%new_ham=haml_val_trial_data(z1%val,z2%val,elecs,1)

        trial_data_array(j)%input_features(1)=trial_data_array(j)%dead
        trial_data_array(j)%input_features(2)=trial_data_array(j)%alive
        trial_data_array(j)%input_features(3)=trial_data_array(j)%old_ham
        trial_data_array(j)%input_features(4)=trial_data_array(j)%orbital
        trial_data_array(j)%input_features(5)=trial_data_array(j)%diag
        trial_data_array(j)%input_features(6)=trial_data_array(j)%rhf
        ind=ind+1
        if(ind.gt.norb)then
            ind=1
        end if
    end do 
  
    call dealloczf(z1)
    call dealloczf(z2)
  
    return 

end subroutine trial_data_creation

function haml_val_trial_data(z1d,z2d,elecs,sm) result(ham_tot)
    implicit none 
    real(wp),dimension(0:),intent(in)::z1d,z2d
    real(wp)::ham_tot
    type(elecintrgl),intent(in)::elecs
    integer,intent(in)::sm
    real(wp),dimension(4,norb)::perts
    real(wp),dimension(norb)::div,bth
    real(wp)::ov,aa,dd,ad,da
    integer::j,k,l
    real(wp),dimension(elecs%num)::ovrlp_vec
    
    ov=1.0d0
    if(sm==0)then
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
    
    ovrlp_vec=1.0d0
    
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
    call val_set(z1)
    return 

end subroutine trial_data_zom

subroutine trial_data_zom_alt(z1,orb)
    implicit none
    type(zombiest)::z1
    integer,intent(in)::orb
    real(wp)::pert
    integer::posneg 
    DOUBLE PRECISION, external::ZBQLU01

    if(ZBQLU01(1).gt.0.4)then
        posneg=-1
    else
        posneg=1
    end if
    pert=(0.01-1.0d-20)*ZBQLU01(1)+1.0d-20
    z1%phi(orb)=z1%phi(orb)+pert*posneg

    z1%val(orb)=sin(z1%phi(orb))
    z1%val(orb+norb)=cos(z1%phi(orb))
    return 
end subroutine trial_data_zom_alt

subroutine layer_intialise(layer)
    implicit none
    type(neural_network_layer),intent(inout)::layer
    DOUBLE PRECISION, external::ZBQLNOR
    integer::j,k 

    do k=1,hidden_size
        do j=1,input_size
            layer%weights_input_hidden(j,k)=ZBQLNOR (0,0.01)
        end do
        layer%weights_hidden_output(k)=ZBQLNOR (0,0.01)
    end do

    layer%biases_hidden = 0
    layer%bias_output = 0

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
real(wp) function sigmoid(x)
    implicit none
    real(wp), intent(in) :: x
    sigmoid = 1.0d0 / (1.0d0 + exp(-x))
end function sigmoid

end module neural_network


