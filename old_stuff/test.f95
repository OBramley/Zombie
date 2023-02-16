Program test
    
    implicit none

    real(kind=8),dimension(:), allocatable::A1,A2,B1,B2,C1,C2,C3,A1_1,A2_1,B1_1,B2_1,c1_1,C2_1!,C3
    real(kind=8),dimension(0:20)::t1,t2
    real(kind=8),dimension(0:10,2)::d1,d2
    integer,dimension(20)::ind
    integer(kind=2),dimension(10)::ind1,ind2,neg
    equivalence (d1,t1)
    equivalence (d2,t2)
    integer::j,k,l,len
    real(kind=8)::val1,val2,val3,start,finish
    integer,dimension(4)::ops

    len=10
    neg=1
    d1(0,:)=0
    ! allocate(d1(len,2))
    ! do j=1,2
    !     do k=1,len 
    !         ! d1(k,j)=k+(j-1)*10
    !         ind(k+(j-1)*10)=k+(j-1)*11
    !     end do
    ! end do
    ! do j=1,10
    !     ind1(j)=j
    !     ind2(j)=j+11
    ! end do
    ! call random_number(d1)
    ! call random_number(d2)

    ! val2=1
    ! call CPU_time(start)
    
    ! do k=1,len 
    !     val2=val2*(d1(k,1)*d2(k,1)+d1(k,2)*d2(k,2))
    ! ! B2(j)=A1(j)*A2(j)
    ! end do
   
    ! call CPU_time(finish)
    ! print*,val2
    ! print*,finish-start

    ! val2=1
    ! call CPU_time(start)
    
    ! do k=1,len 
    !     val2=val2*(t1(k)*t2(k)+t1(k+11)*t2(k+11))
    ! ! B2(j)=A1(j)*A2(j)
    ! end do
   
    ! call CPU_time(finish)
    ! print*,val2
    ! print*,finish-start
   

    ! val2=1
    ! call CPU_time(start)
    
    ! do k=1,len 
    !     val2=val2*(t1(ind1(k))*t2(ind1(k))+t1(ind2(k))*t2(ind2(k)))
    ! ! B2(j)=A1(j)*A2(j)
    ! end do
   
    ! call CPU_time(finish)
    ! print*,val2
    ! print*,finish-start

    ! val2=1
    ! call CPU_time(start)
    
    ! do k=1,len 
    !     val2=val2*(t1(ind1(k))*t2(ind1(k))*neg(k)+t1(ind2(k))*t2(ind2(k)))
    ! end do
   
    ! call CPU_time(finish)
    ! print*,val2
    ! print*,finish-start

    ! ind2(2)=ind1(2)
    ! ind1(2)=0
    ! neg(1)=-1
    ! val2=1
    ! call CPU_time(start)
    
    ! do k=1,len 
    !     val2=val2*(t1(ind1(k))*t2(ind1(k))*neg(k)+t1(ind2(k))*t2(ind2(k)))
    ! end do
   
    ! call CPU_time(finish)
    ! print*,val2
    ! print*,finish-start

    allocate(A1(len),A2(len),B1(len),B2(len),C1(len),c2(len))
    call random_number(A1)
    call random_number(B1)
    call random_number(B2)
    call random_number(A2)
    ops=(/1,3,5,7/)

   
    ! 2.0,1.0,6.0,1.0
    A1 =sin(A1)
    A2 =cos(a2)
    B1 =sin(b1)
    b2=cos(b2)

        b1_1=b1
        b2_1=b2
        val1=b1_1(ops(j))

    C1=a1*b1
    c2=a2*b2
    do j=1,len 

    end do 


    stop
    ! do j=1,20
    !     print*,t1(ind(j))
    ! end do
    ind(13)=ind(2)
    ind(2)=0
    ! do j=1,2
    !     do k=1,len
            ! print*,d1(k,j)
            ! print*,t1((k+(j-1)*11))
    !     end do 
    ! end do


    
    ! do j=1,20
    !     print*,t1(ind(j))
    ! end do
    ! print*,'next'
    ! do j=1,2
    !     do k=1,len
    !         print*,t1(k+(j-1)*10)
    !     end do
    ! end do

    ! allocate(A1(len))
    ! allocate(B1(len))
    ! allocate(A2(len))
    ! allocate(B2(len))
    

    
    ! call random_number(A1)
    ! call random_number(A2)
    ! call random_number(B2)
    ! call random_number(B1)

    ! call CPU_time(start)
    ! val1=product((A1*A2)+(b1*b2))
    ! call CPU_time(finish)
    ! ! print*,val1
    ! print*,finish-start
    ! ! print*,B1(1:len)

    ! val2=1
    ! call CPU_time(start)
    ! do j=1,len
    !     val2=val2*((A1(j)*A2(j))+(B1(j)*B2(j)))
    !     ! B2(j)=A1(j)*A2(j)
    ! end do
    ! call CPU_time(finish)
    ! ! print*,val2
    ! print*,finish-start

    ! val3=1
    ! call CPU_time(start)
    ! allocate(C2(len))
    ! allocate(C1(len))
    ! allocate(C3(len))
    ! do j=1,len
    !     c1(j)=(A1(j)*A2(j))
    !     c2(j)=(B1(j)*B2(j))
    !     c3(j)=(C1(j)+C2(j))
    !     val3=val3*C3(j)
       
    !     ! B2(j)=A1(j)*A2(j)
    ! end do
    ! call CPU_time(finish)
    ! ! print*,val3
    ! print*,finish-start


    ! print*,B2(1:len)
    ! call CPU_time(start)
    ! val1=product((A1(1:4)*B1(1:4))+(A2(1:4)*B2(1:4)))*product((A1(6:11)*B1(6:11))+(A2(6:11)*B2(6:11)))&
    ! *product((A1(13:100)*B1(13:100))+(A2(13:100)*B2(13:100)))*A1(12)*B2(12)*A2(5)*B1(5)
    ! call CPU_time(finish)
    ! print*,finish-start
    ! print*,val1


    ! ! val2=(A1(1)*B1(1)+A2(1)*B2(1))*(A1(2)*B1(2)+A2(2)*B2(2))*(A1(3)*B1(3)+A2(3)*B2(3))*&
    ! ! (A1(4)*B1(4)+A2(4)*B2(4))*(A1(6)*B1(6)+A2(6)*B2(6))*(A1(7)*B1(7)+A2(7)*B2(7))*(A1(8)*B1(8)+A2(8)*B2(8))*&
    ! ! (A1(9)*B1(9)+A2(9)*B2(9))*(A1(10)*B1(10)+A2(10)*B2(10))*A2(5)*B1(5)*(A1(11)*B1(11)+A2(11)*B2(11))*&
    ! ! (A1(12)*B1(12)+A2(12)*B2(12))*(A1(13)*B1(13)+A2(13)*B2(13))*(A1(15)*B1(15)+A2(15)*B2(15))*&
    ! ! (A1(14)*B1(14)+A2(14)*B2(14))*(A1(16)*B1(16)+A2(16)*B2(16))*(A1(17)*B1(17)+A2(17)*B2(17))*(A1(18)*B1(18)+A2(18)*B2(18))&
    ! ! *(A1(19)*B1(19)+A2(19)*B2(19))*(A1(20)*B1(20)+A2(20)*B2(20))
    ! call CPU_time(start)
    ! val2=&
    ! (A1(1)*B1(1)+A2(1)*B2(1))*(A1(2)*B1(2)+A2(2)*B2(2))*(A1(3)*B1(3)+A2(3)*B2(3))*(A1(4)*B1(4)+A2(4)*B2(4))*&
    ! (A1(6)*B1(6)+A2(6)*B2(6))*(A1(7)*B1(7)+A2(7)*B2(7))*(A1(8)*B1(8)+A2(8)*B2(8))*(A1(9)*B1(9)+A2(9)*B2(9))*&
    ! (A1(10)*B1(10)+A2(10)*B2(10))*(A1(11)*B1(11)+A2(11)*B2(11))*&
    ! (A1(13)*B1(13)+A2(13)*B2(13))*(A1(14)*B1(14)+A2(14)*B2(14))*(A1(15)*B1(15)+A2(15)*B2(15))*&
    ! (A1(16)*B1(16)+A2(16)*B2(16))*(A1(17)*B1(17)+A2(17)*B2(17))*(A1(18)*B1(18)+A2(18)*B2(18))*&
    ! (A1(19)*B1(19)+A2(19)*B2(19))*(A1(20)*B1(20)+A2(20)*B2(20))*(A1(21)*B1(21)+A2(21)*B2(21))*&
    ! (A1(22)*B1(22)+A2(22)*B2(22))*(A1(23)*B1(23)+A2(23)*B2(23))*(A1(24)*B1(24)+A2(24)*B2(24))*&
    ! (A1(25)*B1(25)+A2(25)*B2(25))*(A1(26)*B1(26)+A2(26)*B2(26))*(A1(27)*B1(27)+A2(27)*B2(27))*&
    ! (A1(28)*B1(28)+A2(28)*B2(28))*(A1(29)*B1(29)+A2(29)*B2(29))*(A1(30)*B1(30)+A2(30)*B2(30))*&
    ! (A1(31)*B1(31)+A2(31)*B2(31))*(A1(32)*B1(32)+A2(32)*B2(32))*(A1(33)*B1(33)+A2(33)*B2(33))*&
    ! (A1(34)*B1(34)+A2(34)*B2(34))*(A1(35)*B1(35)+A2(35)*B2(35))*(A1(36)*B1(36)+A2(36)*B2(36))*&
    ! (A1(37)*B1(37)+A2(37)*B2(37))*(A1(38)*B1(38)+A2(38)*B2(38))*(A1(39)*B1(39)+A2(39)*B2(39))*&
    ! (A1(40)*B1(40)+A2(40)*B2(40))*(A1(41)*B1(41)+A2(41)*B2(41))*(A1(42)*B1(42)+A2(42)*B2(42))*&
    ! (A1(43)*B1(43)+A2(43)*B2(43))*(A1(44)*B1(44)+A2(44)*B2(44))*(A1(45)*B1(45)+A2(45)*B2(45))*&
    ! (A1(46)*B1(46)+A2(46)*B2(46))*(A1(47)*B1(47)+A2(47)*B2(47))*(A1(48)*B1(48)+A2(48)*B2(48))*&
    ! (A1(49)*B1(49)+A2(49)*B2(49))*(A1(50)*B1(50)+A2(50)*B2(50))*(A1(51)*B1(51)+A2(51)*B2(51))*&
    ! (A1(52)*B1(52)+A2(52)*B2(52))*(A1(53)*B1(53)+A2(53)*B2(53))*(A1(54)*B1(54)+A2(54)*B2(54))*&
    ! (A1(55)*B1(55)+A2(55)*B2(55))*(A1(56)*B1(56)+A2(56)*B2(56))*(A1(57)*B1(57)+A2(57)*B2(57))*&
    ! (A1(58)*B1(58)+A2(58)*B2(58))*(A1(59)*B1(59)+A2(59)*B2(59))*(A1(60)*B1(60)+A2(60)*B2(60))*&
    ! (A1(61)*B1(61)+A2(61)*B2(61))*(A1(62)*B1(62)+A2(62)*B2(62))*(A1(63)*B1(63)+A2(63)*B2(63))*&
    ! (A1(64)*B1(64)+A2(64)*B2(64))*(A1(65)*B1(65)+A2(65)*B2(65))*(A1(66)*B1(66)+A2(66)*B2(66))*&
    ! (A1(67)*B1(67)+A2(67)*B2(67))*(A1(68)*B1(68)+A2(68)*B2(68))*(A1(69)*B1(69)+A2(69)*B2(69))*&
    ! (A1(70)*B1(70)+A2(70)*B2(70))*(A1(71)*B1(71)+A2(71)*B2(71))*&
    ! (A1(72)*B1(72)+A2(72)*B2(72))*(A1(73)*B1(73)+A2(73)*B2(73))*(A1(74)*B1(74)+A2(74)*B2(74))*&
    ! (A1(75)*B1(75)+A2(75)*B2(75))*(A1(76)*B1(76)+A2(76)*B2(76))*(A1(77)*B1(77)+A2(77)*B2(77))*&
    ! (A1(78)*B1(78)+A2(78)*B2(78))*(A1(79)*B1(79)+A2(79)*B2(79))*(A1(80)*B1(80)+A2(80)*B2(80))*&
    ! (A1(81)*B1(81)+A2(81)*B2(81))*(A1(82)*B1(82)+A2(82)*B2(82))*(A1(83)*B1(83)+A2(83)*B2(83))*&
    ! (A1(84)*B1(84)+A2(84)*B2(84))*(A1(85)*B1(85)+A2(85)*B2(85))*(A1(86)*B1(86)+A2(86)*B2(86))*&
    ! (A1(87)*B1(87)+A2(87)*B2(87))*(A1(88)*B1(88)+A2(88)*B2(88))*(A1(89)*B1(89)+A2(89)*B2(89))*&
    ! (A1(90)*B1(90)+A2(90)*B2(90))*(A1(91)*B1(91)+A2(91)*B2(91))*&
    ! (A1(92)*B1(92)+A2(92)*B2(92))*(A1(93)*B1(93)+A2(93)*B2(93))*(A1(94)*B1(94)+A2(94)*B2(94))*&
    ! (A1(95)*B1(95)+A2(95)*B2(95))*(A1(96)*B1(96)+A2(96)*B2(96))*(A1(97)*B1(97)+A2(97)*B2(97))*&
    ! (A1(98)*B1(98)+A2(98)*B2(98))*(A1(99)*B1(99)+A2(99)*B2(99))*(A1(100)*B1(100)+A2(100)*B2(100))*&
    ! A2(5)*B1(5)*(A1(12)*B2(12))
    ! call CPU_time(finish)
    ! print*,finish-start
    ! print*,val2


end program test