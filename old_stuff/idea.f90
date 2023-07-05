Program idea

    implicit none 

    real, target,dimension(2)::a,b
    real, pointer,dimension():: pa1, pa2
    ! real::b

    pa1 =>a
    pa2 => b
    a=4
    print*,pa1
    
    print*,b
    a=2
    pa1=2
    print*,pa1
    print*,b

end program idea