


norb=10
ndet=10
h1count=0
h2count=0

standard="(z(j)%sin({0})*z(k)%sin({0})+z(j)%cos({0})*z(k)%cos({0}))"
negative="(z(j)%cos({0})*z(k)%cos({0})-z(j)%sin({0})*z(k)%sin({0}))"
cr="z(j)%sin({0})*z(k)%cos({0})"
an="z(j)%cos({0})*z(k)%sin({0})"
cran="z(j)%sin({0})*z(k)%sin({0})"

test=standard+negative+"\n"+an

with open('tester.f95','w',encoding="utf-8") as f:
    f.write("MODULE ham2 \n")
    f.write(" \n")
    f.write("   use globvars \n")
    f.write(" \n")
    f.write("   contains \n")
    f.write(" \n")
    f.write(" \n")
    f.write("   subroutine ham_make(ham,z,elecs) \n")
    f.write(" \n")
    f.write("       implicit none\n")
    f.write(" \n")
    f.write("       type(hamiltonian), intent(inout)::ham \n")
    f.write("       type(zombiest),dimension(:),intent(in)::z\n")
    f.write("       type(elecintrgl),intent(in)::elecs \n")
    f.write("       real(kind=8),allocatable,dimension(:)::h1etot,h2etot \n")
    f.write("       integer::j,k,ierr\n")
    f.write(" \n")
    f.write("       if (errorflag .ne. 0) return \n")
    f.write(" \n")
    f.write("       ierr=0\n")
    f.write("       allocate(h1etot("+str(h1count)+"),stat=ierr)\n")
    f.write("       allocate(h2etot("+str(h2count)+"),stat=ierr)\n")
    f.write("       if (ierr/=0) then\n" )
    f.write('           write(0,"(a,i0)") "Error in annihilation and creation array vector allocation. ierr had value ," ierr \n')
    f.write("           errorflag=1\n")
    f.write("           return\n")
    f.write("       end if\n")
    f.write(" \n")
    f.write("       do j=1,ndet\n")
    f.write("           do k=j,ndet\n")
    f.write(" \n")
    f.write(" \n")
    f.write(" \n")
    f.write(" \n")
    f.write(" \n")
    f.write(" \n")
    f.write(" \n")
    f.write(" \n")
    f.write(" \n")
    f.write(" \n")
    f.write("   end subroutine ham_make\n")
    f.write(" \n")
    f.write("END MODULE ham2 \n")


