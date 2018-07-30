program make_mesh_for_gnuplot

character*100 :: input,output
integer :: i,j,k,i1,i2,i3
real*8,allocatable :: meshpoints(:,:)
integer,allocatable :: meshlines(:,:), meshangles(:)

call getarg(1,input)
output=trim(input)//"_triangle"

print *,"outputfile: ", output

open(10,file=trim(output))

open(11,file=trim(input))
	print *,"Opening succeeded"
  read(11,*,err=110,end=110)
  read(11,*,err=110,end=110) i1, i2, i3				!i1 = ncels, i2 = number of areas, i3 = ncels
    allocate(meshpoints(i1,7)); meshpoints=0d0
    allocate(meshlines(i2,4)); meshlines=0
    allocate(meshangles(i2)); meshangles=0
  read(11,*,err=110,end=110)
  print*, i1, i2, i3
  do i=1,i1	!for all cells..
    read(11,*,err=110,end=110) meshpoints(i,:) !..read the coordinates (and other values)
  end do
  print*, "##"
  read(11,*,err=110,end=110)
  do i=1,i2	  				! for all areas
    read(11,*,err=110,end=110) meshangles(i)  ! read the number of edges (Dreieck = 3, Viereck = 4)
    print*, i, meshangles(i)
  end do
110 close(11)
print*, "#"
open(12,file=trim(input)) !reopen the file
  read(12,*,err=111,end=111)
  read(12,*,err=111,end=111) 
  read(12,*,err=111,end=111)
  do i=1,i1   !for all cells
    read(12,*,err=111,end=111) meshpoints(i,:) !read the coordinates (and other values)
  end do
  read(12,*,err=111,end=111)
  do i=1,i2   !for each area
    read(12,*,err=111,end=111) k, meshlines(i,:meshangles(i))  !k = number of edges, meshlines wird gefüllt mit IDs von Eckpunkten
    do j=1,k		!for each edge
      write(10,*) i, j, meshpoints(meshlines(i,j)+1,:)
    end do
    write(10,*) i, j, meshpoints(meshlines(i,1)+1,:)
    write(10,*) ""
  end do
111 close(12)
close(10)

end program
! splot './999out____.off_triangle' us 3:4:5 w linespoints lc "grey", './999out____.off_triangle' us 3:4:5:5 pal pt 7 w p
