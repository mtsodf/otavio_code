program MAIN

external addtwo, addthree


integer X

X = 10

call addtwo(X)

write(*,*) "X = ", X

call addthree(X)

write(*,*) "X = ", X


end program