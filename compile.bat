set PATH=C:\MinGW\bin;%PATH%
gfortran -Og -shared DllTestSub.f90 -o add.dll
gfortran -ffree-line-length-512 MainDllUse.f90 add.dll -o main.exe
main.exe