::This is a windows .bat file which runs the ac_ipl unit tests using the visual studio developer command prompt 

::Invoke the visual studio developer command prompt
call "C:\Program Files (x86)\Microsoft Visual Studio\2022\BuildTools\Common7\Tools\VsDevCmd.bat"

::Set the path for Mgc_home
set MGC_HOME=C:\MGCNoScan\abeemana\sb\sif\ixn\Mgc_home

::All the unit tests
set SOURCES_CPP=rtest_ac_fir_const_coeffs;^
rtest_ac_fir_load_coeffs;^
rtest_ac_fir_prog_coeffs;^
rtest_ac_cic_dec_full;^
rtest_ac_cic_intr_full

::Compile and execute the unit tests
(for %%a in (%SOURCES_CPP%) do ( 
   cl /EHsc /I%MGC_HOME%\shared\include %%a.cpp
   %%a.exe
))
