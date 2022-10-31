l_b = -Wl,--start-group ${MKLROOT}/lib/intel64/libmkl_intel_ilp64.a ${MKLROOT}/lib/intel64/libmkl_intel_thread.a ${MKLROOT}/lib/intel64/libmkl_core.a -Wl,--end-group -liomp5 -lpthread -lm -ldl

program : main.o CSC_Hamiltonian.o Jset.o EIGEN.o
	icpc -o $@ $^ $(l_b)

main.o : main.cpp
	icpc -c $< $(l_b)

CSC_Hamiltonian.o : CSC_Hamiltonian/CSC_Hamiltonian.cpp
	icpc -c $< $(l_b)

Jset.o : Jset/Jset.cpp
	icpc -c $< $(l_b)
	
EIGEN.o : EIGEN/EIGEN.cpp
	icpc -c $< $(l_b)

run : program
	./program

clean:
	rm -f ./program

.PHONY : run clean