FC = gfortran
src_files = app/main.f90

program = build/2d_lambda_phi4
FLAGS = -c

$(program) : build/obj/precision.o build/obj/parameters.o build/obj/pbc.o build/obj/functions.o build/obj/field.o build/obj/observables.o build/obj/dynamics.o build/obj/main.o
	$(FC) $^ -o $@

build/obj/%.o: src/%.f90 build/mod build/obj
	$(FC) $(FLAGS) -J build/mod $< -o $@

build/obj/main.o: app/main.f90 
	$(FC) $(FLAGS) -I build/mod $< -o $@

build:
	mkdir -p $@

build/obj: build
	mkdir -p $@

build/mod: build
	mkdir -p $@

run :
	@{ echo "input/input_parameters.nml"; } | $(program)

clean:
	rm -r build

readme:
	pandoc -f gfm README.md -o doc/README.pdf
