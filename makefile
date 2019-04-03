all:  docs install

docs:
	R -e "devtools::document()"
build:
	(cd ..; R CMD build --no-build-vignettes TrenaProjectLiver)

install:
	(cd ..; R CMD INSTALL TrenaProjectLiver)

check:
	(cd ..; R CMD check `ls -t TrenaProjectLiver) | head -1`)

test:
	for x in inst/unitTests/test_*.R; do echo $$x; R -f $$x; done

